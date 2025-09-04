/*-------------------------------------------------------------------------------------------------*
 *  DFTB+: general package for performing fast atomistic simulations                               *
 *  Copyright (C) 2006 - 2025  DFTB+ developers group                                              *
 *                                                                                                 *
 *  See the LICENSE file for terms of usage and distribution.                                      *
 *-------------------------------------------------------------------------------------------------*/
#include <assert.h>
#include <cuda_runtime.h>
#include <omp.h>
#include <thrust/complex.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <memory>
#include <vector>

#include "kernel.cuh"
#include "slater.cuh"
#include "spharmonics.cuh"
#include "utils.cuh"

// more print statements
constexpr bool debug = false;
// Avoid division by zero
constexpr double INV_R_EPSILON = 1.0e-12;
// amount of shared memory set aside for nEig accumulators
constexpr float SHARED_MEM_FACTOR = 0.95f;
// max output array share of free global memory
constexpr float GLOBAL_MEM_FACTOR = 0.80f;
// Threads per block, multiple of warp size 32
constexpr int block_size = 256;

using complexd = thrust::complex<double>;

// Manages Gpu memory allocation and H2D copy
struct DeviceData {
    // Grid
    DeviceBuffer<double> origin;
    DeviceBuffer<double> gridVecs;

    // System
    DeviceBuffer<double> coords;
    DeviceBuffer<int>    species;
    DeviceBuffer<int>    iStos;

    // Periodic
    DeviceBuffer<double>   latVecs;
    DeviceBuffer<double>   recVecs2pi;
    DeviceBuffer<int>      kIndexes;
    DeviceBuffer<complexd> phases;

    // STO Basis
    DeviceBuffer<int>    sto_angMoms;
    DeviceBuffer<int>    sto_nPows;
    DeviceBuffer<int>    sto_nAlphas;
    DeviceBuffer<double> sto_cutoffsSq;
    DeviceBuffer<double> sto_coeffs;
    DeviceBuffer<double> sto_alphas;
    // Texture for radial LUT
    std::unique_ptr<GpuLutTexture> sto_lut;

    // Eigenvectors
    DeviceBuffer<double>   eigVecsReal;
    DeviceBuffer<complexd> eigVecsCmpl;

    // Output (per-GPU batch buffer)
    DeviceBuffer<complexd> d_valueCmpl_out_batch;
    DeviceBuffer<double>   d_valueReal_out_batch;

    // Constructor handles all H2D allocation and copy
    DeviceData(const GridParams* grid, const SystemParams* system, const PeriodicParams* periodic,
        const StoBasisParams* basis, const CalculationParams* calc, int z_per_batch)
        : origin(grid->origin, 3),
          gridVecs(grid->gridVecs, 9),
          coords(system->coords, (size_t)3 * system->nAtom * system->nCell),
          species(system->species, system->nAtom),
          iStos(system->iStos, system->nSpecies + 1),
          sto_angMoms(basis->angMoms, basis->nStos),
          sto_cutoffsSq(basis->cutoffsSq, basis->nStos) {
        if (basis->useRadialLut) {
            if (debug) printf("Using radial LUT with %d points for %d STOs\n", basis->nLutPoints, basis->nStos);
            sto_lut = std::unique_ptr<GpuLutTexture>(
                new GpuLutTexture(basis->lutGridValues, basis->nLutPoints, basis->nStos));
        } else {
            if (debug) printf("Using direct STO evaluation for %d STOs\n", basis->nStos);
            sto_nPows.assign(basis->nPows, basis->nStos);
            sto_nAlphas.assign(basis->nAlphas, basis->nStos);
            sto_coeffs.assign(basis->coeffs, (size_t)basis->maxNPows * basis->maxNAlphas * basis->nStos);
            sto_alphas.assign(basis->alphas, (size_t)basis->maxNAlphas * basis->nStos);
        }

        if (calc->isRealInput) {
            eigVecsReal.assign(calc->eigVecsReal, (size_t)system->nOrb * calc->nEigIn);
        } else {
            eigVecsCmpl.assign(
                reinterpret_cast<const complexd*>(calc->eigVecsCmpl), (size_t)system->nOrb * calc->nEigIn);
            phases.assign(reinterpret_cast<const complexd*>(periodic->phases), (size_t)system->nCell * calc->nEigIn);
            kIndexes.assign(periodic->kIndexes, calc->nEigIn);
        }
        if (periodic->isPeriodic) {
            latVecs.assign(periodic->latVecs, 9);
            recVecs2pi.assign(periodic->recVecs2pi, 9);
        }

        // Per-GPU batch buffer for the output
        size_t batch_buffer_size_elems = (size_t)grid->nPointsX * grid->nPointsY * z_per_batch * calc->nEigOut;
        if (calc->isRealOutput) {
            d_valueReal_out_batch = DeviceBuffer<double>(batch_buffer_size_elems);
        } else {
            d_valueCmpl_out_batch = DeviceBuffer<complexd>(batch_buffer_size_elems);
        }
    }
};

// Kernel parameters struct to simplify the argument list
struct DeviceKernelParams {
    // Grid
    int           nPointsX, nPointsY, z_per_batch, z_offset_global;
    const double* origin;
    const double* gridVecs;

    // System
    int nAtom, nCell, nOrb;

    const double* coords;
    const int*    species;
    const int*    iStos;

    // Periodic boundary cond.
    bool            isPeriodic;
    const double*   latVecs;
    const double*   recVecs2pi;
    const int*      kIndexes;
    const complexd* phases;

    // STO Basis
    int nStos, maxNPows, maxNAlphas;
    // Texture LUTs
    cudaTextureObject_t lutTex;
    double              inverseLutStep;
    // STO parameters
    const int*    sto_angMoms;
    const int*    sto_nPows;
    const int*    sto_nAlphas;
    const double* sto_cutoffsSq;
    const double* sto_coeffs;
    const double* sto_alphas;

    // Eigenvectors
    int             nEig, nEig_per_pass;
    const double*   eigVecsReal;
    const complexd* eigVecsCmpl;

    // Output (batch pointers)
    double*   valueReal_out_batch;
    complexd* valueCmpl_out_batch;

    // Constructor to initialize the parameters from host data
    // Batch-specific parameters are initialized to zero or nullptr,
    // and need to be set in the loop before kernel launch.
    DeviceKernelParams(DeviceData& data, const GridParams* grid, const SystemParams* system,
        const PeriodicParams* periodic, const StoBasisParams* basis, const CalculationParams* calc,
        int nEig_per_pass_in) {
        // Grid
        origin   = data.origin.get();
        gridVecs = data.gridVecs.get();
        nPointsX = grid->nPointsX;
        nPointsY = grid->nPointsY;

        // System
        nAtom   = system->nAtom;
        nCell   = system->nCell;
        nOrb    = system->nOrb;
        coords  = data.coords.get();
        species = data.species.get();
        iStos   = data.iStos.get();

        // STO Basis
        nStos         = basis->nStos;
        sto_angMoms   = data.sto_angMoms.get();
        sto_cutoffsSq = data.sto_cutoffsSq.get();

        if (basis->useRadialLut) {
            lutTex         = data.sto_lut->get();
            inverseLutStep = basis->inverseLutStep;
        } else {
            maxNPows    = basis->maxNPows;
            maxNAlphas  = basis->maxNAlphas;
            sto_nPows   = data.sto_nPows.get();
            sto_nAlphas = data.sto_nAlphas.get();
            sto_coeffs  = data.sto_coeffs.get();
            sto_alphas  = data.sto_alphas.get();
        }

        // Periodic boundary conditions
        isPeriodic = periodic->isPeriodic;
        latVecs    = data.latVecs.get();
        recVecs2pi = data.recVecs2pi.get();
        kIndexes   = data.kIndexes.get();
        phases     = data.phases.get();

        // Eigenvectors
        nEig        = calc->nEigIn;
        eigVecsReal = data.eigVecsReal.get();
        eigVecsCmpl = data.eigVecsCmpl.get();

        // Output batch buffers
        valueReal_out_batch = data.d_valueReal_out_batch.get();
        valueCmpl_out_batch = data.d_valueCmpl_out_batch.get();

        nEig_per_pass = nEig_per_pass_in;
        // Batch-specific kernel config to be updated in the loop
        z_per_batch     = 0;
        z_offset_global = 0;
    }
};

// ================================================================================================
//  CUDA Kernel.
// ================================================================================================
// To avoid branching (dropped at compile time), we template the kernel 16 ways on boolean flags:
// <isRealInput> decides whether real/complex eigenvectors are used (and adds phases).
// <useRadialLut> decides whether to use texture memory interpolation for the STO radial functions.
// <isPeriodic> enables folding of coords into the unit cell.
// <calcAtomicDensity> squares the basis wavefunction contributions.
//   In this case, the occupation should be passed as the eigenvector.
// <calcTotalChrg> accumulates the density over all states in valueReal_out of shape (x,y,z,1).
//   Here, occupation should be passed by multiplying the eigenvecs with sqrt(occupation).
template <bool isRealInput, bool calcAtomicDensity, bool calcTotalChrg, bool useRadialLut>
__global__ void evaluateKernel(const DeviceKernelParams p) {
    using AccumT = typename std::conditional<(isRealInput), double, complexd>::type;

    // Each thread gets its own private slice of the shared memory buffer for fast accumulation.
    // We have to chunk the eigenstates into nEig_per_pass due to size constraints.
    // (Cuda doesnt allow templating the shared memory type, so we simply recast it.)
    extern __shared__ char shared_workspace[];
    AccumT* point_results_pass = reinterpret_cast<AccumT*>(shared_workspace) + threadIdx.x * p.nEig_per_pass;

    // --- Thread to point mapping ---
    // Map each thread to unique 1d index
    int idx_in_batch          = blockIdx.x * blockDim.x + threadIdx.x;
    int total_points_in_batch = p.nPointsX * p.nPointsY * p.z_per_batch;
    if (idx_in_batch >= total_points_in_batch) return;

    // Map 1d index to point in grid
    int i1        = idx_in_batch % p.nPointsX;
    int i2        = (idx_in_batch / p.nPointsX) % p.nPointsY;
    int i3_batch  = idx_in_batch / (p.nPointsX * p.nPointsY);
    int i3_global = i3_batch + p.z_offset_global;

    // Map point to global coordinates.
    double xyz[3];
    for (int i = 0; i < 3; ++i)
        xyz[i] = p.origin[i] + i1 * p.gridVecs[IDX2F(i, 0, 3)]
                      +        i2 * p.gridVecs[IDX2F(i, 1, 3)]
                      + i3_global * p.gridVecs[IDX2F(i, 2, 3)];

    // If periodic, fold into cell by discarding the non-fractional part in lattice vector
    // multiples.
    if (p.isPeriodic)
        foldCoordsIntoCell(
            xyz, reinterpret_cast<const double(*)[3]>(p.latVecs), reinterpret_cast<const double(*)[3]>(p.recVecs2pi));

    double totChrgAcc = 0.0;
    // --- Loop over eigenstates in chunks that fit in shared memory ---
    for (int eig_base = 0; eig_base < p.nEig; eig_base += p.nEig_per_pass) {
        // Initialize the small, per-pass buffer for this thread
        for (int i = 0; i < p.nEig_per_pass; ++i) {
            point_results_pass[i] = AccumT(0.0);
        }

        // Since we run out of space in point_result_pass[], the spatial calculation
        // is repeated for each chunk of eigenstates.
        // This is to keep the accumulation in fast shared memory.
        for (int iCell = 0; iCell < p.nCell; ++iCell) {
            int orbital_idx_counter = 0;
            for (int iAtom = 0; iAtom < p.nAtom; ++iAtom) {
                int iSpecies = p.species[iAtom] - 1;

                double diff[3];
                for (int i = 0; i < 3; ++i) {
                    diff[i] = xyz[i] - p.coords[IDX3F(i, iAtom, iCell, 3, p.nAtom)];
                }
                double rr = diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2];

                for (int iOrb = p.iStos[iSpecies] - 1; iOrb < p.iStos[iSpecies + 1] - 1; ++iOrb) {
                    int iL = p.sto_angMoms[iOrb];
                    if (rr > p.sto_cutoffsSq[iOrb]) {
                        orbital_idx_counter += 2 * iL + 1;
                        continue;
                    }
                    double r = sqrt(rr);

                    double radialVal;
                    if constexpr (useRadialLut) {
                        double lut_pos = 0.5f + r * p.inverseLutStep; // Add 0.5 to adress texel center (imagine pixels)
                        radialVal      = (double)tex2D<float>(p.lutTex, lut_pos, (float)iOrb + 0.5f);
                    } else {
                        radialVal = getRadialValue(r, iL, iOrb, p.sto_nPows[iOrb], p.sto_nAlphas[iOrb], p.sto_coeffs,
                            p.sto_alphas, p.maxNPows, p.maxNAlphas);
                    }

                    // precompute inverse used across several realTessY calls
                    double inv_r  = (r < INV_R_EPSILON) ? 0.0 : 1.0 / r;

                    for (int iM = -iL; iM <= iL; ++iM) {
                        double val = radialVal * realTessY(iL, iM, diff, inv_r);
                        if constexpr (calcAtomicDensity) val = val * val;

                        // Accumulate into the small shared memory buffer for the current chunk
                        for (int iEig_offset = 0; iEig_offset < p.nEig_per_pass; ++iEig_offset) {
                            int iEig = eig_base + iEig_offset;
                            if (iEig >= p.nEig) break;  // Don't go past the end on the last chunk
                            size_t eig_idx = IDX2F(orbital_idx_counter, iEig, p.nOrb);
                            if constexpr (isRealInput) {
                                point_results_pass[iEig_offset] += val * p.eigVecsReal[eig_idx];
                            } else {
                                point_results_pass[iEig_offset] +=
                                    val * p.phases[IDX2F(iCell, iEig, p.nCell)] * p.eigVecsCmpl[eig_idx];
                            }
                        }
                        orbital_idx_counter++;
                    }
                }
            }
        }

        // Write the complete nEig_per_pass chunk to global memory.
        for (int iEig_offset = 0; iEig_offset < p.nEig_per_pass; ++iEig_offset) {
            int iEig = eig_base + iEig_offset;
            if (iEig >= p.nEig) break;
            size_t out_idx = IDX4F(i1, i2, i3_batch, iEig, p.nPointsX, p.nPointsY, p.z_per_batch);
            if constexpr (isRealInput) {
                if constexpr (calcTotalChrg)
                    totChrgAcc += point_results_pass[iEig_offset] * point_results_pass[iEig_offset];
                else
                    p.valueReal_out_batch[out_idx] = point_results_pass[iEig_offset];

            } else {
                if constexpr (calcTotalChrg)
                    totChrgAcc += thrust::norm(point_results_pass[iEig_offset]);
                else
                    p.valueCmpl_out_batch[out_idx] = point_results_pass[iEig_offset];
            }
        }
    }

    // Density stored in first eig : (x,y,z,1)
    if constexpr (calcTotalChrg) {
        size_t out_idx = IDX4F(i1, i2, i3_batch, 0, p.nPointsX, p.nPointsY, p.z_per_batch);

        p.valueReal_out_batch[out_idx] = totChrgAcc;
    }
}

struct GpuLaunchConfig {
    int    deviceId;
    int    z_start;              // Starting Z-slice for this GPU
    int    z_count;              // Number of Z-slices for this GPU
    int    z_per_batch;          // Number of Z-slices to process per kernel launch
    int    nEig_per_pass;        // Number of eigenstates to accumulate in shared memory
    size_t shared_mem_for_pass;  // Amount of shared memory per block for the accumulators
};

GpuLaunchConfig setup_gpu_config(int deviceId, int numGpus, const GridParams* grid, const CalculationParams* calc) {
    GpuLaunchConfig config;
    config.deviceId = deviceId;

    // Evenly divide work using Z-slices among GPUs
    config.z_count = grid->nPointsZ / numGpus;
    config.z_start = deviceId * config.z_count;
    // Handle uneven Z-slice count: Last GPU takes remaining slices
    if (deviceId == numGpus - 1) config.z_count = grid->nPointsZ - config.z_start;

    // Query available memory sizes with safety margins
    cudaDeviceProp prop;
    size_t         free_global_mem, total_global_mem;
    CHECK_CUDA(cudaMemGetInfo(&free_global_mem, &total_global_mem));
    CHECK_CUDA(cudaGetDeviceProperties(&prop, deviceId));
    size_t available_shared = prop.sharedMemPerBlock * SHARED_MEM_FACTOR;
    size_t available_global = static_cast<size_t>(free_global_mem * GLOBAL_MEM_FACTOR);

    size_t accumulator_number_size = calc->isRealInput ? sizeof(double) : sizeof(complexd);
    size_t output_number_size      = calc->isRealOutput ? sizeof(double) : sizeof(complexd);

    // Determine number of eigenstates that fit into shared memory
    config.nEig_per_pass       = available_shared / (block_size * accumulator_number_size);
    config.nEig_per_pass       = std::min(calc->nEigIn, std::max(1, config.nEig_per_pass));  // clamp to [1, nEigIn]
    config.shared_mem_for_pass = (size_t)config.nEig_per_pass * block_size * accumulator_number_size;

    // Determine max Z-slices that can fit in available (global) memory
    size_t bytes_per_slice = (size_t)grid->nPointsX * grid->nPointsY * calc->nEigOut * output_number_size;
    config.z_per_batch     = std::min(config.z_count, (int)available_global / (int)bytes_per_slice);
    config.z_per_batch     = std::max(1, config.z_per_batch);  // at least 1

    // Debug output
    if (deviceId == 0 && debug) {
        printf("\n--- GPU %d Configuration ---\n", deviceId);
        printf("  Z-slice workload: %d (from index %d to %d)\n", config.z_count, config.z_start,
            config.z_start + config.z_count - 1);
        printf("  Block size: %d threads, %zub shared mem per block, %d eigs of %d per pass\n", block_size,
            config.shared_mem_for_pass, config.nEig_per_pass, calc->nEigIn);
        size_t total_size_valueOut =
            (size_t)grid->nPointsX * grid->nPointsY * grid->nPointsZ * calc->nEigOut * sizeof(double);
        if (!calc->isRealOutput) total_size_valueOut *= 2;
        printf(" (Free device mem: %.2f GB, Grid size: %d x %d x %d (x %d eigs) = %.2f GB)\n", free_global_mem / 1e9,
            grid->nPointsX, grid->nPointsY, grid->nPointsZ, calc->nEigOut, total_size_valueOut / 1e9);
        printf("  Processing Z-slices in batches of %d\n", config.z_per_batch);
    }

    return config;
}
// Since the kernel is templated on the different calculation modes, we cannot simply pass booleans
// at runtime and thus need this dispatch table to call the correct binary.
void dispatchKernel(const DeviceKernelParams* params, bool isRealInput, bool calcAtomicDensity, bool calcTotalChrg,
    bool useRadialLut, int grid_size, size_t shared_mem_for_pass) {
#define CALL_KERNEL(isReal, doAtomic, doChrg, useLut) \
    evaluateKernel<isReal, doAtomic, doChrg, useLut><<<grid_size, block_size, shared_mem_for_pass>>>(*params);

    int idx = (isRealInput ? 1 : 0) + (calcAtomicDensity ? 2 : 0) + (calcTotalChrg ? 4 : 0) + (useRadialLut ? 8 : 0);

    assert(!(calcAtomicDensity && calcTotalChrg));
    assert(!(!isRealInput && calcAtomicDensity));

    switch (idx) {
        case 0:  CALL_KERNEL(false, false, false, false); break;
        case 1:  CALL_KERNEL(true,  false, false, false); break;
        case 3:  CALL_KERNEL(true,  true,  false, false); break;
        case 4:  CALL_KERNEL(false, false, true,  false); break;
        case 5:  CALL_KERNEL(true,  false, true,  false); break;
        case 8:  CALL_KERNEL(false, false, false, true); break;
        case 9:  CALL_KERNEL(true,  false, false, true); break;
        case 11: CALL_KERNEL(true,  true,  false, true); break;
        case 12: CALL_KERNEL(false, false, true,  true); break;
        case 13: CALL_KERNEL(true,  false, true,  true); break;
        default: fprintf(stderr, "Error: invalid kernel configuration index %d\n", idx); exit(EXIT_FAILURE);
    }
#undef CALL_KERNEL
}

// Copy the computed batch back to the correct slice of the final host array.
// The D2H copy will automatically block/ synchronize the kernel for this batch.
// This could be improved by using streams / cudaMemcpyAsync, but currently is not a bottleneck.
// The output array is of fortran shape (x,y,z,nEigOut), thus z-slices are not contiguous.
// We slice on Z instead of nEigOut to save on a little computation in the kernel.
void copyD2H(void* d_src_ptr, void* h_dest_ptr, int nPointsX, int nPointsY, int nPointsZ, int z_per_batch,
    int z_offset_global, const CalculationParams* calc) {
    size_t output_num_size = calc->isRealOutput ? sizeof(double) : sizeof(complexd);
    size_t host_plane_size    = (size_t)nPointsZ * nPointsY * nPointsX * output_num_size;
    size_t device_plane_size  = (size_t)z_per_batch * nPointsY * nPointsX * output_num_size;

    // Memcpy3D could be used to squash this loop.
    for (int iEig = 0; iEig < calc->nEigOut; ++iEig) {
        // From: iEig-th slice of GPU batch buffer
        ptrdiff_t d_offset_bytes = (ptrdiff_t)(iEig * device_plane_size);

        // To: Global Z-position in the iEig-th slice of host buffer
        ptrdiff_t h_offset_bytes = (ptrdiff_t)(iEig * host_plane_size + ((size_t)z_offset_global * nPointsY * nPointsX) * output_num_size);

        CHECK_CUDA(cudaMemcpy((char*)h_dest_ptr + h_offset_bytes, (char*)d_src_ptr + d_offset_bytes, device_plane_size,
            cudaMemcpyDeviceToHost));
    }
}

// =========================================================================
//  C++ Host Interface (callable from C/Fortran)
// =========================================================================
extern "C" void evaluate_on_device_c(const GridParams* grid, const SystemParams* system, const PeriodicParams* periodic,
    const StoBasisParams* basis, const CalculationParams* calc) {
    if (calc->nEigIn * grid->nPointsX * grid->nPointsY * grid->nPointsZ == 0) {
        fprintf(stderr, "Error: Zero-sized dimension in input parameters.\n");
        exit(EXIT_FAILURE);
    }
    if (calc->calcTotalChrg) {
        assert(calc->nEigOut == 1);
    } else {
        assert(calc->nEigOut == calc->nEigIn);
    }
    // We currently assume a hardcoded maximum for the number of powers.
    if (!basis->useRadialLut && basis->maxNPows > STO_MAX_POWS) {
        fprintf(stderr, "Error: maxNPows (%d) exceeds STO_MAX_POWS (%d)\n", basis->maxNPows, STO_MAX_POWS);
        exit(EXIT_FAILURE);
    }

    // Timing events.
    cudaEvent_t startEverything, endEverything, startKernelOnly, endKernelOnly, startCopyOnly, endCopyOnly;

    float totalKernelTime_ms  = 0.0f;
    float totalD2HCopyTime_ms = 0.0f;
    CHECK_CUDA(cudaEventCreate(&startEverything));
    CHECK_CUDA(cudaEventCreate(&endEverything));
    CHECK_CUDA(cudaEventCreate(&startKernelOnly));
    CHECK_CUDA(cudaEventCreate(&endKernelOnly));
    CHECK_CUDA(cudaEventCreate(&startCopyOnly));
    CHECK_CUDA(cudaEventCreate(&endCopyOnly));
    CHECK_CUDA(cudaEventRecord(startEverything));

    // --- Multi-GPU Setup ---
    int numGpus;
    CHECK_CUDA(cudaGetDeviceCount(&numGpus));
    if (numGpus == 0) {
        fprintf(stderr, "No CUDA-enabled GPUs found. Unable to launch Kernel.\n");
        exit(EXIT_FAILURE);
    }
    printf("Libwavegrid: Found %d CUDA-enabled GPUs.\n", numGpus);

#ifndef _OPENMP
    if (numGpus > 1) {
        fprintf(stderr, "\nWARNING: Code not compiled with OpenMP support (-fopenmp). Falling back to single-GPU mode.\n");
        numGpus = 1;
        printf("Running on GPU 0 only.\n");
    }
#endif

// Use OMP to split across available GPUs
// This works irrespective of the number of threads set in OMP_NUM_THREADS.
#pragma omp parallel num_threads(numGpus)
    {
        int deviceId = omp_get_thread_num();
        CHECK_CUDA(cudaSetDevice(deviceId));
        GpuLaunchConfig config = setup_gpu_config(deviceId, numGpus, grid, calc);

        if (config.z_count > 0) {
            // Device allocation and H2D transfer
            DeviceData device_data(grid, system, periodic, basis, calc, config.z_per_batch);

            // Populate Kernel Parameter struct
            DeviceKernelParams deviceParams(device_data, grid, system, periodic, basis, calc, config.nEig_per_pass);

            // --- Per-GPU Kernel Execution Loop ---
            // This loop iterates over the Z-slices assigned to *this* GPU.
            for (int z_offset = 0; z_offset < config.z_count; z_offset += config.z_per_batch) {
                deviceParams.z_per_batch = std::min(config.z_per_batch, config.z_count - z_offset);
                deviceParams.z_offset_global =
                    config.z_start + z_offset;  // required to calculate coordinates in kernel

                int total_points_in_batch = grid->nPointsX * grid->nPointsY * deviceParams.z_per_batch;
                int grid_size             = (total_points_in_batch + block_size - 1) / block_size;

                if (deviceId == 0) CHECK_CUDA(cudaEventRecord(startKernelOnly));

                dispatchKernel(&deviceParams, calc->isRealInput, calc->calcAtomicDensity, calc->calcTotalChrg,
                    basis->useRadialLut, grid_size, config.shared_mem_for_pass);

                if (deviceId == 0) {
                    CHECK_CUDA(cudaEventRecord(endKernelOnly));
                    CHECK_CUDA(cudaEventRecord(startCopyOnly));
                }

                copyD2H(calc->isRealOutput ? (void*)device_data.d_valueReal_out_batch.get()
                                           : (void*)device_data.d_valueCmpl_out_batch.get(),
                    calc->isRealOutput ? (void*)calc->valueReal_out : (void*)calc->valueCmpl_out, grid->nPointsX,
                    grid->nPointsY, grid->nPointsZ, deviceParams.z_per_batch, deviceParams.z_offset_global, calc);

                if (deviceId == 0) {
                    float iterKernel_ms, iterCopy_ms;
                    CHECK_CUDA(cudaEventRecord(endCopyOnly));
                    CHECK_CUDA(cudaEventSynchronize(endCopyOnly));
                    CHECK_CUDA(cudaEventElapsedTime(&iterKernel_ms, startKernelOnly, endKernelOnly));
                    CHECK_CUDA(cudaEventElapsedTime(&iterCopy_ms, endKernelOnly, endCopyOnly));
                    totalKernelTime_ms += iterKernel_ms;
                    totalD2HCopyTime_ms += iterCopy_ms;
                }
            }
        }
    }  // End of omp parallel region

    // Synchronize all devices
    for (int i = 0; i < numGpus; ++i) {
        CHECK_CUDA(cudaSetDevice(i));
        CHECK_CUDA(cudaDeviceSynchronize());
    }

    // Switch back to lead to retrieve timing
    CHECK_CUDA(cudaSetDevice(0));
    CHECK_CUDA(cudaGetLastError());
    CHECK_CUDA(cudaEventRecord(endEverything));
    CHECK_CUDA(cudaEventSynchronize(endEverything));

    float timeEverything;
    CHECK_CUDA(cudaEventElapsedTime(&timeEverything, startEverything, endEverything));

    float overhead = timeEverything - (totalKernelTime_ms + totalD2HCopyTime_ms);
    if (debug) printf("\n--- GPU Timing Results ---\n");
    printf("Total Multi-GPU execution time: %.2f ms\n", timeEverything);
    // Timings are run on device 0 only.
    if (debug) {
        printf("(Lead) Kernel execution: %.2f ms (%.1f%%)\n", totalKernelTime_ms,
            (totalKernelTime_ms / timeEverything) * 100.0);
        printf("(Lead) D2H Copy:         %.2f ms (%.1f%%)\n", totalD2HCopyTime_ms,
            (totalD2HCopyTime_ms / timeEverything) * 100.0);
        printf("(Lead) Overhead:         %.2f ms (%.1f%%)\n", overhead, (overhead / timeEverything) * 100.0);
    }
}
