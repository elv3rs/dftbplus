/*-------------------------------------------------------------------------------------------------*
 *  DFTB+: general package for performing fast atomistic simulations                               *
 *  Copyright (C) 2006 - 2025  DFTB+ developers group                                              *
 *                                                                                                 *
 *  See the LICENSE file for terms of usage and distribution.                                      *
 *-------------------------------------------------------------------------------------------------*/
#include <cuda_runtime.h>
#include <thrust/complex.h>
#include <omp.h>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <assert.h>

#include "kernel.cuh"
#include "utils.cuh"
#include "slater.cuh"

// more print statements
constexpr bool debug = true; 
constexpr double INV_R_EPSILON = 1.0e-12;
// amount of shared memory set aside for nEig accumulators
constexpr float SHARED_MEM_FACTOR = 0.95f;
// max output array share of free global memory
constexpr float GLOBAL_MEM_FACTOR = 0.80f;
// Threads per block, multiple of warp size 32
constexpr int block_size = 256; 

using complexd = thrust::complex<double>;

// Manages Gpu memory allocation
struct DeviceData {
    // Grid
    DeviceBuffer<double> origin;
    DeviceBuffer<double> gridVecs;

    // System
    DeviceBuffer<double> coords;
    DeviceBuffer<int>    species;
    DeviceBuffer<int>    iStos;

    // Periodic
    DeviceBuffer<double> latVecs;
    DeviceBuffer<double> recVecs2pi;
    DeviceBuffer<int>    kIndexes;
    DeviceBuffer<complexd> phases;

    // STO Basis
    DeviceBuffer<int>    sto_angMoms;
    DeviceBuffer<int>    sto_nPows;
    DeviceBuffer<int>    sto_nAlphas;
    DeviceBuffer<double> sto_cutoffsSq;
    DeviceBuffer<double> sto_coeffs;
    DeviceBuffer<double> sto_alphas;
    // Texture for radial LUT
    cudaTextureObject_t sto_lutTex;

    // Eigenvectors
    DeviceBuffer<double> eigVecsReal;
    DeviceBuffer<complexd> eigVecsCmpl;

    // Constructor handles all H2D allocation and copy
    DeviceData(const GridParams* grid, const SystemParams* system, const PeriodicParams* periodic, const StoBasisParams* basis, const CalculationParams* calc)
        : origin(grid->origin, 3),
          gridVecs(grid->gridVecs, 9),
          coords(system->coords, (size_t)3 * system->nAtom * system->nCell),
          species(system->species, system->nAtom),
          iStos(system->iStos, system->nSpecies + 1),
    {
        if (basis->useRadialLut) {
            // Convert the Fortran passed doubles to floats
            size_t totalLutValues = (size_t)basis->nStos * basis->nLutPoints;
            std::vector<float> lutGridValuesFloat(totalLutValues);
            for (size_t i = 0; i < totalLutValues; ++i) {
                lutGridValuesFloat[i] = (float)(basis->lutGridValues[i]);
            }

            // Allocate
            cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
            cudaExtent shape = make_cudaExtent(basis->nLutPoints, basis->nStos, 0);
            cudaArray_t lutArray;
            CHECK_CUDA(cudaMalloc3DArray(&lutArray, &channelDesc, shape, cudaArrayLayered));

            // Copy data to array
            cudaMemcpy2DToArray(
                lutArray,                              // dst array
                0, 0,                                  // no offset in dst
                lutGridValuesFloat.data(),             // src pointer
                basis->nLutPoints * sizeof(float),     // src pitch (for alignment, bytes tp next row)
                basis->nLutPoints * sizeof(float),     // width in bytes
                basis->nStos,                          // height (number of cached stos)
                cudaMemcpyHostToDevice
            );

            // Prepare texture object properties
            cudaResourceDesc resDesc{};
            resDesc.resType = cudaResourceTypeArray;
            resDesc.res.array.array = lutArray;

            cudaTextureDesc texDesc{};
            texDesc.addressMode[0] = cudaAddressModeClamp;
            texDesc.addressMode[1] = cudaAddressModeClamp;
            texDesc.filterMode = cudaFilterModeLinear; 
            texDesc.readMode = cudaReadModeElementType; // dont normalize our lut values
            texDesc.normalizedCoords = 0; // access using [0, N-1]
            // Create texture object
            CHECK_CUDA(cudaCreateTextureObject(&sto_lutTex, &resDesc, &texDesc, nullptr));

        } else {
            sto_angMoms.assign(basis->sto_angMoms, basis->nStos);
            sto_nPows.assign(basis->sto_nPows, basis->nStos);
            sto_nAlphas.assign(basis->sto_nAlphas, basis->nStos);
            sto_cutoffsSq.assign(basis->sto_cutoffsSq, basis->nStos);
            sto_coeffs.assign(basis->sto_coeffs, (size_t)basis->maxNPows * basis->maxNAlphas * basis->nStos);
            sto_alphas.assign(basis->sto_alphas, (size_t)basis->maxNAlphas * basis->nStos);
        }

        if (calc->isRealInput) {
            eigVecsReal.assign(calc->eigVecsReal, (size_t)system->nOrb * calc->nEigIn);
        } else {
            eigVecsCmpl.assign(reinterpret_cast<const complexd*>(calc->eigVecsCmpl), (size_t)system->nOrb * calc->nEigIn);
            phases.assign(reinterpret_cast<const complexd*>(periodic->phases), (size_t)system->nCell * calc->nEigIn);
            kIndexes.assign(periodic->kIndexes, calc->nEigIn);
        }
        if (periodic->isPeriodic) {
            latVecs.assign(periodic->latVecs, 9);
            recVecs2pi.assign(periodic->recVecs2pi, 9);
        }
    }
};

// Kernel parameters
struct DeviceKernelParams {
    // Grid
    int nPointsX, nPointsY, nPointsZ_batch, z_offset_global;
    const double* origin;
    const double* gridVecs;

    // System
    int nAtom, nCell, nOrb;
    const double* coords;
    const int*    species;
    const int*    iStos;

    // Periodic boundary cond.
    bool isPeriodic;
    const double* latVecs;
    const double* recVecs2pi;
    const int*    kIndexes;
    const complexd* phases;

    // STO Basis
    int nStos, maxNPows, maxNAlphas;
    // Texture LUTs
    cudaTextureObject* lutTex;
    double inverseLutStep;
    // STO parameters
    const int*    sto_angMoms;
    const int*    sto_nPows;
    const int*    sto_nAlphas;
    const double* sto_cutoffsSq;
    const double* sto_coeffs;
    const double* sto_alphas;

    // Eigenvectors
    int nEig, nEig_per_pass;
    const double* eigVecsReal;
    const complexd* eigVecsCmpl;

    // Output (batch pointers), varied in batch loop
    double* valueReal_out_batch;
    complexd* valueCmpl_out_batch;

    // Constructor to initialize the parameters from host data
    // Batch-specific parameters are initialized to zero or nullptr,
    // and need to be set in the loop before kernel launch.
    DeviceKernelParams(
        const DeviceData& data,
        const GridParams* grid,
        const SystemParams* system,
        const PeriodicParams* periodic,
        const StoBasisParams* basis,
        const CalculationParams* calc
    ) {
        // Grid
        origin = data.origin.get();
        gridVecs = data.gridVecs.get();
        nPointsX = grid->nPointsX;
        nPointsY = grid->nPointsY;

        // System
        nAtom = system->nAtom;
        nCell = system->nCell;
        nOrb = system->nOrb;
        coords = data.coords.get();
        species = data.species.get();
        iStos = data.iStos.get();

        // STO Basis
        nStos = basis->nStos;
        if (basis->useRadialLut) {
            lutTex = data.sto_lutTex;
            inverseLutStep = basis->inverseLutStep;
        } else {
            maxNPows = basis->maxNPows;
            maxNAlphas = basis->maxNAlphas;
            sto_angMoms = data.sto_angMoms.get();
            sto_nPows = data.sto_nPows.get();
            sto_nAlphas = data.sto_nAlphas.get();
            sto_cutoffsSq = data.sto_cutoffsSq.get();
            sto_coeffs = data.sto_coeffs.get();
            sto_alphas = data.sto_alphas.get();
        }


        // Periodic boundary conditions
        isPeriodic = periodic->isPeriodic;
        latVecs = data.latVecs.get();
        recVecs2pi = data.recVecs2pi.get();
        kIndexes = data.kIndexes.get();
        phases = data.phases.get();

        // Eigenvectors
        nEig = calc->nEigIn;
        eigVecsReal = data.eigVecsReal.get();
        eigVecsCmpl = data.eigVecsCmpl.get();

        // Batch-specific kernel config to be updated in the loop
        nPointsZ_batch = 0;
        z_offset_global = 0;
        nEig_per_pass = 0;
        valueReal_out_batch = nullptr;
        valueCmpl_out_batch = nullptr;
    }
};



// =========================================================================
//  CUDA Kernel.
// =========================================================================
// To avoid branching (dropped at compile time), we template the kernel 16 ways on (isRealInput, calcDensity, calcTotalChrg, useLUT).
// useLUT decides whether to use texture memory interpolation for STO radial functions.
// isPeriodic decides whether to fold coords into unit cell.
// isRealInput decides whether to use real/complex eigenvectors (and adds phases)
// calcAtomicDensity squares the basis wavefunction contributions, result in valueReal_out of shape (x,y,z,n)
// calcTotalChrg accumulates the density over all states, leading to valueReal_out of shape (x,y,z,1).
// User is responsible for providing eigenvec multiplied with sqrt(occupation) if needed.
template <bool isRealInput, bool calcAtomicDensity, bool calcTotalChrg, bool useLUT>
__global__ void evaluateKernel(const DeviceKernelParams p)
{
    using AccumT = typename std::conditional<(isRealInput), double, complexd>::type;
    
    // Each thread gets its own private slice of the shared memory buffer for fast accumulation.
    // We have to chunk the eigenstates into nEig_per_pass due to size constraints.
    // (Cuda doesnt allow templating the shared memory type, so we simply recast it.)
    extern __shared__ char shared_workspace[];
    AccumT* point_results_pass = reinterpret_cast<AccumT*>(shared_workspace) + threadIdx.x * p.nEig_per_pass;


    // --- Thread to point mapping ---
    // Map each thread to unique 1d index
    int idx_in_batch = blockIdx.x * blockDim.x + threadIdx.x;
    int total_points_in_batch = p.nPointsX * p.nPointsY * p.nPointsZ_batch;
    if (idx_in_batch >= total_points_in_batch) return;

    // Map 1d index to point in grid
    int i1 = idx_in_batch % p.nPointsX;
    int i2 = (idx_in_batch / p.nPointsX) % p.nPointsY;
    int i3_batch = idx_in_batch / (p.nPointsX * p.nPointsY);
    int i3_global = i3_batch + p.z_offset_global; 

    // Map point to global coordinates.
    double xyz[3];
    for (int i = 0; i < 3; ++i) 
        xyz[i] = p.origin[i] + i1 * p.gridVecs[IDX2F(i, 0, 3)]
                             + i2 * p.gridVecs[IDX2F(i, 1, 3)]
                             + i3_global * p.gridVecs[IDX2F(i, 2, 3)];
    
    // If periodic, fold into cell by discarding the non-fractional part in lattice vector multiples.
    if (p.isPeriodic) 
        foldCoordsIntoCell(xyz, reinterpret_cast<const double (*)[3]>(p.latVecs), reinterpret_cast<const double (*)[3]>(p.recVecs2pi));
    

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

                    if constexpr (useLUT) {
                        double lut_pos = 0.5f + r * p.inverseLutStep;
                        double radialVal = static_cast<double>(tex1DLayered<float>(p.lutTex, lut_pos, iOrb)[0]);
                    } else {
                        double radialVal = getRadialValue(
                            r, iL, iOrb, p.sto_nPows[iOrb], p.sto_nAlphas[iOrb],
                            p.sto_coeffs, p.sto_alphas, p.maxNPows, p.maxNAlphas);
                    } 


                    // precompute inverse used across several realTessY calls
                    double inv_r = (r < INV_R_EPSILON) ? 0.0 : 1.0 / r;
                    double inv_r2 = inv_r * inv_r;

                    for (int iM = -iL; iM <= iL; ++iM) {
                        double val = radialVal * realTessY(iL, iM, diff, inv_r, inv_r2);
                        
                        // Accumulate into the small shared memory buffer for the current chunk
                        for (int iEig_offset = 0; iEig_offset < p.nEig_per_pass; ++iEig_offset) {
                            int iEig = eig_base + iEig_offset;
                            if (iEig >= p.nEig) break; // Don't go past the end on the last chunk
                            size_t eig_idx = IDX2F(orbital_idx_counter, iEig, p.nOrb);
                            if constexpr (isRealInput) {
                                point_results_pass[iEig_offset] += val * p.eigVecsReal[eig_idx];
                            } else {
                                point_results_pass[iEig_offset] += val * p.phases[IDX2F(iCell, iEig, p.nCell)] * p.eigVecsCmpl[eig_idx];
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
            size_t out_idx = IDX4F(i1, i2, i3_batch, iEig, p.nPointsX, p.nPointsY, p.nPointsZ_batch);
            if constexpr (isRealInput) {
                if constexpr (calcTotalChrg) {
                    totChrgAcc += point_results_pass[iEig_offset] * point_results_pass[iEig_offset];
                } else if (calcAtomicDensity) {
                    p.valueReal_out_batch[out_idx] = point_results_pass[iEig_offset] * point_results_pass[iEig_offset];
                } else {
                    p.valueReal_out_batch[out_idx] = point_results_pass[iEig_offset];
                }
            } else {
                if constexpr (calcTotalChrg) {
                    totChrgAcc += thrust::norm(point_results_pass[iEig_offset]);
                } else if constexpr (calcAtomicDensity) {
                    p.valueReal_out_batch[out_idx] = thrust::norm(point_results_pass[iEig_offset]);
                } else {
                    p.valueCmpl_out_batch[out_idx] = point_results_pass[iEig_offset];
                }
            }
        }
    }

    // Density stored in first eig : (x,y,z, 1)
    if constexpr (calcTotalChrg) {
        size_t out_idx = IDX4F(i1, i2, i3_batch, 0, p.nPointsX, p.nPointsY, p.nPointsZ_batch);
        p.valueReal_out_batch[out_idx] = totChrgAcc; 
    }
}



// =========================================================================
//  C++ Host Interface (callable from C/Fortran)
// =========================================================================
extern "C" void evaluate_on_device_c(
    const GridParams* grid,
    const SystemParams* system,
    const PeriodicParams* periodic,
    const StoBasisParams* basis,
    const CalculationParams* calc
){
    // Since we use these often, derefence them and add to namespace
    int nPointsX = grid->nPointsX;
    int nPointsY = grid->nPointsY;
    int nPointsZ = grid->nPointsZ;
    bool isRealOutput = calc->isRealInput || calc->calcAtomicDensity;

    if (calc->nEigIn == 0 || nPointsZ == 0) return; // Nothing to do
    if (calc->calcTotalChrg) {
        assert(calc->nEigOut == 1);
    } else {
        assert(calc->nEigOut == calc->nEigIn);
    }
    
    
    // We currently assume a hardcoded maximum for the number of powers.
    if (basis->maxNPows > STO_MAX_POWS) {
        fprintf(stderr, "Error: maxNPows (%d) exceeds STO_MAX_POWS (%d)\n", basis->maxNPows, STO_MAX_POWS);
        exit(EXIT_FAILURE);
    }
    
    // Timing events.
    cudaEvent_t startEverything, endEverything;
    cudaEvent_t startKernelOnly, endKernelOnly, startCopyOnly, endCopyOnly;
    float totalKernelTime_ms = 0.0f;
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
    printf("Found %d GPUs.", numGpus);

#ifndef _OPENMP
    if (numGpus > 1) {
    printf("\nWARNING: Code not compiled with OpenMP support (-fopenmp). Falling back to single-GPU mode.\n");
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

        // --- Work Distribution: Divide Z-slices among GPUs ---
        int z_slices_per_gpu = nPointsZ / numGpus;
        int z_start_for_device = deviceId * z_slices_per_gpu;
        // Handle uneven Z-slice count
        int z_count_for_device = (deviceId == numGpus - 1) ? (nPointsZ - z_start_for_device) : z_slices_per_gpu;

        if (z_count_for_device > 0) {
            // --- Allocate on and copy data to gpu---
            DeviceData device_data(grid, system, periodic, basis, calc);

            // --- Per-GPU Kernel Configuration ---
            cudaDeviceProp prop;
            CHECK_CUDA(cudaGetDeviceProperties(&prop, deviceId));
            
            // Determine available shared memory for nEig_per_pass
            size_t available_shared = prop.sharedMemPerBlock * SHARED_MEM_FACTOR;
            size_t number_size = isRealOutput ? sizeof(double) : sizeof(complexd);
            int nEig_per_pass = available_shared / (block_size * number_size);
            if (nEig_per_pass == 0) nEig_per_pass = 1;
            if (nEig_per_pass > calc->nEigIn) nEig_per_pass = calc->nEigIn;
            size_t shared_mem_for_pass = (size_t)nEig_per_pass * block_size * number_size;

            
            // Determine the number of Z-slices to process in a single batch
            size_t free_mem, total_mem;
            CHECK_CUDA(cudaMemGetInfo(&free_mem, &total_mem));
            size_t available_for_batch = static_cast<size_t>(free_mem * GLOBAL_MEM_FACTOR);
            size_t z_slice_size_bytes = (size_t)nPointsX * nPointsY * calc->nEigOut * number_size;
            
            // Determine max Z-slices that can fit in available (global) memory
            int z_batch_size = z_count_for_device; 
            if (z_slice_size_bytes > 0 && ((size_t)z_count_for_device * z_slice_size_bytes) > available_for_batch) {
                z_batch_size = available_for_batch / z_slice_size_bytes;
                if (z_batch_size == 0) z_batch_size = 1;
            }


            // Per-GPU batch buffer for the output
            size_t batch_buffer_size_elems = (size_t)nPointsX * nPointsY * std::min(z_count_for_device, z_batch_size) * calc->nEigOut;
            // todo: move this to deviceParams
            DeviceBuffer<complexd> d_valueCmpl_out_batch;
            DeviceBuffer<double> d_valueReal_out_batch;
            if (calc->isRealInput || calc->calcAtomicDensity) {
                d_valueReal_out_batch = DeviceBuffer<double>(batch_buffer_size_elems);
            } else {
                d_valueCmpl_out_batch = DeviceBuffer<complexd>(batch_buffer_size_elems);
            }

            // Debug output
            #pragma omp critical
            if (deviceId == 0 && debug) {
                printf("\n--- GPU %d (Lead) Configuration ---\n", deviceId);
                printf("  Z-slice workload: %d (from index %d to %d)\n", z_count_for_device, z_start_for_device, z_start_for_device + z_count_for_device - 1);
                printf("  Block size: %d threads, %zub shared mem per block, %d eigs per pass\n",
                    block_size, shared_mem_for_pass, nEig_per_pass);
                size_t total_size_valueOut = (size_t)nPointsX * nPointsY * nPointsZ * calc->nEigOut * sizeof(double);
                if (!calc->isRealInput && !calc->calcAtomicDensity) total_size_valueOut *= 2; 
                printf(" (Free device mem: %.2f GB, Grid size: %d x %d x %d (x %d eigs) = %.2f GB)\n",
                    free_mem / 1e9, nPointsX, nPointsY, nPointsZ, calc->nEigOut,
                    total_size_valueOut / 1e9);
                printf("  Processing Z-slices in batches of %d\n", z_batch_size);

            }

            // --- Populate Kernel Parameter struct ---
            DeviceKernelParams deviceParams(device_data, grid, system, periodic, basis, calc);
            deviceParams.nEig_per_pass = nEig_per_pass; // Set remaining params
            deviceParams.valueReal_out_batch = d_valueReal_out_batch.get();
            deviceParams.valueCmpl_out_batch = d_valueCmpl_out_batch.get();


            // --- Per-GPU Kernel Execution Loop ---
            // This loop iterates over the Z-slices assigned to *this* GPU.
            for (int z_offset_in_device_chunk = 0; z_offset_in_device_chunk < z_count_for_device; z_offset_in_device_chunk += z_batch_size) {
                deviceParams.nPointsZ_batch = std::min(z_batch_size, z_count_for_device - z_offset_in_device_chunk);

                int total_points_in_batch = nPointsX * nPointsY * deviceParams.nPointsZ_batch;
                if (total_points_in_batch == 0) continue;

                // The global z_offset is what the kernel needs to calculate correct coordinates
                deviceParams.z_offset_global = z_start_for_device + z_offset_in_device_chunk;
                int grid_size = (total_points_in_batch + block_size - 1) / block_size;
                if(deviceId == 0) {
                    CHECK_CUDA(cudaEventRecord(startKernelOnly));
                }
                 
                #define CALL_KERNEL(isReal, doAtomic, doChrg, useLut) \
                    evaluateKernel<isReal, doAtomic, doChrg, useLut> \
                        <<<grid_size, block_size, shared_mem_for_pass>>>(deviceParams);

                int idx = (calc->isRealInput     ? 1 : 0)
                        + (calc->calcAtomicDensity ? 2 : 0)
                        + (calc->calcTotalChrg     ? 4 : 0);
                        + (basis->useRadialLut     ? 8 : 0);

                switch (idx) {
                    case 0:  CALL_KERNEL(false, false, false, false); break;
                    case 1:  CALL_KERNEL(true,  false, false, false); break;
                    case 2:  CALL_KERNEL(false, true,  false, false); break;
                    case 3:  CALL_KERNEL(true,  true,  false, false); break;
                    case 4:  CALL_KERNEL(false, false, true,  false); break;
                    case 5:  CALL_KERNEL(true,  false, true,  false); break;
                    case 6:  CALL_KERNEL(false, true,  true,  false); break;
                    case 7:  CALL_KERNEL(true,  true,  true,  false); break;
                    case 8:  CALL_KERNEL(false, false, false, true);  break;
                    case 9:  CALL_KERNEL(true,  false, false, true);  break;
                    case 10: CALL_KERNEL(false, true,  false, true);  break;
                    case 11: CALL_KERNEL(true,  true,  false, true);  break;
                    case 12: CALL_KERNEL(false, false, true,  true);  break;
                    case 13: CALL_KERNEL(true,  false, true,  true);  break;
                    case 14: CALL_KERNEL(false, true,  true,  true);  break;
                    case 15: CALL_KERNEL(true,  true,  true,  true);  break;
                }
                                    


                if(deviceId == 0) {
                    CHECK_CUDA(cudaEventRecord(endKernelOnly));
                    CHECK_CUDA(cudaEventRecord(startCopyOnly));
                }

                // --- D2H Copy ---
                // Copy the computed batch back to the correct slice of the final host array.
                // The D2H copy will automatically block/ synchronize the kernel for this batch.
                // This could be improved by using streams / cudaMemcpyAsync.
                void *d_src_ptr = (isRealOutput ? (void*)d_valueReal_out_batch.get() : (void*)d_valueCmpl_out_batch.get());
                void* h_dest_ptr = (isRealOutput ? (void*)calc->valueReal_out : (void*)calc->valueCmpl_out);

                size_t host_plane_size = (size_t)nPointsZ * nPointsY * nPointsX * number_size;
                size_t device_plane_size = (size_t)deviceParams.nPointsZ_batch * nPointsY * nPointsX * number_size;

                for(int iEig = 0; iEig < calc->nEigOut; ++iEig) {
                    // From: iEig-th slice of GPU batch buffer
                    ptrdiff_t d_offset_bytes = (ptrdiff_t)(iEig * device_plane_size);
                    
                    // To: Global Z-position in the iEig-th slice of host buffer
                    ptrdiff_t h_offset_bytes = (ptrdiff_t)(iEig * host_plane_size + ( (size_t)deviceParams.z_offset_global * nPointsY * nPointsX) * number_size);

                    CHECK_CUDA(cudaMemcpy((char*)h_dest_ptr + h_offset_bytes, (char*)d_src_ptr + d_offset_bytes, device_plane_size, cudaMemcpyDeviceToHost));
                }


                if(deviceId == 0) {
                    CHECK_CUDA(cudaEventRecord(endCopyOnly));
                    CHECK_CUDA(cudaEventSynchronize(endCopyOnly));

                    float iterKernel_ms, iterCopy_ms;
                    CHECK_CUDA(cudaEventElapsedTime(&iterKernel_ms, startKernelOnly, endKernelOnly));
                    CHECK_CUDA(cudaEventElapsedTime(&iterCopy_ms, endKernelOnly, endCopyOnly));

                    totalKernelTime_ms += iterKernel_ms;
                    totalD2HCopyTime_ms += iterCopy_ms;
                }
            }
        }
    } // End of omp parallel region

    // Synchronize all devices
    for(int i = 0; i < numGpus; ++i) {
        CHECK_CUDA(cudaSetDevice(i));
        CHECK_CUDA(cudaDeviceSynchronize());
    }

    // Switch back to lead to retrieve timing
    CHECK_CUDA(cudaSetDevice(0));
    CHECK_CUDA(cudaGetLastError());

    // --- Final Timing ---
    CHECK_CUDA(cudaEventRecord(endEverything));
    CHECK_CUDA(cudaEventSynchronize(endEverything));

    float timeEverything;
    CHECK_CUDA(cudaEventElapsedTime(&timeEverything, startEverything, endEverything));


    float overhead = timeEverything - (totalKernelTime_ms + totalD2HCopyTime_ms);
    if(debug)
    printf("\n--- GPU Timing Results ---\n");
    printf("Total Multi-GPU execution time: %.2f ms\n", timeEverything);
    if(debug){
    printf("(Lead) Kernel execution: %.2f ms (%.1f%%)\n", totalKernelTime_ms, (totalKernelTime_ms / timeEverything) * 100.0);
    printf("(Lead) D2H Copy:         %.2f ms (%.1f%%)\n", totalD2HCopyTime_ms, (totalD2HCopyTime_ms / timeEverything) * 100.0);
    }
}
