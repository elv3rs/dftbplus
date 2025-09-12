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
#include "host_logic.cuh"
#include "kernel_sig.cuh"
#include "device_params.cuh"
#include "slater.cuh"
#include "spharmonics.cuh"
#include "utils.cuh"

// Avoid division by zero
constexpr double INV_R_EPSILON = 1.0e-12;

using complexd = thrust::complex<double>;

// ================================================================================================
//  Main MolOrb CUDA Kernel
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

    // If periodic, fold into cell by discarding the non-fractional part in lattice vector multiples.
    if (p.isPeriodic) foldCoordsIntoCell(xyz, p.latVecs, p.recVecs2pi);

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


// C++ Host Interface (extern "C", then called from Fortran)
// Handles the high-level flow of querying devices, and then
// concurrently preparing data, launching kernels, and copying back results.
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

    // Debug kernel timing
    GpuTimer everything_timer(true), kernel_timer, d2h_timer;


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
                deviceParams.z_offset_global = config.z_start + z_offset;  // required to calculate coordinates in kernel

                int total_points_in_batch = grid->nPointsX * grid->nPointsY * deviceParams.z_per_batch;
                int grid_size             = (total_points_in_batch + block_size - 1) / block_size;

                if (deviceId == 0) kernel_timer.start();

                dispatchKernel(&deviceParams, calc->isRealInput, calc->calcAtomicDensity, calc->calcTotalChrg,
                    basis->useRadialLut, grid_size, config.shared_mem_for_pass);

                if (deviceId == 0) {
                    kernel_timer.stop();
                    d2h_timer.start();
                }

                copyD2H(calc->isRealOutput ? (void*)device_data.d_valueReal_out_batch.get()
                                           : (void*)device_data.d_valueCmpl_out_batch.get(),
                    calc->isRealOutput ? (void*)calc->valueReal_out : (void*)calc->valueCmpl_out, grid->nPointsX,
                    grid->nPointsY, grid->nPointsZ, deviceParams.z_per_batch, deviceParams.z_offset_global, calc);

                if (deviceId == 0) d2h_timer.stop();
            }
        }
    }  // End of omp parallel region

    // Synchronize all devices
    for (int i = 0; i < numGpus; ++i) {
        CHECK_CUDA(cudaSetDevice(i));
        CHECK_CUDA(cudaDeviceSynchronize());
    }
    CHECK_CUDA(cudaSetDevice(0));
    CHECK_CUDA(cudaGetLastError());
    everything_timer.stop();

    if (debug) {
        printf("\nGPU 0 execution time: %.1f ms\n", everything_timer.elapsed_ms());
        float kernel_share = kernel_timer.elapsed_ms() / everything_timer.elapsed_ms();
        printf(" -> Kernel:   %.1f ms (%.1f%%)\n", kernel_timer.elapsed_ms(), kernel_share * 100.0f);
        printf(" -> D2H copy: %.1f ms (%.1f%%)\n", d2h_timer.elapsed_ms(), (1.0f - kernel_share) * 100.0f);
    }
}
