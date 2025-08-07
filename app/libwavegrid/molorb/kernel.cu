#include <cuda_runtime.h>
#include <omp.h>
#include <cstdio>
#include <cmath>
#include <algorithm>

#include "kernel.cuh"
#include "utils.cuh"
#include "slater.cuh"



// =========================================================================
//  CUDA Kernel.
//  We might want to separate the arguments into structs for better maintainability.
// =========================================================================
__global__ void evaluateKernel(
    const int  nPointsX, const int nPointsY, int nPointsZ_batch, int z_offset, int nEig, int nOrb, int nStos,
    int maxNPows, int maxNAlphas, int nAtom, int nCell, int nEig_per_pass,
    const double* __restrict__ origin, const double* __restrict__ gridVecs, const double* __restrict__ eigVecsReal,
    const double* __restrict__ coords, const int* __restrict__ species, const int* __restrict__ iStos,
    const int* __restrict__ sto_angMoms, const int* __restrict__ sto_nPows, const int* __restrict__ sto_nAlphas,
    const double* __restrict__ sto_cutoffsSq, const double* __restrict__ sto_coeffs, const double* __restrict__ sto_alphas,
    double* valueReal_out_batch)

{
    extern __shared__ double shared_workspace[];

    // Each thread gets its own private slice of the shared memory buffer for fast accumulation.
    // We have to chunk the eigenstates into nEig_per_pass due to size constraints.
    double* point_results_pass = &shared_workspace[threadIdx.x * nEig_per_pass];

    // --- Thread to point mapping ---
    // Map each thread to unique 1d index
    int idx_in_batch = blockIdx.x * blockDim.x + threadIdx.x;
    int total_points_in_batch = nPointsX * nPointsY * nPointsZ_batch;
    if (idx_in_batch >= total_points_in_batch) return;

    // Map 1d index to point in grid
    int i1 = idx_in_batch % nPointsX;
    int i2 = (idx_in_batch / nPointsX) % nPointsY;
    int i3_batch = idx_in_batch / (nPointsX * nPointsY);
    int i3_global = i3_batch + z_offset;

    // Map point to global coordinates
    double xyz[3];
    xyz[0] = origin[0] + i1 * gridVecs[IDX2F(0, 0, 3)] + i2 * gridVecs[IDX2F(0, 1, 3)] + i3_global * gridVecs[IDX2F(0, 2, 3)];
    xyz[1] = origin[1] + i1 * gridVecs[IDX2F(1, 0, 3)] + i2 * gridVecs[IDX2F(1, 1, 3)] + i3_global * gridVecs[IDX2F(1, 2, 3)];
    xyz[2] = origin[2] + i1 * gridVecs[IDX2F(2, 0, 3)] + i2 * gridVecs[IDX2F(2, 1, 3)] + i3_global * gridVecs[IDX2F(2, 2, 3)];

    // --- Loop over eigenstates in chunks that fit in shared memory ---
    for (int eig_base = 0; eig_base < nEig; eig_base += nEig_per_pass) {
        
        // Initialize the small, per-pass buffer for this thread
        for (int i = 0; i < nEig_per_pass; ++i) {
            point_results_pass[i] = 0.0;
        }

        // The spatial calculation is repeated for each chunk of eigenstates.
        // This is to keep the accumulation in fast shared memory.
        for (int iCell = 0; iCell < nCell; ++iCell) {
            int orbital_idx_counter = 0; 
            for (int iAtom = 0; iAtom < nAtom; ++iAtom) {
                int iSpecies = species[iAtom] - 1;
                double diff[3];
                diff[0] = xyz[0] - coords[IDX3F(0, iAtom, iCell, 3, nAtom)];
                diff[1] = xyz[1] - coords[IDX3F(1, iAtom, iCell, 3, nAtom)];
                diff[2] = xyz[2] - coords[IDX3F(2, iAtom, iCell, 3, nAtom)];
                double r_sq = diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2];

                for (int iOrb = iStos[iSpecies] - 1; iOrb < iStos[iSpecies + 1] - 1; ++iOrb) {
                    int iL = sto_angMoms[iOrb];
                    if (r_sq > sto_cutoffsSq[iOrb]) {
                        orbital_idx_counter += 2 * iL + 1;
                        continue;
                    }
                    double r = sqrt(r_sq);
                    
                    double radialVal = getRadialValue(
                        r, iL, iOrb, sto_nPows[iOrb], sto_nAlphas[iOrb],
                        sto_coeffs, sto_alphas, maxNPows, maxNAlphas);


                    // precompute inverse used across several realTessY calls
                    double inv_r = (r < 1.e-12) ? 0.0 : 1.0 / r;
                    double inv_r2 = inv_r * inv_r;

                    for (int iM = -iL; iM <= iL; ++iM) {
                        double val = radialVal * realTessY(iL, iM, diff, inv_r, inv_r2);
                        
                        // Accumulate into the small shared memory buffer for the current chunk
                        for (int iEig_offset = 0; iEig_offset < nEig_per_pass; ++iEig_offset) {
                            int iEig = eig_base + iEig_offset;
                            if (iEig >= nEig) break; // Don't go past the end on the last chunk
                            size_t eig_idx = IDX2F(orbital_idx_counter, iEig, nOrb);
                            point_results_pass[iEig_offset] += val * eigVecsReal[eig_idx];
                        }
                        orbital_idx_counter++;
                    }
                }
            }
        }

        // Write the complete nEig_per_pass chunk to global memory.
        for (int iEig_offset = 0; iEig_offset < nEig_per_pass; ++iEig_offset) {
            int iEig = eig_base + iEig_offset;
            if (iEig >= nEig) break;
            size_t out_idx = IDX4F(i1, i2, i3_batch, iEig, nPointsX, nPointsY, nPointsZ_batch);
            valueReal_out_batch[out_idx] = point_results_pass[iEig_offset];
        }
    }
}



// =========================================================================
//  C++ Host Interface (callable from C/Fortran)
// =========================================================================
extern "C" void evaluate_on_device_c(
    const int nPointsX, const int nPointsY, const int nPointsZ, const int nEig,
    const int nOrb, const int nStos, const int maxNPows, const int maxNAlphas,
    const int nAtom, const int nCell, const int nSpecies,
    const double* h_origin, const double* h_gridVecs, const double* h_eigVecsReal,
    const double* h_coords, const int* h_species, const int* h_iStos,
    const int* h_sto_angMoms, const int* h_sto_nPows, const int* h_sto_nAlphas,
    const double* h_sto_cutoffsSq, const double* h_sto_coeffs, const double* h_sto_alphas,
    double* h_valueReal_out)
{
    if (nEig == 0 || nPointsZ == 0) return; // Nothing to do
    
    // We currently assume a hardcoded maximum for the number of powers.
    if (maxNPows > STO_MAX_POWS) {
        fprintf(stderr, "Error: maxNPows (%d) exceeds STO_MAX_POWS (%d)\n", maxNPows, STO_MAX_POWS);
        exit(EXIT_FAILURE);
    }

    size_t total_size_valueReal = (size_t)nPointsX * nPointsY * nPointsZ * nEig * sizeof(double);
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
        fprintf(stderr, "Error: No CUDA-enabled GPUs found.\n");
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
    #pragma omp parallel num_threads(numGpus) 
    {
        int deviceId = omp_get_thread_num();
        CHECK_CUDA(cudaSetDevice(deviceId));

        // --- Work Distribution: Divide Z-slices among GPUs ---
        int z_slices_per_gpu = nPointsZ / numGpus;
        int z_start_for_device = deviceId * z_slices_per_gpu;
        int z_count_for_device = (deviceId == numGpus - 1) ? (nPointsZ - z_start_for_device) : z_slices_per_gpu;

        if (z_count_for_device > 0) {
            // --- Per-GPU Data Allocation and H2D Copy ---
            // Each thread allocates data on its own assigned GPU.
            DeviceBuffer<double> d_origin(h_origin, 3);
            DeviceBuffer<double> d_gridVecs(h_gridVecs, 9);
            DeviceBuffer<double> d_eigVecsReal(h_eigVecsReal, (size_t)nOrb * nEig);
            DeviceBuffer<double> d_coords(h_coords, (size_t)3 * nAtom * nCell);
            DeviceBuffer<int>    d_species(h_species, nAtom);
            DeviceBuffer<int>    d_iStos(h_iStos, nSpecies + 1);
            DeviceBuffer<int>    d_sto_angMoms(h_sto_angMoms, nStos);
            DeviceBuffer<int>    d_sto_nPows(h_sto_nPows, nStos);
            DeviceBuffer<int>    d_sto_nAlphas(h_sto_nAlphas, nStos);
            DeviceBuffer<double> d_sto_cutoffsSq(h_sto_cutoffsSq, nStos);
            DeviceBuffer<double> d_sto_coeffs(h_sto_coeffs, (size_t)maxNPows * maxNAlphas * nStos);
            DeviceBuffer<double> d_sto_alphas(h_sto_alphas, (size_t)maxNAlphas * nStos);

            // --- Per-GPU Kernel Configuration ---
            int block_size = 256;
            cudaDeviceProp prop;
            CHECK_CUDA(cudaGetDeviceProperties(&prop, deviceId));
            
            // Determine available shared memory for nEig_per_pass
            size_t available_shared = prop.sharedMemPerBlock * 0.95;
            int nEig_per_pass = available_shared / (block_size * sizeof(double));
            if (nEig_per_pass == 0) nEig_per_pass = 1;
            if (nEig_per_pass > nEig) nEig_per_pass = nEig;
            size_t shared_mem_for_pass = (size_t)nEig_per_pass * block_size * sizeof(double);
            
            // Determine the number of Z-slices to process in a single batch
            size_t free_mem, total_mem;
            CHECK_CUDA(cudaMemGetInfo(&free_mem, &total_mem));
            size_t available_for_batch = static_cast<size_t>(free_mem * 0.8);
            size_t z_slice_size_bytes = (size_t)nPointsX * nPointsY * nEig * sizeof(double);
            
            // Determine max Z-slices that can fit in available (global) memory
            int z_batch_size = z_count_for_device; 
            if (z_slice_size_bytes > 0 && ((size_t)z_count_for_device * z_slice_size_bytes) > available_for_batch) {
                z_batch_size = available_for_batch / z_slice_size_bytes;
                if (z_batch_size == 0) z_batch_size = 1;
            }

            #pragma omp critical
            if (deviceId == 0) {
                printf("\n--- GPU %d (Lead) Configuration ---\n", deviceId);
                printf("  Z-slice workload: %d (from index %d to %d)\n", z_count_for_device, z_start_for_device, z_start_for_device + z_count_for_device - 1);
                printf("  Block size: %d threads, %zub shared mem per block, %d eigs per pass\n",
                    block_size, shared_mem_for_pass, nEig_per_pass);
                printf(" (Free device mem: %.2f GB, Grid size: %d x %d x %d (x %d eigs) = %.2f GB)\n",
                    free_mem / 1e9, nPointsX, nPointsY, nPointsZ, nEig,
                    total_size_valueReal / 1e9);
                printf("  Processing Z-slices in batches of %d\n", z_batch_size);

            }

            // Per-GPU batch buffer for the output
            size_t batch_buffer_size_elems = (size_t)nPointsX * nPointsY * std::min(z_count_for_device, z_batch_size) * nEig;
            DeviceBuffer<double> d_valueReal_out_batch(batch_buffer_size_elems);


            // --- Per-GPU Kernel Execution Loop ---
            // This loop iterates over the Z-slices assigned to *this* GPU.
            for (int z_offset_in_device_chunk = 0; z_offset_in_device_chunk < z_count_for_device; z_offset_in_device_chunk += z_batch_size) {
                int current_nPointsZ_batch = std::min(z_batch_size, z_count_for_device - z_offset_in_device_chunk);
                int total_points_in_batch = nPointsX * nPointsY * current_nPointsZ_batch;
                if (total_points_in_batch == 0) continue;

                // The global z_offset is what the kernel needs to calculate correct coordinates
                int z_offset_global = z_start_for_device + z_offset_in_device_chunk;
                int grid_size = (total_points_in_batch + block_size - 1) / block_size;
                if(deviceId == 0) {
                    CHECK_CUDA(cudaEventRecord(startKernelOnly));
                }

                evaluateKernel<<<grid_size, block_size, shared_mem_for_pass>>>(
                    nPointsX, nPointsY, current_nPointsZ_batch, z_offset_global, nEig, nOrb, nStos,
                    maxNPows, maxNAlphas, nAtom, nCell, nEig_per_pass,
                    d_origin.get(), d_gridVecs.get(), d_eigVecsReal.get(),
                    d_coords.get(), d_species.get(), d_iStos.get(),
                    d_sto_angMoms.get(), d_sto_nPows.get(), d_sto_nAlphas.get(),
                    d_sto_cutoffsSq.get(), d_sto_coeffs.get(), d_sto_alphas.get(),
                    d_valueReal_out_batch.get()
                );

                if(deviceId == 0) {
                    CHECK_CUDA(cudaEventRecord(endKernelOnly));
                    CHECK_CUDA(cudaEventRecord(startCopyOnly));
                }

                // --- Per-GPU D2H Copy ---
                // Copy the computed batch back to the correct slice of the final host array.
                // The D2H copy will synchronize the kernel for this batch.
                for(int iEig = 0; iEig < nEig; ++iEig) {
                    size_t plane_size_bytes = (size_t)current_nPointsZ_batch * nPointsY * nPointsX * sizeof(double);

                    // Source pointer in this GPU's batch buffer
                    const double* d_src_ptr_eig = d_valueReal_out_batch.get() + (size_t)iEig * current_nPointsZ_batch * nPointsY * nPointsX;
                    
                    // Destination pointer in the final large host output array.
                    // The offset is calculated using the GLOBAL Z-offset.
                    double* h_dest_ptr_eig = h_valueReal_out + (size_t)iEig * nPointsZ * nPointsY * nPointsX + (size_t)z_offset_global * nPointsY * nPointsX;
                    
                    CHECK_CUDA(cudaMemcpy(h_dest_ptr_eig, d_src_ptr_eig, plane_size_bytes, cudaMemcpyDeviceToHost));
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
    printf("\n--- GPU Timing Results ---\n");
    printf("Total Multi-GPU execution time: %.2f ms\n", timeEverything);
    printf("(Lead) Kernel execution: %.2f ms (%.1f%%)\n", totalKernelTime_ms, (totalKernelTime_ms / timeEverything) * 100.0);
    printf("(Lead) D2H Copy:         %.2f ms (%.1f%%)\n", totalD2HCopyTime_ms, (totalD2HCopyTime_ms / timeEverything) * 100.0);
}
