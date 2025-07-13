#include <cuda_runtime.h>
#include <omp.h>
#include <cstdio>
#include <cmath>
#include <algorithm>

#include "kernel.cuh"
#include "utils.cuh"
#include "slater.cuh"



// =========================================================================
//  CUDA Kernel
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

    int idx_in_batch = blockIdx.x * blockDim.x + threadIdx.x;
    int total_points_in_batch = nPointsX * nPointsY * nPointsZ_batch;
    if (idx_in_batch >= total_points_in_batch) return;

    // Each thread gets its own private slice of the shared memory buffer.
    // This buffer is small enough to fit, but only holds results for nEig_per_pass eigenstates.
    double* point_results_pass = &shared_workspace[threadIdx.x * nEig_per_pass];

    int i1 = idx_in_batch % nPointsX;
    int i2 = (idx_in_batch / nPointsX) % nPointsY;
    int i3_batch = idx_in_batch / (nPointsX * nPointsY);
    int i3_global = i3_batch + z_offset;

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


                    // Timer calculate inverse once 
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

        // After all spatial contributions are summed for this chunk, write results to global memory.
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
    if (nEig == 0) return; // Nothing to do
    if (maxNPows > STO_MAX_POWS) {
        fprintf(stderr, "Error: maxNPows (%d) exceeds STO_MAX_POWS (%d)\n", maxNPows, STO_MAX_POWS);
        exit(EXIT_FAILURE);
    }

    // Timing events
    cudaEvent_t startKernelTimer, endKernelTimer, startCopyTimer, endCopyTimer, 
                startEverything, endEverything;
    CHECK_CUDA(cudaEventCreate(&startEverything));
    CHECK_CUDA(cudaEventCreate(&endEverything));
    CHECK_CUDA(cudaEventCreate(&startKernelTimer));
    CHECK_CUDA(cudaEventCreate(&endKernelTimer));
    CHECK_CUDA(cudaEventCreate(&startCopyTimer));
    CHECK_CUDA(cudaEventCreate(&endCopyTimer));

    // Copy Mem to device
    CHECK_CUDA(cudaEventRecord(startEverything));
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


    size_t total_size_valueReal = (size_t)nPointsX * nPointsY * nPointsZ * nEig * sizeof(double);

    // BATCHING & DYNAMIC KERNEL CONFIGURATION
    int block_size = 256;
    int nEig_per_pass;
    size_t shared_mem_for_pass;

    int deviceId;
    CHECK_CUDA(cudaGetDevice(&deviceId));
    cudaDeviceProp prop;
    CHECK_CUDA(cudaGetDeviceProperties(&prop, deviceId));

    size_t available_shared = prop.sharedMemPerBlock * 0.95;
    nEig_per_pass = available_shared / (block_size * sizeof(double));
    if (nEig_per_pass == 0) nEig_per_pass = 1;
    if (nEig_per_pass > nEig) nEig_per_pass = nEig;
    shared_mem_for_pass = (size_t)nEig_per_pass * block_size * sizeof(double);

    printf("Kernel Configuration:\n");
    printf("  Block size: %d threads\n", block_size);
    printf("  Eigenstates per pass: %d (out of %d total)\n", nEig_per_pass, nEig);
    printf("  Shared memory per block: %zu bytes (Device max: %zu bytes)\n", shared_mem_for_pass, prop.sharedMemPerBlock);

    size_t free_mem, total_mem;
    CHECK_CUDA(cudaMemGetInfo(&free_mem, &total_mem));
    size_t available_for_batch = static_cast<size_t>(free_mem * 0.8);
    size_t z_slice_size_bytes = (size_t)nPointsX * nPointsY * nEig * sizeof(double);

    int z_batch_size = nPointsZ;
    if (z_slice_size_bytes > 0 && total_size_valueReal > available_for_batch) {
        z_batch_size = available_for_batch / z_slice_size_bytes;
        if (z_batch_size == 0) z_batch_size = 1;
    }
    printf("Grid processing: Z-slices will be processed in batches of %d\n", z_batch_size);
    printf(" (Free device mem: %.2f GB, Grid size: %d x %d x %d (x %d eigs) = %.2f GB)\n",
           free_mem / 1e9, nPointsX, nPointsY, nPointsZ, nEig,
           total_size_valueReal / 1e9);


    // Batch buffer
    size_t batch_buffer_size_elems = (size_t)nPointsX * nPointsY * std::min(nPointsZ, z_batch_size) * nEig;
    DeviceBuffer<double> d_valueReal_out_batch(batch_buffer_size_elems);

    
    float totalKernelTime_ms = 0.0f;
    float totalD2HCopyTime_ms = 0.0f;

    for (int z_offset = 0; z_offset < nPointsZ; z_offset += z_batch_size) {
        int current_nPointsZ = std::min(z_batch_size, nPointsZ - z_offset);
        int total_points_in_batch = nPointsX * nPointsY * current_nPointsZ;
        if (total_points_in_batch == 0) continue;
        int grid_size = (total_points_in_batch + block_size - 1) / block_size;

        CHECK_CUDA(cudaEventRecord(startKernelTimer));

        evaluateKernel<<<grid_size, block_size, shared_mem_for_pass>>>(
            nPointsX, nPointsY, current_nPointsZ, z_offset, nEig, nOrb, nStos,
            maxNPows, maxNAlphas, nAtom, nCell, nEig_per_pass,
            d_origin.get(), d_gridVecs.get(), d_eigVecsReal.get(),
            d_coords.get(), d_species.get(), d_iStos.get(),
            d_sto_angMoms.get(), d_sto_nPows.get(), d_sto_nAlphas.get(),
            d_sto_cutoffsSq.get(), d_sto_coeffs.get(), d_sto_alphas.get(),
            d_valueReal_out_batch.get()
        );

        CHECK_CUDA(cudaEventRecord(endKernelTimer));
        CHECK_CUDA(cudaEventRecord(startCopyTimer));
        
        // The D2H copy will automatically block.
        for(int iEig=0; iEig < nEig; ++iEig) {
            size_t plane_size_bytes = (size_t)current_nPointsZ * nPointsY * nPointsX * sizeof(double);
            const double* d_src_ptr_eig = d_valueReal_out_batch.get() + (size_t)iEig * current_nPointsZ * nPointsY * nPointsX;
            double* h_dest_ptr_eig = h_valueReal_out + (size_t)iEig * nPointsZ * nPointsY * nPointsX + (size_t)z_offset * nPointsY * nPointsX;
            CHECK_CUDA(cudaMemcpy(h_dest_ptr_eig, d_src_ptr_eig, plane_size_bytes, cudaMemcpyDeviceToHost));
        }

        CHECK_CUDA(cudaEventRecord(endCopyTimer));
        CHECK_CUDA(cudaEventSynchronize(endCopyTimer));
        
        float iterKernel_ms, iterCopy_ms;
        CHECK_CUDA(cudaEventElapsedTime(&iterKernel_ms, startKernelTimer, endKernelTimer));
        CHECK_CUDA(cudaEventElapsedTime(&iterCopy_ms, endKernelTimer, endCopyTimer)); 
        
        totalKernelTime_ms += iterKernel_ms;
        totalD2HCopyTime_ms += iterCopy_ms;
    }
    
    CHECK_CUDA(cudaGetLastError());

    // Timing
    CHECK_CUDA(cudaEventRecord(endEverything));
    CHECK_CUDA(cudaEventSynchronize(endEverything));

    float timeEverything;
    CHECK_CUDA(cudaEventElapsedTime(&timeEverything, startEverything, endEverything));
    float overhead = timeEverything- totalKernelTime_ms - totalD2HCopyTime_ms;

    printf("\n--- Timing Results ---\n");
    printf("Kernel execution: %.2f ms (%.1f%%)\n", totalKernelTime_ms, (totalKernelTime_ms / timeEverything) * 100.0);
    printf("D2H Copy:         %.2f ms (%.1f%%)\n", totalD2HCopyTime_ms, (totalD2HCopyTime_ms / timeEverything) * 100.0);
    printf("Other:            %.2f ms (%.1f%%)\n", overhead, (overhead / timeEverything) * 100.0);
}
