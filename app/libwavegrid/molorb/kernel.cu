#include <cuda_runtime.h>
#include <omp.h>
#include <cstdio>
#include <cmath>
#include <algorithm> // For std::min

#include "kernel.cuh"

// Helper macro for robust CUDA calls
#define CHECK_CUDA(call) do { \
    cudaError_t err = call; \
    if (err != cudaSuccess) { \
        fprintf(stderr, "CUDA Error in %s at line %d: %s\n", __FILE__, __LINE__, cudaGetErrorString(err)); \
        exit(EXIT_FAILURE); \
    } \
} while (0)


// Helper macros for column-major (Fortran-style) indexing.
#define IDX2F(i, j, lda) ((j) * (size_t)(lda) + (i))
#define IDX3F(i, j, k, lda, ldb) (((k) * (size_t)(ldb) + (j)) * (size_t)(lda) + (i))
#define IDX4F(i, j, k, l, lda, ldb, ldc) ((((l) * (size_t)(ldc) + (k)) * (size_t)(ldb) + (j)) * (size_t)(lda) + (i))



// =========================================================================
//  CUDA Device Functions
// =========================================================================
__device__ __forceinline__ double realTessY_device_opt(int ll, int mm, const double* diff, double inv_r, double inv_r2) {
    const double x = diff[0];
    const double y = diff[1];
    const double z = diff[2];

    const double x_r = x * inv_r;
    const double y_r = y * inv_r;
    const double z_r = z * inv_r;

    switch (ll) {
        case 0: return 0.2820947917738782; // 1/sqrt(4*PI)
        case 1:
            switch (mm) {
                case -1: return 0.4886025119029198 * y_r;
                case  0: return 0.4886025119029198 * z_r;
                case  1: return 0.4886025119029198 * x_r;
            }
            break;
        case 2:
            {
                const double xx_r2 = x * x * inv_r2;
                const double yy_r2 = y * y * inv_r2;
                const double zz_r2 = z * z * inv_r2;
                const double xy_r2 = x * y * inv_r2;
                const double yz_r2 = y * z * inv_r2;
                const double xz_r2 = x * z * inv_r2;
                switch (mm) {
                    case -2: return 1.092548430592079 * xy_r2;
                    case -1: return 1.092548430592079 * yz_r2;
                    case  0: return -0.3153915652525200 * (-2.0 * zz_r2 + xx_r2 + yy_r2);
                    case  1: return 1.092548430592079 * xz_r2;
                    case  2: return 0.5462742152960395 * (xx_r2 - yy_r2);
                }
            }
            break;
        case 3:
            {
                const double x_r3 = x_r * inv_r2;
                const double y_r3 = y_r * inv_r2;
                const double z_r3 = z_r * inv_r2;
                switch (mm) {
                    case -3: return 0.5900435899266435 * y_r3 * (3.0 * x * x - y * y);
                    case -2: return 2.890611442640554  * x_r3 * y * z;
                    case -1: return -0.4570457994644658 * y_r3 * (-4.0 * z * z + x * x + y * y);
                    case  0: return -0.3731763325901155 * z_r3 * (-2.0 * z * z + 3.0 * x * x + 3.0 * y * y);
                    case  1: return -0.4570457994644658 * x_r3 * (-4.0 * z * z + x * x + y * y);
                    case  2: return 1.445305721320277  * z_r3 * (x * x - y * y);
                    case  3: return 0.5900435899266435 * x_r3 * (x * x - 3.0 * y * y);
                }
            }
            break;
    }
    return 0.0;
}

__device__ __forceinline__ double getRadialValue(
    double r, int iL, int nPows, int nAlphas,
    const double* coeffs, const double* alphas)
{
    constexpr int STO_TMP_POWS_SIZE = 16;
    double sto_tmp_pows[STO_TMP_POWS_SIZE];
    double sto_tmp_rexp = (iL == 0 && r < 1.0e-12) ? 1.0 : pow(r, iL);
    for (int ii = 0; ii < nPows; ++ii) {
        sto_tmp_pows[ii] = sto_tmp_rexp;
        sto_tmp_rexp *= r;
    }
    double radialVal = 0.0;
    for (int ii = 0; ii < nAlphas; ++ii) {
        double term = 0.0;
        for (int jj = 0; jj < nPows; ++jj) {
            term += coeffs[IDX3F(jj, ii, iL, nPows, nAlphas)] * sto_tmp_pows[jj];
        }
        radialVal += term * exp(alphas[IDX2F(ii, iL, nAlphas)] * r);
    }
    return radialVal;
}



// =========================================================================
//  CUDA Kernel
// =========================================================================

__global__ void evaluateKernel(
    int nPointsX, int nPointsY, int nPointsZ_batch, int z_offset, int nEig, int nOrb, int nStos,
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
        // This is a trade-off to keep the accumulation in fast shared memory.
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
                        r, iL, sto_nPows[iOrb], sto_nAlphas[iOrb],
                        &sto_coeffs[IDX2F(iOrb * maxNPows, 0, maxNPows)],
                        &sto_alphas[IDX2F(0, iOrb, maxNAlphas)]
                    );

                    // Only calculate inverse once 
                    double inv_r = (r < 1.e-12) ? 0.0 : 1.0 / r;
                    double inv_r2 = inv_r * inv_r;

                    for (int iM = -iL; iM <= iL; ++iM) {
                        double val = radialVal * realTessY_device_opt(iL, iM, diff, inv_r, inv_r2);
                        
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

    // Allocation of constant data on device
    double *d_origin, *d_gridVecs, *d_eigVecsReal, *d_coords, *d_sto_cutoffsSq, *d_sto_coeffs, *d_sto_alphas;
    int *d_species, *d_iStos, *d_sto_angMoms, *d_sto_nPows, *d_sto_nAlphas;

    size_t size_origin = 3 * sizeof(double);
    size_t size_gridVecs = 9 * sizeof(double);
    size_t size_eigVecsReal = (size_t)nOrb * nEig * sizeof(double);
    size_t size_coords = 3 * (size_t)nAtom * nCell * sizeof(double);
    size_t size_species = nAtom * sizeof(int);
    size_t size_iStos = (nSpecies + 1) * sizeof(int);
    size_t size_sto_angMoms = nStos * sizeof(int);
    size_t size_sto_nPows = nStos * sizeof(int);
    size_t size_sto_nAlphas = nStos * sizeof(int);
    size_t size_sto_cutoffsSq = nStos * sizeof(double);
    size_t size_sto_coeffs = (size_t)maxNPows * maxNAlphas * nStos * sizeof(double);
    size_t size_sto_alphas = (size_t)maxNAlphas * nStos * sizeof(double);
    size_t total_size_valueReal = (size_t)nPointsX * nPointsY * nPointsZ * nEig * sizeof(double);

    cudaEvent_t startInit, endInit, startKernel, endKernel, startFinalise, endFinalise;
    CHECK_CUDA(cudaEventCreate(&startInit));
    CHECK_CUDA(cudaEventCreate(&endInit));
    CHECK_CUDA(cudaEventCreate(&startKernel));
    CHECK_CUDA(cudaEventCreate(&endKernel));
    CHECK_CUDA(cudaEventCreate(&startFinalise));
    CHECK_CUDA(cudaEventCreate(&endFinalise));
    
    CHECK_CUDA(cudaEventRecord(startInit));
    CHECK_CUDA(cudaMalloc(&d_origin, size_origin));
    CHECK_CUDA(cudaMalloc(&d_gridVecs, size_gridVecs));
    CHECK_CUDA(cudaMalloc(&d_eigVecsReal, size_eigVecsReal));
    CHECK_CUDA(cudaMalloc(&d_coords, size_coords));
    CHECK_CUDA(cudaMalloc(&d_species, size_species));
    CHECK_CUDA(cudaMalloc(&d_iStos, size_iStos));
    CHECK_CUDA(cudaMalloc(&d_sto_angMoms, size_sto_angMoms));
    CHECK_CUDA(cudaMalloc(&d_sto_nPows, size_sto_nPows));
    CHECK_CUDA(cudaMalloc(&d_sto_nAlphas, size_sto_nAlphas));
    CHECK_CUDA(cudaMalloc(&d_sto_cutoffsSq, size_sto_cutoffsSq));
    CHECK_CUDA(cudaMalloc(&d_sto_coeffs, size_sto_coeffs));
    CHECK_CUDA(cudaMalloc(&d_sto_alphas, size_sto_alphas));
    CHECK_CUDA(cudaMemcpy(d_origin, h_origin, size_origin, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_gridVecs, h_gridVecs, size_gridVecs, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_eigVecsReal, h_eigVecsReal, size_eigVecsReal, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_coords, h_coords, size_coords, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_species, h_species, size_species, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_iStos, h_iStos, size_iStos, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_sto_angMoms, h_sto_angMoms, size_sto_angMoms, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_sto_nPows, h_sto_nPows, size_sto_nPows, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_sto_nAlphas, h_sto_nAlphas, size_sto_nAlphas, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_sto_cutoffsSq, h_sto_cutoffsSq, size_sto_cutoffsSq, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_sto_coeffs, h_sto_coeffs, size_sto_coeffs, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_sto_alphas, h_sto_alphas, size_sto_alphas, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaEventRecord(endInit));
    CHECK_CUDA(cudaEventSynchronize(endInit));


    // BATCHING & DYNAMIC KERNEL CONFIGURATION
    CHECK_CUDA(cudaEventRecord(startKernel));

    int block_size = 256;
    int nEig_per_pass;
    size_t shared_mem_for_pass;

    int deviceId;
    CHECK_CUDA(cudaGetDevice(&deviceId));
    cudaDeviceProp prop;
    CHECK_CUDA(cudaGetDeviceProperties(&prop, deviceId));

    // Leave a 5% buffer for compiler's internal use of shared memory
    size_t available_shared = prop.sharedMemPerBlock * 0.95;
    nEig_per_pass = available_shared / (block_size * sizeof(double));
    if (nEig_per_pass == 0) nEig_per_pass = 1; // Must process at least one at a time
    if (nEig_per_pass > nEig) nEig_per_pass = nEig;

    shared_mem_for_pass = (size_t)nEig_per_pass * block_size * sizeof(double);

    printf("Kernel Configuration:\n");
    printf("  Block size: %d threads\n", block_size);
    printf("  Eigenstates per pass: %d (out of %d total)\n", nEig_per_pass, nEig);
    printf("  Shared memory per block: %zu bytes (Device max: %zu bytes)\n", shared_mem_for_pass, prop.sharedMemPerBlock);


    // Batching logic for Z-slices (to handle very large grids)
    size_t free_mem, total_mem;
    CHECK_CUDA(cudaMemGetInfo(&free_mem, &total_mem));
    size_t available_for_batch = static_cast<size_t>(free_mem * 0.8);
    size_t z_slice_size_bytes = (size_t)nPointsX * nPointsY * nEig * sizeof(double);

    int z_batch_size = nPointsZ; // Default to processing all at once
    if (z_slice_size_bytes > 0 && total_size_valueReal > available_for_batch) {
        z_batch_size = available_for_batch / z_slice_size_bytes;
        if (z_batch_size == 0) z_batch_size = 1; 
    }
    printf("Grid processing: Z-slices will be processed in batches of %d\n", z_batch_size);

    double* d_valueReal_out_batch;
    size_t batch_buffer_size_bytes = (size_t)nPointsX * nPointsY * std::min(nPointsZ, z_batch_size) * nEig * sizeof(double);
    CHECK_CUDA(cudaMalloc(&d_valueReal_out_batch, batch_buffer_size_bytes));
    
    for (int z_offset = 0; z_offset < nPointsZ; z_offset += z_batch_size) {
        int current_nPointsZ = std::min(z_batch_size, nPointsZ - z_offset);
        int total_points_in_batch = nPointsX * nPointsY * current_nPointsZ;
        if (total_points_in_batch == 0) continue;
        int grid_size = (total_points_in_batch + block_size - 1) / block_size;

        evaluateKernel<<<grid_size, block_size, shared_mem_for_pass>>>(
            nPointsX, nPointsY, current_nPointsZ, z_offset, nEig, nOrb, nStos,
            maxNPows, maxNAlphas, nAtom, nCell, nEig_per_pass,
            d_origin, d_gridVecs, d_eigVecsReal,
            d_coords, d_species, d_iStos,
            d_sto_angMoms, d_sto_nPows, d_sto_nAlphas,
            d_sto_cutoffsSq, d_sto_coeffs, d_sto_alphas,
            d_valueReal_out_batch
        );

        size_t width_bytes = (size_t)nPointsX * nPointsY * current_nPointsZ * sizeof(double);
        
        // We now copy one full eigenstate plane at a time due to the complex memory layouts.
        for(int iEig=0; iEig < nEig; ++iEig) {
            const double* d_src_ptr_eig = d_valueReal_out_batch + iEig * current_nPointsZ * nPointsY * nPointsX;
            double* h_dest_ptr_eig = h_valueReal_out + iEig * nPointsZ * nPointsY * nPointsX + (size_t)z_offset * nPointsY * nPointsX;
            CHECK_CUDA(cudaMemcpy(h_dest_ptr_eig, d_src_ptr_eig, width_bytes, cudaMemcpyDeviceToHost));
        }
    }
    
    CHECK_CUDA(cudaGetLastError());
    CHECK_CUDA(cudaDeviceSynchronize());
    CHECK_CUDA(cudaEventRecord(endKernel));
    CHECK_CUDA(cudaEventSynchronize(endKernel));

    // Free Device Memory 
    CHECK_CUDA(cudaEventRecord(startFinalise));
    CHECK_CUDA(cudaFree(d_valueReal_out_batch));
    CHECK_CUDA(cudaFree(d_origin));
    CHECK_CUDA(cudaFree(d_gridVecs));
    CHECK_CUDA(cudaFree(d_eigVecsReal));
    CHECK_CUDA(cudaFree(d_coords));
    CHECK_CUDA(cudaFree(d_species));
    CHECK_CUDA(cudaFree(d_iStos));
    CHECK_CUDA(cudaFree(d_sto_angMoms));
    CHECK_CUDA(cudaFree(d_sto_nPows));
    CHECK_CUDA(cudaFree(d_sto_nAlphas));
    CHECK_CUDA(cudaFree(d_sto_cutoffsSq));
    CHECK_CUDA(cudaFree(d_sto_coeffs));
    CHECK_CUDA(cudaFree(d_sto_alphas));
    CHECK_CUDA(cudaEventRecord(endFinalise));
    CHECK_CUDA(cudaEventSynchronize(endFinalise));

    float timeInit, timeKernel, timeFinalise;
    CHECK_CUDA(cudaEventElapsedTime(&timeInit, startInit, endInit));
    CHECK_CUDA(cudaEventElapsedTime(&timeKernel, startKernel, endKernel));
    CHECK_CUDA(cudaEventElapsedTime(&timeFinalise, startFinalise, endFinalise));

    printf("Initialization time: %.2f ms\n", timeInit);
    printf("Kernel execution & D2H copy time: %.2f ms\n", timeKernel);
    printf("Finalization time: %.2f ms\n", timeFinalise);
}
