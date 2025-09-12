/*-------------------------------------------------------------------------------------------------*
 *  DFTB+: general package for performing fast atomistic simulations                               *
 *  Copyright (C) 2006 - 2025  DFTB+ developers group                                              *
 *                                                                                                 *
 *  See the LICENSE file for terms of usage and distribution.                                      *
 *-------------------------------------------------------------------------------------------------*/
#include "host_logic.cuh"
#include <assert.h>
#include <omp.h>
#include <algorithm>
#include <cstdio>


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

    // Strided Memcpy3D could be used to squash this loop.
    for (int iEig = 0; iEig < calc->nEigOut; ++iEig) {
        // From: iEig-th slice of GPU batch buffer
        ptrdiff_t d_offset_bytes = (ptrdiff_t)(iEig * device_plane_size);

        // To: Global Z-position in the iEig-th slice of host buffer
        ptrdiff_t h_offset_bytes = (ptrdiff_t)(iEig * host_plane_size + ((size_t)z_offset_global * nPointsY * nPointsX) * output_num_size);

        CHECK_CUDA(cudaMemcpy((char*)h_dest_ptr + h_offset_bytes, (char*)d_src_ptr + d_offset_bytes, device_plane_size,
            cudaMemcpyDeviceToHost));
    }
}


