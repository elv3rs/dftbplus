/*-------------------------------------------------------------------------------------------------*
 *  DFTB+: general package for performing fast atomistic simulations                               *
 *  Copyright (C) 2006 - 2025  DFTB+ developers group                                              *
 *                                                                                                 *
 *  See the LICENSE file for terms of usage and distribution.                                      *
 *-------------------------------------------------------------------------------------------------*/
// This file contains the D2H copy of the computed batches back to the final host array.
#include "host_logic.cuh"


/** 
 * @brief Copies a batch of computed eigenstates from device to host memory.
 *
 * The D2H copy will automatically block/ synchronize the kernel for this batch.
 * This could be improved by using streams / cudaMemcpyAsync, but currently is not a bottleneck.
 * The output array is of fortran shape (x,y,z,nEigOut), thus z-slices are not contiguous.
 * We slice on Z instead of nEigOut to save on a little computation in the kernel.
 *
 * @param d_src_ptr        Pointer to the source data on the device (GPU).
 * @param h_dest_ptr       Pointer to the destination data on the host (CPU).
 * @param nPointsX         Number of grid points in the X dimension.
 * @param nPointsY         Number of grid points in the Y dimension.
 * @param nPointsZ         Total number of grid points in the Z dimension.
 * @param z_per_batch      Number of Z-slices processed in the current batch.
 * @param z_offset_global  Global Z-offset indicating where this batch fits in the full grid.
 * @param calc             Pointer to the CalculationParams structure containing calculation flags.
 */
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


