#include <cuda_runtime.h>
#include <cstdio>
#include <cmath>

#include "kernel.cuh"

// Helper macro for robust CUDA calls
#define CHECK_CUDA(call) do { \
    cudaError_t err = call; \
    if (err != cudaSuccess) { \
        fprintf(stderr, "CUDA Error in %s at line %d: %s\n", __FILE__, __LINE__, cudaGetErrorString(err)); \
        exit(EXIT_FAILURE); \
    } \
} while (0)

// =========================================================================
//  CUDA Device Functions
// =========================================================================

// Port of the realTessY function for the GPU
__device__ double realTessY_device(int ll, int mm, const double* diff, double rr) {
    constexpr double epsilon = 1.0e-12;

    if (rr < epsilon && ll != 0) {
        return 0.0;
    }

    const double xx = diff[0];
    const double yy = diff[1];
    const double zz = diff[2];

    switch (ll) {
        case 0: return 0.2820947917738782;
        case 1:
            switch (mm) {
                case -1: return 0.4886025119029198 * yy / rr;
                case 0:  return 0.4886025119029198 * zz / rr;
                case 1:  return 0.4886025119029198 * xx / rr;
            }
            break;
        case 2:
            switch (mm) {
                case -2: return 1.092548430592079 * xx * yy / (rr * rr);
                case -1: return 1.092548430592079 * yy * zz / (rr * rr);
                case 0:  return -0.3153915652525200 * (-2.0 * zz * zz + xx * xx + yy * yy) / (rr * rr);
                case 1:  return 1.092548430592079 * xx * zz / (rr * rr);
                case 2:  return 0.5462742152960395 * (xx * xx - yy * yy) / (rr * rr);
            }
            break;
        case 3:
            switch (mm) {
                case -3: return 0.5900435899266435 * yy * (3.0 * xx * xx - yy * yy) / (rr * rr * rr);
                case -2: return 2.890611442640554 * xx * yy * zz / (rr * rr * rr);
                case -1: return -0.4570457994644658 * yy * (-4.0 * zz * zz + xx * xx + yy * yy) / (rr * rr * rr);
                case 0:  return -0.3731763325901155 * zz * (-2.0 * zz * zz + 3.0 * xx * xx + 3.0 * yy * yy) / (rr * rr * rr);
                case 1:  return -0.4570457994644658 * xx * (-4.0 * zz * zz + xx * xx + yy * yy) / (rr * rr * rr);
                case 2:  return 1.445305721320277 * zz * (xx * xx - yy * yy) / (rr * rr * rr);
                case 3:  return 0.5900435899266435 * xx * (xx * xx - 3.0 * yy * yy) / (rr * rr * rr);
            }
            break;
    }
    return 0.0; // Should not be reached for valid ll, mm
}


// =========================================================================
//  CUDA Kernel
// =========================================================================

// Helper macros for column-major (Fortran-style) indexing.
#define IDX2F(i, j, lda) ((j) * (size_t)(lda) + (i))
#define IDX3F(i, j, k, lda, ldb) (((k) * (size_t)(ldb) + (j)) * (size_t)(lda) + (i))
#define IDX4F(i, j, k, l, lda, ldb, ldc) ((((l) * (size_t)(ldc) + (k)) * (size_t)(ldb) + (j)) * (size_t)(lda) + (i))

__global__ void evaluateKernel(
    int nPointsX, int nPointsY, int nPointsZ, int nEig, int nOrb, int nStos,
    int maxNPows, int maxNAlphas, int nAtom, int nCell,
    const double* origin, const double* gridVecs, const double* eigVecsReal,
    const double* coords, const int* species, const int* iStos,
    const int* sto_angMoms, const int* sto_nPows, const int* sto_nAlphas,
    const double* sto_cutoffs, const double* sto_coeffs, const double* sto_alphas,
    double* valueReal_out)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int total_points = nPointsX * nPointsY * nPointsZ;
    if (idx >= total_points) return;

    int i1 = idx % nPointsX;
    int i2 = (idx / nPointsX) % nPointsY;
    int i3 = idx / (nPointsX * nPointsY);

    double xyz[3];
    xyz[0] = origin[0] + i1 * gridVecs[IDX2F(0, 0, 3)] + i2 * gridVecs[IDX2F(0, 1, 3)] + i3 * gridVecs[IDX2F(0, 2, 3)];
    xyz[1] = origin[1] + i1 * gridVecs[IDX2F(1, 0, 3)] + i2 * gridVecs[IDX2F(1, 1, 3)] + i3 * gridVecs[IDX2F(1, 2, 3)];
    xyz[2] = origin[2] + i1 * gridVecs[IDX2F(2, 0, 3)] + i2 * gridVecs[IDX2F(2, 1, 3)] + i3 * gridVecs[IDX2F(2, 2, 3)];


    for (int i = 0; i < nEig; ++i) {
        valueReal_out[IDX4F(i1, i2, i3, i, nPointsX, nPointsY, nPointsZ)] = 0.0;
    
    }

    for (int iCell = 0; iCell < nCell; ++iCell) {
        int orbital_idx_counter = 0; // 'ind' in Fortran is reset for each cell
        for (int iAtom = 0; iAtom < nAtom; ++iAtom) {
            // Note: In Fortran, iStos is 1-based. C-arrays are 0-based.
            int iSpecies = species[iAtom] - 1;

            double diff[3];
            diff[0] = xyz[0] - coords[IDX3F(0, iAtom, iCell, 3, nAtom)];
            diff[1] = xyz[1] - coords[IDX3F(1, iAtom, iCell, 3, nAtom)];
            diff[2] = xyz[2] - coords[IDX3F(2, iAtom, iCell, 3, nAtom)];
            double r = sqrt(diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]);

            for (int iOrb = iStos[iSpecies] - 1; iOrb < iStos[iSpecies + 1] - 1; ++iOrb) {
                int iL = sto_angMoms[iOrb];
                if (r > sto_cutoffs[iOrb]) {
                    orbital_idx_counter += 2 * iL + 1;
                    continue;
                }

                double sto_tmp_rexp;
                if (iL == 0 && r < 1.0e-12) {
                    sto_tmp_rexp = 1.0;
                } else {
                    sto_tmp_rexp = pow(r, iL);
                }
                // Guaranteed to be large enough
                constexpr int STO_TMP_POWS_SIZE = 16;
                double sto_tmp_pows[STO_TMP_POWS_SIZE];
                int current_sto_nPows = sto_nPows[iOrb];
                for (int ii = 0; ii < current_sto_nPows; ++ii) {
                    sto_tmp_pows[ii] = sto_tmp_rexp;
                    sto_tmp_rexp *= r;
                }

                double radialVal = 0.0;
                int current_sto_nAlphas = sto_nAlphas[iOrb];
                for (int ii = 0; ii < current_sto_nAlphas; ++ii) {
                    double term = 0.0;
                    for (int jj = 0; jj < current_sto_nPows; ++jj) {
                        term += sto_coeffs[IDX3F(jj, ii, iOrb, maxNPows, maxNAlphas)] * sto_tmp_pows[jj];
                    }
                    radialVal += term * exp(sto_alphas[IDX2F(ii, iOrb, maxNAlphas)] * r);
                }


                for (int iM = -iL; iM <= iL; ++iM) {
                    double val = radialVal * realTessY_device(iL, iM, diff, r);


                    for (int iEig = 0; iEig < nEig; ++iEig) {
                        size_t out_idx = IDX4F(i1, i2, i3, iEig, nPointsX, nPointsY, nPointsZ);
                        size_t eig_idx = IDX2F(orbital_idx_counter, iEig, nOrb);
                        valueReal_out[out_idx] += val * eigVecsReal[eig_idx];
                    }


                    orbital_idx_counter++;
                }
            }
        }
    }


}


// =========================================================================
//  C++ Host Interface (callable from C/Fortran)
// =========================================================================
extern "C" void evaluate_on_device_c(
    const int* nPointsX, const int* nPointsY, const int* nPointsZ, const int* nEig,
    const int* nOrb, const int* nStos, const int* maxNPows, const int* maxNAlphas,
    const int* nAtom, const int* nCell, const int* nSpecies,
    const double* h_origin, const double* h_gridVecs, const double* h_eigVecsReal,
    const double* h_coords, const int* h_species, const int* h_iStos,
    const int* h_sto_angMoms, const int* h_sto_nPows, const int* h_sto_nAlphas,
    const double* h_sto_cutoffs, const double* h_sto_coeffs, const double* h_sto_alphas,
    double* h_valueReal_out)
{
    // Dereference pointers from Fortran
    int _nPointsX = *nPointsX, _nPointsY = *nPointsY, _nPointsZ = *nPointsZ, _nEig = *nEig;
    int _nOrb = *nOrb, _nSpecies = *nSpecies, _nStos = *nStos, _maxNPows = *maxNPows, _maxNAlphas = *maxNAlphas;
    int _nAtom = *nAtom, _nCell = *nCell;

    // Allocate Device Memory
    double *d_origin, *d_gridVecs, *d_eigVecsReal, *d_coords, *d_sto_cutoffs, *d_sto_coeffs, *d_sto_alphas, *d_valueReal_out;
    int *d_species, *d_iStos, *d_sto_angMoms, *d_sto_nPows, *d_sto_nAlphas;

    size_t size_origin = 3 * sizeof(double);
    size_t size_gridVecs = 9 * sizeof(double);
    size_t size_eigVecsReal = (size_t)_nOrb * _nEig * sizeof(double);
    size_t size_coords = 3 * (size_t)_nAtom * _nCell * sizeof(double);
    size_t size_species = _nAtom * sizeof(int);
    size_t size_iStos = (_nSpecies + 1) * sizeof(int);
    size_t size_sto_angMoms = _nStos * sizeof(int);
    size_t size_sto_nPows = _nStos * sizeof(int);
    size_t size_sto_nAlphas = _nStos * sizeof(int);
    size_t size_sto_cutoffs = _nStos * sizeof(double);
    size_t size_sto_coeffs = (size_t)_maxNPows * _maxNAlphas * _nStos * sizeof(double);
    size_t size_sto_alphas = (size_t)_maxNAlphas * _nStos * sizeof(double);
    size_t size_valueReal = (size_t)_nPointsX * _nPointsY * _nPointsZ * _nEig * sizeof(double);

    CHECK_CUDA(cudaMalloc(&d_origin, size_origin));
    CHECK_CUDA(cudaMalloc(&d_gridVecs, size_gridVecs));
    CHECK_CUDA(cudaMalloc(&d_eigVecsReal, size_eigVecsReal));
    CHECK_CUDA(cudaMalloc(&d_coords, size_coords));
    CHECK_CUDA(cudaMalloc(&d_species, size_species));
    CHECK_CUDA(cudaMalloc(&d_iStos, size_iStos));
    CHECK_CUDA(cudaMalloc(&d_sto_angMoms, size_sto_angMoms));
    CHECK_CUDA(cudaMalloc(&d_sto_nPows, size_sto_nPows));
    CHECK_CUDA(cudaMalloc(&d_sto_nAlphas, size_sto_nAlphas));
    CHECK_CUDA(cudaMalloc(&d_sto_cutoffs, size_sto_cutoffs));
    CHECK_CUDA(cudaMalloc(&d_sto_coeffs, size_sto_coeffs));
    CHECK_CUDA(cudaMalloc(&d_sto_alphas, size_sto_alphas));
    CHECK_CUDA(cudaMalloc(&d_valueReal_out, size_valueReal));

    // Copy Data from Host to Device
    CHECK_CUDA(cudaMemcpy(d_origin, h_origin, size_origin, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_gridVecs, h_gridVecs, size_gridVecs, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_eigVecsReal, h_eigVecsReal, size_eigVecsReal, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_coords, h_coords, size_coords, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_species, h_species, size_species, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_iStos, h_iStos, size_iStos, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_sto_angMoms, h_sto_angMoms, size_sto_angMoms, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_sto_nPows, h_sto_nPows, size_sto_nPows, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_sto_nAlphas, h_sto_nAlphas, size_sto_nAlphas, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_sto_cutoffs, h_sto_cutoffs, size_sto_cutoffs, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_sto_coeffs, h_sto_coeffs, size_sto_coeffs, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_sto_alphas, h_sto_alphas, size_sto_alphas, cudaMemcpyHostToDevice));
    
    // Kernel Launch
    int total_points = _nPointsX * _nPointsY * _nPointsZ;
    int block_size = 256;
    int grid_size = (total_points + block_size - 1) / block_size;
    
    evaluateKernel<<<grid_size, block_size>>>(
        _nPointsX, _nPointsY, _nPointsZ, _nEig, _nOrb, _nStos,
        _maxNPows, _maxNAlphas, _nAtom, _nCell,
        d_origin, d_gridVecs, d_eigVecsReal,
        d_coords, d_species, d_iStos,
        d_sto_angMoms, d_sto_nPows, d_sto_nAlphas,
        d_sto_cutoffs, d_sto_coeffs, d_sto_alphas,
        d_valueReal_out
    );
    
    CHECK_CUDA(cudaGetLastError());
    CHECK_CUDA(cudaDeviceSynchronize());
    
    // Copy Result from Device to Host
    CHECK_CUDA(cudaMemcpy(h_valueReal_out, d_valueReal_out, size_valueReal, cudaMemcpyDeviceToHost));
    
    // Free Device Memory
    CHECK_CUDA(cudaFree(d_origin));
    CHECK_CUDA(cudaFree(d_gridVecs));
    CHECK_CUDA(cudaFree(d_eigVecsReal));
    CHECK_CUDA(cudaFree(d_coords));
    CHECK_CUDA(cudaFree(d_species));
    CHECK_CUDA(cudaFree(d_iStos));
    CHECK_CUDA(cudaFree(d_sto_angMoms));
    CHECK_CUDA(cudaFree(d_sto_nPows));
    CHECK_CUDA(cudaFree(d_sto_nAlphas));
    CHECK_CUDA(cudaFree(d_sto_cutoffs));
    CHECK_CUDA(cudaFree(d_sto_coeffs));
    CHECK_CUDA(cudaFree(d_sto_alphas));
    CHECK_CUDA(cudaFree(d_valueReal_out));
}
