#ifndef KERNEL_CUH_
#define KERNEL_CUH_

// Use extern "C" to make this function callable from Fortran
#ifdef __cplusplus
extern "C" {
#endif

void evaluate_on_device_c(
    const int* nPointsX, const int* nPointsY, const int* nPointsZ, const int* nEig,
    const int* nOrb, const int* nStos, const int* maxNPows, const int* maxNAlphas,
    const int* nAtom, const int* nCell, const int* nSpecies,
    const double* origin, const double* gridVecs, const double* eigVecsReal,
    const double* coords, const int* species, const int* iStos,
    const int* sto_angMoms, const int* sto_nPows, const int* sto_nAlphas,
    const double* sto_cutoffs, const double* sto_coeffs, const double* sto_alphas,
    double* valueReal_out);

#ifdef __cplusplus
}
#endif

#endif // KERNEL_CUH_
