#ifndef KERNEL_CUH_
#define KERNEL_CUH_

// Use extern "C" to make this function callable from Fortran.
// This interface remains unchanged from the original.
#ifdef __cplusplus
extern "C" {
#endif


void evaluate_on_device_c(
    const int nPointsX, const int nPointsY, const int nPointsZ,
    const int nEigIn, const int nEigOut, const int nOrb, const int nStos,
    const int maxNPows, const int maxNAlphas,
    const int nAtom, const int nCell, const int nSpecies,
    const int isReal, const int isPeriodic, const int isDensityCalc,
    const double* origin,               // [3]
    const double* gridVecs,             // [3][3]
    const double* eigVecsReal,          // [nOrb][nEig]
    const cuDoubleComplex* eigVecsCmpl, // [nOrb][nEig]
    const double* coords,               // [3][nAtom][nCell]
    const int* species,                 // [nAtom]
    const int* iStos,                   // [nSpecies+1]
    const double* latVecs,              // [3][3]
    const double* recVecs2p,            // [3][3]
    const int* kIndexes,                // [nEig]
    const cuDoubleComplex* phases,      // [nCell][nEig]
    const int* sto_angMoms,             // [nStos]
    const int* sto_nPows,               // [nStos]
    const int* sto_nAlphas,             // [nStos]
    const double* sto_cutoffsSq,        // [nStos]
    const double* sto_coeffs,           // [maxNPows][maxNAlphas][nStos]
    const double* sto_alphas,           // [maxNAlphas][nStos]
    double* valueReal_out,              // [nPointsX][nPointsY][nPointsZ][nEig]
    cuDoubleComplex* valueCmpl_out      // [nPointsX][nPointsY][nPointsZ][nEig]
);




#ifdef __cplusplus
}
#endif

#endif // KERNEL_CUH_
