#ifndef KERNEL_CUH_
#define KERNEL_CUH_

#include <cuComplex.h>

#ifdef __cplusplus
extern "C" {
#endif

// Calculation grid
typedef struct {
    int nPointsX;
    int nPointsY;
    int nPointsZ;
    const double* origin;   // [3]
    const double* gridVecs; // [3][3]
} GridParams;

typedef struct {
    int nAtom;
    int nCell;
    int nSpecies;
    int nOrb;
    const double* coords;   // [3][nAtom][nCell]
    const int* species;     // [nAtom]
    const int* iStos;       // [nSpecies+1]
} SystemParams;

// Additional System information describing periodic boundary conditions
typedef struct {
    int isPeriodic;
    const double* latVecs;        // [3][3]
    const double* recVecs2p;      // [3][3]
    const int* kIndexes;          // [nEig]
    const cuDoubleComplex* phases;// [nCell][nEigIn]
} PeriodicParams;

// Basis parameters describing orbitals
typedef struct {
    int nStos;
    int maxNPows;
    int maxNAlphas;
    const int* sto_angMoms;       // [nStos]
    const int* sto_nPows;         // [nStos]
    const int* sto_nAlphas;       // [nStos]
    const double* sto_cutoffsSq;  // [nStos]
    const double* sto_coeffs;     // [maxNPows][maxNAlphas][nStos]
    const double* sto_alphas;     // [maxNAlphas][nStos]
} StoBasisParams;


// Coefficient Input and Control Flags
typedef struct {
    int nEigIn;
    int nEigOut; // accDensity ? 1 : nEigIn
    int isRealInput; 
    int isDensityCalc;
    int accDensity;
    const double* eigVecsReal;          // [nOrb][nEigIn]
    const cuDoubleComplex* eigVecsCmpl; // [nOrb][nEigIn]
    double* valueReal_out;              // [nPointsX][nPointsY][nPointsZ][nEigOut]
    cuDoubleComplex* valueCmpl_out    ; // [nPointsX][nPointsY][nPointsZ][nEigOut]
} CalculationParams;


void evaluate_on_device_c(
    const GridParams* grid,
    const SystemParams* system,
    const PeriodicParams* periodic,
    const StoBasisParams* basis,
    const CalculationParams* calc
);


#ifdef __cplusplus
}
#endif

#endif // KERNEL_CUH_
