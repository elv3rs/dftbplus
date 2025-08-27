/*-------------------------------------------------------------------------------------------------*
 *  DFTB+: general package for performing fast atomistic simulations                               *
 *  Copyright (C) 2006 - 2025  DFTB+ developers group                                              *
 *                                                                                                 *
 *  See the LICENSE file for terms of usage and distribution.                                      *
 *-------------------------------------------------------------------------------------------------*/
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
    const double* coords; // [3][nAtom][nCell]
    const int* species;   // [nAtom]
    const int* iStos;     // [nSpecies+1]
} SystemParams;

// Additional System information describing periodic boundary conditions
typedef struct {
    bool isPeriodic;
    const double* latVecs;        // [3][3]
    const double* recVecs2pi;     // [3][3]
    const int* kIndexes;          // [nEig]
    const cuDoubleComplex* phases;// [nCell][nEigIn]
} PeriodicParams;

// Basis parameters describing orbitals
typedef struct {
    bool useRadialLut;
    int nStos;
    int nLutPoints;

    double inverseLutStep;
    const double* lutGridValues;  // [nStos][nLutPoints]

    int maxNPows;
    int maxNAlphas;
    const int* angMoms;      // [nStos]
    const int* nPows;        // [nStos]
    const int* nAlphas;      // [nStos]
    const double* cutoffsSq; // [nStos]
    const double* coeffs;    // [maxNPows][maxNAlphas][nStos]
    const double* alphas;    // [maxNAlphas][nStos]
} StoBasisParams;


// Coefficient Input and Control Flags
typedef struct {
    bool isRealInput; 
    bool isRealOutput;
    bool calcAtomicDensity;
    bool calcTotalChrg;
    int nEigIn;
    int nEigOut;                        // calcTotalChrg ? 1 : nEigIn
    const double* eigVecsReal;          // [nOrb][nEigIn]
    const cuDoubleComplex* eigVecsCmpl; // [nOrb][nEigIn]
    double* valueReal_out;              // [nPointsX][nPointsY][nPointsZ][nEigOut]
    cuDoubleComplex* valueCmpl_out;     // [nPointsX][nPointsY][nPointsZ][nEigOut]
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
