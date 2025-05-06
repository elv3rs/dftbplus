#include <vector>
#include <cmath>
#include <numeric> // For std::fill
#include <limits>
#include <cassert>
// #include <iostream> // Only if needed for debugging

// Constants
// const double PI = 3.14159265358979323846; // Original, ensure no redefinition if included elsewhere

// Anonymous namespace for internal constants and helper functions not part of the API
namespace {

// Maximum number of polynomial terms for STO radial part.
// User mentioned "I havent seen any with more than 15". Set to 16 as a safe upper bound.
const int MAX_NPOW_CONST = 16;

// Helper Functions (kept inline as they are small)
inline double norm2_sq(double x, double y, double z) {
    return x * x + y * y + z * z;
}

inline double norm2(double x, double y, double z) {
    return std::sqrt(norm2_sq(x, y, z));
}

inline void matvec_mult33_colmajor(const double* A, const double* v, double* out) {
    out[0] = A[0] * v[0] + A[3] * v[1] + A[6] * v[2];
    out[1] = A[1] * v[0] + A[4] * v[1] + A[7] * v[2];
    out[2] = A[2] * v[0] + A[5] * v[1] + A[8] * v[2];
}

// Slater Orbital Data Structure (already POD-like)
struct SlaterOrbitalData {
    int angMom;
    double cutoff;
    int nPow; 
    int nAlpha;
    const double* aa; 
    const double* alpha;
};

double getRadialDirect(double rr, const SlaterOrbitalData& stoData) {
    assert(stoData.nPow <= MAX_NPOW_CONST && "nPow exceeds MAX_NPOW_CONST");
    if (stoData.nPow == 0) return 0.0; 

    if (rr < std::numeric_limits<double>::epsilon()) {
        if (stoData.angMom > 0) {
            return 0.0; 
        }
        rr = 0.0; 
    }

    double pows[MAX_NPOW_CONST]; 

    if (rr == 0.0) { 
        assert(stoData.angMom == 0);
        pows[0] = 1.0; 
        for (int i = 1; i < stoData.nPow; ++i) {
            pows[i] = 0.0; 
        }
    } else { 
        double r_val = rr; 
        double r_to_angMom;

        switch (stoData.angMom) {
            case 0: r_to_angMom = 1.0; break;
            case 1: r_to_angMom = r_val; break;
            case 2: r_to_angMom = r_val * r_val; break;
            case 3: r_to_angMom = r_val * r_val * r_val; break; 
            default: r_to_angMom = std::pow(r_val, stoData.angMom); 
        }
        
        pows[0] = r_to_angMom;
        for (int i = 1; i < stoData.nPow; ++i) {
            pows[i] = pows[i-1] * r_val;
        }
    }

    double sto_val = 0.0;
    const double* current_aa_coeffs_start = stoData.aa;
    for (int i = 0; i < stoData.nAlpha; ++i) { 
        const double* aa_for_this_alpha = current_aa_coeffs_start + (long long)i * stoData.nPow;
        double poly_sum = 0.0;
        for (int j = 0; j < stoData.nPow; ++j) { 
            poly_sum += aa_for_this_alpha[j] * pows[j];
        }
        sto_val += poly_sum * std::exp(stoData.alpha[i] * rr); 
    }
    return sto_val;
}

// Constants for Real Tesseral Spherical Harmonics Y_lm (l<=3)
const double C00 = 0.2820947917738782;  // 1/sqrt(4pi)
const double C1  = 0.4886025119029198;  // sqrt(3/4pi)
const double C20 = 0.3153915652525200;  // sqrt(5/16pi) * (3z2-r2)/r2 form
const double C21 = 1.092548430592079;  // sqrt(15/4pi) for dxz, dyz, dxy (Nobel/ACES convention for dxy,dyz,dxz)
const double C22 = 0.5462742152960395;  // sqrt(15/16pi) for dx2-y2
const double C30 = 0.3731763325901155;  // sqrt(7/16pi)
const double C31 = 0.4570457994644658;  // sqrt(21/32pi)*sqrt(2) (Nobel/ACES convention for fxyz-like terms)
const double C32 = 1.445305721320277;   // sqrt(105/16pi)
const double C33 = 0.5900435899266435;  // sqrt(35/32pi)*sqrt(2) (Nobel/ACES convention for fx3-like terms)

double realTessY(int ll, int mm, double x, double y, double z, double rr) {
    assert(ll >= 0 && ll <= 3); 
    assert(std::abs(mm) <= ll);

    if (ll == 0) {
        return C00;
    }
    if (rr < std::numeric_limits<double>::epsilon()) { // For ll > 0
        return 0.0;
    }

    double x_r = x / rr;
    double y_r = y / rr;
    double z_r = z / rr;

    // Original code used divisions by rr, rr2, rr3. Switched to x_r, y_r, z_r for clarity.
    // Ensure expressions match. Example: (2zz-xx-yy)/rr2 = (3zz-rr2)/rr2 = 3*z_r*z_r - 1
    // The C20 constant is for (2zz-xx-yy)/rr2 form.
    const double rr2 = rr*rr; // Retain for forms that are simpler with it.

    switch (ll) {
        case 1:
            switch (mm) {
                case -1: return C1 * y_r; 
                case  0: return C1 * z_r; 
                case  1: return C1 * x_r; 
            } break;
        case 2:
            switch (mm) {
                case -2: return C21 * x_r * y_r;                       
                case -1: return C21 * y_r * z_r;                       
                case  0: return C20 * (2.0 * z*z - x*x - y*y) / rr2; 
                case  1: return C21 * x_r * z_r;                       
                case  2: return C22 * (x*x - y*y) / rr2;         
            } break;
        case 3:
             // Using x_r, y_r, z_r for f-orbitals to match common normalized forms
            switch (mm) {
                case -3: return C33 * y_r * (3.0 * x_r*x_r - y_r*y_r);
                case -2: return C32 * x_r * y_r * z_r;
                case -1: return C31 * y_r * (4.0 * z_r*z_r - x_r*x_r - y_r*y_r); // y(5z^2-r^2)/r^3 form
                case  0: return C30 * z_r * (2.0 * z_r*z_r - 3.0 * x_r*x_r - 3.0 * y_r*y_r); // z(5z^2-3r^2)/r^3 form
                case  1: return C31 * x_r * (4.0 * z_r*z_r - x_r*x_r - y_r*y_r); // x(5z^2-r^2)/r^3 form
                case  2: return C32 * z_r * (x_r*x_r - y_r*y_r);
                case  3: return C33 * x_r * (x_r*x_r - 3.0 * y_r*y_r);
            } break;
    }
    assert(false && "Fell through realTessY switch statements"); 
    return 0.0;
}

} // end anonymous namespace


// Main Evaluation Function (C Interface)
extern "C" {

void evaluatePointwise_cpp_real(
// Grid definition
const double* origin,         
const double* gridVecs,       
int nGridX, int nGridY, int nGridZ, 

// Eigenvectors
const double* eigVecsReal,    
int nEig,                     

// System Geometry & Basis
int nAtom,
int nOrb,
const double* coords,         
const int* species,           
const int* iStos,             
int nSpecies,                 

// STO Data
const int* sto_angMom,        
const double* sto_cutoff,     
const int* sto_nPow,          
const int* sto_nAlpha,        
const double* sto_aa_flat,    
const double* sto_alpha_flat, 
const int* sto_aa_offsets,    
const int* sto_alpha_offsets, 
int nTotalStos,               

// Periodicity & Calculation Type
bool tPeriodic,
const double* latVecs,        
const double* recVecs2p,      
int nCell,                    
const double* cellVec,        // Unused if coords contains all cell images explicitly
bool tAddDensities,           

// Output Arrays
double* valueReal            
) {

// --- Temporary Storage ---
// atomAllOrbVal removed
std::vector<double> atomOrbValReal(nOrb); 
std::vector<int> nonZeroIndices;
nonZeroIndices.reserve(nOrb);


// --- Grid Point Iteration ---
long long nGridPointsTotal = (long long)nGridX * nGridY * nGridZ; 
for (int i3 = 0; i3 < nGridZ; ++i3) {
    double p3[3] = { gridVecs[0 + 3*2] * i3, gridVecs[1 + 3*2] * i3, gridVecs[2 + 3*2] * i3 };
    for (int i2 = 0; i2 < nGridY; ++i2) {
        double p2[3] = { gridVecs[0 + 3*1] * i2, gridVecs[1 + 3*1] * i2, gridVecs[2 + 3*1] * i2 };
        for (int i1 = 0; i1 < nGridX; ++i1) {
            double p1[3] = { gridVecs[0 + 3*0] * i1, gridVecs[1 + 3*0] * i1, gridVecs[2 + 3*0] * i1 };

            double xyz[3] = {
                p1[0] + p2[0] + p3[0] + origin[0],
                p1[1] + p2[1] + p3[1] + origin[1],
                p1[2] + p2[2] + p3[2] + origin[2]
            };

            if (tPeriodic) {
                double xyz_orig[3] = {xyz[0], xyz[1], xyz[2]};
                double frac[3];
                matvec_mult33_colmajor(recVecs2p, xyz_orig, frac);
                frac[0] -= std::floor(frac[0]);
                frac[1] -= std::floor(frac[1]);
                frac[2] -= std::floor(frac[2]);
            }

            // --- Calculate AO contributions for this point ---
            std::fill(atomOrbValReal.begin(), atomOrbValReal.end(), 0.0);

            int currentOrbIndex = 0;
            for (int iAtom = 0; iAtom < nAtom; ++iAtom) {
                int sp = species[iAtom] - 1; 
                assert(sp >= 0 && sp < nSpecies);
                int stoStartIdx = iStos[sp] - 1;
                int stoEndIdx = iStos[sp + 1] - 1;

                for (int iSto = stoStartIdx; iSto < stoEndIdx; ++iSto) {
                    SlaterOrbitalData currentStoData;
                    currentStoData.angMom = sto_angMom[iSto];
                    currentStoData.cutoff = sto_cutoff[iSto];
                    currentStoData.nPow = sto_nPow[iSto];
                    currentStoData.nAlpha = sto_nAlpha[iSto];
                    currentStoData.aa = sto_aa_flat + sto_aa_offsets[iSto];
                    currentStoData.alpha = sto_alpha_flat + sto_alpha_offsets[iSto];
                    int ll = currentStoData.angMom;

                    if (currentStoData.nPow == 0 && ll > 0) { // Optimization: if nPow=0 for l>0, STO is 0.
                         for (int mm_dummy = -ll; mm_dummy <= ll; ++mm_dummy) { // Still need to advance currentOrbIndex
                            atomOrbValReal[currentOrbIndex++] = 0.0; // Explicitly set to 0.0
                        }
                        continue; // Skip to next STO
                    }


                    for (int mm = -ll; mm <= ll; ++mm) {
                        double orbitalSumOverCells = 0.0;
                        for (int iCell = 0; iCell < nCell; ++iCell) {
                            const double* atomCoord = coords + 3 * (iAtom + (long long)nAtom * iCell);
                            double diff[3] = { xyz[0] - atomCoord[0],
                                               xyz[1] - atomCoord[1],
                                               xyz[2] - atomCoord[2] };
                            double rr = norm2(diff[0], diff[1], diff[2]);
                            
                            if (rr <= currentStoData.cutoff) {
                                if (currentStoData.nPow == 0) { // Should only be reachable if ll=0 for nPow=0 to be non-zero
                                    // If angMom=0, nPow=0, getRadialDirect returns 0.
                                    // This case effectively means this STO contributes 0.
                                    // No need to add to orbitalSumOverCells.
                                } else {
                                    double radialVal = getRadialDirect(rr, currentStoData);
                                    double angularVal = realTessY(ll, mm, diff[0], diff[1], diff[2], rr);
                                    orbitalSumOverCells += radialVal * angularVal;
                                }
                            }
                        } // end cell loop

                        if (tAddDensities) {
                            atomOrbValReal[currentOrbIndex] = orbitalSumOverCells * orbitalSumOverCells;
                        } else {
                            atomOrbValReal[currentOrbIndex] = orbitalSumOverCells;
                        }
                        currentOrbIndex++;
                    } // end mm loop
                } // end iSto loop
            } // end iAtom loop
            assert(currentOrbIndex == nOrb);

            // --- Sum contributions and calculate MO values ---
            long long currentGridPointOffset = (long long)i1 + (long long)nGridX * (i2 + (long long)nGridY * i3); 

            bool allCalculatedValuesZero = true;
            for(int iOrb = 0; iOrb < nOrb; ++iOrb) {
                if (std::abs(atomOrbValReal[iOrb]) > 1e-15) { 
                    allCalculatedValuesZero = false;
                    break;
                }
            }

            if (allCalculatedValuesZero) {
               for(int iE = 0; iE < nEig; ++iE) {
                   valueReal[currentGridPointOffset + nGridPointsTotal * iE] = 0.0;
               }
               continue; 
            }

            nonZeroIndices.clear();
            for(int iOrb = 0; iOrb < nOrb; ++iOrb) {
                if (std::abs(atomOrbValReal[iOrb]) > 1e-15) { 
                    nonZeroIndices.push_back(iOrb);
                }
            }
            int nNonZero = nonZeroIndices.size();
            assert(nNonZero > 0);


            // Calculate final MO value (dot product with eigenvectors)
            for (int iE = 0; iE < nEig; ++iE) {
                double finalValue = 0.0;
                const double* eigVec = eigVecsReal + (long long)nOrb * iE; 
                for (int idx = 0; idx < nNonZero; ++idx) {
                     int iOrb = nonZeroIndices[idx];
                     finalValue += atomOrbValReal[iOrb] * eigVec[iOrb];
                }
                valueReal[currentGridPointOffset + nGridPointsTotal * iE] = finalValue;
            }
        } // end i1 loop
    } // end i2 loop
} // end i3 loop

} // end extern "C" evaluatePointwise_cpp_real

} // end extern "C"
