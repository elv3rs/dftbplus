#include "utils.cuh"
// Maximum number of powers in STOs, for static array sizing
constexpr int STO_MAX_POWS = 16;

/**
 * @brief Namespace for real spherical harmonics constants up to l=3.
 *
 *  C_{l,m} = sqrt((2l+1)/(4*PI) * (l-m)!/(l+m)!)
 */
namespace RealTessYConsts {
static const double C00 = 0.2820947917738782;  // 1 / sqrt(4*PI)
static const double C11 = 0.4886025119029198;  // sqrt(3 / (4*PI))
static const double C20 = 0.3153915652525200;  // sqrt(5 / (16*PI))
static const double C21 = 1.092548430592079;   // sqrt(15 / (4*PI)) / 2
static const double C22 = 0.5462742152960395;  // sqrt(15 / (16*PI))
static const double C30 = 0.3731763325901155;  // sqrt(7 / (16*PI))
static const double C31 = 0.4570457994644658;  // sqrt(42 / (64*PI))
static const double C32 = 1.445305721320277;   // sqrt(105 / (16*PI)) / 2
static const double C33 = 0.5900435899266435;  // sqrt(70 / (64*PI))
}  // namespace RealTessYConsts

// =========================================================================
//  CUDA Device Functions
// =========================================================================

/**
 * @brief Computes real tesseral spherical harmonics Y_lm(r) up to l=3.
 * @param ll      Orbital quantum number
 * @param mm      Magnetic quantum number
 * @param diff    Pointer to (x, y, z) vector
 * @param inv_r   Pre-calculated 1/r
 * @param inv_r2  Pre-calculated 1/r^2
 * @return The value of the real spherical harmonic.
 */
__device__ __forceinline__ double realTessY(int ll, int mm, const double* diff, double inv_r, double inv_r2) {
    const double x = diff[0];
    const double y = diff[1];
    const double z = diff[2];

    const double x_r = x * inv_r;
    const double y_r = y * inv_r;
    const double z_r = z * inv_r;

    using namespace RealTessYConsts;

    switch (ll) {
        case 0: return C00;
        case 1:
            switch (mm) {
                case -1: return C11 * y_r;
                case 0: return C11 * z_r;
                case 1: return C11 * x_r;
            }
            break;
        case 2: {
            const double xx_r2 = x_r * x_r;
            const double yy_r2 = y_r * y_r;
            const double zz_r2 = z_r * z_r;
            const double xy_r2 = x_r * y_r;
            const double yz_r2 = y_r * z_r;
            const double xz_r2 = x_r * x_r;
            switch (mm) {
                case -2: return C21 * xy_r2;
                case -1: return C21 * yz_r2;
                case 0: return C20 * (3.0 * zz_r2 - 1.0);
                case 1: return C21 * xz_r2;
                case 2: return C22 * (xx_r2 - yy_r2);
            }
        } break;
        case 3: {
            const double z_r2 = z_r * z_r;
            switch (mm) {
                case -3: return C33 * y_r * (3.0 * x_r * x_r - y_r * y_r);
                case -2: return C32 * x_r * y_r * z_r;
                case -1: return C31 * y_r * (5.0 * z_r2 - 1.0);
                case 0: return C30 * z_r * (5.0 * z_r2 - 3.0);
                case 1: return C31 * x_r * (5.0 * z_r2 - 1.0);
                case 2: return C32 * z_r * (x_r * x_r - y_r * y_r);
                case 3: return C33 * x_r * (x_r * x_r - 3.0 * y_r * y_r);
            }
        } break;
    }
    return 0.0;
}

/**
 * @brief Computes the radial part of a Slater-type orbital (STO).
 * @param r         Distance from center
 * @param iL        Angular momentum
 * @param iOrb      Index of the orbital
 * @param nPows     Number of polynomial powers in the STO
 * @param nAlphas   Number of exponential terms in the STO
 * @param coeffs    Coefficients for the polynomial terms
 * @param alphas    Exponential decay constants
 * @param maxNPows  Maximum number of polynomial powers across all orbitals
 * @param maxNAlphas Maximum number of alphas across all orbitals
 * @return The radial value of the STO at distance r.
 */
__device__ __forceinline__ double getRadialValue(double r, int iL, int iOrb, int nPows, int nAlphas,
    const double* coeffs, const double* alphas, int maxNPows, int maxNAlphas) {
    double sto_tmp_pows[STO_MAX_POWS];

    double sto_tmp_rexp = 1.0;

    if (iL > 0 || r > 1.0e-12) {
        for (int p = 0; p < iL; ++p) {
            sto_tmp_rexp *= r;
        }
    }
    for (int ii = 0; ii < nPows; ++ii) {
        sto_tmp_pows[ii] = sto_tmp_rexp;
        sto_tmp_rexp *= r;
    }

    double radialVal = 0.0;
    for (int ii = 0; ii < nAlphas; ++ii) {
        double term = 0.0;
        for (int jj = 0; jj < nPows; ++jj) {
            term += coeffs[IDX3F(jj, ii, iOrb, maxNPows, maxNAlphas)] * sto_tmp_pows[jj];
        }
        radialVal += term * exp(alphas[IDX2F(ii, iOrb, maxNAlphas)] * r);
    }

    return radialVal;
}

__device__ __forceinline__ void matmul3x3_vec(const double mat[3][3], const double vec[3], double result[3]) {
    for (int i = 0; i < 3; i++) {
        result[i] = mat[i][0] * vec[0] + mat[i][1] * vec[1] + mat[i][2] * vec[2];
    }
}

__device__ __forceinline__ void foldCoordsIntoCell(
    double xyz[3], const double latVecs[3][3], const double recVecs2p[3][3]) {
    double frac[3];
    matmul3x3_vec(recVecs2p, xyz, frac);

    for (int i = 0; i < 3; i++) {
        frac[i] -= floor(frac[i]);
    }

    matmul3x3_vec(latVecs, frac, xyz);
}
