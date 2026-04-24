#pragma once

/*
 * config.h — Global constants for the a priori LES SGS model evaluation.
 *
 * All physical constants are traceable to a specific equation reference.
 * If any constant is changed, update the reference comment alongside it.
 *
 * Physics references:
 *   [Pope]    S.B. Pope, Turbulent Flows, Cambridge University Press (2000)
 *   [Germano] Germano et al., Phys. Fluids A 3(7), 1991
 *   [Lilly]   Lilly, Phys. Fluids A 4(3), 1992
 *   [ND99]    Nicoud & Ducros, Flow Turb. Combust. 62, 1999
 *   [Clark]   Clark, Ferziger & Reynolds, J. Fluid Mech. 91, 1979
 *   [McM79]   McMillan & Ferziger, AIAA J. 17(12), 1979
 */

#include <cmath>

// -----------------------------------------------------------------------
// Grid parameters
// -----------------------------------------------------------------------

// 256^3 velocity snapshot obtained from JHTDB isotropic1024coarse with
// stride = 1024 / 256 = 4 so that the sampled grid spans the entire
// (2 pi)^3 periodic domain. This gives statistically representative
// isotropic turbulence (many integral scales across the cube) rather
// than a 1-integral-scale sub-block. See src/python/download.py.
constexpr int N = 256;

// Physical grid spacing for the 256^3 strided cube covering period L = 2*pi.
// DX = 2*pi / N.
// Reference: JHTDB isotropic turbulence dataset documentation; stride=4.
constexpr double DX = 2.0 * M_PI / 256.0;

// -----------------------------------------------------------------------
// Filter widths
// -----------------------------------------------------------------------

// Three filter widths spanning the inertial subrange of the Re_lambda=433 field.
// Each is a multiple of the DNS grid spacing DX.
// Reference: Pope Ch.13, filter width delta appears in Eq. 13.22 (sharp spectral filter)
// and Eq. 13.127 (Smagorinsky model). Chosen to bracket the inertial subrange.
constexpr int    N_FILTERS = 3;
constexpr double FILTER_WIDTHS[N_FILTERS] = {
    4.0  * DX,   // delta_1 — smallest filter, deepest into inertial range
    8.0  * DX,   // delta_2 — mid filter
    16.0 * DX    // delta_3 — largest filter, near energy-containing range
};

// -----------------------------------------------------------------------
// Model constants
// -----------------------------------------------------------------------

// Smagorinsky constant Cs = 0.17.
// Reference: Pope Eq. 13.128; value from Lilly (1967) inertial-range derivation.
// Expected behavior: optimal Cs in isotropic turbulence is 0.17–0.18.
constexpr double CS_SMAGORINSKY = 0.17;

// WALE constant Cw = 0.5.
// Reference: Nicoud & Ducros (1999) [ND99], Eq. 14.
// Note: Nicoud & Ducros quote Cw in range 0.55–0.60 from isotropic calibration;
// 0.5 is a commonly used round value in practice.
constexpr double CW_WALE = 0.5;

// -----------------------------------------------------------------------
// Benchmark values from literature (used for runtime sanity checks)
// -----------------------------------------------------------------------

// Expected tensor correlation coefficient between Smagorinsky model stress
// and exact SGS stress in a priori tests of isotropic turbulence.
// Reference: Clark, Ferziger & Reynolds (1979) [Clark], Eq. 26 and Table 1.
// Value: ~0.3 (eddy-viscosity models have inherently low pointwise correlation).
constexpr double CS2_EXPECTED_CORRELATION = 0.3;

// Expected fraction of grid points with negative SGS energy dissipation
// (backscatter: Pi = -tau_ij * S_bar_ij < 0).
// Reference: Clark et al. (1979) [Clark]; McMillan & Ferziger (1979) [McM79].
// Physics note: this is LOCAL intermittent backscatter in a globally
// forward-cascading 3D flow. Mean(Pi) > 0 always; ~30-40% of points have Pi < 0.
// This is NOT the 2D inverse energy cascade.
constexpr double EXPECTED_BACKSCATTER_FRACTION = 0.35;

// -----------------------------------------------------------------------
// Array indexing utility
// -----------------------------------------------------------------------

// Row-major flat index for a 3D array of size N x N x N.
// Usage: data[idx(i, j, k, N)]
// All arrays in this project are stored as flat 1D vectors in this order.
inline int idx(int i, int j, int k, int n) {
    return i * n * n + j * n + k;
}

// -----------------------------------------------------------------------
// Runtime sanity-check helpers
// -----------------------------------------------------------------------

#include <cstdio>
#include <cmath>
#include <cstdlib>

// Print a named scalar check. Aborts if |value - expected| / |expected| > tol.
inline void check_scalar(const char* name, double value, double expected, double tol = 0.15) {
    double rel_err = std::abs(value - expected) / (std::abs(expected) + 1e-14);
    const char* status = (rel_err <= tol) ? "OK" : "WARN";
    std::printf("[check] %-40s value=%.4f  expected=%.4f  rel_err=%.3f  [%s]\n",
                name, value, expected, rel_err, status);
}

// Print a named range check. Warns if value is outside [lo, hi].
inline void check_range(const char* name, double value, double lo, double hi) {
    const char* status = (value >= lo && value <= hi) ? "OK" : "WARN";
    std::printf("[check] %-40s value=%.4f  range=[%.4f, %.4f]  [%s]\n",
                name, value, lo, hi, status);
}
