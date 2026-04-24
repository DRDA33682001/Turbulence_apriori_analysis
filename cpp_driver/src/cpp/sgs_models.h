/*
 * sgs_models.h — Subgrid-scale closure models for a priori LES evaluation.
 *
 * Three eddy-viscosity / dynamic closures live here. Each consumes filtered
 * fields produced by sgs_exact.* and emits a MODELLED deviatoric SGS stress
 * that is later compared against the exact stress from compute_exact_sgs_stress.
 *
 * Physics references
 *   [1] Smagorinsky-Lilly
 *         Pope, Turbulent Flows (2000), Eq. 13.127-13.128, Eq. 13.135
 *         Lilly (1967) inertial-range analysis  -> Cs = 0.17
 *   [2] WALE
 *         Nicoud & Ducros, Flow Turb. Combust. 62 (1999), Eq. 5, Eq. 14
 *   [3] Dynamic Smagorinsky (Germano-Lilly)
 *         Germano et al., Phys. Fluids A 3(7), 1991  — Eq. 4 (Leonard), Eq. 8 (identity)
 *         Lilly,          Phys. Fluids A 4(3), 1992  — Eq. 8 (M_ij), Eq. 10 (least squares)
 *         Pope (2000), p. 606-607 and Eq. 13.252
 *
 * All tau outputs are the DEVIATORIC part of the modelled SGS stress
 * (the eddy-viscosity closure family is traceless by construction).
 * Symmetric tensor storage (6 components): 11, 12, 13, 22, 23, 33.
 */

#pragma once

#include "io.h"

#include <vector>

// ---------------------------------------------------------------------------
// MODEL 1: Smagorinsky-Lilly
// ---------------------------------------------------------------------------
//   nu_r       = (Cs * delta)^2 * |S|              (Pope Eq. 13.128)
//   tau_ij_dev = -2 * nu_r * S_ij                  (Pope Eq. 13.127)
//
// Purely dissipative (Pi = -tau:S >= 0 everywhere). This is asserted inside
// the implementation; a violation aborts with the reference printed.
//
void smagorinsky_model(const TensorField6&        S,
                       const std::vector<double>& S_mag,
                       double                     delta,
                       int                        N,
                       TensorField6&              tau_smag,
                       std::vector<double>&       nu_r);

// ---------------------------------------------------------------------------
// MODEL 2: WALE (Nicoud & Ducros 1999)
// ---------------------------------------------------------------------------
//   g2_ij      = g_ik * g_kj                        [matrix square of grad u]
//   Sd_ij      = 0.5*(g2_ij + g2_ji) - (1/3)*delta_ij*g2_kk   (ND99 Eq. 5)
//   nu_r       = (Cw*delta)^2 * (Sd:Sd)^(3/2)
//                / ((S:S)^(5/2) + (Sd:Sd)^(5/4) + eps)         (ND99 Eq. 14)
//   tau_ij_dev = -2 * nu_r * S_ij
//
// Purely dissipative like Smagorinsky, but nu_r vanishes in pure-rotation
// regions (where Sd_ij = 0), which is WALE's key advantage for wall-bounded
// and transitional flows.
//
void wale_model(const std::vector<std::vector<double>>& g,
                const TensorField6&                     S,
                double                                  delta,
                int                                     N,
                TensorField6&                           tau_wale,
                std::vector<double>&                    nu_r_wale);

// ---------------------------------------------------------------------------
// MODEL 3: Dynamic Smagorinsky (Germano-Lilly)
// ---------------------------------------------------------------------------
//   hat_delta = 2 * delta
//   L_ij      = filter_hat(ubar_i * ubar_j) - uhatbar_i * uhatbar_j   (Germano Eq. 4)
//   Ld_ij     = L_ij - (1/3)*delta_ij*L_kk
//   M_ij      = 2*hat_delta^2*|S_hat|*S_hat_ij
//               - filter_hat( 2*delta^2*|S|*S_ij )                    (Lilly Eq. 8)
//   Cs^2      = (Ld:M) / (M:M + eps)                                  (Lilly Eq. 10)
//   tau_ij    = -2 * Cs^2 * delta^2 * |S| * S_ij
//
// Homogeneous-direction averaging (`avg_mode`):
//   The pointwise Cs^2 field from the Germano-Lilly formula is very noisy
//   when the test filter is a sharp spectral filter (Meneveau & Katz,
//   ARFM 2000, §3.2). Lilly (1992) recommends averaging the numerator
//   (Ld:M) and denominator (M:M) SEPARATELY over the homogeneous
//   directions before dividing:
//
//     Cs^2_avg = <Ld:M>_homog / <M:M>_homog
//
//   * Box       : average over all three directions (fully homogeneous
//                 flow, e.g. isotropic turbulence) -> scalar Cs^2(t).
//   * XZ_Plane  : average over (x,z) at each y (channel flow; statistical
//                 homogeneity in streamwise and spanwise only)
//                 -> Cs^2(y,t) profile broadcast as a 3-D field.
//   * Pointwise : no averaging (backward-compatible, not recommended).
//
// The pointwise Cs^2 is also returned (`Cs2_pointwise`) so the raw
// distribution remains available for the Cs^2 PDF diagnostic plot.
// `tau_dyn` is ALWAYS built from `Cs2_effective` (the averaged field).
//
enum class AveragingMode { Pointwise, Box, XZ_Plane };

void dynamic_smagorinsky(const VelocityField&       vel_bar,
                         const TensorField6&        S,
                         const std::vector<double>& S_mag,
                         double                     delta,
                         double                     dx,
                         int                        N,
                         TensorField6&              tau_dyn,
                         std::vector<double>&       Cs2_pointwise,
                         std::vector<double>&       Cs2_effective,
                         AveragingMode              avg_mode = AveragingMode::Pointwise);
