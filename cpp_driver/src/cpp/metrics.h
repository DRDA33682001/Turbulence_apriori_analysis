/*
 * metrics.h — A priori evaluation metrics for SGS model assessment.
 *
 * All functions compare a modelled SGS stress tensor against the exact
 * (DNS-derived) stress produced by compute_exact_sgs_stress in sgs_exact.h.
 * No modelling assumptions live here — these are diagnostic quantities only.
 *
 * Physics references
 *   [Clark79]  Clark, Ferziger & Reynolds, J. Fluid Mech. 91 (1979)
 *                — defines tensor correlation coefficient, backscatter PDF
 *   [McM79]    McMillan & Ferziger, AIAA J. 17(12) (1979)
 *                — direct testing methodology, alignment angle
 *   [Pope]     Pope, Turbulent Flows (2000)
 *                — Pi definition Eq. 13.123, strain-rate Eq. 13.73–13.74
 *   [Liu94]    Liu, Meneveau & Katz, J. Fluid Mech. 275 (1994)
 *                — alignment PDF measurements in physical experiments
 *
 * Tensor storage convention (6 components in order): 11, 12, 13, 22, 23, 33
 * This matches TensorField6::components and all callers in this project.
 *
 * Index mapping: for component index c in [0,5]
 *   c=0 -> (i=0,j=0)=11,  c=1 -> (i=0,j=1)=12,  c=2 -> (i=0,j=2)=13
 *   c=3 -> (i=1,j=1)=22,  c=4 -> (i=1,j=2)=23,  c=5 -> (i=2,j=2)=33
 * Off-diagonal components (12, 13, 23) appear once in the 6-component list
 * but represent TWO entries of the full 3x3 tensor. The double-contraction
 * A_ij * B_ij for symmetric tensors therefore uses weight 2 for off-diagonal
 * components. This is made explicit in every function below.
 */

#pragma once

#include "io.h"

#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// METRIC 1: Tensor correlation coefficient
// ---------------------------------------------------------------------------
//
// Physics
//   Pearson correlation between the model stress tensor tau_model and the
//   exact SGS stress tau_exact. All 6 unique components at all N^3 grid
//   points are stacked into a single flat vector for each field; the scalar
//   Pearson rho is computed over this combined sample.
//
//   This is the standard a priori metric introduced in Clark et al. (1979)
//   and is the most-cited single number for evaluating eddy-viscosity models.
//   It measures POINTWISE alignment in the tensor sense — a high correlation
//   means the model not only captures the right sign and magnitude but also
//   the correct spatial organisation of the stress field.
//
// Formula
//   rho = sum_p (A_p - <A>) * (B_p - <B>)
//         / sqrt( sum_p(A_p-<A>)^2 * sum_p(B_p-<B>)^2 )
//   where A = flattened tau_model, B = flattened tau_exact (both length 6*N^3)
//
// Expected values (from a priori studies of isotropic turbulence)
//   Smagorinsky (Cs=0.17)    : rho ~ 0.2 – 0.4   [Clark et al. 1979, Table 1]
//   WALE                     : rho ~ 0.2 – 0.4   [similar, eddy-viscosity class]
//   Dynamic Smagorinsky      : rho ~ 0.4 – 0.6   [improved by local Cs adjustment]
//   Perfect model            : rho = 1.0
//
// Warning thresholds
//   rho > 0.6 for Smagorinsky : suspect — eddy-viscosity models cannot achieve
//     this in isotropic turbulence; check for sign errors or double-counting
//   rho < 0.1 for any model   : implementation error likely
//
double correlation_coefficient(const TensorField6& tau_model,
                               const TensorField6& tau_exact,
                               int                 N);

// ---------------------------------------------------------------------------
// METRIC 2: SGS energy dissipation field
// ---------------------------------------------------------------------------
//
// Physics
//   Pi = -tau_ij * S_ij    (Pope Eq. 13.123, double contraction over i,j)
//
//   Pi > 0 : forward scatter — energy transferred from resolved to sub-grid
//   Pi < 0 : backscatter     — energy transferred from sub-grid to resolved
//
//   In 3D isotropic turbulence the MEAN Pi is always positive (net forward
//   cascade) but a significant FRACTION of grid points exhibit Pi < 0
//   (local backscatter), which is physically well-established.
//
// Tracelessness note
//   For the full tau (not tau_dev), the contraction with S gives the same
//   result as contracting tau_dev with S, because S is traceless
//   (incompressibility: S_kk = du_k/dx_k = 0 for each k, and summing gives
//   S_11+S_22+S_33 = div(u_bar) = 0). Therefore:
//       tau_ij * S_ij = tau_dev_ij * S_ij  (S_kk = 0 makes the trace term vanish)
//   Callers may pass either tau or tau_dev; the result is identical. The
//   function documentation uses tau for generality.
//
// Double-contraction storage weight
//   Because tau and S are symmetric and stored as 6 unique components, the
//   full double sum tau_ij * S_ij must double-count the off-diagonal entries:
//       tau:S = tau_11*S_11 + tau_22*S_22 + tau_33*S_33
//             + 2*(tau_12*S_12 + tau_13*S_13 + tau_23*S_23)
//
// Output
//   Pi : flat N^3 vector, Pi[p] = -( tau:S )[p] at each grid point
//
void compute_dissipation_field(const TensorField6&  tau,
                               const TensorField6&  S,
                               int                  N,
                               std::vector<double>& Pi);

// ---------------------------------------------------------------------------
// METRIC 3: Backscatter fraction
// ---------------------------------------------------------------------------
//
// Physics
//   Returns the fraction of grid points where Pi < 0 (local backscatter).
//   This is a fundamental test of whether a model can represent the
//   intermittent reverse energy transfer observed in DNS data.
//
// Expected values
//   DNS exact stress          : 0.30 – 0.40   [Clark et al. 1979, Fig. 4]
//   Smagorinsky               : 0.000          [purely dissipative by construction;
//                                               assert failure triggers if > 1e-12]
//   WALE                      : 0.000          [purely dissipative by construction;
//                                               assert failure triggers if > 1e-12]
//   Dynamic Smagorinsky       : 0.15 – 0.30   [partial recovery via negative Cs^2]
//
// Assertion
//   The function takes an optional model_name argument. If model_name is
//   "smagorinsky" or "wale" (case-insensitive), it asserts that the returned
//   fraction is < 1e-12 and aborts with an explanatory message if violated.
//
double backscatter_fraction(const std::vector<double>& Pi,
                            int                        N,
                            const std::string&         model_name = "");

// ---------------------------------------------------------------------------
// METRIC 4: Alignment angle between tau_dev and S
// ---------------------------------------------------------------------------
//
// Physics
//   cos_theta_p = (tau_dev_ij * S_ij)_p / ( |tau_dev|_p * |S|_p )
//
//   where the tensor norms use the double-contraction inner product:
//     |A|^2 = A_ij * A_ij  (with weight 2 for off-diagonals as above)
//
//   Under the eddy-viscosity hypothesis tau_dev_ij = -2 * nu_r * S_ij,
//   so tau_dev is exactly anti-parallel to S, giving cos_theta = -1 everywhere.
//   A broad distribution of cos_theta is therefore evidence that the
//   eddy-viscosity hypothesis is poorly satisfied locally.
//
// Expected values (from a priori and experimental studies)
//   Perfect eddy-viscosity alignment : cos_theta = -1.000 everywhere
//   DNS exact vs S                   : mean(cos_theta) ~ -0.3 to -0.5
//                                      (broad distribution, NOT peaked at -1)
//                                      [Liu, Meneveau & Katz 1994, JFM 275]
//   Smagorinsky / WALE               : cos_theta = -1 EXACTLY
//                                      (they are proportional to S by construction)
//   Dynamic Smagorinsky              : cos_theta = -1 EXACTLY
//                                      (same structural form, only Cs^2 varies)
//
// Regularization
//   If |tau_dev|_p < 1e-15 or |S|_p < 1e-15 at a grid point, the angle is
//   undefined. cos_theta[p] is set to 0 (uncorrelated convention) and that
//   point is excluded from any downstream statistics.
//
// Output
//   cos_theta : flat N^3 vector in [-1, 1] (or 0 for regularised points)
//
void compute_alignment_angle(const TensorField6&  tau_dev,
                             const TensorField6&  S,
                             int                  N,
                             std::vector<double>& cos_theta);

// ---------------------------------------------------------------------------
// METRIC 5: Aggregate driver — compute all metrics, print summary, write HDF5
// ---------------------------------------------------------------------------
//
// Calls correlation_coefficient, compute_dissipation_field, backscatter_fraction,
// and compute_alignment_angle in sequence, then writes all scalar results to
// HDF5 via write_metrics_hdf5 and prints a human-readable physics summary.
//
// Parameters
//   tau_exact_dev    : exact deviatoric SGS stress (from sgs_exact.h)
//   tau_model        : modelled deviatoric SGS stress (from sgs_models.h)
//   S                : filtered strain-rate tensor (from sgs_exact.h)
//   nu_r             : eddy-viscosity field from model (N^3 scalar)
//   N                : grid size per axis
//   delta            : filter width used for this pass (for printed header)
//   model_name       : label used in HDF5 group name and printed header
//   results_filepath : HDF5 file path for write_metrics_hdf5
//
// Printed output format
//   === [model_name] at delta=[delta] ===
//   Correlation coefficient rho  : X.XXX  (literature: ~0.3 Smag)
//   Backscatter fraction         : XX.X%  (DNS truth: ~35%)
//   Mean eddy viscosity          : X.XXe-XX
//   Mean SGS dissipation         : X.XXe-XX
//   Mean cos_theta               : X.XXX  (perfect EV = -1.000)
//
//   Physics assessments:
//   [WARNING / NOTE / OK lines based on thresholds below]
//
// Physics assessment thresholds
//   rho < 0.1                              -> WARNING: correlation very low, check implementation
//   backscatter > 0 for smag/wale          -> ERROR: (asserted in backscatter_fraction)
//   mean cos_theta > -0.1                  -> NOTE: eddy-viscosity hypothesis poorly satisfied
//
void compute_all_metrics(const TensorField6&        tau_exact_dev,
                         const TensorField6&        tau_model,
                         const TensorField6&        S,
                         const std::vector<double>& nu_r,
                         int                        N,
                         double                     delta,
                         const std::string&         model_name,
                         const std::string&         results_filepath);
