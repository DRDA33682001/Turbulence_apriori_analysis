/*
 * metrics.cpp — Implementation of a priori SGS model evaluation metrics.
 *
 * See metrics.h for full physics documentation and expected value ranges.
 *
 * Implementation notes
 *   - Double contraction of symmetric tensors stored as 6 components:
 *       A:B = A_11*B_11 + A_22*B_22 + A_33*B_33
 *           + 2*(A_12*B_12 + A_13*B_13 + A_23*B_23)
 *     The factor-of-2 on off-diagonals is consistently applied everywhere.
 *
 *   - All loops iterate over N^3 grid points. The inner loop body is a
 *     small, branch-free arithmetic expression suitable for auto-vectorisation.
 *
 *   - No heap allocation inside the hot loops — all output vectors are
 *     pre-sized by the caller or resized once at function entry.
 *
 *   - Kahan summation is NOT used; for N=256, N^3=16.7M, double-precision
 *     round-off on a sum of O(1) values is O(N^3 * eps) ~ 3e-9, acceptable
 *     for the metrics here (all normalised to O(1)).
 */

#include "metrics.h"
#include "io.h"       // TensorField6, ScalarField, write_metrics_hdf5

#include <algorithm>  // std::transform
#include <cassert>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <numeric>    // std::inner_product
#include <stdexcept>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

// Lower-case copy of s, used for case-insensitive model name checks.
static std::string to_lower(const std::string& s) {
    std::string out(s.size(), '\0');
    std::transform(s.begin(), s.end(), out.begin(),
                   [](unsigned char c){ return std::tolower(c); });
    return out;
}

// Double contraction A:B for symmetric tensors in 6-component storage.
// Component order: 11(c=0), 12(c=1), 13(c=2), 22(c=3), 23(c=4), 33(c=5).
// Off-diagonal components are counted twice (symmetry of both tensors).
// Returns A_ij * B_ij summed over all i,j in {1,2,3}.
static inline double double_contract_6(const double* __restrict__ a,
                                       const double* __restrict__ b) {
    return a[0]*b[0]                   // 11*11
         + a[3]*b[3]                   // 22*22
         + a[5]*b[5]                   // 33*33
         + 2.0*(a[1]*b[1]             // 2 * 12*12
              + a[2]*b[2]             // 2 * 13*13
              + a[4]*b[4]);           // 2 * 23*23
}

// Squared Frobenius norm of a symmetric tensor: |A|^2 = A_ij * A_ij.
static inline double frobenius_sq_6(const double* __restrict__ a) {
    return double_contract_6(a, a);
}

// ---------------------------------------------------------------------------
// METRIC 1: correlation_coefficient
// ---------------------------------------------------------------------------

double correlation_coefficient(const TensorField6& tau_model,
                               const TensorField6& tau_exact,
                               int                 N) {
    if (tau_model.N != N || tau_exact.N != N)
        throw std::runtime_error("correlation_coefficient: N mismatch");
    if (static_cast<int>(tau_model.components.size()) != 6 ||
        static_cast<int>(tau_exact.components.size())  != 6)
        throw std::runtime_error("correlation_coefficient: expected 6 components");

    // Frobenius tensor correlation coefficient (Clark, Ferziger & Reynolds 1979):
    //
    //     rho = <A_ij B_ij> / ( sqrt(<A_ij A_ij>) * sqrt(<B_ij B_ij>) )
    //
    // where <.> is the volume average over all N^3 grid points and repeated
    // indices are summed (full double contraction over i,j in {1,2,3}).
    // This is an UNCENTERED inner-product correlation: means are NOT subtracted.
    // Off-diagonal components (12,13,23) are counted twice because their
    // symmetric counterparts (21,31,32) also contribute — handled by
    // double_contract_6.
    const long long Np3 = (long long)N * N * N;

    const double* a[6];
    const double* b[6];
    for (int c = 0; c < 6; ++c) {
        a[c] = tau_model.components[c].data();
        b[c] = tau_exact.components[c].data();
    }

    double sum_AB = 0.0, sum_AA = 0.0, sum_BB = 0.0;
    for (long long p = 0; p < Np3; ++p) {
        const double ap[6] = { a[0][p], a[1][p], a[2][p], a[3][p], a[4][p], a[5][p] };
        const double bp[6] = { b[0][p], b[1][p], b[2][p], b[3][p], b[4][p], b[5][p] };
        sum_AB += double_contract_6(ap, bp);
        sum_AA += double_contract_6(ap, ap);
        sum_BB += double_contract_6(bp, bp);
    }

    const double inv_Np3 = 1.0 / static_cast<double>(Np3);
    const double mean_AB = sum_AB * inv_Np3;
    const double mean_AA = sum_AA * inv_Np3;
    const double mean_BB = sum_BB * inv_Np3;

    // Guard against degenerate (zero-norm) fields.
    const double denom = std::sqrt(mean_AA * mean_BB);
    if (denom < 1e-30) {
        std::printf("[metrics] correlation_coefficient: denominator near zero "
                    "(degenerate field?) — returning 0\n");
        return 0.0;
    }

    return mean_AB / denom;
}

// ---------------------------------------------------------------------------
// METRIC 2: compute_dissipation_field
// ---------------------------------------------------------------------------

void compute_dissipation_field(const TensorField6&  tau,
                               const TensorField6&  S,
                               int                  N,
                               std::vector<double>& Pi) {
    if (tau.N != N || S.N != N)
        throw std::runtime_error("compute_dissipation_field: N mismatch");
    if (static_cast<int>(tau.components.size()) != 6 ||
        static_cast<int>(S.components.size())   != 6)
        throw std::runtime_error("compute_dissipation_field: expected 6 components");

    const long long Np3 = (long long)N * N * N;
    Pi.resize(Np3);

    // Pointers for each of the 6 symmetric tensor components.
    // Order: 11(0), 12(1), 13(2), 22(3), 23(4), 33(5).
    const double* t[6], *s[6];
    for (int c = 0; c < 6; ++c) {
        t[c] = tau.components[c].data();
        s[c] = S.components[c].data();
    }

    // Pi_p = -tau_ij * S_ij at each point p.
    // tau:S using symmetric storage (off-diagonals counted twice):
    //   tau:S = t11*s11 + t22*s22 + t33*s33
    //         + 2*(t12*s12 + t13*s13 + t23*s23)
    //
    // Sign convention: Pi > 0 is FORWARD scatter (energy removed from resolved
    // scales), Pi < 0 is BACKSCATTER. The minus sign is physical: when the
    // model stress aligns with -S (dissipative), the dot product -tau:S > 0.
    //
    // Note on tracelessness: for incompressible flow, S is traceless
    // (S_11+S_22+S_33 = div(u_bar) = 0), so tau:S == tau_dev:S exactly.
    // Callers may pass either tau or tau_dev — the result is the same.
    for (long long p = 0; p < Np3; ++p) {
        // Component scratch arrays at point p for double_contract_6:
        double tp[6] = { t[0][p], t[1][p], t[2][p], t[3][p], t[4][p], t[5][p] };
        double sp[6] = { s[0][p], s[1][p], s[2][p], s[3][p], s[4][p], s[5][p] };
        Pi[p] = -double_contract_6(tp, sp);
    }
}

// ---------------------------------------------------------------------------
// METRIC 3: backscatter_fraction
// ---------------------------------------------------------------------------

double backscatter_fraction(const std::vector<double>& Pi,
                            int                        N,
                            const std::string&         model_name) {
    const long long Np3 = (long long)N * N * N;
    if (static_cast<long long>(Pi.size()) != Np3)
        throw std::runtime_error("backscatter_fraction: Pi length != N^3");

    long long count_neg = 0;
    for (long long p = 0; p < Np3; ++p)
        if (Pi[p] < 0.0) ++count_neg;

    const double frac = static_cast<double>(count_neg) / static_cast<double>(Np3);

    // Purely dissipative models (Smagorinsky, WALE) must have zero backscatter
    // by construction: nu_r >= 0 everywhere, so Pi = -tau:S = 2*nu_r*|S|^2 >= 0.
    // Any nonzero backscatter fraction indicates an implementation error.
    const std::string lname = to_lower(model_name);
    const bool is_purely_dissipative =
        (lname.find("smagorinsky") != std::string::npos && lname.find("dynamic") == std::string::npos)
        || lname.find("wale") != std::string::npos;

    if (is_purely_dissipative && frac > 1e-12) {
        std::fprintf(stderr,
            "[metrics] ERROR: backscatter_fraction for model '%s' = %.6e > 0.\n"
            "  Smagorinsky and WALE are purely dissipative (nu_r >= 0 by construction),\n"
            "  so Pi = 2*nu_r*|S|^2 >= 0 at every grid point.\n"
            "  A nonzero backscatter fraction indicates a sign error in tau or Pi.\n"
            "  Reference: Pope Eq. 13.127 (tau_dev = -2*nu_r*S) and Eq. 13.123 (Pi).\n",
            model_name.c_str(), frac);
        std::abort();
    }

    return frac;
}

// ---------------------------------------------------------------------------
// METRIC 4: compute_alignment_angle
// ---------------------------------------------------------------------------

void compute_alignment_angle(const TensorField6&  tau_dev,
                             const TensorField6&  S,
                             int                  N,
                             std::vector<double>& cos_theta) {
    if (tau_dev.N != N || S.N != N)
        throw std::runtime_error("compute_alignment_angle: N mismatch");
    if (static_cast<int>(tau_dev.components.size()) != 6 ||
        static_cast<int>(S.components.size())       != 6)
        throw std::runtime_error("compute_alignment_angle: expected 6 components");

    const long long Np3 = (long long)N * N * N;
    cos_theta.resize(Np3);

    const double* td[6], *sv[6];
    for (int c = 0; c < 6; ++c) {
        td[c] = tau_dev.components[c].data();
        sv[c] = S.components[c].data();
    }

    // cos_theta_p = (tau_dev : S)_p / ( |tau_dev|_p * |S|_p )
    //
    // The eddy-viscosity hypothesis (EVH) predicts tau_dev_ij = -2*nu_r*S_ij,
    // which gives tau_dev : S = -2*nu_r*(S:S) < 0 and
    // |tau_dev| = 2*nu_r*|S|, so cos_theta = -1 everywhere when EVH holds.
    //
    // DNS a priori studies (Clark 1979, Liu et al. 1994) find a broad
    // distribution peaked somewhere between -1 and 0, reflecting that the
    // exact SGS stress is NOT well described by an eddy-viscosity form.
    //
    // Regularization: if either norm is below 1e-15, the angle is undefined
    // (zero-stress or zero-strain-rate point). Set cos_theta = 0 at those
    // points to avoid division by zero and flag them as "uncorrelated".
    constexpr double EPSILON = 1e-15;

    for (long long p = 0; p < Np3; ++p) {
        double tp[6] = { td[0][p], td[1][p], td[2][p], td[3][p], td[4][p], td[5][p] };
        double sp[6] = { sv[0][p], sv[1][p], sv[2][p], sv[3][p], sv[4][p], sv[5][p] };

        const double dot  = double_contract_6(tp, sp);
        const double n_td = std::sqrt(frobenius_sq_6(tp));
        const double n_s  = std::sqrt(frobenius_sq_6(sp));

        if (n_td < EPSILON || n_s < EPSILON) {
            cos_theta[p] = 0.0;  // undefined — regularised to uncorrelated
        } else {
            // Clamp to [-1, 1] to guard against floating-point rounding.
            cos_theta[p] = std::max(-1.0, std::min(1.0, dot / (n_td * n_s)));
        }
    }
}

// ---------------------------------------------------------------------------
// METRIC 5: compute_all_metrics
// ---------------------------------------------------------------------------

void compute_all_metrics(const TensorField6&        tau_exact_dev,
                         const TensorField6&        tau_model,
                         const TensorField6&        S,
                         const std::vector<double>& nu_r,
                         int                        N,
                         double                     delta,
                         const std::string&         model_name,
                         const std::string&         results_filepath) {
    const long long Np3 = (long long)N * N * N;

    std::printf("\n");
    std::printf("=== %s at delta=%.5f ===\n", model_name.c_str(), delta);

    // ------------------------------------------------------------------
    // 1. Correlation coefficient (model vs exact stress)
    // ------------------------------------------------------------------
    const double rho = correlation_coefficient(tau_model, tau_exact_dev, N);
    std::printf("Correlation coefficient rho  : %6.3f  (literature: ~0.3 Smag)\n", rho);

    // ------------------------------------------------------------------
    // 2. SGS energy dissipation field from MODEL stress
    // ------------------------------------------------------------------
    std::vector<double> Pi_model;
    compute_dissipation_field(tau_model, S, N, Pi_model);

    // Mean model SGS dissipation.
    double sum_Pi_model = 0.0;
    for (long long p = 0; p < Np3; ++p) sum_Pi_model += Pi_model[p];
    const double mean_Pi_model = sum_Pi_model / static_cast<double>(Np3);

    // ------------------------------------------------------------------
    // 3. Backscatter fraction (from model dissipation field)
    // ------------------------------------------------------------------
    const double bs_frac = backscatter_fraction(Pi_model, N, model_name);
    std::printf("Backscatter fraction         : %5.1f%%  (DNS truth: ~35%%)\n",
                bs_frac * 100.0);

    // ------------------------------------------------------------------
    // 4. Mean eddy viscosity
    // ------------------------------------------------------------------
    double sum_nu = 0.0;
    for (long long p = 0; p < static_cast<long long>(nu_r.size()); ++p)
        sum_nu += nu_r[p];
    const double mean_nu = (nu_r.empty()) ? 0.0 : sum_nu / static_cast<double>(nu_r.size());
    std::printf("Mean eddy viscosity          : %+.3e\n", mean_nu);

    // ------------------------------------------------------------------
    // 5. Mean SGS dissipation from the EXACT stress (DNS truth reference)
    // ------------------------------------------------------------------
    std::vector<double> Pi_exact;
    compute_dissipation_field(tau_exact_dev, S, N, Pi_exact);

    double sum_Pi_exact = 0.0;
    for (long long p = 0; p < Np3; ++p) sum_Pi_exact += Pi_exact[p];
    const double mean_Pi_exact = sum_Pi_exact / static_cast<double>(Np3);

    std::printf("Mean SGS dissipation (model) : %+.3e\n", mean_Pi_model);
    std::printf("Mean SGS dissipation (DNS)   : %+.3e\n", mean_Pi_exact);

    // DNS backscatter fraction (for reference; no assertion applied here).
    const double bs_frac_dns = backscatter_fraction(Pi_exact, N, /* no assertion */ "exact");
    std::printf("Fraction Pi < 0 (DNS exact)  : %5.1f%%\n", bs_frac_dns * 100.0);

    // ------------------------------------------------------------------
    // 6. Alignment angle between tau_dev (exact) and S
    // ------------------------------------------------------------------
    // We compare the EXACT deviatoric stress against S so the alignment
    // metric characterises how well the DNS stress satisfies the EVH —
    // i.e. it is a property of the physics, not the model. When comparing
    // tau_model vs S, eddy-viscosity models always give cos_theta = -1
    // by construction (tau_model_dev = -2*nu_r*S), so the interesting
    // quantity is always the exact stress.
    std::vector<double> cos_theta;
    compute_alignment_angle(tau_exact_dev, S, N, cos_theta);

    double sum_ct = 0.0;
    for (long long p = 0; p < Np3; ++p) sum_ct += cos_theta[p];
    const double mean_cos_theta = sum_ct / static_cast<double>(Np3);
    std::printf("Mean cos_theta (exact vs S)  : %+6.3f  (perfect EV = -1.000)\n",
                mean_cos_theta);

    // ------------------------------------------------------------------
    // Physics assessments
    // ------------------------------------------------------------------
    std::printf("\nPhysics assessments:\n");

    bool all_ok = true;

    if (rho < 0.1) {
        std::printf("  WARNING: correlation rho=%.3f is very low — check implementation "
                    "(typical Smagorinsky: 0.2–0.4)\n", rho);
        all_ok = false;
    }
    // NOTE: Smagorinsky/WALE backscatter already asserted inside backscatter_fraction.
    // Here we additionally warn if Dynamic gives implausibly high backscatter.
    {
        const std::string lname = to_lower(model_name);
        if (lname.find("dynamic") != std::string::npos && bs_frac > 0.40) {
            std::printf("  WARNING: dynamic Smagorinsky backscatter=%.1f%% > 40%% — "
                        "global clipping of Cs^2 may be needed for a posteriori use\n",
                        bs_frac * 100.0);
            all_ok = false;
        }
    }
    if (mean_cos_theta > -0.1) {
        std::printf("  NOTE: mean cos_theta=%.3f > -0.1 — "
                    "eddy-viscosity hypothesis poorly satisfied for this stress field\n",
                    mean_cos_theta);
        all_ok = false;
    }
    if (all_ok) {
        std::printf("  OK: all metrics within expected ranges\n");
    }

    // ------------------------------------------------------------------
    // 7. Write scalar metrics to HDF5
    // ------------------------------------------------------------------
    write_metrics_hdf5(results_filepath,
                       model_name,
                       rho,
                       bs_frac,
                       mean_nu,
                       mean_Pi_model,
                       mean_cos_theta);

    std::printf("\n");
}
