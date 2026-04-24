/*
 * sgs_models.cpp — Implementation of the three SGS closure models declared
 * in sgs_models.h. See the header for equation references and conventions.
 *
 * Notes on conventions used throughout this file:
 *
 *   - Symmetric 3x3 tensors use 6-component storage in the order
 *     (11, 12, 13, 22, 23, 33) = indices (0,1,2,3,4,5). This matches
 *     TensorField6 in io.h and the outputs of sgs_exact.cpp.
 *
 *   - The strain-rate magnitude supplied from sgs_exact.cpp is defined as
 *         |S| = sqrt( 2 * S_ij * S_ij )
 *     so the double contraction S:S == S_ij*S_ij == |S|^2 / 2. We use this
 *     identity to express the Smagorinsky dissipation compactly.
 *
 *   - The velocity gradient tensor g from compute_velocity_gradient_tensor
 *     is stored as 9 fields with g[3*i + j] = d(u_bar_i)/d(x_j).
 */

#include "sgs_models.h"
#include "sgs_exact.h"
#include "filter.h"
#include "config.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

namespace {

// Double contraction T:U = sum_{i,j} T_ij * U_ij for two symmetric 6-component
// tensors stored in the (11,12,13,22,23,33) convention. Off-diagonal terms
// contribute twice because the symmetry of both T and U means e.g. T_12*U_12
// appears as T_12*U_12 plus T_21*U_21.
inline double contract_sym6(const TensorField6& T, const TensorField6& U, size_t q) {
    return  T.components[0][q] * U.components[0][q]
          + T.components[3][q] * U.components[3][q]
          + T.components[5][q] * U.components[5][q]
          + 2.0 * ( T.components[1][q] * U.components[1][q]
                  + T.components[2][q] * U.components[2][q]
                  + T.components[4][q] * U.components[4][q] );
}

} // anonymous namespace

// ===========================================================================
// MODEL 1: Smagorinsky-Lilly
// ===========================================================================
// Reference equations
//   Pope (2000) Eq. 13.127:  tau_ij^dev = -2 * nu_r * S_ij
//   Pope (2000) Eq. 13.128:  nu_r = (Cs * delta)^2 * |S|
//   Pope (2000) Eq. 13.135:  Cs = 0.17  (Lilly 1967 inertial-range)
//
// Physical constraints (asserted below):
//   nu_r >= 0 everywhere                                      [|S| >= 0]
//   Pi = -tau_ij * S_ij = nu_r * |S|^2 >= 0 everywhere        [Pope p.589]
//   tau_smag is symmetric by construction                     [tau_ij = c*S_ij]
//
void smagorinsky_model(const TensorField6&        S,
                       const std::vector<double>& S_mag,
                       double                     delta,
                       int                        N,
                       TensorField6&              tau_smag,
                       std::vector<double>&       nu_r)
{
    const size_t n          = (size_t)N * N * N;
    const double Cs_delta   = CS_SMAGORINSKY * delta;
    const double Cs_delta_sq = Cs_delta * Cs_delta;

    std::printf("[smagorinsky] delta=%.6f  Cs=%.3f  (Cs*delta)^2=%.6e\n",
                delta, CS_SMAGORINSKY, Cs_delta_sq);

    nu_r.assign(n, 0.0);
    tau_smag.N    = N;
    tau_smag.name = "tau_smag";
    tau_smag.components.assign(6, std::vector<double>(n, 0.0));

    double   min_nu  = 1.0e300;
    double   max_nu  = 0.0;
    double   mean_nu = 0.0;

    // -----------------------------------------------------------------------
    // Main loop: nu_r = (Cs*delta)^2 * |S|  (Pope Eq. 13.128)
    //            tau  = -2 * nu_r * S_ij    (Pope Eq. 13.127)
    // -----------------------------------------------------------------------
    for (size_t q = 0; q < n; ++q) {
        const double nu = Cs_delta_sq * S_mag[q];

        // Assertion: nu_r >= 0 everywhere (Pope Eq. 13.128; |S| >= 0 by def.).
        if (nu < 0.0) {
            std::fprintf(stderr,
                "[smagorinsky] FATAL: nu_r<0 at index %zu (nu=%.6e). "
                "Violates Pope Eq. 13.128 (|S|>=0). Abort.\n", q, nu);
            std::abort();
        }

        nu_r[q]  = nu;
        mean_nu += nu;
        if (nu < min_nu) min_nu = nu;
        if (nu > max_nu) max_nu = nu;

        // tau_ij_dev = -2 * nu_r * S_ij (Pope Eq. 13.127). The closure is
        // symmetric because S_ij is — we just scale all 6 unique components.
        const double c = -2.0 * nu;
        tau_smag.components[0][q] = c * S.components[0][q];   // tau_11
        tau_smag.components[1][q] = c * S.components[1][q];   // tau_12
        tau_smag.components[2][q] = c * S.components[2][q];   // tau_13
        tau_smag.components[3][q] = c * S.components[3][q];   // tau_22
        tau_smag.components[4][q] = c * S.components[4][q];   // tau_23
        tau_smag.components[5][q] = c * S.components[5][q];   // tau_33
    }
    mean_nu /= (double)n;

    std::printf("[smagorinsky]   nu_r   min=%.6e  mean=%.6e  max=%.6e\n",
                min_nu, mean_nu, max_nu);
    std::printf("[smagorinsky]   (mean nu_r should scale as delta^(4/3) in "
                "the inertial subrange)\n");

    // -----------------------------------------------------------------------
    // Dissipation check: Pi = -tau:S.
    // Using tau_ij = -2*nu*S_ij, the double contraction tau:S = -2*nu*(S:S).
    // |S|^2 = 2*S_ij*S_ij  =>  S:S = |S|^2 / 2, hence
    //    Pi = -tau:S = 2*nu*(S:S) = nu * |S|^2  >=  0.
    // Reference: Pope p. 589 — Smagorinsky is purely forward-cascading.
    // -----------------------------------------------------------------------
    long long n_backscatter = 0;
    double    min_Pi = 1.0e300;
    double    mean_Pi = 0.0;
    for (size_t q = 0; q < n; ++q) {
        const double Pi = nu_r[q] * S_mag[q] * S_mag[q];
        if (Pi < 0.0)        ++n_backscatter;
        if (Pi < min_Pi)     min_Pi = Pi;
        mean_Pi += Pi;
    }
    mean_Pi /= (double)n;

    if (n_backscatter > 0) {
        std::fprintf(stderr,
            "[smagorinsky] FATAL: %lld points with Pi<0. Smagorinsky must be "
            "purely dissipative (Pope p.589). Abort.\n",
            n_backscatter);
        std::abort();
    }
    std::printf("[smagorinsky]   Pi     min=%.6e  mean=%.6e  [OK: no backscatter]\n",
                min_Pi, mean_Pi);
}

// ===========================================================================
// MODEL 2: WALE
// ===========================================================================
// Reference equations
//   Nicoud & Ducros (1999), Eq. 5:
//       Sd_ij = 0.5*(g2_ij + g2_ji) - (1/3)*delta_ij*g2_kk
//     with g2_ij = g_ik * g_kj  (matrix square of the velocity-gradient tensor)
//   Nicoud & Ducros (1999), Eq. 14:
//       nu_r = (Cw*delta)^2 * (Sd:Sd)^(3/2)
//              / ((S:S)^(5/2) + (Sd:Sd)^(5/4))
//   Modelled stress (eddy-viscosity closure):
//       tau_ij^dev = -2 * nu_r * S_ij
//
void wale_model(const std::vector<std::vector<double>>& g,
                const TensorField6&                     S,
                double                                  delta,
                int                                     N,
                TensorField6&                           tau_wale,
                std::vector<double>&                    nu_r_wale)
{
    const size_t n           = (size_t)N * N * N;
    const double Cw_delta    = CW_WALE * delta;
    const double Cw_delta_sq = Cw_delta * Cw_delta;
    constexpr double EPS     = 1.0e-20;   // stabiliser in WALE denominator

    std::printf("[wale] delta=%.6f  Cw=%.3f  (Cw*delta)^2=%.6e\n",
                delta, CW_WALE, Cw_delta_sq);

    nu_r_wale.assign(n, 0.0);
    tau_wale.N    = N;
    tau_wale.name = "tau_wale";
    tau_wale.components.assign(6, std::vector<double>(n, 0.0));

    // g[3*i + j] = d(u_bar_i)/d(x_j). Alias each of the 9 fields for clarity.
    const double* g00 = g[0].data();  // du1/dx1
    const double* g01 = g[1].data();  // du1/dx2
    const double* g02 = g[2].data();  // du1/dx3
    const double* g10 = g[3].data();  // du2/dx1
    const double* g11 = g[4].data();  // du2/dx2
    const double* g12 = g[5].data();  // du2/dx3
    const double* g20 = g[6].data();  // du3/dx1
    const double* g21 = g[7].data();  // du3/dx2
    const double* g22 = g[8].data();  // du3/dx3

    double    mean_nu    = 0.0;
    double    max_nu     = 0.0;
    long long n_negative = 0;

    for (size_t q = 0; q < n; ++q) {
        // -------------------------------------------------------------------
        // Step (a): matrix square of the gradient tensor, g2_ij = g_ik * g_kj.
        // Not symmetric in general — compute all 9 entries.
        // -------------------------------------------------------------------
        const double G00 = g00[q]*g00[q] + g01[q]*g10[q] + g02[q]*g20[q];
        const double G01 = g00[q]*g01[q] + g01[q]*g11[q] + g02[q]*g21[q];
        const double G02 = g00[q]*g02[q] + g01[q]*g12[q] + g02[q]*g22[q];
        const double G10 = g10[q]*g00[q] + g11[q]*g10[q] + g12[q]*g20[q];
        const double G11 = g10[q]*g01[q] + g11[q]*g11[q] + g12[q]*g21[q];
        const double G12 = g10[q]*g02[q] + g11[q]*g12[q] + g12[q]*g22[q];
        const double G20 = g20[q]*g00[q] + g21[q]*g10[q] + g22[q]*g20[q];
        const double G21 = g20[q]*g01[q] + g21[q]*g11[q] + g22[q]*g21[q];
        const double G22 = g20[q]*g02[q] + g21[q]*g12[q] + g22[q]*g22[q];

        // -------------------------------------------------------------------
        // Step (b): trace of g^2.
        // -------------------------------------------------------------------
        const double g2_kk = G00 + G11 + G22;

        // -------------------------------------------------------------------
        // Step (c): Sd_ij — traceless symmetric part of g^2 (Nicoud & Ducros
        // 1999, Eq. 5). We only need the 6 unique components.
        // -------------------------------------------------------------------
        const double third_tr = g2_kk / 3.0;
        const double Sd11 = G00 - third_tr;
        const double Sd22 = G11 - third_tr;
        const double Sd33 = G22 - third_tr;
        const double Sd12 = 0.5 * (G01 + G10);
        const double Sd13 = 0.5 * (G02 + G20);
        const double Sd23 = 0.5 * (G12 + G21);

        // -------------------------------------------------------------------
        // Step (d): scalar invariants (double contractions over all i,j).
        //   Sd:Sd = Sd_ij * Sd_ij  (symmetric -> diag + 2*off-diag)
        //   S:S   = S_ij  * S_ij
        // -------------------------------------------------------------------
        const double Sd_sq = Sd11*Sd11 + Sd22*Sd22 + Sd33*Sd33
                           + 2.0 * (Sd12*Sd12 + Sd13*Sd13 + Sd23*Sd23);

        const double S11 = S.components[0][q];
        const double S12 = S.components[1][q];
        const double S13 = S.components[2][q];
        const double S22 = S.components[3][q];
        const double S23 = S.components[4][q];
        const double S33 = S.components[5][q];
        const double S_sq = S11*S11 + S22*S22 + S33*S33
                          + 2.0 * (S12*S12 + S13*S13 + S23*S23);

        // -------------------------------------------------------------------
        // Step (e): WALE eddy viscosity, Nicoud & Ducros (1999) Eq. 14.
        //   nu_r = (Cw*delta)^2 * (Sd:Sd)^{3/2}
        //                        / ((S:S)^{5/2} + (Sd:Sd)^{5/4} + eps)
        // Both numerator and denominator are non-negative; EPS regularises
        // the fully laminar case where both S:S and Sd:Sd vanish.
        // -------------------------------------------------------------------
        const double numer = std::pow(Sd_sq, 1.5);
        const double denom = std::pow(S_sq, 2.5) + std::pow(Sd_sq, 1.25) + EPS;
        const double nu    = Cw_delta_sq * numer / denom;

        if (nu < 0.0) ++n_negative;

        nu_r_wale[q] = nu;
        mean_nu     += nu;
        if (nu > max_nu) max_nu = nu;

        // -------------------------------------------------------------------
        // Step (f): modelled stress, tau_wale_ij = -2 * nu_r * S_ij.
        // -------------------------------------------------------------------
        const double c = -2.0 * nu;
        tau_wale.components[0][q] = c * S11;
        tau_wale.components[1][q] = c * S12;
        tau_wale.components[2][q] = c * S13;
        tau_wale.components[3][q] = c * S22;
        tau_wale.components[4][q] = c * S23;
        tau_wale.components[5][q] = c * S33;
    }
    mean_nu /= (double)n;

    // Assertion: nu_r_wale >= 0 everywhere (WALE is purely dissipative).
    if (n_negative > 0) {
        std::fprintf(stderr,
            "[wale] FATAL: %lld points with nu_r<0. WALE must be purely "
            "dissipative (Nicoud & Ducros 1999, Eq. 14). Abort.\n",
            n_negative);
        std::abort();
    }

    // Diagnostic: near-zero nu_r_wale in low-strain / rotation-dominated
    // regions. WALE's advantage over Smagorinsky is that it vanishes there,
    // whereas Smagorinsky's nu_r only vanishes where |S| itself is zero.
    const double threshold = 0.01 * max_nu;   // 1% of peak
    long long n_near_zero = 0;
    for (size_t q = 0; q < n; ++q)
        if (nu_r_wale[q] < threshold) ++n_near_zero;

    std::printf("[wale]   nu_r   mean=%.6e  max=%.6e\n", mean_nu, max_nu);
    std::printf("[wale]   near-zero nu_r (< 1%% of max) = %lld / %zu  (%.2f%%)\n",
                n_near_zero, n, 100.0 * (double)n_near_zero / (double)n);
    std::printf("[wale]   (Smagorinsky vanishes only where |S|=0; WALE also "
                "vanishes in pure-rotation regions.)\n");
}

// ===========================================================================
// MODEL 3: Dynamic Smagorinsky (Germano-Lilly)
// ===========================================================================
// Reference equations
//   Germano (1991) Eq. 4 / Pope (2000) Eq. 13.252:
//       L_ij = filter_hat(ubar_i ubar_j) - uhatbar_i * uhatbar_j
//   Germano (1991) Eq. 8:  scale-similarity / Germano identity
//   Lilly (1992) Eq. 8 / Pope (2000) §13.6:
//       M_ij = -2*( hat_delta^2*|S_hat|*S_hat_ij
//                   - filter_hat( delta^2*|S|*S_ij ) )
//   Lilly (1992) Eq. 10:
//       Cs^2 = (Ld_ij * M_ij) / (M_ij * M_ij)    [least squares, double sum]
//
// Cs^2 is allowed to be negative pointwise (local backscatter). Clark et al.
// (1979) report 30-40% of isotropic-turbulence points having Cs^2 < 0. We
// warn if the observed fraction is inconsistent with that.
//
void dynamic_smagorinsky(const VelocityField&       vel_bar,
                         const TensorField6&        S,
                         const std::vector<double>& S_mag,
                         double                     delta,
                         double                     dx,
                         int                        N,
                         TensorField6&              tau_dyn,
                         std::vector<double>&       Cs2_pointwise,
                         std::vector<double>&       Cs2_effective,
                         AveragingMode              avg_mode)
{
    const size_t n         = (size_t)N * N * N;
    // Test-filter ratio alpha = 2 (Germano 1991, Lilly 1992 convention).
    const double hat_delta = 2.0 * delta;
    // Test filter implemented as a physical-space TOP-HAT (box) of width
    // hat_delta.  The primary filter remains sharp spectral; the test
    // filter is top-hat because the sharp spectral filter is ill-posed in
    // the Germano identity (Meneveau & Katz, ARFM 2000, sec. 3.2) and
    // produces wrong-signed mean Cs^2 when used as the test filter.
    // Box radius in grid points: half the box width.
    int test_box_r = static_cast<int>(std::lround(hat_delta / (2.0 * dx)));
    if (test_box_r < 1) test_box_r = 1;
    constexpr double EPS   = 1.0e-20;

    const char* mode_str = (avg_mode == AveragingMode::Pointwise) ? "Pointwise"
                         : (avg_mode == AveragingMode::Box)        ? "Box (isotropic)"
                                                                   : "XZ-plane (channel)";
    std::printf("[dyn_smag] delta=%.6f  hat_delta=%.6f (= 2*delta)  "
                "test_filter=top-hat (r=%d, window=%d*dx)  avg=%s\n",
                delta, hat_delta, test_box_r, 2 * test_box_r + 1, mode_str);

    // -----------------------------------------------------------------------
    // (b) Germano Leonard stress L_ij (Germano 1991 Eq. 4, Pope Eq. 13.252)
    //       L_ij = filter_hat(ubar_i * ubar_j) - uhatbar_i * uhatbar_j
    //     L_ij represents the stress contribution from scales between delta
    //     and hat_delta and is computable directly from the filtered field.
    // -----------------------------------------------------------------------
    VelocityField vbar_hat;
    vbar_hat.N  = N;
    vbar_hat.dx = dx;
    vbar_hat.u1.assign(n, 0.0);
    vbar_hat.u2.assign(n, 0.0);
    vbar_hat.u3.assign(n, 0.0);
    // Test filter = 3-D top-hat (physical-space box) of radius `test_box_r`.
    box_filter_3d(vel_bar.u1.data(), vbar_hat.u1.data(), N, test_box_r);
    box_filter_3d(vel_bar.u2.data(), vbar_hat.u2.data(), N, test_box_r);
    box_filter_3d(vel_bar.u3.data(), vbar_hat.u3.data(), N, test_box_r);

    TensorField6 L;
    L.N    = N;
    L.name = "L_dyn";
    L.components.assign(6, std::vector<double>(n, 0.0));

    const int pair_i[6] = {0, 0, 0, 1, 1, 2};
    const int pair_j[6] = {0, 1, 2, 1, 2, 2};
    const std::vector<double>* ub [3] = { &vel_bar.u1,  &vel_bar.u2,  &vel_bar.u3 };
    const std::vector<double>* ubh[3] = { &vbar_hat.u1, &vbar_hat.u2, &vbar_hat.u3 };

    for (int p = 0; p < 6; ++p) {
        const int a = pair_i[p];
        const int b = pair_j[p];

        // filter_hat(ubar_a * ubar_b) — product formed first, then filtered.
        // This must filter the product (not filter each factor separately);
        // the distinction is the same subtlety Pope flags at p.558.
        // Uses top-hat test filter to match the velocity test filter above.
        box_filter_product_3d(ub[a]->data(), ub[b]->data(),
                              L.components[p].data(), N, test_box_r);

        // Subtract uhatbar_a * uhatbar_b pointwise (Germano Eq. 4).
        const double* ah = ubh[a]->data();
        const double* bh = ubh[b]->data();
        double*       out = L.components[p].data();
        for (size_t q = 0; q < n; ++q)
            out[q] -= ah[q] * bh[q];
    }

    // -----------------------------------------------------------------------
    // (c) Deviatoric Leonard stress  Ld_ij = L_ij - (1/3)*delta_ij*L_kk
    // -----------------------------------------------------------------------
    TensorField6 Ld;
    Ld.N    = N;
    Ld.name = "Ld_dyn";
    Ld.components.assign(6, std::vector<double>(n, 0.0));
    for (size_t q = 0; q < n; ++q) {
        const double L11 = L.components[0][q];
        const double L22 = L.components[3][q];
        const double L33 = L.components[5][q];
        const double third_tr = (L11 + L22 + L33) / 3.0;
        Ld.components[0][q] = L11 - third_tr;
        Ld.components[1][q] = L.components[1][q];   // off-diagonals unchanged
        Ld.components[2][q] = L.components[2][q];
        Ld.components[3][q] = L22 - third_tr;
        Ld.components[4][q] = L.components[4][q];
        Ld.components[5][q] = L33 - third_tr;
    }

    // -----------------------------------------------------------------------
    // (d) S_hat and |S_hat| from the test-filtered velocity uhatbar.
    //     We reuse the spectral-derivative machinery from sgs_exact.cpp
    //     rather than differentiating S itself — directly differentiating
    //     ubar_hat is the cleanest way to ensure the result is consistent
    //     with the Germano identity.
    // -----------------------------------------------------------------------
    std::vector<std::vector<double>> g_hat;
    compute_velocity_gradient_tensor(vbar_hat, dx, g_hat);

    TensorField6        S_hat;
    std::vector<double> S_hat_mag;
    compute_strain_rate_tensor(g_hat, N, S_hat, S_hat_mag);

    // -----------------------------------------------------------------------
    // (e) M_ij (Lilly 1992 Eq. 8 / Pope 2000 §13.6)
    //
    //   From the Germano identity: L_ij^dev = C_s^2 * M_ij, where
    //       M_ij = -2*( hat_delta^2 * |S_hat| * S_hat_ij
    //                   - filter_hat( delta^2 * |S| * S_ij ) )
    //
    //   The leading NEGATIVE sign is required so that L_ij = C_s^2 * M_ij
    //   with C_s^2 > 0 for forward-cascading flow (Pope Eq. 13.253-13.254).
    //   Using +2 flips the sign of every C_s^2_eff (bug: negative eddy
    //   viscosity everywhere).
    //
    // CRITICAL: filter the COMBINED field alpha_ij = 2*delta^2*|S|*S_ij,
    // not S_ij alone. Filtering factors separately is the most common
    // implementation bug for the dynamic model. Reference: Pope p.607.
    // -----------------------------------------------------------------------
    TensorField6 M;
    M.N    = N;
    M.name = "M_dyn";
    M.components.assign(6, std::vector<double>(n, 0.0));

    const double two_hd2 = 2.0 * hat_delta * hat_delta;
    const double two_d2  = 2.0 * delta     * delta;

    std::vector<double> alpha     (n, 0.0);
    std::vector<double> filt_alpha(n, 0.0);

    for (int p = 0; p < 6; ++p) {
        // alpha_ij = 2*delta^2 * |S| * S_ij — assembled first on the grid.
        for (size_t q = 0; q < n; ++q)
            alpha[q] = two_d2 * S_mag[q] * S.components[p][q];

        // filter_hat(alpha_ij) using the same top-hat test filter as velocity.
        box_filter_3d(alpha.data(), filt_alpha.data(), N, test_box_r);

        // M_ij = -2*hat_delta^2*|S_hat|*S_hat_ij + filt_alpha_ij
        // Note: sign is negative on the first term (Lilly 1992 / Pope §13.6).
        double* out = M.components[p].data();
        for (size_t q = 0; q < n; ++q)
            out[q] = filt_alpha[q]
                   - two_hd2 * S_hat_mag[q] * S_hat.components[p][q];
    }

    // -----------------------------------------------------------------------
    // (f) Numerator and denominator fields for Lilly least-squares.
    //     num[q] = (Ld:M)[q],  den[q] = (M:M)[q].
    // -----------------------------------------------------------------------
    std::vector<double> num_LdM(n, 0.0);
    std::vector<double> den_MM (n, 0.0);
    for (size_t q = 0; q < n; ++q) {
        num_LdM[q] = contract_sym6(Ld, M, q);
        den_MM [q] = contract_sym6(M,  M, q);
    }

    // Pointwise Germano-Lilly Cs^2 field (diagnostic; fat-tailed with a sharp
    // spectral test filter). Kept for the Cs^2 PDF plot.
    Cs2_pointwise.assign(n, 0.0);
    long long n_neg_point = 0;
    double    mean_Cs2_pt = 0.0;
    double    min_Cs2_pt  =  1.0e300;
    double    max_Cs2_pt  = -1.0e300;
    for (size_t q = 0; q < n; ++q) {
        const double Cs2 = num_LdM[q] / (den_MM[q] + EPS);
        Cs2_pointwise[q] = Cs2;
        mean_Cs2_pt += Cs2;
        if (Cs2 < 0.0)        ++n_neg_point;
        if (Cs2 < min_Cs2_pt) min_Cs2_pt = Cs2;
        if (Cs2 > max_Cs2_pt) max_Cs2_pt = Cs2;
    }
    mean_Cs2_pt /= (double)n;

    // -----------------------------------------------------------------------
    // (f.2) Effective Cs^2 field with homogeneous-direction averaging applied
    //       to the numerator and denominator SEPARATELY (Lilly 1992).
    //       This is what `tau_dyn` is actually built from.
    // -----------------------------------------------------------------------
    Cs2_effective.assign(n, 0.0);

    if (avg_mode == AveragingMode::Pointwise) {
        // No averaging — same as pointwise field.
        for (size_t q = 0; q < n; ++q)
            Cs2_effective[q] = num_LdM[q] / (den_MM[q] + EPS);
    }
    else if (avg_mode == AveragingMode::Box) {
        // Full-box average (assumes full statistical homogeneity in all 3
        // directions, appropriate for forced isotropic turbulence).
        long double sum_num = 0.0L, sum_den = 0.0L;
        for (size_t q = 0; q < n; ++q) {
            sum_num += num_LdM[q];
            sum_den += den_MM [q];
        }
        const double Cs2_scalar =
            static_cast<double>(sum_num / (sum_den + (long double)EPS));
        for (size_t q = 0; q < n; ++q)
            Cs2_effective[q] = Cs2_scalar;
        std::printf("[dyn_smag]   box-averaged scalar Cs^2 = %+.6e  "
                    "(Lilly 1992 target: %.6e)\n",
                    Cs2_scalar, CS_SMAGORINSKY * CS_SMAGORINSKY);
    }
    else {
        // XZ-plane averaging (channel flow: homogeneous in x and z, NOT y).
        // Row-major layout: idx(i,j,k) = i*N*N + j*N + k where i is x-index
        // (fastest-varying block), j is y-index, k is z-index. Average at
        // each fixed y over all (i,k) pairs.
        std::vector<long double> plane_num(N, 0.0L);
        std::vector<long double> plane_den(N, 0.0L);
        std::vector<double>      Cs2_y    (N, 0.0);
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                for (int k = 0; k < N; ++k) {
                    const size_t q = (size_t)i * N * N + (size_t)j * N + k;
                    plane_num[j] += num_LdM[q];
                    plane_den[j] += den_MM [q];
                }
        for (int j = 0; j < N; ++j) {
            Cs2_y[j] = static_cast<double>(
                plane_num[j] / (plane_den[j] + (long double)EPS));
        }
        // Broadcast back to a 3-D field uniform in (x,z).
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                for (int k = 0; k < N; ++k) {
                    const size_t q = (size_t)i * N * N + (size_t)j * N + k;
                    Cs2_effective[q] = Cs2_y[j];
                }
        // Summary: Cs^2(y) range and mid-plane value.
        double Cs2_min_y =  1.0e300, Cs2_max_y = -1.0e300, Cs2_mean_y = 0.0;
        for (int j = 0; j < N; ++j) {
            if (Cs2_y[j] < Cs2_min_y) Cs2_min_y = Cs2_y[j];
            if (Cs2_y[j] > Cs2_max_y) Cs2_max_y = Cs2_y[j];
            Cs2_mean_y += Cs2_y[j];
        }
        Cs2_mean_y /= (double)N;
        std::printf("[dyn_smag]   plane-averaged Cs^2(y): "
                    "min=%+.4e  mean=%+.4e  max=%+.4e  "
                    "(Lilly 1992 target: %.4e)\n",
                    Cs2_min_y, Cs2_mean_y, Cs2_max_y,
                    CS_SMAGORINSKY * CS_SMAGORINSKY);
    }

    // Recompute summary stats on the effective field for the physics check
    // further below.
    long long n_neg_Cs2 = 0;
    double    mean_Cs2  = 0.0;
    double    min_Cs2   =  1.0e300;
    double    max_Cs2   = -1.0e300;
    for (size_t q = 0; q < n; ++q) {
        const double c = Cs2_effective[q];
        mean_Cs2 += c;
        if (c < 0.0)     ++n_neg_Cs2;
        if (c < min_Cs2) min_Cs2 = c;
        if (c > max_Cs2) max_Cs2 = c;
    }
    mean_Cs2 /= (double)n;

    // -----------------------------------------------------------------------
    // (g) Modelled stress  tau_dyn_ij = -2 * Cs^2_eff * delta^2 * |S| * S_ij
    //     Uses the HOMOGENEOUSLY-AVERAGED Cs^2_effective field (Lilly 1992).
    //     Cs^2_eff can still be negative locally — physical backscatter.
    // -----------------------------------------------------------------------
    tau_dyn.N    = N;
    tau_dyn.name = "tau_dyn";
    tau_dyn.components.assign(6, std::vector<double>(n, 0.0));

    const double d2 = delta * delta;
    for (size_t q = 0; q < n; ++q) {
        const double coef = -2.0 * Cs2_effective[q] * d2 * S_mag[q];
        tau_dyn.components[0][q] = coef * S.components[0][q];
        tau_dyn.components[1][q] = coef * S.components[1][q];
        tau_dyn.components[2][q] = coef * S.components[2][q];
        tau_dyn.components[3][q] = coef * S.components[3][q];
        tau_dyn.components[4][q] = coef * S.components[4][q];
        tau_dyn.components[5][q] = coef * S.components[5][q];
    }

    // Print a diagnostic line comparing pointwise vs effective Cs^2 stats so
    // the reader can see the ill-conditioning of the pointwise procedure.
    std::printf("[dyn_smag]   Cs^2 pointwise: min=%+.4e mean=%+.4e max=%+.4e  "
                "(negative frac %.1f%%)\n",
                min_Cs2_pt, mean_Cs2_pt, max_Cs2_pt,
                100.0 * (double)n_neg_point / (double)n);

    // -----------------------------------------------------------------------
    // Diagnostics on the EFFECTIVE (averaged) Cs^2 field.
    // -----------------------------------------------------------------------
    const double frac_neg = (double)n_neg_Cs2 / (double)n;
    const double Cs2_target = CS_SMAGORINSKY * CS_SMAGORINSKY;
    std::printf("[dyn_smag]   Cs2 effective: min=%+.4e mean=%+.4e max=%+.4e  "
                "(fraction Cs^2_eff < 0: %.2f%%)\n",
                min_Cs2, mean_Cs2, max_Cs2, 100.0 * frac_neg);
    std::printf("[dyn_smag]   mean Cs2_eff = %+.5f  "
                "(Lilly 1992 target for homogeneous flow: %.5f)\n",
                mean_Cs2, Cs2_target);

    // Pointwise-mode only: the fraction of negative Cs^2 should reflect the
    // Clark (1979) ~30-40% backscatter statistic. With averaging, the
    // effective field is highly structured (a constant for Box or a y-profile
    // for XZ_Plane), so this fraction is not a meaningful backscatter proxy.
    if (avg_mode == AveragingMode::Pointwise) {
        if (frac_neg < 0.01) {
            std::fprintf(stderr,
                "[dyn_smag] WARN: only %.2f%% negative Cs2 (pointwise) — "
                "Germano identity likely broken. Check Pope Eq. 13.252.\n",
                100.0 * frac_neg);
        }
        if (frac_neg > 0.60) {
            std::fprintf(stderr,
                "[dyn_smag] WARN: %.2f%% negative Cs2 (pointwise) > 60%% — "
                "likely sign error in M_ij. Check Lilly 1992 Eq. 8.\n",
                100.0 * frac_neg);
        }
    }

    // Mean Cs^2_eff should recover Cs=0.17 (=> Cs^2=0.0289) for fully
    // homogeneous flow. A 50% deviation is a warning, not fatal — coarse
    // averaging over a moderate-Re volume can legitimately be off.
    if (std::abs(mean_Cs2 - Cs2_target) > 0.5 * Cs2_target) {
        std::fprintf(stderr,
            "[dyn_smag] WARN: mean Cs2_eff=%+.5f deviates strongly from "
            "expected %.5f (Lilly 1992). Acceptable for a moderate-Re "
            "sub-volume; investigate if persistent across filter widths.\n",
            mean_Cs2, Cs2_target);
    }
}
