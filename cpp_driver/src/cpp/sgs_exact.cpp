/*
 * sgs_exact.cpp — Exact SGS stress, velocity gradient, strain-rate tensor.
 *
 * See sgs_exact.h for the mathematical definitions and references. This file
 * is purely an implementation of those formulae. No modelling happens here —
 * everything is derived directly from the DNS field by filtering and spectral
 * differentiation.
 *
 * Spectral derivatives
 *   We maintain a small FFTW plan cache local to this translation unit, with
 *   the same structure as the one in filter.cpp. FFTW wisdom is global state
 *   inside the FFTW library, so if filter.cpp has already loaded wisdom from
 *   disk, our plan builds here are effectively free. The plans themselves
 *   are intentionally leaked at program exit (same policy as filter.cpp).
 *
 * Nyquist handling
 *   The derivative of a real signal at the Nyquist frequency cannot be stored
 *   in an r2c layout without breaking Hermitian symmetry, because i * k *
 *   (real number) is purely imaginary and there is no negative-k partner to
 *   complex-conjugate it with. We therefore zero the Nyquist mode on the axis
 *   being differentiated before multiplying by i*k_j. This matches the
 *   convention used by references/spectralDNS when computing derivatives on
 *   a triply-periodic box.
 */

#include "sgs_exact.h"
#include "config.h"
#include "filter.h"

#include <fftw3.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

// ---------------------------------------------------------------------------
// Spectral-derivative plan cache (file-local)
// ---------------------------------------------------------------------------

namespace {

struct SpectralDerivContext {
    int            N             = 0;
    double*        real_buf      = nullptr;   // N^3 doubles
    fftw_complex*  freq_buf      = nullptr;   // N*N*(N/2+1) complex — working
    fftw_complex*  saved_coeffs  = nullptr;   // N*N*(N/2+1) complex — fwd result
    fftw_plan      plan_fwd      = nullptr;
    fftw_plan      plan_bwd      = nullptr;
};

SpectralDerivContext& get_deriv_context(int N) {
    static SpectralDerivContext ctx;
    if (ctx.N == N && ctx.plan_fwd && ctx.plan_bwd)
        return ctx;

    // Rebuild on first use or if N changed.
    if (ctx.plan_fwd)     { fftw_destroy_plan(ctx.plan_fwd);   ctx.plan_fwd = nullptr; }
    if (ctx.plan_bwd)     { fftw_destroy_plan(ctx.plan_bwd);   ctx.plan_bwd = nullptr; }
    if (ctx.real_buf)     { fftw_free(ctx.real_buf);           ctx.real_buf = nullptr; }
    if (ctx.freq_buf)     { fftw_free(ctx.freq_buf);           ctx.freq_buf = nullptr; }
    if (ctx.saved_coeffs) { fftw_free(ctx.saved_coeffs);       ctx.saved_coeffs = nullptr; }

    ctx.N = N;
    const size_t n_real = (size_t)N * N * N;
    const size_t n_comp = (size_t)N * N * (N / 2 + 1);

    ctx.real_buf     = (double*)       fftw_malloc(sizeof(double)       * n_real);
    ctx.freq_buf     = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n_comp);
    ctx.saved_coeffs = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n_comp);
    if (!ctx.real_buf || !ctx.freq_buf || !ctx.saved_coeffs) {
        std::fprintf(stderr, "[sgs_exact] FATAL: fftw_malloc failed for N=%d\n", N);
        std::abort();
    }

    ctx.plan_fwd = fftw_plan_dft_r2c_3d(N, N, N,
                                        ctx.real_buf, ctx.freq_buf,
                                        FFTW_MEASURE);
    ctx.plan_bwd = fftw_plan_dft_c2r_3d(N, N, N,
                                        ctx.freq_buf, ctx.real_buf,
                                        FFTW_MEASURE);
    if (!ctx.plan_fwd || !ctx.plan_bwd) {
        std::fprintf(stderr, "[sgs_exact] FATAL: fftw_plan creation failed for N=%d\n", N);
        std::abort();
    }

    std::printf("[sgs_exact] deriv plan ready: N=%d  r2c+c2r  FFTW_MEASURE\n", N);
    return ctx;
}

inline int signed_wavenumber(int i, int N) {
    return (i <= N / 2) ? i : i - N;
}

// Multiply the r2c coefficients in `freq` by (i * k_axis) in place. The
// Nyquist mode on the differentiation axis is explicitly zeroed so that the
// result respects Hermitian symmetry and inverse-transforms cleanly to a real
// field.
//
//   freq         : coefficients in FFTW r2c layout, N x N x (N/2+1)
//   N            : grid size
//   axis         : 0, 1, or 2 — which spatial direction to differentiate by
//   two_pi_over_L: 2*pi / (N*dx), conversion from signed index to physical k
//
void multiply_by_ik(fftw_complex* freq, int N, int axis, double two_pi_over_L) {
    const int    Nzh     = N / 2 + 1;
    const bool   even_N  = (N % 2 == 0);
    const int    nyq_ax  = N / 2;

    for (int i = 0; i < N; ++i) {
        const int nx = signed_wavenumber(i, N);
        for (int j = 0; j < N; ++j) {
            const int ny = signed_wavenumber(j, N);
            const size_t row = ((size_t)i * N + j) * Nzh;
            for (int k = 0; k < Nzh; ++k) {
                const size_t ci = row + (size_t)k;

                // Nyquist detection for the differentiation axis.
                bool is_nyquist = false;
                if (even_N) {
                    if (axis == 0 && i == nyq_ax) is_nyquist = true;
                    if (axis == 1 && j == nyq_ax) is_nyquist = true;
                    if (axis == 2 && k == nyq_ax) is_nyquist = true;
                }
                if (is_nyquist) {
                    freq[ci][0] = 0.0;
                    freq[ci][1] = 0.0;
                    continue;
                }

                double kj;
                if      (axis == 0) kj = two_pi_over_L * nx;
                else if (axis == 1) kj = two_pi_over_L * ny;
                else                kj = two_pi_over_L * k;   // axis 2: n_z = k

                // (i * kj) * (re + i*im) = -kj*im + i*(kj*re)
                const double re = freq[ci][0];
                const double im = freq[ci][1];
                freq[ci][0] = -kj * im;
                freq[ci][1] =  kj * re;
            }
        }
    }
}

} // anonymous namespace

// ===========================================================================
// compute_exact_sgs_stress  —  Pope Eq. 13.93
// ===========================================================================
//
// The order of operations matters:
//   filter(u_i u_j)  requires computing the product on the FULL DNS grid
//                    then filtering the result.
//   u_bar_i * u_bar_j means filter u_i and u_j separately first.
// These two operations give different results — their difference IS tau_ij.
// Reference: Pope (2000) p. 558, discussion after Eq. 13.93.

void compute_exact_sgs_stress(const VelocityField& vel,
                              double               delta,
                              double               dx,
                              TensorField6&        tau,
                              TensorField6&        tau_dev,
                              VelocityField&       vel_bar)
{
    const int    N      = vel.N;
    const size_t n_real = (size_t)N * N * N;

    std::printf("[sgs_exact] compute_exact_sgs_stress  N=%d  delta=%.6f\n", N, delta);

    // -----------------------------------------------------------------------
    // Step (a): filter each velocity component to obtain u_bar_i.
    // -----------------------------------------------------------------------
    vel_bar.N  = N;
    vel_bar.dx = dx;
    vel_bar.u1.assign(n_real, 0.0);
    vel_bar.u2.assign(n_real, 0.0);
    vel_bar.u3.assign(n_real, 0.0);

    spectral_filter_3d(vel.u1.data(), vel_bar.u1.data(), N, delta, dx);
    spectral_filter_3d(vel.u2.data(), vel_bar.u2.data(), N, delta, dx);
    spectral_filter_3d(vel.u3.data(), vel_bar.u3.data(), N, delta, dx);

    // -----------------------------------------------------------------------
    // Step (b): for each symmetric pair (i,j), compute
    //     tau_ij = filter(u_i u_j) - u_bar_i * u_bar_j
    // Storage order (matches TensorField6): 11, 12, 13, 22, 23, 33.
    // -----------------------------------------------------------------------
    tau.N    = N;
    tau.name = "tau_exact";
    tau.components.assign(6, std::vector<double>(n_real, 0.0));

    const int pair_i[6] = {0, 0, 0, 1, 1, 2};
    const int pair_j[6] = {0, 1, 2, 1, 2, 2};

    const std::vector<double>* u_raw[3]    = { &vel.u1,     &vel.u2,     &vel.u3     };
    const std::vector<double>* u_filt[3]   = { &vel_bar.u1, &vel_bar.u2, &vel_bar.u3 };

    for (int p_idx = 0; p_idx < 6; ++p_idx) {
        const int a = pair_i[p_idx];
        const int b = pair_j[p_idx];

        // filter(u_a * u_b) — product formed on DNS grid, then filtered.
        filter_product(u_raw[a]->data(), u_raw[b]->data(),
                       tau.components[p_idx].data(), N, delta, dx);

        // Subtract u_bar_a * u_bar_b elementwise.
        const double* ua_f = u_filt[a]->data();
        const double* ub_f = u_filt[b]->data();
        double*       out  = tau.components[p_idx].data();
        for (size_t q = 0; q < n_real; ++q)
            out[q] -= ua_f[q] * ub_f[q];
    }

    // -----------------------------------------------------------------------
    // Step (c): deviatoric part tau_dev_ij = tau_ij - (1/3) tr(tau) delta_ij.
    // -----------------------------------------------------------------------
    tau_dev.N    = N;
    tau_dev.name = "tau_exact_dev";
    tau_dev.components.assign(6, std::vector<double>(n_real, 0.0));

    for (size_t q = 0; q < n_real; ++q) {
        const double t11 = tau.components[0][q];
        const double t22 = tau.components[3][q];
        const double t33 = tau.components[5][q];
        const double third_trace = (t11 + t22 + t33) / 3.0;

        tau_dev.components[0][q] = t11 - third_trace;
        tau_dev.components[1][q] = tau.components[1][q];   // off-diagonal
        tau_dev.components[2][q] = tau.components[2][q];   // off-diagonal
        tau_dev.components[3][q] = t22 - third_trace;
        tau_dev.components[4][q] = tau.components[4][q];   // off-diagonal
        tau_dev.components[5][q] = t33 - third_trace;
    }

    // -----------------------------------------------------------------------
    // Physical verification printouts.
    //
    // Expected for isotropic turbulence:
    //   mean(tau_12)  ~ 0         (off-diagonal stress has no preferred sign)
    //   rms(tau_11)   > 0         (non-trivial normal stress)
    //   mean(tr tau)  > 0         (= 2 * k_sgs, the SGS turbulent KE)
    //
    // Note on pointwise trace:
    //   For a positivity-preserving filter (Gaussian, box) the trace is
    //   non-negative everywhere: tr(tau) = 2 k_sgs >= 0 pointwise.
    //   For the SHARP SPECTRAL filter used here the bound is only guaranteed
    //   IN THE MEAN — pointwise negative values are allowed because the
    //   sharp-cutoff kernel is not positive in physical space. We therefore
    //   print the minimum and the fraction of negative points as a diagnostic
    //   rather than aborting on a negative value.
    // -----------------------------------------------------------------------
    double   mean_tau12 = 0.0;
    double   sq_tau11   = 0.0;
    double   mean_trace = 0.0;
    double   min_trace  = 1.0e300;
    long long n_neg_trace = 0;

    for (size_t q = 0; q < n_real; ++q) {
        const double t11   = tau.components[0][q];
        const double t12   = tau.components[1][q];
        const double t22   = tau.components[3][q];
        const double t33   = tau.components[5][q];
        const double trace = t11 + t22 + t33;

        mean_tau12 += t12;
        sq_tau11   += t11 * t11;
        mean_trace += trace;
        if (trace < min_trace)  min_trace = trace;
        if (trace < 0.0)        ++n_neg_trace;
    }
    mean_tau12 /= (double)n_real;
    sq_tau11   /= (double)n_real;
    mean_trace /= (double)n_real;
    const double rms_tau11 = std::sqrt(sq_tau11);

    std::printf("[sgs_exact]   mean(tau_12)       = %+.6e   (should be ~0)\n", mean_tau12);
    std::printf("[sgs_exact]   rms(tau_11)        = %.6e    (should be > 0)\n", rms_tau11);
    std::printf("[sgs_exact]   mean(tr tau)       = %+.6e   (= 2 k_sgs, should be > 0)\n",
                mean_trace);
    std::printf("[sgs_exact]   min(tr tau)        = %+.6e   (sharp-filter: may be < 0 locally)\n",
                min_trace);
    std::printf("[sgs_exact]   points with tr<0   = %lld / %zu  (%.2f%%)\n",
                n_neg_trace, n_real, 100.0 * (double)n_neg_trace / (double)n_real);
}

// ===========================================================================
// compute_velocity_gradient_tensor  —  spectral differentiation
// ===========================================================================

void compute_velocity_gradient_tensor(const VelocityField&              vel_bar,
                                      double                            dx,
                                      std::vector<std::vector<double>>& g)
{
    const int    N      = vel_bar.N;
    const size_t n_real = (size_t)N * N * N;
    const size_t n_comp = (size_t)N * N * (N / 2 + 1);
    const double L      = (double)N * dx;
    const double two_pi_over_L = 2.0 * M_PI / L;
    const double norm   = 1.0 / (double)n_real;

    std::printf("[sgs_exact] compute_velocity_gradient_tensor  N=%d  L=%.6f\n", N, L);

    SpectralDerivContext& ctx = get_deriv_context(N);

    // Output: 9 fields, each N^3.
    g.assign(9, std::vector<double>(n_real, 0.0));

    const std::vector<double>* u_ptrs[3] = {
        &vel_bar.u1, &vel_bar.u2, &vel_bar.u3
    };

    // For each velocity component i: one forward FFT, then three i*k multiplies
    // and inverse FFTs (one per spatial derivative direction j).
    for (int i = 0; i < 3; ++i) {
        std::memcpy(ctx.real_buf, u_ptrs[i]->data(), n_real * sizeof(double));
        fftw_execute(ctx.plan_fwd);

        // Cache the forward coefficients so we can reuse them for all three
        // derivative directions without re-running the forward transform.
        std::memcpy(ctx.saved_coeffs, ctx.freq_buf, n_comp * sizeof(fftw_complex));

        for (int j = 0; j < 3; ++j) {
            // Restore cached forward coefficients into the working buffer.
            std::memcpy(ctx.freq_buf, ctx.saved_coeffs, n_comp * sizeof(fftw_complex));

            // Apply i * k_j in place (Nyquist on axis j is zeroed inside).
            multiply_by_ik(ctx.freq_buf, N, j, two_pi_over_L);

            // Inverse transform and normalise.
            fftw_execute(ctx.plan_bwd);
            double* out = g[3 * i + j].data();
            for (size_t q = 0; q < n_real; ++q)
                out[q] = ctx.real_buf[q] * norm;
        }
    }

    // -----------------------------------------------------------------------
    // Verification: on a triply periodic box, the mean of any derivative
    // component must be zero (the constant mode is killed by multiplication
    // by i*k with k=0 giving 0). Print the worst offender as a sanity check.
    // -----------------------------------------------------------------------
    double worst_abs_mean = 0.0;
    int    worst_comp     = -1;
    for (int c = 0; c < 9; ++c) {
        double sum = 0.0;
        for (size_t q = 0; q < n_real; ++q) sum += g[c][q];
        const double mean = sum / (double)n_real;
        if (std::abs(mean) > worst_abs_mean) {
            worst_abs_mean = std::abs(mean);
            worst_comp     = c;
        }
    }
    std::printf("[sgs_exact]   max|mean(g)| over 9 components = %.3e  (component %d; should be ~0)\n",
                worst_abs_mean, worst_comp);
}

// ===========================================================================
// compute_strain_rate_tensor  —  Pope Eq. 13.73-13.74
// ===========================================================================

void compute_strain_rate_tensor(const std::vector<std::vector<double>>& g,
                                int                                     N,
                                TensorField6&                           S,
                                std::vector<double>&                    S_mag)
{
    const size_t n_real = (size_t)N * N * N;

    S.N    = N;
    S.name = "S_bar";
    S.components.assign(6, std::vector<double>(n_real, 0.0));
    S_mag.assign(n_real, 0.0);

    // g[3*i + j] = du_i / dx_j, 0-based. Extract as scalar fields for clarity.
    const double* g00 = g[0].data();  // du1/dx1
    const double* g01 = g[1].data();  // du1/dx2
    const double* g02 = g[2].data();  // du1/dx3
    const double* g10 = g[3].data();  // du2/dx1
    const double* g11 = g[4].data();  // du2/dx2
    const double* g12 = g[5].data();  // du2/dx3
    const double* g20 = g[6].data();  // du3/dx1
    const double* g21 = g[7].data();  // du3/dx2
    const double* g22 = g[8].data();  // du3/dx3

    double* S11_out = S.components[0].data();
    double* S12_out = S.components[1].data();
    double* S13_out = S.components[2].data();
    double* S22_out = S.components[3].data();
    double* S23_out = S.components[4].data();
    double* S33_out = S.components[5].data();

    for (size_t q = 0; q < n_real; ++q) {
        const double S11 = g00[q];
        const double S22 = g11[q];
        const double S33 = g22[q];
        const double S12 = 0.5 * (g01[q] + g10[q]);
        const double S13 = 0.5 * (g02[q] + g20[q]);
        const double S23 = 0.5 * (g12[q] + g21[q]);

        S11_out[q] = S11;
        S12_out[q] = S12;
        S13_out[q] = S13;
        S22_out[q] = S22;
        S23_out[q] = S23;
        S33_out[q] = S33;

        // |S|^2 = 2 * S_ij * S_ij (double contraction over all i, j)
        //       = 2 * [ S11^2 + S22^2 + S33^2 + 2*(S12^2 + S13^2 + S23^2) ]
        const double s2 = 2.0 * ( S11*S11 + S22*S22 + S33*S33
                                + 2.0 * (S12*S12 + S13*S13 + S23*S23) );
        S_mag[q] = std::sqrt(s2);
    }

    // -----------------------------------------------------------------------
    // Verification: |S| is non-negative by construction. Print min/mean/max
    // of |S| and the trace of S (which should be ~0 because incompressibility
    // gives div(u_bar) = 0 = tr(S) for the filtered field).
    // -----------------------------------------------------------------------
    double min_Smag = 1.0e300;
    double max_Smag = 0.0;
    double mean_Smag = 0.0;
    double mean_trS  = 0.0;
    double max_abs_trS = 0.0;

    for (size_t q = 0; q < n_real; ++q) {
        const double s = S_mag[q];
        if (s < min_Smag) min_Smag = s;
        if (s > max_Smag) max_Smag = s;
        mean_Smag += s;

        const double trS = S11_out[q] + S22_out[q] + S33_out[q];
        mean_trS += trS;
        const double a = std::abs(trS);
        if (a > max_abs_trS) max_abs_trS = a;
    }
    mean_Smag /= (double)n_real;
    mean_trS  /= (double)n_real;

    std::printf("[sgs_exact] strain-rate stats\n");
    std::printf("[sgs_exact]   min|S|             = %.6e   (must be >= 0)\n", min_Smag);
    std::printf("[sgs_exact]   max|S|             = %.6e\n", max_Smag);
    std::printf("[sgs_exact]   mean|S|            = %.6e   (O(1/eta); eta_K ~ 0.003)\n",
                mean_Smag);
    std::printf("[sgs_exact]   mean(tr S)         = %+.3e  (should be ~0 by incompressibility)\n",
                mean_trS);
    std::printf("[sgs_exact]   max|tr S|          = %.3e   (incompressibility round-off)\n",
                max_abs_trS);
}
