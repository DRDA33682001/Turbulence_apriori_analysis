/*
 * filter.cpp — Sharp spectral cutoff filter, FFTW-backed.
 *
 * Physics: Pope Ch. 13, Eq. 13.22. The sharp spectral filter is the
 * projection onto the subspace of Fourier modes with |k| <= kc = pi/delta.
 * In Fourier space it is a unit indicator function, so applying it can only
 * remove energy, never add it — we assert this after every call.
 *
 * Parseval's theorem guarantees that computing the "energy" sum_i u_i^2 in
 * physical space is equivalent to summing |c_k|^2 in spectral space, up to
 * a fixed normalization. We therefore measure energy in physical space,
 * which is cheaper and avoids any ambiguity about the r2c storage weights.
 *
 * FFTW wavenumber layout (r2c of shape N x N x N → N x N x (N/2+1) complex):
 *   axis 0 (i=0..N-1):   n_x = (i <= N/2) ? i : i - N
 *   axis 1 (j=0..N-1):   n_y = (j <= N/2) ? j : j - N
 *   axis 2 (k=0..N/2):   n_z = k               [r2c: non-negative only]
 * Physical wavenumber on each axis: k_phys = 2*pi * n / L, with L = N*dx.
 *
 * FFTW reference:
 *   https://www.fftw.org/fftw3_doc/Real_002ddata-DFTs.html
 *   https://www.fftw.org/fftw3_doc/Multi_002dDimensional-DFTs-of-Real-Data.html
 *   https://www.fftw.org/fftw3_doc/Wisdom.html
 *
 * Cross-check against references/spectralDNS: their solvers build K[0], K[1]
 * via fftfreq and K[2] via rfftfreq (multiplied by 2*pi/L), then zero out
 * coefficients whose |K| exceeds the cutoff. This file follows the same
 * convention exactly.
 */

#include "filter.h"
#include "config.h"

#include <fftw3.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

// ---------------------------------------------------------------------------
// Plan cache
// ---------------------------------------------------------------------------
//
// FFTW plan creation with FFTW_MEASURE takes noticeable time (tens of ms to
// seconds at N=256). We cache plans + scratch buffers in function-local
// static state, keyed by N, so that the first call in a process pays the
// measurement cost and every subsequent call is a pure execute.
//
// Persistence across program runs is handled via fftw_import/export_wisdom:
// with wisdom loaded, even the first plan build in a new process is
// near-instant.
//
// Lifetime note: these resources are intentionally leaked at program exit.
// FFTW does not require plan_destroy for correct results, and wiring an
// atexit hook would add noise for no observable benefit.

namespace {

struct FFTFilterContext {
    int            N        = 0;
    double*        real_buf = nullptr;   // N^3 doubles
    fftw_complex*  freq_buf = nullptr;   // N * N * (N/2+1) complex
    fftw_plan      plan_fwd = nullptr;
    fftw_plan      plan_bwd = nullptr;
};

constexpr const char* WISDOM_PATH = "../results/fftw_wisdom.dat";

FFTFilterContext& get_context(int N) {
    static FFTFilterContext ctx;
    if (ctx.N == N && ctx.plan_fwd && ctx.plan_bwd)
        return ctx;

    // A different N was requested, or first use. Tear down any previous
    // state before rebuilding.
    if (ctx.plan_fwd) { fftw_destroy_plan(ctx.plan_fwd); ctx.plan_fwd = nullptr; }
    if (ctx.plan_bwd) { fftw_destroy_plan(ctx.plan_bwd); ctx.plan_bwd = nullptr; }
    if (ctx.real_buf) { fftw_free(ctx.real_buf); ctx.real_buf = nullptr; }
    if (ctx.freq_buf) { fftw_free(ctx.freq_buf); ctx.freq_buf = nullptr; }

    // Try to load cached wisdom (may be absent on first-ever run — harmless).
    const int got_wisdom = fftw_import_wisdom_from_filename(WISDOM_PATH);
    std::printf("[filter] fftw wisdom %s '%s'\n",
                got_wisdom ? "loaded from" : "not found at", WISDOM_PATH);

    ctx.N = N;
    const size_t n_real = (size_t)N * N * N;
    const size_t n_comp = (size_t)N * N * (N / 2 + 1);

    ctx.real_buf = (double*)       fftw_malloc(sizeof(double)       * n_real);
    ctx.freq_buf = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n_comp);
    if (!ctx.real_buf || !ctx.freq_buf) {
        std::fprintf(stderr, "[filter] FATAL: fftw_malloc failed for N=%d\n", N);
        std::abort();
    }

    // FFTW_MEASURE is allowed to clobber the arrays during planning, which is
    // fine — the buffers are uninitialised at this point. We plan both
    // directions before feeding data in.
    ctx.plan_fwd = fftw_plan_dft_r2c_3d(N, N, N,
                                        ctx.real_buf, ctx.freq_buf,
                                        FFTW_MEASURE);
    ctx.plan_bwd = fftw_plan_dft_c2r_3d(N, N, N,
                                        ctx.freq_buf, ctx.real_buf,
                                        FFTW_MEASURE);
    if (!ctx.plan_fwd || !ctx.plan_bwd) {
        std::fprintf(stderr, "[filter] FATAL: fftw_plan creation failed for N=%d\n", N);
        std::abort();
    }

    // Persist wisdom so the next run starts fast.
    const int saved = fftw_export_wisdom_to_filename(WISDOM_PATH);
    if (!saved) {
        std::fprintf(stderr, "[filter] WARN: could not save fftw wisdom to '%s'\n",
                     WISDOM_PATH);
    } else {
        std::printf("[filter] fftw wisdom saved to '%s'\n", WISDOM_PATH);
    }

    std::printf("[filter] plan ready: N=%d  r2c+c2r  FFTW_MEASURE\n", N);
    return ctx;
}

// Signed integer wavenumber for FFTW's full-axis layout.
inline int signed_wavenumber(int i, int N) {
    return (i <= N / 2) ? i : i - N;
}

} // anonymous namespace

// ---------------------------------------------------------------------------
// spectral_filter_3d
// ---------------------------------------------------------------------------

void spectral_filter_3d(const double* u_in,
                        double*       u_out,
                        int           N,
                        double        delta,
                        double        dx)
{
    FFTFilterContext& ctx = get_context(N);

    const size_t n_real = (size_t)N * N * N;
    const int    Nzh    = N / 2 + 1;

    // Energy before (physical space — Parseval-equivalent to spectral).
    double energy_in = 0.0;
    for (size_t p = 0; p < n_real; ++p)
        energy_in += u_in[p] * u_in[p];

    // Load input into FFTW buffer.
    std::memcpy(ctx.real_buf, u_in, n_real * sizeof(double));

    // Forward r2c: real_buf -> freq_buf.
    fftw_execute(ctx.plan_fwd);

    // Apply sharp cutoff in physical wavenumbers.
    const double L           = static_cast<double>(N) * dx;
    const double two_pi_over_L = 2.0 * M_PI / L;
    const double kc          = M_PI / delta;
    const double kc2         = kc * kc;

    long long n_killed = 0;
    long long n_total_modes = (long long)N * N * Nzh;

    for (int i = 0; i < N; ++i) {
        const double kx = two_pi_over_L * signed_wavenumber(i, N);
        const double kx2 = kx * kx;
        for (int j = 0; j < N; ++j) {
            const double ky = two_pi_over_L * signed_wavenumber(j, N);
            const double ky2 = ky * ky;
            // Row base index into freq_buf for this (i, j) slab.
            const size_t row = ((size_t)i * N + j) * Nzh;
            for (int k = 0; k < Nzh; ++k) {
                // r2c third axis: only non-negative wavenumbers stored.
                const double kz = two_pi_over_L * k;
                const double k2 = kx2 + ky2 + kz * kz;
                if (k2 > kc2) {
                    const size_t idx_c = row + (size_t)k;
                    ctx.freq_buf[idx_c][0] = 0.0;
                    ctx.freq_buf[idx_c][1] = 0.0;
                    ++n_killed;
                }
            }
        }
    }

    // Backward c2r: freq_buf -> real_buf. (c2r destroys freq_buf, which is
    // fine — we don't reuse it outside this call.)
    fftw_execute(ctx.plan_bwd);

    // FFTW transforms are unnormalised: fwd.bwd = N^3 * identity.
    const double norm = 1.0 / static_cast<double>(n_real);
    for (size_t p = 0; p < n_real; ++p)
        u_out[p] = ctx.real_buf[p] * norm;

    // Energy after.
    double energy_out = 0.0;
    for (size_t p = 0; p < n_real; ++p)
        energy_out += u_out[p] * u_out[p];

    const double ratio = energy_out / (energy_in + 1e-300);
    const double kept_frac = 1.0 - (double)n_killed / (double)n_total_modes;

    std::printf("[filter] spectral_filter_3d  N=%d  delta=%.6f  kc=%.4f\n",
                N, delta, kc);
    std::printf("[filter]   modes kept : %lld / %lld  (%.2f%%)\n",
                n_total_modes - n_killed, n_total_modes, 100.0 * kept_frac);
    std::printf("[filter]   energy     : in=%.6e  out=%.6e  ratio=%.4f\n",
                energy_in, energy_out, ratio);

    // Sharp cutoff is a projection — it cannot increase the L2 norm of the
    // field. If it does, something is deeply wrong (wavenumber layout,
    // normalisation, aliasing). Fail loudly so the bug cannot propagate.
    //
    // Tolerance: allow a few ULPs of round-off above 1.0.
    if (ratio > 1.0 + 1e-10) {
        std::fprintf(stderr,
            "\n[filter] *** FATAL *** filtered energy exceeds input:\n"
            "          energy_in  = %.12e\n"
            "          energy_out = %.12e\n"
            "          ratio      = %.12e  (must be <= 1)\n"
            "        Sharp spectral cutoff is a projection; this means the\n"
            "        implementation is buggy (wavenumber layout, plan type,\n"
            "        normalisation, or aliasing). Aborting.\n\n",
            energy_in, energy_out, ratio);
        std::abort();
    }
}

// ---------------------------------------------------------------------------
// filter_product
// ---------------------------------------------------------------------------
//
// tau_ij in a priori tests requires filter(ui*uj). This is NOT the same as
// filter(ui) * filter(uj) — the difference is exactly the Leonard/SGS stress.
// We therefore form the pointwise product on the DNS grid first, then feed
// it through the sharp-cutoff filter.

void filter_product(const double* ui,
                    const double* uj,
                    double*       result,
                    int           N,
                    double        delta,
                    double        dx)
{
    const size_t n_real = (size_t)N * N * N;

    // Scratch for the pointwise product. A heap allocation is acceptable:
    // this is called O(6 * n_filters) times per run, not per-grid-point.
    std::vector<double> product(n_real);
    for (size_t p = 0; p < n_real; ++p)
        product[p] = ui[p] * uj[p];

    spectral_filter_3d(product.data(), result, N, delta, dx);
}

// ============================================================================
// TOP-HAT (BOX) FILTER — classical test filter for the dynamic procedure
// ============================================================================
//
// A 3-D periodic top-hat filter that averages each point over a cube of
// (2r+1)^3 neighbours. Implemented as three separable 1-D box filters, each
// using a periodic prefix-sum so the cost is O(N^3) independent of r.
//
// Reference
//   Classical top-hat / box filter described in Germano (1991) and Lilly
//   (1992). Meneveau & Katz (2000) ARFM sec. 3.2 recommends it over the
//   sharp spectral filter as a TEST filter for the dynamic procedure.

namespace {

// Apply a periodic 1-D box filter along one axis of a 3-D row-major array
// of size N^3. The axis is selected by `stride` and `n_outer`. For each
// 1-D line of length N, output[i] = mean(input[i-r .. i+r]) with periodic
// wrap.  Cost: O(N) per line, O(N^3) total.
//
//   src, dst: flat arrays of length N^3 (must NOT alias each other)
//   N       : grid size per axis
//   r       : filter half-width in grid points
//   stride  : stride between successive 1-D samples within a line
//             (= N*N for axis 0, N for axis 1, 1 for axis 2)
//   n_lines : number of 1-D lines = N^3 / N = N^2 (used to iterate)
static void box_filter_1d_axis(const double* src,
                               double*       dst,
                               int           N,
                               int           r,
                               ptrdiff_t     stride)
{
    const int window = 2 * r + 1;
    const double inv_window = 1.0 / static_cast<double>(window);

    // Number of 1-D lines to filter, equal to N^2.
    const ptrdiff_t n_lines = static_cast<ptrdiff_t>(N) * N;

    // Enumerate starting offsets of each 1-D line. A line is a set of N
    // elements separated by `stride`. For a row-major (i,j,k) array:
    //   - axis 2 (stride=1):    lines start at (i*N + j) * N, one per (i,j)
    //   - axis 1 (stride=N):    lines start at i*N^2 + k, one per (i,k)
    //   - axis 0 (stride=N^2):  lines start at j*N + k, one per (j,k)
    // For each case, the line-start offsets cover [0, N^2) packed in a
    // predictable way.  We generate them by stepping through (a, b) with
    // the axis whose stride is NOT selected.

    // Because the three cases differ only in which two axes are "outer",
    // we handle them uniformly by picking the outer axes from the stride.
    // Small helper: returns the axis (0,1,2) corresponding to the stride.
    auto axis_from_stride = [N](ptrdiff_t s) -> int {
        if (s == 1)           return 2;
        if (s == N)           return 1;
        if (s == (ptrdiff_t)N * N) return 0;
        return -1;
    };
    const int filter_axis = axis_from_stride(stride);
    (void)filter_axis;  // not used beyond documenting intent

    // For each 1-D line:
    //   - build a periodic prefix sum of length N+1 (cheap)
    //   - compute output[i] = (prefix[i+r+1] - prefix[i-r]) / window
    //     with indices taken modulo N via pre-extended indexing.
    // To avoid the modulo branch inside the inner loop, we extend the
    // line by replicating its front and back on a scratch buffer of size
    // N + 2*r.  Cost per line: O(N + r) with a tight inner loop.

    std::vector<double> ext(static_cast<size_t>(N) + 2 * r);

    // Generate line starts in canonical order.  We iterate over the two
    // axes that AREN'T `filter_axis` in lexicographic order (outer axis
    // varies slowest).  For each (a, b) pair we compute the base offset
    // as a * outer_stride_a + b * outer_stride_b.
    ptrdiff_t outer_stride_a = 0, outer_stride_b = 0, N_a = 0, N_b = 0;
    if (filter_axis == 0) {            // lines along axis 0, outer = (1,2)
        outer_stride_a = N;             // axis 1
        outer_stride_b = 1;             // axis 2
        N_a = N;
        N_b = N;
    } else if (filter_axis == 1) {     // lines along axis 1, outer = (0,2)
        outer_stride_a = (ptrdiff_t)N * N;
        outer_stride_b = 1;
        N_a = N;
        N_b = N;
    } else {                            // lines along axis 2, outer = (0,1)
        outer_stride_a = (ptrdiff_t)N * N;
        outer_stride_b = N;
        N_a = N;
        N_b = N;
    }

    for (ptrdiff_t a = 0; a < N_a; ++a) {
        for (ptrdiff_t b = 0; b < N_b; ++b) {
            const ptrdiff_t base = a * outer_stride_a + b * outer_stride_b;

            // Load the 1-D line from the sparse src array into `ext`
            // extended with periodic wrap on both ends.
            for (int i = 0; i < N; ++i)
                ext[static_cast<size_t>(i + r)] = src[base + (ptrdiff_t)i * stride];
            for (int i = 0; i < r; ++i) {
                ext[static_cast<size_t>(i)]         = src[base + (ptrdiff_t)(N - r + i) * stride];  // front wrap
                ext[static_cast<size_t>(N + r + i)] = src[base + (ptrdiff_t)i * stride];            // back wrap
            }

            // Rolling sum over window = 2r+1 points of `ext`.
            double running = 0.0;
            for (int k = 0; k < window; ++k) running += ext[static_cast<size_t>(k)];
            dst[base + 0 * stride] = running * inv_window;
            for (int i = 1; i < N; ++i) {
                running += ext[static_cast<size_t>(i + 2 * r)];
                running -= ext[static_cast<size_t>(i - 1)];
                dst[base + (ptrdiff_t)i * stride] = running * inv_window;
            }
        }
    }
}

} // anonymous namespace

void box_filter_3d(const double* u_in,
                   double*       u_out,
                   int           N,
                   int           r)
{
    // Clamp the radius so the window cannot exceed the periodic domain.
    if (r < 0) r = 0;
    if (r > N / 2) r = N / 2;
    if (r == 0) {
        // Trivial copy.
        std::copy_n(u_in, (size_t)N * N * N, u_out);
        return;
    }

    const size_t n = (size_t)N * N * N;
    std::vector<double> scratch(n);

    // Apply in sequence along axis 2 (innermost), axis 1, axis 0.
    // Source and destination alternate between scratch and u_out so the
    // final result lands in u_out.
    box_filter_1d_axis(u_in,           scratch.data(), N, r, /*stride=*/1);
    box_filter_1d_axis(scratch.data(), u_out,          N, r, /*stride=*/N);
    box_filter_1d_axis(u_out,          scratch.data(), N, r, /*stride=*/(ptrdiff_t)N * N);
    std::copy(scratch.begin(), scratch.end(), u_out);
}

void box_filter_product_3d(const double* ui,
                           const double* uj,
                           double*       result,
                           int           N,
                           int           r)
{
    const size_t n = (size_t)N * N * N;
    std::vector<double> product(n);
    for (size_t p = 0; p < n; ++p)
        product[p] = ui[p] * uj[p];
    box_filter_3d(product.data(), result, N, r);
}
