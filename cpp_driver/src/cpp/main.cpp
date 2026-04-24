/*
 * main.cpp — Master driver for the a priori LES SGS model evaluation.
 *
 * Usage
 *   ./les_apriori [iso|channel] [delta_index 0|1|2|all] [--debug]
 *                 [--input <path/to/file.h5>] [--tag <label>]
 *
 * Pipeline (per filter width delta):
 *   1. compute_exact_sgs_stress   -> vel_bar, tau, tau_dev
 *   2. compute_velocity_gradient_tensor(vel_bar) -> grad
 *   3. compute_strain_rate_tensor(grad)          -> S, |S|
 *   4. smagorinsky_model                         -> tau_smag, nu_r_smag
 *   5. wale_model                                -> tau_wale, nu_r_wale
 *   6. dynamic_smagorinsky                       -> tau_dyn,  Cs2_field
 *   7. compute_all_metrics for each of the three models
 *   8. Write tau_exact, tau_smag, tau_wale, tau_dyn and per-point fields
 *      (nu_r, Cs2, Pi, cos_theta) into results/<tag>_delta<N>_<model>.h5
 *
 * Debug checks (--debug, aborts with exit code 1 on first failure):
 *   (a) Energy ratio |ubar|^2 / |u|^2 < 1                         [Pope Eq. 13.22]
 *   (b) Smagorinsky backscatter fraction = 0                      [Pope Eq. 13.127]
 *   (c) Dynamic Smagorinsky Cs^2 < 0 fraction in [0.30, 0.40]     [Germano 1991]
 *   (d) Correlation(tau_exact, tau_exact) = 1                     [definition]
 *   (e) Mean(Pi_exact) > 0                                         [Pope p.589]
 */

#include "config.h"
#include "io.h"
#include "sgs_exact.h"
#include "sgs_models.h"
#include "metrics.h"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// Timing helper
// ---------------------------------------------------------------------------

using Clock      = std::chrono::steady_clock;
using TimePoint  = std::chrono::time_point<Clock>;

static double seconds_since(const TimePoint& t0) {
    const auto now = Clock::now();
    return std::chrono::duration<double>(now - t0).count();
}

// RAII timer that prints the elapsed wall time on destruction.
struct ScopedTimer {
    const char* label;
    TimePoint   t0;
    ScopedTimer(const char* l) : label(l), t0(Clock::now()) {}
    ~ScopedTimer() {
        std::printf("[time] %-38s %7.3f s\n", label, seconds_since(t0));
    }
};

// ---------------------------------------------------------------------------
// Debug-check plumbing
// ---------------------------------------------------------------------------
//
// A failed check is non-recoverable: it means the numerical or physical
// foundation of the study is broken, so the program exits with code 1 and
// prints a reference to the literature value that was expected.
static void fail_check(const char* name,
                       const char* expected,
                       const char* got,
                       const char* reference) {
    std::fprintf(stderr, "\nPHYSICS CHECK FAILED: %s\n",  name);
    std::fprintf(stderr, "Expected: %s\n",                expected);
    std::fprintf(stderr, "Got: %s\n",                     got);
    std::fprintf(stderr, "Reference: %s\n",               reference);
    std::exit(1);
}

// Fraction of a double array strictly less than zero.
static double frac_negative(const std::vector<double>& x) {
    if (x.empty()) return 0.0;
    long long neg = 0;
    for (double v : x) if (v < 0.0) ++neg;
    return static_cast<double>(neg) / static_cast<double>(x.size());
}

// Mean of a double array.
static double mean_vec(const std::vector<double>& x) {
    if (x.empty()) return 0.0;
    long double s = 0.0L;
    for (double v : x) s += v;
    return static_cast<double>(s / static_cast<long double>(x.size()));
}

// Sum of squares of a vector (for energy ratio).
static double sumsq(const std::vector<double>& x) {
    long double s = 0.0L;
    for (double v : x) s += static_cast<long double>(v) * v;
    return static_cast<double>(s);
}

// ---------------------------------------------------------------------------
// Per-model summary row (for the final table)
// ---------------------------------------------------------------------------
struct ModelRow {
    std::string model;
    double      delta         = 0.0;
    double      rho           = 0.0;
    double      backscatter   = 0.0;
    double      mean_cos      = 0.0;
    double      mean_nu       = 0.0;
    double      mean_Pi_model = 0.0;
};

// ---------------------------------------------------------------------------
// CLI parsing
// ---------------------------------------------------------------------------
struct CliOpts {
    std::string dataset    = "iso";   // "iso" or "channel"
    int         delta_idx  = -1;      // -1 = all filters, else 0..2
    bool        debug      = false;
    std::string input_path = "";      // --input: override HDF5 input path
    std::string tag        = "";      // --tag: override output filename prefix
};

static void print_usage_and_exit() {
    std::fprintf(stderr,
        "Usage: ./les_apriori [iso|channel] [0|1|2|all] [--debug]"
        " [--input <path>] [--tag <label>]\n"
        "  dataset     : 'iso' (isotropic_256.h5) or 'channel' (channel_256.h5)\n"
        "  delta index : 0=%.5f  1=%.5f  2=%.5f  or 'all' for every filter\n"
        "  --debug     : enable physics sanity checks, exit 1 on failure\n"
        "  --input     : override input HDF5 path (default: ../data/<dataset>_256.h5)\n"
        "  --tag       : override output filename prefix (default: 'iso' or 'channel')\n",
        FILTER_WIDTHS[0], FILTER_WIDTHS[1], FILTER_WIDTHS[2]);
    std::exit(2);
}

static CliOpts parse_cli(int argc, char* argv[]) {
    CliOpts opts;
    int positional = 0;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--debug") {
            opts.debug = true;
        } else if (arg == "-h" || arg == "--help") {
            print_usage_and_exit();
        } else if (positional == 0) {
            if (arg != "iso" && arg != "channel") {
                std::fprintf(stderr, "error: dataset must be 'iso' or 'channel' (got '%s')\n",
                             arg.c_str());
                print_usage_and_exit();
            }
            opts.dataset = arg;
            ++positional;
        } else if (positional == 1) {
            if (arg == "all") {
                opts.delta_idx = -1;
            } else {
                char* end = nullptr;
                long v = std::strtol(arg.c_str(), &end, 10);
                if (!end || *end != '\0' || v < 0 || v >= N_FILTERS) {
                    std::fprintf(stderr, "error: delta index must be 0..%d or 'all' (got '%s')\n",
                                 N_FILTERS - 1, arg.c_str());
                    print_usage_and_exit();
                }
                opts.delta_idx = static_cast<int>(v);
            }
            ++positional;
        } else if (arg == "--input") {
            if (++i >= argc) {
                std::fprintf(stderr, "error: --input requires a path argument\n");
                print_usage_and_exit();
            }
            opts.input_path = argv[i];
        } else if (arg == "--tag") {
            if (++i >= argc) {
                std::fprintf(stderr, "error: --tag requires a label argument\n");
                print_usage_and_exit();
            }
            opts.tag = argv[i];
        } else {
            std::fprintf(stderr, "error: unexpected argument '%s'\n", arg.c_str());
            print_usage_and_exit();
        }
    }
    return opts;
}

// ---------------------------------------------------------------------------
// Per-filter pipeline
// ---------------------------------------------------------------------------
//
// Runs the full SGS evaluation at a single filter width and populates three
// ModelRow entries (one per model) into `rows`. If `debug` is true, runs the
// extra physics checks and calls fail_check on any violation.
//
static void run_one_filter(const VelocityField&      vf,
                           double                    delta,
                           int                       delta_idx,
                           const std::string&        dataset_tag,
                           bool                      debug,
                           std::vector<ModelRow>&    rows) {
    const int    N  = vf.N;
    const double dx = vf.dx;

    std::printf("\n");
    std::printf("=============================================================\n");
    std::printf("  Filter pass: delta_index=%d  delta=%.6f  (%.1f * dx)\n",
                delta_idx, delta, delta / dx);
    std::printf("=============================================================\n");

    const TimePoint t_pass_start = Clock::now();

    // -----------------------------------------------------------------------
    // 1. Exact SGS stress + filtered velocity
    // -----------------------------------------------------------------------
    TensorField6  tau_exact, tau_exact_dev;
    VelocityField vel_bar;
    {
        ScopedTimer T("compute_exact_sgs_stress");
        compute_exact_sgs_stress(vf, delta, dx, tau_exact, tau_exact_dev, vel_bar);
    }

    // DEBUG (a): energy ratio must be <= 1 (sharp cutoff removes energy).
    // We use u1 alone as a representative check; any component would do.
    if (debug) {
        const double e_in  = sumsq(vf.u1);
        const double e_out = sumsq(vel_bar.u1);
        const double ratio = e_out / (e_in + 1e-300);
        std::printf("[debug] energy ratio E(ubar)/E(u) on u1 = %.6f\n", ratio);
        if (!(ratio < 1.0 + 1e-9)) {
            char got[64]; std::snprintf(got, sizeof got, "ratio = %.6f", ratio);
            fail_check("filter energy ratio",
                       "ratio < 1 (sharp spectral cutoff cannot add energy)",
                       got,
                       "Pope, Turbulent Flows (2000), Eq. 13.22");
        }
    }

    // -----------------------------------------------------------------------
    // 2. Velocity gradient + strain-rate tensor
    // -----------------------------------------------------------------------
    std::vector<std::vector<double>> grad;
    TensorField6                     S;
    std::vector<double>              S_mag;
    {
        ScopedTimer T("compute_velocity_gradient_tensor");
        compute_velocity_gradient_tensor(vel_bar, dx, grad);
    }
    {
        ScopedTimer T("compute_strain_rate_tensor");
        compute_strain_rate_tensor(grad, N, S, S_mag);
    }

    // -----------------------------------------------------------------------
    // 3. Three SGS models
    // -----------------------------------------------------------------------
    TensorField6        tau_smag;
    std::vector<double> nu_r_smag;
    {
        ScopedTimer T("smagorinsky_model");
        smagorinsky_model(S, S_mag, delta, N, tau_smag, nu_r_smag);
    }

    TensorField6        tau_wale;
    std::vector<double> nu_r_wale;
    {
        ScopedTimer T("wale_model");
        wale_model(grad, S, delta, N, tau_wale, nu_r_wale);
    }

    TensorField6        tau_dyn;
    std::vector<double> Cs2_pointwise;   // raw Germano-Lilly pointwise (for PDF diagnostic)
    std::vector<double> Cs2_effective;   // homogeneously averaged (used for tau_dyn)
    {
        ScopedTimer T("dynamic_smagorinsky");
        // Choose homogeneous-direction averaging based on dataset topology:
        //   iso      -> Box       (statistically homogeneous in all 3 directions)
        //   channel  -> XZ_Plane  (homogeneous in x,z only; Cs^2 varies with y)
        const AveragingMode avg_mode =
            (dataset_tag.rfind("channel", 0) == 0) ? AveragingMode::XZ_Plane
                                                   : AveragingMode::Box;
        dynamic_smagorinsky(vel_bar, S, S_mag, delta, dx, N,
                            tau_dyn, Cs2_pointwise, Cs2_effective, avg_mode);
    }

    // -----------------------------------------------------------------------
    // DEBUG (d) and (e): exact-stress self-consistency checks
    // -----------------------------------------------------------------------
    // Compute exact diagnostics once; reuse in checks and per-model HDF5 writes.
    // Alignment-angle definition (SGS stress vs strain-rate):
    //   cos(theta) = (tau'_ij S_ij) / (|tau'| |S|).
    std::vector<double> Pi_exact;
    std::vector<double> cos_theta_exact;
    compute_dissipation_field(tau_exact_dev, S, N, Pi_exact);
    compute_alignment_angle(tau_exact_dev, S, N, cos_theta_exact);

    if (debug) {
        // (d) correlation of exact stress with itself must be 1.
        const double rho_self = correlation_coefficient(tau_exact_dev, tau_exact_dev, N);
        std::printf("[debug] correlation(tau_exact, tau_exact) = %.6f\n", rho_self);
        if (std::fabs(rho_self - 1.0) > 1e-9) {
            char got[64]; std::snprintf(got, sizeof got, "rho = %.6f", rho_self);
            fail_check("self-correlation of exact stress",
                       "rho = 1.000000",
                       got,
                       "Pearson correlation definition (rho(X,X) = 1)");
        }

        // (e) mean Pi_exact must be positive (net forward cascade in 3D).
        const double mean_Pi = mean_vec(Pi_exact);
        std::printf("[debug] mean Pi_exact = %+.6e\n", mean_Pi);
        if (!(mean_Pi > 0.0)) {
            char got[64]; std::snprintf(got, sizeof got, "mean_Pi = %+.6e", mean_Pi);
            fail_check("mean exact SGS dissipation",
                       "mean_Pi > 0 (net forward cascade in 3D turbulence)",
                       got,
                       "Pope, Turbulent Flows (2000), p.589 and Eq. 13.123");
        }
    }

    // -----------------------------------------------------------------------
    // 4. Per-model metrics + HDF5 output
    // -----------------------------------------------------------------------
    // Results file name: results/<tag>_delta<idx>_<model>.h5
    auto results_path = [&](const char* model) {
        char buf[256];
        std::snprintf(buf, sizeof buf, "../results/%s_delta%d_%s.h5",
                      dataset_tag.c_str(), delta_idx, model);
        return std::string(buf);
    };

    // We need Pi and cos_theta per model for the scalar outputs and summary.
    auto process_model = [&](const std::string&   model_name,
                             const TensorField6&  tau_mod,
                             std::vector<double>& nu_field) {
        const std::string fpath = results_path(model_name.c_str());

        // Run the aggregate metrics driver (writes scalars + prints summary).
        {
            ScopedTimer T(("compute_all_metrics " + model_name).c_str());
            compute_all_metrics(tau_exact_dev, tau_mod, S, nu_field,
                                N, delta, model_name, fpath);
        }

        // Per-point scalar fields for plotting.
        std::vector<double> Pi_mod, cos_theta_mod;
        compute_dissipation_field(tau_mod, S, N, Pi_mod);
        compute_alignment_angle(tau_mod, S, N, cos_theta_mod);

        write_tensor_hdf5(tau_mod,       fpath, "tau_" + model_name);
        write_tensor_hdf5(tau_exact_dev, fpath, "tau_exact_dev");

        // Eddy-viscosity field: nu_r for smag/wale, Cs^2 for dynamic.
        ScalarField nu_s{nu_field, N, "nu_r"};
        write_scalar_hdf5(nu_s, fpath, (model_name == "dynamic") ? "Cs2" : "nu_r");

        ScalarField Pi_s{Pi_mod, N, "Pi"};
        write_scalar_hdf5(Pi_s, fpath, "Pi");

        ScalarField cs_s{cos_theta_mod, N, "cos_theta"};
        write_scalar_hdf5(cs_s, fpath, "cos_theta");

        // Also save the DNS-exact dissipation alongside the model for easy
        // comparison in downstream plots.
        ScalarField pie_s{Pi_exact, N, "Pi_exact"};
        write_scalar_hdf5(pie_s, fpath, "Pi_exact");
        ScalarField cte{cos_theta_exact, N, "cos_theta_exact"};
        write_scalar_hdf5(cte, fpath, "cos_theta_exact");

        // Collect row for the final cross-filter summary table.
        ModelRow row;
        row.model         = model_name;
        row.delta         = delta;
        row.rho           = correlation_coefficient(tau_mod, tau_exact_dev, N);
        row.backscatter   = frac_negative(Pi_mod);
        row.mean_cos      = mean_vec(cos_theta_mod);
        row.mean_nu       = mean_vec(nu_field);
        row.mean_Pi_model = mean_vec(Pi_mod);
        rows.push_back(row);

        return Pi_mod;
    };

    std::vector<double> Pi_smag = process_model("smagorinsky", tau_smag, nu_r_smag);
    std::vector<double> Pi_wale = process_model("wale",        tau_wale, nu_r_wale);
    // For dynamic model: pass the POINTWISE field as nu_field so:
    //   - /Cs2 in the HDF5 (written by process_model) contains the raw
    //     Germano-Lilly pointwise distribution (preserves the Cs^2 PDF plot)
    //   - mean_nu_r attribute reflects the pointwise mean (matches existing CSV)
    // The physically meaningful EFFECTIVE field is additionally written as
    // /Cs2_eff immediately below, alongside its mean as an HDF5 attribute.
    std::vector<double> Pi_dyn  = process_model("dynamic",     tau_dyn,  Cs2_pointwise);

    // Additional write: the homogeneously-averaged Cs^2 field (what tau_dyn
    // was actually built from). Exposed so downstream analysis can inspect
    // the scalar (Box) or y-profile (XZ_Plane) Cs^2 used by the model.
    {
        const std::string dyn_path = results_path("dynamic");
        ScalarField cs_eff{Cs2_effective, N, "Cs2_eff"};
        write_scalar_hdf5(cs_eff, dyn_path, "Cs2_eff");
    }

    // -----------------------------------------------------------------------
    // DEBUG (b) and (c): model-specific structural checks
    // -----------------------------------------------------------------------
    if (debug) {
        // (b) Smagorinsky backscatter fraction must be exactly zero.
        // Note: backscatter_fraction already asserts this and aborts, but we
        // run an explicit check here so the failure message uses the standard
        // fail_check format with the literature reference.
        const double bs_smag = frac_negative(Pi_smag);
        std::printf("[debug] Smagorinsky backscatter fraction = %.3e\n", bs_smag);
        if (bs_smag > 1e-12) {
            char got[64]; std::snprintf(got, sizeof got, "fraction = %.6e", bs_smag);
            fail_check("Smagorinsky zero-backscatter",
                       "fraction = 0 (nu_r >= 0 and Pi = 2*nu_r*|S|^2 >= 0 everywhere)",
                       got,
                       "Pope, Turbulent Flows (2000), Eq. 13.127");
        }

        // (c) Dynamic Smagorinsky: fraction of grid points with Cs^2 < 0
        // in the POINTWISE field should fall in [0.30, 0.40] for the
        // Germano-Lilly procedure applied to 3D isotropic turbulence.
        // (The effective field after Box/XZ_Plane averaging is spatially
        // structured, so this check is against the raw pointwise field.)
        const double neg_Cs2 = frac_negative(Cs2_pointwise);
        std::printf("[debug] Dynamic Smag  Cs^2 < 0 fraction   = %.3f\n", neg_Cs2);
        if (neg_Cs2 < 0.30 || neg_Cs2 > 0.40) {
            char got[64]; std::snprintf(got, sizeof got, "fraction = %.3f", neg_Cs2);
            fail_check("Dynamic Smagorinsky backscatter range",
                       "fraction of Cs^2 < 0 in [0.30, 0.40]",
                       got,
                       "Germano et al., Phys. Fluids A 3(7) 1991, Eq. 8-10");
        }
    }

    std::printf("[time] %-38s %7.3f s  (total for this filter pass)\n",
                "PIPELINE TOTAL", seconds_since(t_pass_start));
}

// ---------------------------------------------------------------------------
// Final summary table across all filter widths
// ---------------------------------------------------------------------------
static void print_summary_table(const std::vector<ModelRow>& rows) {
    if (rows.empty()) return;

    std::printf("\n");
    std::printf("=============================================================\n");
    std::printf("  Final summary: a priori metrics across all filter widths\n");
    std::printf("=============================================================\n");
    std::printf("%-14s %10s %8s %10s %12s %12s\n",
                "model", "delta", "rho", "bs_frac", "mean_cos", "mean_nu_r");
    std::printf("%-14s %10s %8s %10s %12s %12s\n",
                "-----", "-----", "---", "-------", "--------", "---------");
    for (const ModelRow& r : rows) {
        std::printf("%-14s %10.5f %8.3f %9.1f%% %12.3f %+12.3e\n",
                    r.model.c_str(),
                    r.delta,
                    r.rho,
                    r.backscatter * 100.0,
                    r.mean_cos,
                    r.mean_nu);
    }
    std::printf("\n");
    std::printf("Reference expectations (isotropic, Re_lambda=433):\n");
    std::printf("  Smagorinsky   rho ~ 0.2-0.4   bs_frac = 0.0%%\n");
    std::printf("  WALE          rho ~ 0.2-0.4   bs_frac = 0.0%%\n");
    std::printf("  Dynamic Smag  rho ~ 0.4-0.6   bs_frac may be non-zero\n");
    std::printf("  Exact field   mean_cos_theta ~ -0.3 to -0.5 (Liu et al. 1994)\n");
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------
int main(int argc, char* argv[]) {
    const CliOpts opts = parse_cli(argc, argv);

    std::printf("les_apriori: a priori LES SGS model evaluation\n");
    std::printf("  Grid size        : %d^3\n",           N);
    std::printf("  DNS grid spacing : %.6f\n",           DX);
    std::printf("  Filter widths    : %.5f  %.5f  %.5f\n",
                FILTER_WIDTHS[0], FILTER_WIDTHS[1], FILTER_WIDTHS[2]);
    std::printf("  Cs (Smagorinsky) : %.2f\n",           CS_SMAGORINSKY);
    std::printf("  Cw (WALE)        : %.2f\n",           CW_WALE);
    std::printf("  Dataset          : %s\n",             opts.dataset.c_str());
    if (opts.delta_idx < 0)
        std::printf("  Delta index      : all (0, 1, 2)\n");
    else
        std::printf("  Delta index      : %d\n",         opts.delta_idx);
    std::printf("  Debug checks     : %s\n",             opts.debug ? "ON" : "off");
    std::printf("\n");

    // -----------------------------------------------------------------------
    // Step 1: read velocity field
    // -----------------------------------------------------------------------
    const std::string data_path =
        opts.input_path.empty()
        ? ((opts.dataset == "iso") ? "../data/isotropic_256.h5"
                                   : "../data/channel_256.h5")
        : opts.input_path;
    const std::string dataset_tag =
        opts.tag.empty()
        ? ((opts.dataset == "iso") ? "iso" : "channel")
        : opts.tag;

    std::printf("  Input file       : %s\n",             data_path.c_str());
    std::printf("  Output tag       : %s\n",             dataset_tag.c_str());
    std::printf("\n");

    const TimePoint t_program_start = Clock::now();

    VelocityField vf;
    {
        ScopedTimer T("read_velocity_hdf5");
        vf = read_velocity_hdf5(data_path);
    }

    if (vf.N != N) {
        std::fprintf(stderr,
            "[FATAL] file N=%d differs from compile-time N=%d. "
            "Rebuild with matching config.h.\n", vf.N, N);
        return 1;
    }

    // -----------------------------------------------------------------------
    // Step 2: field statistics sanity check (already printed by read_velocity_hdf5,
    //         add per-component min/max here as the spec says "sanity check").
    // -----------------------------------------------------------------------
    {
        const long long total = (long long)vf.N * vf.N * vf.N;
        double min_u[3] = { vf.u1[0], vf.u2[0], vf.u3[0] };
        double max_u[3] = { vf.u1[0], vf.u2[0], vf.u3[0] };
        const std::vector<double>* vecs[3] = { &vf.u1, &vf.u2, &vf.u3 };
        for (int c = 0; c < 3; ++c)
            for (double v : *vecs[c]) {
                if (v < min_u[c]) min_u[c] = v;
                if (v > max_u[c]) max_u[c] = v;
            }
        std::printf("[io] per-component range (sanity check):\n");
        std::printf("[io]   u1  min=%+10.4e  max=%+10.4e\n", min_u[0], max_u[0]);
        std::printf("[io]   u2  min=%+10.4e  max=%+10.4e\n", min_u[1], max_u[1]);
        std::printf("[io]   u3  min=%+10.4e  max=%+10.4e\n", min_u[2], max_u[2]);
        std::printf("[io]   total points = %lld\n", total);
    }

    // -----------------------------------------------------------------------
    // Step 3: run pipeline at one or all filter widths
    // -----------------------------------------------------------------------
    std::vector<ModelRow> rows;
    if (opts.delta_idx < 0) {
        for (int k = 0; k < N_FILTERS; ++k)
            run_one_filter(vf, FILTER_WIDTHS[k], k, dataset_tag, opts.debug, rows);
    } else {
        run_one_filter(vf, FILTER_WIDTHS[opts.delta_idx], opts.delta_idx,
                       dataset_tag, opts.debug, rows);
    }

    // -----------------------------------------------------------------------
    // Step 4: final summary
    // -----------------------------------------------------------------------
    print_summary_table(rows);

    std::printf("[time] %-38s %7.3f s  (wall clock, whole program)\n",
                "TOTAL", seconds_since(t_program_start));

    return 0;
}
