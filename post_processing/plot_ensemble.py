#!/usr/bin/env python3
"""
plot_ensemble.py — Ensemble-level plots, tables, and stability report.

Consolidates what used to live in three separate scripts:
  * ensemble_visualize.py       (5 publication figures in figures/ensemble/)
  * ensemble_stats.py           (ensemble summary CSV + AIAA LaTeX table)
  * advanced_plots_ensemble.py  (per-snapshot correlation convergence table
                                 and stability_report.csv)

The ensemble here is 8 snapshots × 3 filter widths × 3 models = 72 HDF5
files produced by cpp_driver/run_full_ensemble.sh.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Optional

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
SNAPSHOTS     = ["t1", "t2", "t3", "t4"]
DATASETS      = ["iso", "channel"]
MODELS        = ["smagorinsky", "wale", "dynamic"]
DELTA_INDICES = [0, 1, 2]
FILTER_RATIOS = [4, 8, 16]
STABILITY_THRESHOLD = 0.01

MODEL_COLORS = {
    "smagorinsky": "#2166ac",
    "wale":        "#d6604d",
    "dynamic":     "#1a9641",
}
MODEL_LABELS = {
    "smagorinsky": "Smagorinsky",
    "wale":        "WALE",
    "dynamic":     "Dynamic Smag.",
}

sns.set_style("whitegrid")
sns.set_context("paper", font_scale=1.35)
plt.rcParams.update({
    "font.family":       "DejaVu Sans",
    "lines.linewidth":   1.8,
    "axes.grid":         True,
    "grid.alpha":        0.30,
    "axes.edgecolor":    "#444444",
    "axes.linewidth":    0.9,
    "legend.framealpha": 0.92,
})


# ---------------------------------------------------------------------------
# IO helpers
# ---------------------------------------------------------------------------

def _rfile(results: Path, ds: str, t: str, di: int, model: str) -> Path:
    return results / f"{ds}_{t}_delta{di}_{model}.h5"


def _load_attr(h5_path: Path, model: str, attr: str) -> Optional[float]:
    if not h5_path.exists():
        return None
    try:
        with h5py.File(h5_path, "r") as f:
            return float(f[model].attrs[attr])
    except Exception as exc:
        print(f"  [warn] {h5_path.name}: {exc}")
        return None


def _load_flat(h5_path: Path, dset: str) -> Optional[np.ndarray]:
    if not h5_path.exists():
        return None
    try:
        with h5py.File(h5_path, "r") as f:
            return f[dset][()].astype(np.float64).ravel()
    except Exception as exc:
        print(f"  [warn] {h5_path.name}: {exc}")
        return None


def _histogram(arr: np.ndarray, bins: np.ndarray) -> np.ndarray:
    counts, _ = np.histogram(arr, bins=bins, density=True)
    return counts


def _collect_snapshots(results: Path, dataset: str, di: int,
                       model: str, attr: str) -> List[float]:
    vals: List[float] = []
    for t in SNAPSHOTS:
        v = _load_attr(_rfile(results, dataset, t, di, model), model, attr)
        if v is not None:
            vals.append(v)
    return vals


# ===========================================================================
# Figure 1 — Correlation vs Filter Width (per dataset, with snapshot scatter)
# ===========================================================================
_FIG1_YLIM = (-0.17, 0.34)


def figure1_correlation(results: Path, figures: Path, dataset: str) -> None:
    fig, ax = plt.subplots(figsize=(7.5, 5.3))
    rng = np.random.RandomState(42)

    for model in MODELS:
        means, stds = [], []
        for di in DELTA_INDICES:
            vals = _collect_snapshots(results, dataset, di, model, "correlation")
            means.append(float(np.mean(vals)) if vals else np.nan)
            stds.append(float(np.std(vals, ddof=0)) if vals else np.nan)
            if vals:
                x_jit = FILTER_RATIOS[di] * (1 + rng.uniform(-0.04, 0.04, len(vals)))
                ax.scatter(x_jit, vals, s=26, alpha=0.55,
                           color=MODEL_COLORS[model], zorder=3,
                           edgecolor="white", linewidth=0.5)

        ax.errorbar(FILTER_RATIOS, means, yerr=stds,
                    fmt="o-", color=MODEL_COLORS[model], label=MODEL_LABELS[model],
                    capsize=5, lw=2.2, markersize=8.5,
                    markeredgecolor="white", markeredgewidth=1.0, zorder=5)

    ax.axhline(0.30, color="#888888", linestyle=(0, (4, 4)),
               lw=1.1, alpha=0.55, zorder=1)
    ax.text(4, 0.305,
            r"McMillan & Ferziger (1979): $\rho \approx 0.30$ at $\Delta = 2\Delta x$",
            ha="left", va="bottom", fontsize=9.5, color="#555555", style="italic")
    ax.axhline(0.0, color="k", lw=0.6, alpha=0.45, zorder=1)

    ax.set_xscale("log", base=2)
    ax.set_xticks(FILTER_RATIOS)
    ax.set_xticklabels([str(r) for r in FILTER_RATIOS])
    ax.set_xlim(3.5, 18.5)
    ax.set_ylim(*_FIG1_YLIM)
    ax.set_xlabel(r"Filter ratio $\Delta/\Delta x$")
    ax.set_ylabel(r"Correlation coefficient $\rho$")
    ax.legend(loc="upper right", ncol=1)

    out = figures / f"ensemble_correlation_{dataset}.png"
    fig.tight_layout()
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig1] {out.name}")


# ===========================================================================
# Figure 2 — Backscatter grouped bar chart
# ===========================================================================

def figure2_backscatter(results: Path, figures: Path, dataset: str) -> None:
    fig, ax = plt.subplots(figsize=(8.2, 5.3))
    rng = np.random.RandomState(42)

    bar_w = 0.26
    offsets = {"smagorinsky": -bar_w, "wale": 0.0, "dynamic": bar_w}
    x = np.arange(len(FILTER_RATIOS))

    ax.axhspan(30, 40, color="#ffd966", alpha=0.22, zorder=1,
               label="DNS a priori: ~30–40%")

    for model in MODELS:
        means, stds, per_di = [], [], []
        for di in DELTA_INDICES:
            vals = [100.0 * v for v in _collect_snapshots(
                results, dataset, di, model, "backscatter_frac")]
            per_di.append(vals)
            means.append(float(np.mean(vals)) if vals else 0.0)
            stds.append(float(np.std(vals, ddof=0)) if vals else 0.0)

        ax.bar(x + offsets[model], means, bar_w,
               yerr=stds, color=MODEL_COLORS[model], label=MODEL_LABELS[model],
               alpha=0.88, capsize=4, edgecolor="black", linewidth=0.8, zorder=3,
               error_kw={"lw": 1.4, "zorder": 4})

        for di, vals in zip(DELTA_INDICES, per_di):
            if vals:
                x_jit = (x[di] + offsets[model]) + rng.uniform(-0.04, 0.04, len(vals))
                ax.scatter(x_jit, vals, color="black", s=14, zorder=5, alpha=0.80,
                           edgecolor="white", linewidth=0.5)

    ax.set_xticks(x)
    ax.set_xticklabels([str(r) for r in FILTER_RATIOS])
    ax.set_xlabel(r"Filter ratio $\Delta/\Delta x$")
    ax.set_ylabel("Backscatter fraction (%)")
    ax.set_ylim(0, 60)
    ax.axhline(0, color="k", lw=0.7)
    ax.legend(loc="upper left", ncol=1, framealpha=0.95)

    out = figures / f"ensemble_backscatter_{dataset}.png"
    fig.tight_layout()
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig2] {out.name}")


# ===========================================================================
# Figure 3 — SGS dissipation PDF (min-max band + mean line, 4 sources)
# ===========================================================================

def figure3_dissipation_pdf(results: Path, figures: Path) -> None:
    n_bins = 250

    lo_list, hi_list = [], []
    for t in SNAPSHOTS:
        for src, dset in [("smagorinsky", "Pi_exact"),
                          ("smagorinsky", "Pi"),
                          ("wale",        "Pi"),
                          ("dynamic",     "Pi")]:
            arr = _load_flat(_rfile(results, "iso", t, 0, src), dset)
            if arr is not None:
                lo_list.append(float(np.percentile(arr, 0.3)))
                hi_list.append(float(np.percentile(arr, 99.7)))
                del arr

    if not lo_list:
        print("[fig3] no data — skipping")
        return

    lo, hi = min(lo_list), max(hi_list)
    bin_edges = np.linspace(lo, hi, n_bins + 1)
    centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    def collect(model, dset):
        pdfs = []
        for t in SNAPSHOTS:
            arr = _load_flat(_rfile(results, "iso", t, 0, model), dset)
            if arr is not None:
                pdfs.append(_histogram(arr, bin_edges))
                del arr
        return np.asarray(pdfs) if pdfs else None

    exact = collect("smagorinsky", "Pi_exact")
    smag  = collect("smagorinsky", "Pi")
    wale  = collect("wale",        "Pi")
    dyn   = collect("dynamic",     "Pi")

    fig, ax = plt.subplots(figsize=(9.2, 6.0))
    ax.axvspan(lo, 0.0, color="#fff2cc", alpha=0.65, zorder=0)

    def plot_band(pdfs, color, label, linestyle):
        if pdfs is None:
            return
        mean = pdfs.mean(axis=0)
        lo_b = pdfs.min(axis=0)
        hi_b = pdfs.max(axis=0)
        mband = (lo_b > 0) & (hi_b > 0)
        if mband.any():
            ax.fill_between(centers[mband], lo_b[mband], hi_b[mband],
                            color=color, alpha=0.20, zorder=2, linewidth=0)
        mline = mean > 0
        ax.semilogy(centers[mline], mean[mline],
                    color=color, lw=2.4, linestyle=linestyle, label=label, zorder=4)

    plot_band(exact, "black",                     "Exact DNS",     "-")
    plot_band(smag,  MODEL_COLORS["smagorinsky"], "Smagorinsky",   "-")
    plot_band(wale,  MODEL_COLORS["wale"],        "WALE",          "--")
    plot_band(dyn,   MODEL_COLORS["dynamic"],     "Dynamic Smag.", "-")

    ax.axvline(0.0, color="k", linestyle=":", lw=1.3, zorder=3)
    ax.text(0.02, 0.60, "BACKSCATTER  " + r"($\Pi < 0$)",
            transform=ax.transAxes, ha="left", va="top",
            fontsize=10, fontweight="bold", color="#9c6500", alpha=0.92, zorder=6)
    ax.text(0.98, 0.60, "FORWARD SCATTER  " + r"($\Pi > 0$)",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=10, fontweight="bold", color="#444444", alpha=0.70, zorder=6)

    ax.set_xlabel(r"SGS Dissipation  $\Pi = -\tau_{ij}^{dev}\, S_{ij}$   (code units)")
    ax.set_ylabel("PDF (log scale)")
    ax.legend(loc="upper right", fontsize=10, ncol=2)

    out = figures / "ensemble_dissipation_pdf_iso.png"
    fig.tight_layout()
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig3] {out.name}")


# ===========================================================================
# Figure 4 — Alignment-angle PDF, all three Δ overlaid
# ===========================================================================

def figure4_alignment_pdf(results: Path, figures: Path) -> None:
    bins = np.linspace(-1.0, 1.0, 151)
    centers = 0.5 * (bins[:-1] + bins[1:])

    fig, ax = plt.subplots(figsize=(8.2, 5.5))
    delta_palette = {0: "#2c7bb6", 1: "#fdae61", 2: "#d7191c"}

    for di, ratio in zip(DELTA_INDICES, FILTER_RATIOS):
        pdfs = []
        for t in SNAPSHOTS:
            arr = _load_flat(
                _rfile(results, "iso", t, di, "smagorinsky"), "cos_theta_exact")
            if arr is not None:
                pdfs.append(_histogram(arr, bins))
                del arr
        if not pdfs:
            continue
        pdfs = np.asarray(pdfs)
        mean = pdfs.mean(axis=0)
        lo = pdfs.min(axis=0)
        hi = pdfs.max(axis=0)
        mean_cos = float(np.sum(mean * centers) / np.sum(mean))

        color = delta_palette[di]
        ax.fill_between(centers, lo, hi, color=color,
                        alpha=0.22, zorder=2, linewidth=0)
        ax.plot(centers, mean, color=color, lw=2.3, zorder=4,
                label=rf"$\Delta = {ratio}\Delta x$   "
                      rf"$\langle\cos\theta\rangle = {mean_cos:+.2f}$")

    ymax_data = ax.get_ylim()[1]
    arrow_base_y = 0.55 * ymax_data
    arrow_tip_y = 0.06 * ymax_data
    ax.annotate("", xy=(-1.0, arrow_tip_y), xytext=(-0.52, arrow_base_y),
                arrowprops=dict(arrowstyle="-|>", color="#6a3d9a",
                                lw=2.2, shrinkA=0, shrinkB=0), zorder=5)
    ax.text(-0.50, arrow_base_y,
            "Eddy-viscosity models\n"
            r"predict $\cos\theta \equiv -1$" + "\n"
            r"($\delta$-function at $-1$)",
            color="#6a3d9a", fontsize=10.5, ha="left", va="center",
            bbox=dict(boxstyle="round,pad=0.35", fc="white",
                      ec="#6a3d9a", alpha=0.95), zorder=6)
    ax.axvline(0.0, color="gray", linestyle=":", lw=1.0, alpha=0.6)

    ax.set_xlim(-1.0, 1.0)
    cur_ylim = ax.get_ylim()
    ax.set_ylim(cur_ylim[0], cur_ylim[1] * 1.25)
    ax.set_xlabel(r"$\cos\theta$   (alignment of exact $\tau^{dev}$ with resolved strain $S$)")
    ax.set_ylabel("PDF")
    ax.legend(loc="upper right", fontsize=10, ncol=1)

    out = figures / "ensemble_alignment_pdf_iso.png"
    fig.tight_layout()
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig4] {out.name}")


# ===========================================================================
# Figure 5 — Stability summary (iso + channel, shared y-axis)
# ===========================================================================

def figure5_stability(results: Path, figures: Path) -> None:
    data: Dict[str, Dict[str, List[float]]] = {}
    for dataset in DATASETS:
        data[dataset] = {
            model: _collect_snapshots(results, dataset, 1, model, "correlation")
            for model in MODELS
        }

    all_vals = [v for ds in data.values() for lst in ds.values() for v in lst]
    if not all_vals:
        print("[fig5] no data — skipping")
        return
    ymin = min(all_vals) - 0.03
    ymax = max(all_vals) + 0.035

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.2), sharey=True)
    rng = np.random.RandomState(42)

    for ax, dataset in zip(axes, DATASETS):
        x = np.arange(len(MODELS))
        for i, model in enumerate(MODELS):
            vals = data[dataset][model]
            m = float(np.mean(vals)) if vals else 0.0
            s = float(np.std(vals, ddof=0)) if vals else 0.0
            ax.bar(i, m, width=0.55, color=MODEL_COLORS[model], alpha=0.72,
                   edgecolor="black", linewidth=0.8, zorder=2)
            ax.errorbar(i, m, yerr=s, fmt="none", color="black",
                        capsize=6, lw=1.5, zorder=4)
            if vals:
                x_jit = i + rng.uniform(-0.15, 0.15, len(vals))
                ax.scatter(x_jit, vals, color="black", s=32, zorder=5,
                           edgecolor="white", linewidth=0.7, alpha=0.95)
        ax.axhline(0.0, color="k", lw=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels([MODEL_LABELS[m] for m in MODELS])
        ax.set_xlabel("Isotropic" if dataset == "iso" else "Channel Flow")
        ax.grid(True, axis="y", alpha=0.3)

    axes[0].set_ylim(ymin, ymax)
    axes[0].set_ylabel(r"Correlation coefficient $\rho$")

    out = figures / "ensemble_stability_summary.png"
    fig.tight_layout()
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig5] {out.name}")


# ===========================================================================
# Ensemble-summary tables (CSV + LaTeX) and stability report
# ===========================================================================

def _load_metrics_for_stats(h5_path: Path, model: str) -> Dict[str, float]:
    out: Dict[str, float] = {}
    try:
        with h5py.File(h5_path, "r") as f:
            if model in f:
                grp = f[model]
                for key in ("correlation", "backscatter_frac", "mean_nu_r",
                            "mean_Pi", "mean_cos_theta"):
                    if key in grp.attrs:
                        out[key] = float(grp.attrs[key])
            if "Pi_exact" in f:
                Pi_e = f["Pi_exact"][:]
                out["mean_Pi_exact"] = float(Pi_e.mean())
                out["backscatter_frac_exact"] = float((Pi_e < 0).mean())
    except Exception as exc:
        print(f"  [warn] {h5_path.name}: {exc}")
    return out


def aggregate_stats(results: Path) -> pd.DataFrame:
    rows = []
    keys = ["correlation", "backscatter_frac", "mean_nu_r", "mean_Pi",
            "mean_cos_theta", "mean_Pi_exact", "backscatter_frac_exact"]

    for ds in DATASETS:
        for model in MODELS:
            for delta in DELTA_INDICES:
                snap = {k: [] for k in keys}
                for t in SNAPSHOTS:
                    fp = _rfile(results, ds, t, delta, model)
                    if fp.exists():
                        m = _load_metrics_for_stats(fp, model)
                        for k in keys:
                            if k in m:
                                snap[k].append(m[k])
                if snap["correlation"]:
                    row = {"dataset": ds, "model": model, "delta": delta}
                    for k, v in snap.items():
                        row[f"{k}_mean"] = np.mean(v) if v else np.nan
                        row[f"{k}_std"] = np.std(v) if v else np.nan
                    rows.append(row)
    return pd.DataFrame(rows)


def format_latex_table(df: pd.DataFrame) -> str:
    df_mid = df[df["delta"] == 1].copy()
    df_mid["model"] = df_mid["model"].map(MODEL_LABELS)

    latex = (
        "\\begin{table}[h]\n\\centering\n"
        "\\caption{Ensemble-averaged metrics across 4 snapshots ($\\Delta = 8\\Delta x$).}\n"
        "\\begin{tabular}{l c c c c}\n\\hline\\hline\n"
        "Dataset & Model & Correlation ($\\rho$) & Backscatter \\% & Mean $\\nu_r$ \\\\\n"
        "\\hline\n"
    )
    for ds in DATASETS:
        ds_label = "Isotropic" if ds == "iso" else "Channel"
        for _, r in df_mid[df_mid["dataset"] == ds].iterrows():
            latex += (f"{ds_label} & {r['model']} & "
                      f"{r['correlation_mean']:.3f} $\\pm$ {r['correlation_std']:.3f} & "
                      f"{r['backscatter_frac_mean']*100:.1f}\\% $\\pm$ "
                      f"{r['backscatter_frac_std']*100:.2f}\\% & "
                      f"{r['mean_nu_r_mean']:.3e} \\\\\n")
        latex += "\\hline\n"
    latex += ("\\hline\\hline\n\\end{tabular}\n"
              "\\label{tab:ensemble_stats}\n\\end{table}")
    return latex


def build_stability_report(results: Path) -> pd.DataFrame:
    """Per-snapshot correlation table (delta=1) + STABLE/UNSTABLE flag."""
    rows = []
    for ds in DATASETS:
        for model in MODELS:
            snap_vals = {t: _load_attr(_rfile(results, ds, t, 1, model),
                                       model, "correlation")
                         for t in SNAPSHOTS}
            valid = [v for v in snap_vals.values() if v is not None]
            mean = float(np.mean(valid)) if valid else float("nan")
            std = float(np.std(valid, ddof=0)) if valid else float("nan")
            stable = ("STABLE" if (not np.isnan(std) and std < STABILITY_THRESHOLD)
                      else "UNSTABLE")
            rows.append({
                "dataset": ds,
                "model": MODEL_LABELS[model],
                **{f"corr_{t}": snap_vals[t] for t in SNAPSHOTS},
                "mean": mean,
                "std": std,
                "stable": stable,
            })
    return pd.DataFrame(rows)


def print_stability_table(df: pd.DataFrame) -> None:
    W = 92
    print()
    print("=" * W)
    print("  SNAPSHOT CONVERGENCE TABLE — Pearson ρ at Δ = 8Δx")
    print("=" * W)
    print(f"{'Dataset':<10} {'Model':<20}"
          f"{'t1':>10} {'t2':>10} {'t3':>10} {'t4':>10}"
          f"  {'Mean±Std':>18}  {'Flag':<10}")
    print("-" * W)
    for _, row in df.iterrows():
        ds_label = "ISO" if row["dataset"] == "iso" else "CHANNEL"
        t_strs = []
        for t in SNAPSHOTS:
            v = row[f"corr_{t}"]
            t_strs.append("     N/A  " if v is None or
                          (isinstance(v, float) and np.isnan(v))
                          else f"{v:+.5f}")
        mean_std = (f"{row['mean']:+.5f}±{row['std']:.5f}"
                    if not np.isnan(row["mean"]) else "N/A")
        print(f"{ds_label:<10} {row['model']:<20}"
              f"{t_strs[0]:>10} {t_strs[1]:>10} {t_strs[2]:>10} {t_strs[3]:>10}"
              f"  {mean_std:>18}  {row['stable']:<10}")
    print("=" * W)
    print()


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def _default_paths():
    here = Path(__file__).resolve().parent
    root = here.parent
    return (root / "cpp_driver" / "results",
            root / "figures" / "ensemble")


def main() -> int:
    default_results, default_figures = _default_paths()
    ap = argparse.ArgumentParser(
        description="Generate ensemble figures + summary tables."
    )
    ap.add_argument("--results", type=Path, default=default_results,
                    help="Directory of C++ HDF5 result files.")
    ap.add_argument("--figures", type=Path, default=default_figures,
                    help="Output directory for ensemble PNG figures.")
    ap.add_argument("--tables-out", type=Path, default=None,
                    help="Directory for CSV/LaTeX tables "
                         "(default: --results).")
    args = ap.parse_args()

    args.figures.mkdir(parents=True, exist_ok=True)
    tables_out = args.tables_out or args.results
    tables_out.mkdir(parents=True, exist_ok=True)

    print(f"  results     : {args.results}")
    print(f"  figures     : {args.figures}")
    print(f"  tables out  : {tables_out}")
    print()

    # 5 publication figures
    for ds in DATASETS:
        figure1_correlation(args.results, args.figures, ds)
        figure2_backscatter(args.results, args.figures, ds)
    figure3_dissipation_pdf(args.results, args.figures)
    figure4_alignment_pdf(args.results, args.figures)
    figure5_stability(args.results, args.figures)

    # Summary CSV + LaTeX
    stats = aggregate_stats(args.results)
    if stats.empty:
        print("\n[stats] no result files found — skipping tables")
    else:
        csv_out = tables_out / "ensemble_summary.csv"
        stats.to_csv(csv_out, index=False)
        print(f"\n[stats] CSV → {csv_out}")

        tex = format_latex_table(stats)
        tex_out = tables_out / "table_stats.tex"
        tex_out.write_text(tex)
        print(f"[stats] LaTeX → {tex_out}")

    # Snapshot-convergence / stability report
    conv = build_stability_report(args.results)
    if conv.empty:
        print("\n[stability] no data — skipping")
    else:
        print_stability_table(conv)
        conv_csv = tables_out / "stability_report.csv"
        conv.to_csv(conv_csv, index=False)
        print(f"[stability] CSV → {conv_csv}")

    print(f"\nDone — {args.figures}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
