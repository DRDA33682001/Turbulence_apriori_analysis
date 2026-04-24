#!/usr/bin/env python3
"""
plot_figures.py — Single-snapshot and diagnostic plots for the a priori
LES SGS model evaluation.

Consolidates what used to live in five separate scripts:
  * advanced_plots.py             (scatter cloud, fat tails, corr/dissip
                                   profiles, backscatter maps, structure vs
                                   stress co-location, dynamic Cs^2 PDF)
  * project_completeness_plots.py (channel Reynolds stresses, exact vs
                                   model backscatter, dynamic Cs^2 effective)
  * validate_channel.py           (channel law-of-the-wall validation)
  * make_band_diagram.py          (filter-band schematic over E(k))
  * visualize.py                  (E(k) energy spectrum; other entries in
                                   the old file were stale and have been
                                   dropped)

Physics literature basis:
    Clark, Ferziger & Reynolds (1979)   — scatter cloud, fat tails,
                                          exact vs model backscatter
    Piomelli, Cabot, Moin & Lee (1988)  — wall-normal profiles
    Meneveau & Katz (2000)              — backscatter co-location
    Moser, Kim & Mansour (1999)         — Reynolds-stress benchmark
    Millikan (1938)                     — log-law of the wall

Produces every non-ensemble PNG/PDF in figures/ — see the `main()`
function at the bottom for the full list.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Optional, Tuple

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np


# ---------------------------------------------------------------------------
# Global constants (mirror C++ config.h)
# ---------------------------------------------------------------------------
N: int = 256
DX: float = 2.0 * np.pi / N
FILTER_RATIOS = [4, 8, 16]                   # Δ/dx for the three filter widths
FILTER_WIDTHS = np.array([r * DX for r in FILTER_RATIOS])
SNAPSHOTS = ["t1", "t2", "t3", "t4"]
MODELS = ["smagorinsky", "wale", "dynamic"]
DELTAS = [0, 1, 2]
DELTA_IDX_MID = 1    # mid-filter (Δ = 8·dx) — default for single-delta figures

MODEL_COLORS = {
    "smagorinsky": "#2166ac",
    "wale":        "#d6604d",
    "dynamic":     "#1a9641",
    "exact":       "#333333",
}
MODEL_LABELS = {
    "smagorinsky": "Smagorinsky",
    "wale":        "WALE",
    "dynamic":     "Dynamic Smag.",
    "exact":       "Exact DNS",
}

plt.rcParams.update({
    "font.size":        12,
    "axes.labelsize":   13,
    "axes.titlesize":   13,
    "legend.fontsize":  11,
    "figure.dpi":       150,
    "lines.linewidth":  1.8,
    "axes.grid":        True,
    "grid.alpha":       0.3,
})


# ---------------------------------------------------------------------------
# Path / IO helpers
# ---------------------------------------------------------------------------

def _rfile(results: Path, ds: str, t: str, di: int, model: str) -> Path:
    """Path to a C++ results HDF5 for one (dataset, snapshot, Δ, model)."""
    return results / f"{ds}_{t}_delta{di}_{model}.h5"


def _vfile(data: Path, ds: str, t: str) -> Path:
    """Path to the velocity HDF5 for one (dataset, snapshot)."""
    name = "isotropic_256" if ds == "iso" else "channel_256"
    return data / f"{name}_{t}.h5"


def _load_component(fname: Path, group: str, comp: str) -> np.ndarray:
    with h5py.File(fname, "r") as f:
        return f[group][comp][:]


def _load_scalar(fname: Path, dset: str) -> np.ndarray:
    with h5py.File(fname, "r") as f:
        return f[dset][:]


def _pearson(a: np.ndarray, b: np.ndarray) -> float:
    ac = a - a.mean()
    bc = b - b.mean()
    denom = np.sqrt((ac @ ac) * (bc @ bc))
    return float((ac @ bc) / (denom + 1e-30))


# ---------------------------------------------------------------------------
# 1. Band diagram — SGS filter bands over idealised E(k)
# ---------------------------------------------------------------------------

def plot_band_diagram(figures: Path) -> None:
    """Schematic E(k) with the three SGS filter bands shaded."""
    k = np.logspace(np.log10(1.5), np.log10(400), 600)
    k0, k_diss, k_max = 4.0, 230.0, 128.0

    def E(k_arr):
        out = np.empty_like(k_arr)
        for i, ki in enumerate(k_arr):
            if ki < k0:
                out[i] = ki ** 2
            elif ki <= k_diss:
                out[i] = ki ** (-5 / 3)
            else:
                out[i] = k_diss ** (-5 / 3) * np.exp(-0.5 * (ki / k_diss - 1) ** 2 / 0.08)
        return out / (k0 ** (-5 / 3))

    Ek = E(k)

    dx = 2 * np.pi / 256
    k_c = {
        r"$\Delta = 4\Delta x$":  np.pi / (4 * dx),
        r"$\Delta = 8\Delta x$":  np.pi / (8 * dx),
        r"$\Delta = 16\Delta x$": np.pi / (16 * dx),
    }
    band_colors = ["#3A86FF", "#FF6B35", "#6ABF69"]

    fig, ax = plt.subplots(figsize=(5.5, 3.6))
    ax.loglog(k, Ek, "k-", lw=1.8, zorder=5)

    k_ref = np.array([5.0, 60.0])
    ax.loglog(k_ref, 0.8 * k_ref ** (-5 / 3) / k_ref[0] ** (-5 / 3),
              "k--", lw=0.9, alpha=0.55, zorder=3)
    ax.text(18, 0.065, r"$k^{-5/3}$", fontsize=8, color="0.35")

    for (label, kc), color in zip(k_c.items(), band_colors):
        k_band = k[(k >= kc) & (k <= k_max)]
        E_band = Ek[(k >= kc) & (k <= k_max)]
        if len(k_band) < 2:
            continue
        ax.fill_between(k_band, 1e-5, E_band, alpha=0.22, color=color, zorder=2)
        E_at_kc = np.interp(kc, k, Ek)
        ax.axvline(kc, color=color, lw=1.1, ls=":", zorder=4)
        ax.text(kc * 1.05, E_at_kc * 1.5, f"$k_c={kc:.0f}$",
                fontsize=7, color=color, va="bottom")

    E_at_kmax = np.interp(k_max, k, Ek)
    ax.axvline(k_max, color="0.35", lw=1.5, ls="--", zorder=4)
    ax.text(k_max * 1.06, E_at_kmax * 3.5,
            r"$k_{\max}=128$" "\n(JHTDB cutoff)",
            fontsize=7.5, color="0.25", va="top", ha="left")

    ax.text(2.2, 1.6, "Energy-\ncontaining", fontsize=7, color="0.4",
            ha="center", va="bottom")
    ax.annotate("", xy=(k0 * 0.9, 1.2), xytext=(k0 * 2.5, 1.2),
                arrowprops=dict(arrowstyle="->", color="0.45", lw=0.8))
    ax.text(55, 3.5, "Inertial\nsubrange", fontsize=7, color="0.4",
            ha="center", va="bottom")
    ax.text(230, 0.3, "Diss.\nrange", fontsize=7, color="0.4",
            ha="center", va="top", style="italic")

    ax.set_xlabel(r"Wavenumber $k$", fontsize=9)
    ax.set_ylabel(r"$E(k)$  [arbitrary units]", fontsize=9)
    ax.set_xlim(1.5, 380)
    ax.set_ylim(5e-4, 8)
    ax.tick_params(labelsize=8)

    handles = [mpatches.Patch(color="k", label=r"$E(k)\propto k^{-5/3}$")]
    for (label, _), color in zip(k_c.items(), band_colors):
        handles.append(mpatches.Patch(facecolor=color, alpha=0.35,
                                      label=f"SGS eval. band ({label})"))
    ax.legend(handles=handles, fontsize=7, loc="lower left",
              framealpha=0.8, borderpad=0.6)

    fig.tight_layout()
    for ext in ("pdf", "png"):
        out = figures / f"band_diagram.{ext}"
        fig.savefig(out, dpi=200, bbox_inches="tight")
        print(f"[band_diagram] {out.name}")
    plt.close(fig)


# ---------------------------------------------------------------------------
# 2. Isotropic 3-D energy spectrum E(k) from DNS velocity
# ---------------------------------------------------------------------------

def plot_energy_spectrum(data: Path, figures: Path, snapshot: str = "t1") -> None:
    """3-D radially-averaged E(k) from the isotropic velocity file."""
    h5_path = _vfile(data, "iso", snapshot)
    with h5py.File(h5_path, "r") as f:
        u1 = f["u1"][()].astype(np.float64)
        u2 = f["u2"][()].astype(np.float64)
        u3 = f["u3"][()].astype(np.float64)

    N_ = u1.shape[0]
    norm = 1.0 / N_ ** 3
    U1 = np.fft.fftn(u1) * norm
    U2 = np.fft.fftn(u2) * norm
    U3 = np.fft.fftn(u3) * norm

    E_mode = 0.5 * (np.abs(U1) ** 2 + np.abs(U2) ** 2 + np.abs(U3) ** 2)

    ki = np.fft.fftfreq(N_, d=1.0 / N_).astype(int)
    kx, ky, kz = ki[:, None, None], ki[None, :, None], ki[None, None, :]
    k_shell = np.round(
        np.sqrt((kx ** 2 + ky ** 2 + kz ** 2).astype(np.float64))
    ).astype(int)

    k_max = N_ // 2
    E_k = np.bincount(k_shell.ravel(), weights=E_mode.ravel(),
                      minlength=k_max + 1)
    k_vals = np.arange(1, k_max + 1)
    E_vals = E_k[1:k_max + 1]

    kc_vals = [N_ / (2.0 * n) for n in FILTER_RATIOS]
    kc_colors = [MODEL_COLORS["smagorinsky"], MODEL_COLORS["wale"],
                 MODEL_COLORS["dynamic"]]

    fig, ax = plt.subplots(figsize=(6, 5))
    ax.loglog(k_vals, E_vals, color="k", lw=1.5, label=r"$E(k)$")

    ref_idx = np.searchsorted(k_vals, 10)
    k_ref = np.array([5.0, k_max * 0.7])
    y_ref = E_vals[ref_idx] * (k_ref / k_vals[ref_idx]) ** (-5.0 / 3.0)
    ax.loglog(k_ref, y_ref, "k--", lw=1.0, alpha=0.55, label=r"$k^{-5/3}$")

    for i, (kci, col) in enumerate(zip(kc_vals, kc_colors)):
        ax.axvline(kci, color=col, lw=1.2, ls=":",
                   label=rf"$k_c\,(\Delta_{{{i+1}}},{FILTER_RATIOS[i]}\Delta x)$")

    ax.set_xlabel(r"Wavenumber $k$")
    ax.set_ylabel(r"$E(k)$")
    ax.legend(fontsize=9, loc="lower left")
    fig.tight_layout()

    out = figures / "iso_energy_spectrum.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"[energy_spectrum] {out.name}")


# ---------------------------------------------------------------------------
# 3. Scatter cloud: τ₁₂_model vs τ₁₂_exact (three-panel)
# ---------------------------------------------------------------------------

def plot_scatter_cloud(results: Path, figures: Path, snapshot: str = "t1") -> None:
    """Hexbin joint density of model τ₁₂ vs exact τ₁₂ for all three models."""
    models_info = [
        ("smagorinsky", "tau_smagorinsky", "Smagorinsky"),
        ("wale",        "tau_wale",        "WALE"),
        ("dynamic",     "tau_dynamic",     "Dynamic Smag."),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(14.5, 5.8), constrained_layout=True)
    rng = np.random.default_rng(seed=42)

    ref_file = _rfile(results, "iso", snapshot, DELTA_IDX_MID, "smagorinsky")
    exact_all = _load_component(ref_file, "tau_exact_dev", "T_12").ravel()
    lim = float(np.percentile(np.abs(exact_all), 99.5))

    last_hb = None
    for ax, (model, group, label) in zip(axes, models_info):
        fname = _rfile(results, "iso", snapshot, DELTA_IDX_MID, model)
        exact = _load_component(fname, "tau_exact_dev", "T_12").ravel()
        mod = _load_component(fname, group, "T_12").ravel()
        r = _pearson(exact, mod)

        idx = rng.choice(len(exact), size=min(200_000, len(exact)), replace=False)
        last_hb = ax.hexbin(exact[idx], mod[idx], gridsize=70, cmap="Blues",
                            mincnt=1, norm=mcolors.LogNorm(),
                            extent=(-lim, lim, -lim, lim))
        ax.plot([-lim, lim], [-lim, lim], "r--", lw=1.3, label=r"$y = x$")
        ax.axhline(0, color="k", lw=0.5, alpha=0.4)
        ax.axvline(0, color="k", lw=0.5, alpha=0.4)
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_aspect("equal")
        ax.set_xlabel(rf"$\tau_{{12}}^{{\rm exact}}$   ({label})")
        ax.set_ylabel(rf"$\tau_{{12}}^{{\rm {label.split()[0]}}}$")
        ax.text(0.04, 0.96, rf"$\rho_{{\tau_{{12}}}} = {r:+.3f}$",
                transform=ax.transAxes, va="top", ha="left", fontsize=11,
                bbox=dict(boxstyle="round,pad=0.3", fc="white",
                          ec=MODEL_COLORS[model], alpha=0.90))

    cb = fig.colorbar(last_hb, ax=axes, shrink=0.88, pad=0.02, aspect=30)
    cb.set_label("Point count (log scale)")

    out = figures / "iso_scatter_cloud_all_models.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[scatter_cloud] {out.name}")


# ---------------------------------------------------------------------------
# 4. Fat-tail PDF of τ₁₂ on log-y scale
# ---------------------------------------------------------------------------

def plot_fat_tails(results: Path, figures: Path, snapshot: str = "t1") -> None:
    fname_smag = _rfile(results, "iso", snapshot, DELTA_IDX_MID, "smagorinsky")
    fname_wale = _rfile(results, "iso", snapshot, DELTA_IDX_MID, "wale")
    fname_dyn  = _rfile(results, "iso", snapshot, DELTA_IDX_MID, "dynamic")

    exact   = _load_component(fname_smag, "tau_exact_dev",   "T_12").ravel()
    smag    = _load_component(fname_smag, "tau_smagorinsky", "T_12").ravel()
    wale    = _load_component(fname_wale, "tau_wale",        "T_12").ravel()
    dynamic = _load_component(fname_dyn,  "tau_dynamic",     "T_12").ravel()

    sig_e = float(exact.std())
    bins = np.linspace(-6, 6, 240)

    def _flatness(data: np.ndarray) -> float:
        d = data - data.mean()
        v = float((d * d).mean())
        return float(((d ** 4).mean()) / (v * v + 1e-30))

    fig, ax = plt.subplots(figsize=(8.8, 5.0))

    def _semilogy_pdf(data, label, color, ls="-"):
        h, edges = np.histogram(data / sig_e, bins=bins, density=True)
        xc = 0.5 * (edges[:-1] + edges[1:])
        mask = h > 0
        ax.semilogy(xc[mask], h[mask], ls, color=color, lw=2.0,
                    label=f"{label}   (flatness = {_flatness(data):.2f})")

    _semilogy_pdf(exact,   "Exact DNS",           "#222222")
    _semilogy_pdf(smag,    "Smagorinsky",         MODEL_COLORS["smagorinsky"], "--")
    _semilogy_pdf(wale,    "WALE",                MODEL_COLORS["wale"],        "--")
    _semilogy_pdf(dynamic, "Dynamic Smagorinsky", MODEL_COLORS["dynamic"],     "-.")

    x_ref = np.linspace(-6, 6, 500)
    ax.semilogy(x_ref, np.exp(-0.5 * x_ref ** 2) / np.sqrt(2.0 * np.pi),
                "k:", lw=1.3, label="Gaussian reference (flatness = 3.00)")

    ax.set_xlim(-6, 6)
    ax.set_ylim(1e-6, 1.5)
    ax.set_xlabel(r"$\tau_{12}\,/\,\sigma_{\tau_{12}^{\rm exact}}$")
    ax.set_ylabel("PDF (log scale)")
    ax.legend(loc="lower center", fontsize=9)
    fig.tight_layout()

    out = figures / "iso_fat_tails_stress_pdf.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"[fat_tails] {out.name}")


# ---------------------------------------------------------------------------
# 5. Channel correlation profile ρ(y⁺)
# ---------------------------------------------------------------------------

def plot_channel_correlation_profile(results: Path, data: Path, figures: Path,
                                     snapshot: str = "t1") -> None:
    with h5py.File(_vfile(data, "channel", snapshot), "r") as f:
        y_plus = f["y_plus"][:]
    half = int(np.argmax(y_plus)) + 1
    yp_half = y_plus[:half]

    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    for model, ls, col, label in [
        ("smagorinsky", "-",  "C0", "Smagorinsky"),
        ("wale",        "--", "C1", "WALE"),
        ("dynamic",     "-.", "C2", "Dynamic"),
    ]:
        fname = _rfile(results, "channel", snapshot, DELTA_IDX_MID, model)
        exact = _load_component(fname, "tau_exact_dev", "T_12")
        mod = _load_component(fname, f"tau_{model}", "T_12")
        rho = np.empty(half)
        for j in range(half):
            a = exact[:, j, :].ravel()
            b = mod[:, j, :].ravel()
            rho[j] = _pearson(a, b)
        ax.plot(yp_half, rho, ls, color=col, label=label)

    ax.axhline(0, color="k", lw=0.8, ls=":")
    ax.set_xscale("log")
    ax.set_xlim(max(yp_half[1], 0.5), yp_half[-1])
    ax.set_xlabel(r"$y^+$")
    ax.set_ylabel(r"$\rho\!\left(\tau_{12}^{\rm model},\,\tau_{12}^{\rm exact}\right)$")
    ax.legend()
    fig.tight_layout()

    out = figures / "channel_correlation_profile.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"[channel_correlation_profile] {out.name}")


# ---------------------------------------------------------------------------
# 6. Channel SGS dissipation profile ⟨Π⟩(y⁺)
# ---------------------------------------------------------------------------

def plot_channel_dissipation_profile(results: Path, data: Path, figures: Path,
                                     snapshot: str = "t1") -> None:
    with h5py.File(_vfile(data, "channel", snapshot), "r") as f:
        y_plus = f["y_plus"][:]
    half = int(np.argmax(y_plus)) + 1
    yp_half = y_plus[:half]

    fig, ax = plt.subplots(figsize=(6.5, 4.5))

    fname_smag = _rfile(results, "channel", snapshot, DELTA_IDX_MID, "smagorinsky")
    pi_exact = _load_scalar(fname_smag, "Pi_exact")
    ax.plot(yp_half, pi_exact.mean(axis=(0, 2))[:half],
            "k-", lw=2.2, label="Exact DNS", zorder=5)
    del pi_exact

    for model, ls, col, label in [
        ("smagorinsky", "-",  "C0", "Smagorinsky"),
        ("wale",        "--", "C1", "WALE"),
        ("dynamic",     "-.", "C2", "Dynamic"),
    ]:
        pi_mod = _load_scalar(
            _rfile(results, "channel", snapshot, DELTA_IDX_MID, model), "Pi")
        ax.plot(yp_half, pi_mod.mean(axis=(0, 2))[:half],
                ls, color=col, label=label)
        del pi_mod

    ax.axhline(0, color="k", lw=0.8, ls=":")
    ax.set_xscale("log")
    ax.set_xlim(max(yp_half[1], 0.5), yp_half[-1])
    ax.set_xlabel(r"$y^+$")
    ax.set_ylabel(r"Mean SGS dissipation $\langle\Pi\rangle$")
    ax.legend()
    fig.tight_layout()

    out = figures / "channel_dissipation_profile.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"[channel_dissipation_profile] {out.name}")


# ---------------------------------------------------------------------------
# 7. Backscatter maps (exact, dynamic, combined) — z-slice
# ---------------------------------------------------------------------------

def plot_backscatter_maps(results: Path, figures: Path, snapshot: str = "t1") -> None:
    k_slice = 128
    fname_dyn = _rfile(results, "iso", snapshot, DELTA_IDX_MID, "dynamic")
    with h5py.File(fname_dyn, "r") as f:
        pi_exact = f["Pi_exact"][:, :, k_slice]
        pi_dyn = f["Pi"][:, :, k_slice]

    vmax = float(np.percentile(np.abs(pi_exact), 98))
    kw = dict(origin="lower", cmap="RdBu_r", vmin=-vmax, vmax=vmax,
              extent=[0, 2 * np.pi, 0, 2 * np.pi], aspect="equal")

    for field, cbarlabel, name in [
        (pi_exact, r"$\Pi^{\rm exact}$", "iso_backscatter_exact_zslice.png"),
        (pi_dyn,   r"$\Pi^{\rm dyn}$",   "iso_backscatter_dynamic_zslice.png"),
    ]:
        fig, ax = plt.subplots(figsize=(5.2, 4.6))
        im = ax.imshow(field.T, **kw)
        cb = fig.colorbar(im, ax=ax, shrink=0.85, pad=0.02)
        cb.set_label(cbarlabel, rotation=0, labelpad=10)
        ax.set_xlabel(r"$x_1$")
        ax.set_ylabel(r"$x_2$")
        fig.tight_layout()
        fig.savefig(figures / name, dpi=150)
        plt.close(fig)
        print(f"[backscatter] {name}")

    fig, axes = plt.subplots(1, 2, figsize=(11, 5.6), constrained_layout=True)
    for ax, field, panel, title in [
        (axes[0], pi_exact, "a",
         r"Exact DNS  $\Pi^{\rm exact} = -\tau^{\rm exact}_{ij} S_{ij}$"),
        (axes[1], pi_dyn,   "b",
         r"Dynamic Smag.  $\Pi^{\rm dyn} = -\tau^{\rm dyn}_{ij} S_{ij}$"),
    ]:
        im = ax.imshow(field.T, **kw)
        ax.set_xlabel(r"$x_1$")
        ax.set_ylabel(r"$x_2$")
        bs_frac = float((field < 0).mean()) * 100.0
        ax.text(0.03, 0.97, f"({panel})  backscatter area: {bs_frac:.1f}%",
                transform=ax.transAxes, va="top", ha="left", fontsize=10,
                bbox=dict(boxstyle="round,pad=0.3", fc="white",
                          ec="gray", alpha=0.88))
    cb = fig.colorbar(im, ax=axes, shrink=0.88, pad=0.02, aspect=30)
    cb.set_label(r"$\Pi$  (shared range, clipped at DNS 98th pct.)")

    out = figures / "iso_backscatter_comparison_zslice.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[backscatter] {out.name}")


# ---------------------------------------------------------------------------
# 8. Structure vs stress co-location (vorticity + τ₁₂ slices)
# ---------------------------------------------------------------------------

def plot_structure_vs_stress(results: Path, data: Path, figures: Path,
                             snapshot: str = "t1") -> None:
    k_slice = 128
    dx = 2.0 * np.pi / 256

    with h5py.File(_vfile(data, "iso", snapshot), "r") as f:
        u1 = f["u1"][:, :, k_slice].astype(np.float64)
        u2 = f["u2"][:, :, k_slice].astype(np.float64)
    omega3 = np.gradient(u2, dx, axis=0) - np.gradient(u1, dx, axis=1)

    fname_smag = _rfile(results, "iso", snapshot, DELTA_IDX_MID, "smagorinsky")
    tau12 = _load_component(fname_smag, "tau_exact_dev", "T_12")[:, :, k_slice]

    ext = [0, 2 * np.pi, 0, 2 * np.pi]
    base = dict(origin="lower", extent=ext, aspect="equal")

    for field, cbarlabel, name in [
        (omega3, r"$\omega_3$",              "iso_z_vorticity_slice.png"),
        (tau12,  r"$\tau_{12}^{\rm exact}$", "iso_sgs_stress_tau12_slice.png"),
    ]:
        vmax = float(np.percentile(np.abs(field), 97))
        fig, ax = plt.subplots(figsize=(6.0, 5.2))
        im = ax.imshow(field.T, cmap="RdBu_r", vmin=-vmax, vmax=vmax, **base)
        cb = fig.colorbar(im, ax=ax, shrink=0.85, pad=0.02)
        cb.set_label(cbarlabel, rotation=0, labelpad=12)
        ax.set_xlabel(r"$x_1$")
        ax.set_ylabel(r"$x_2$")
        fig.tight_layout()
        fig.savefig(figures / name, dpi=150)
        plt.close(fig)
        print(f"[structure_vs_stress] {name}")

    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.8), constrained_layout=True)
    fig.get_layout_engine().set(w_pad=0.08, wspace=0.15)
    for ax, field, panel, cbarlabel in [
        (axes[0], omega3, "a", r"$\omega_3$"),
        (axes[1], tau12,  "b", r"$\tau_{12}^{\rm exact}$"),
    ]:
        vmax = float(np.percentile(np.abs(field), 97))
        im = ax.imshow(field.T, cmap="RdBu_r", vmin=-vmax, vmax=vmax, **base)
        cb = fig.colorbar(im, ax=ax, shrink=0.85, pad=0.04)
        cb.set_label(cbarlabel, rotation=0, labelpad=14)
        ax.set_xlabel(r"$x_1$")
        ax.set_ylabel(r"$x_2$")
        ax.text(0.03, 0.97, f"({panel})", transform=ax.transAxes,
                va="top", ha="left", fontsize=11, fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.25", fc="white",
                          ec="gray", alpha=0.85))

    out = figures / "iso_structure_vs_stress_zslice.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[structure_vs_stress] {out.name}")


# ---------------------------------------------------------------------------
# 9. Dynamic Cs² pointwise PDF
# ---------------------------------------------------------------------------

def plot_dynamic_cs2_pdf(results: Path, figures: Path, snapshot: str = "t1") -> None:
    fname = _rfile(results, "iso", snapshot, DELTA_IDX_MID, "dynamic")
    cs2 = _load_scalar(fname, "Cs2").ravel()

    lo = float(np.percentile(cs2, 0.3))
    hi = float(np.percentile(cs2, 99.7))
    lo = min(lo, -abs(hi) * 0.5)
    bins = np.linspace(lo, hi, 220)
    centers = 0.5 * (bins[:-1] + bins[1:])
    width = bins[1] - bins[0]
    counts, _ = np.histogram(cs2, bins=bins, density=True)

    fig, ax = plt.subplots(figsize=(7.8, 5.0))
    ax.axvspan(lo, 0.0, color="#fff2cc", alpha=0.55, zorder=0)

    neg = centers < 0.0
    ax.bar(centers[neg], counts[neg], width=width, color="#f4a582", alpha=0.85,
           edgecolor="none", label=r"$C_s^2 < 0$  (backscatter: $\nu_r<0$)")
    ax.bar(centers[~neg], counts[~neg], width=width, color="#1a9641", alpha=0.55,
           edgecolor="none", label=r"$C_s^2 \geq 0$  (forward scatter)")
    ax.plot(centers, counts, color="#1a9641", lw=1.1, zorder=5)

    ax.axvline(0.0, color="k", lw=1.3, ls="--", zorder=6, label=r"$C_s^2 = 0$")
    ax.axvline(0.0289, color="#555555", lw=1.2, ls=":", zorder=6,
               label=r"Lilly (1967): $C_s^2 = 0.17^2 \approx 0.029$")

    ax.set_xlabel(r"Dynamic Smagorinsky coefficient  $C_s^2$")
    ax.set_ylabel("PDF")
    ax.legend(loc="upper right", fontsize=9)
    fig.tight_layout()

    out = figures / "iso_dynamic_cs2_pdf.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"[dynamic_cs2_pdf] {out.name}")


# ---------------------------------------------------------------------------
# 10. Channel law-of-the-wall validation
# ---------------------------------------------------------------------------

def plot_channel_law_of_wall(data: Path, figures: Path, snapshot: str = "t1") -> None:
    RE_TAU = 1_000.0
    with h5py.File(_vfile(data, "channel", snapshot), "r") as f:
        u1 = f["u1"][:]
        y_plus = f["y_plus"][:]

    U_mean = u1.mean(axis=(0, 2))
    del u1

    i_ctr = int(np.argmax(y_plus))
    half = i_ctr + 1
    yp = y_plus[:half]
    U_h = U_mean[:half]

    mask_visc = (yp >= 1.0) & (yp <= 8.0)
    u_tau = float(np.linalg.lstsq(
        yp[mask_visc].reshape(-1, 1), U_h[mask_visc], rcond=None)[0][0])
    U_plus = U_h / u_tau

    yp_ref = np.geomspace(0.15, 1_300.0, 1_000)
    U_visc = yp_ref
    U_log = (1.0 / 0.41) * np.log(yp_ref) + 5.2

    fig, ax = plt.subplots(figsize=(8, 5.5))
    ax.semilogx(yp_ref, U_visc, "k--", lw=1.4,
                label=r"Viscous sublayer: $U^+ = y^+$")
    ax.semilogx(yp_ref, U_log, "k-.", lw=1.4,
                label=r"Log law: $U^+ = \frac{1}{0.41}\ln y^+ + 5.2$")

    mask_plot = yp >= 1.0
    ax.semilogx(yp[mask_plot], U_plus[mask_plot],
                "C0-o", ms=3.5, mew=0, lw=2.0,
                label=fr"JHTDB channel (Re_τ = {int(RE_TAU)}, 256³ sub-sample)")

    for xv, txt in [(5, r"$y^+ = 5$"), (30, r"$y^+ = 30$")]:
        ax.axvline(xv, color="silver", lw=0.9, ls=":")
        ax.text(xv * 1.08, 0.7, txt, fontsize=9, color="gray", va="bottom")

    ax.text(1.5, 2, "Viscous\nsublayer", fontsize=9, ha="center", color="gray")
    ax.text(13, 8, "Buffer\nlayer", fontsize=9, ha="center", color="gray")
    ax.text(100, 13, "Log-law region", fontsize=9, ha="center", color="gray")
    ax.text(700, 20, "Wake region", fontsize=9, ha="center", color="gray")

    ax.set_xscale("log")
    ax.set_xlim(0.8, 1_300)
    ax.set_ylim(0, 28)
    ax.set_xlabel(r"$y^+$  (log scale)")
    ax.set_ylabel(r"$U^+  =  \langle U \rangle / u_\tau$")
    ax.legend(loc="upper left")
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()

    out = figures / "channel_law_of_the_wall_validation.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"[law_of_the_wall] {out.name}  (u_τ ≈ {u_tau:.5f})")


# ---------------------------------------------------------------------------
# 11. Channel Reynolds-stress profiles (ensemble-averaged)
# ---------------------------------------------------------------------------

def plot_channel_reynolds_stresses(data: Path, figures: Path) -> None:
    uu_list, vv_list, ww_list, uv_list = [], [], [], []
    u_tau_list = []
    y_plus = None

    for t in SNAPSHOTS:
        fpath = _vfile(data, "channel", t)
        if not fpath.exists():
            print(f"  [channel_reynolds] missing {fpath.name}, skipping")
            continue
        with h5py.File(fpath, "r") as f:
            u = f["u1"][:]
            v = f["u2"][:]
            w = f["u3"][:]
            yp_cur = f["y_plus"][:]
        if y_plus is None:
            y_plus = yp_cur

        U = u.mean(axis=(0, 2))
        V = v.mean(axis=(0, 2))
        W = w.mean(axis=(0, 2))
        up = u - U[None, :, None]
        vp = v - V[None, :, None]
        wp = w - W[None, :, None]

        uu_list.append((up ** 2).mean(axis=(0, 2)))
        vv_list.append((vp ** 2).mean(axis=(0, 2)))
        ww_list.append((wp ** 2).mean(axis=(0, 2)))
        uv_list.append((up * vp).mean(axis=(0, 2)))

        half_i = int(np.argmax(yp_cur)) + 1
        yp_h = yp_cur[:half_i]
        U_h = U[:half_i]
        mask = (yp_h >= 1.0) & (yp_h <= 8.0)
        u_tau_list.append(float(np.linalg.lstsq(
            yp_h[mask].reshape(-1, 1), U_h[mask], rcond=None)[0][0]))

    if not uu_list:
        print("  [channel_reynolds] no data — skipping")
        return

    uu_m = np.mean(uu_list, axis=0)
    vv_m = np.mean(vv_list, axis=0)
    ww_m = np.mean(ww_list, axis=0)
    uv_m = np.mean(uv_list, axis=0)
    u_tau = float(np.mean(u_tau_list))

    half = int(np.argmax(y_plus)) + 1
    yp = y_plus[:half]
    uu_p = uu_m[:half] / u_tau ** 2
    vv_p = vv_m[:half] / u_tau ** 2
    ww_p = ww_m[:half] / u_tau ** 2
    uv_p = uv_m[:half] / u_tau ** 2

    # Moser 1999 Re_tau~590 benchmark overlay (visual reference)
    moser_yp = np.array([5, 10, 15, 20, 30, 50, 100, 200, 300, 500])
    moser_uu = np.array([2.2, 6.3, 7.6, 7.3, 6.0, 4.9, 4.3, 3.6, 3.2, 2.7])
    moser_vv = np.array([0.01, 0.08, 0.19, 0.31, 0.58, 0.87, 1.0, 1.0, 0.98, 0.85])
    moser_ww = np.array([0.28, 1.0, 1.6, 1.9, 2.0, 2.1, 1.9, 1.7, 1.5, 1.2])
    moser_uv = np.array([-0.04, -0.22, -0.40, -0.55, -0.70, -0.82,
                         -0.86, -0.77, -0.65, -0.43])

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2))
    ax = axes[0]
    ax.semilogx(yp, uu_p, "-",  lw=2.0, color="#d62728",
                label=r"$\langle u'u' \rangle^+$  (this study)")
    ax.semilogx(yp, vv_p, "--", lw=2.0, color="#1f77b4",
                label=r"$\langle v'v' \rangle^+$")
    ax.semilogx(yp, ww_p, "-.", lw=2.0, color="#2ca02c",
                label=r"$\langle w'w' \rangle^+$")
    ax.semilogx(moser_yp, moser_uu, "o", color="#d62728", ms=5, alpha=0.5,
                mec="white", mew=0.7, label="Moser 1999 (Re_τ≈590)")
    ax.semilogx(moser_yp, moser_vv, "s", color="#1f77b4", ms=5, alpha=0.5,
                mec="white", mew=0.7)
    ax.semilogx(moser_yp, moser_ww, "^", color="#2ca02c", ms=5, alpha=0.5,
                mec="white", mew=0.7)
    ax.set_xlabel(r"$y^+$   (diagonal Reynolds stresses)")
    ax.set_ylabel(r"Reynolds normal stress $\langle u_i' u_i' \rangle / u_\tau^2$")
    ax.set_xlim(0.5, 1200)
    ax.set_ylim(0, max(12.0, 1.1 * uu_p.max()))
    ax.legend(loc="upper right", fontsize=9)

    ax = axes[1]
    ax.semilogx(yp, -uv_p, "-", lw=2.0, color="#8c564b",
                label=r"$-\langle u'v' \rangle^+$  (this study)")
    ax.semilogx(moser_yp, -moser_uv, "D", color="#8c564b", ms=5, alpha=0.5,
                mec="white", mew=0.7, label="Moser 1999 (Re_τ≈590)")
    ax.axhline(1.0, color="gray", lw=0.8, ls=":", alpha=0.6)
    ax.text(1, 1.02, "total-stress limit", color="gray", fontsize=8.5)
    ax.set_xlabel(r"$y^+$   (Reynolds shear stress)")
    ax.set_ylabel(r"Reynolds shear stress $-\langle u'v' \rangle / u_\tau^2$")
    ax.set_xlim(0.5, 1200)
    ax.set_ylim(-0.05, 1.15)
    ax.legend(loc="upper right", fontsize=9)

    fig.tight_layout()
    out = figures / "channel_reynolds_stresses.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[channel_reynolds] {out.name}")


# ---------------------------------------------------------------------------
# 12. Exact vs model backscatter (bar chart, ensemble-averaged)
# ---------------------------------------------------------------------------

def plot_exact_vs_model_backscatter(results: Path, figures: Path) -> None:
    def _collect(ds, di):
        bs_exact, bs_model = [], {m: [] for m in MODELS}
        for t in SNAPSHOTS:
            fp_ref = _rfile(results, ds, t, di, "smagorinsky")
            if fp_ref.exists():
                with h5py.File(fp_ref, "r") as f:
                    Pi_e = f["Pi_exact"][:]
                    bs_exact.append(float((Pi_e < 0).mean()))
            for m in MODELS:
                fp = _rfile(results, ds, t, di, m)
                if fp.exists():
                    with h5py.File(fp, "r") as f:
                        if m in f and "backscatter_frac" in f[m].attrs:
                            bs_model[m].append(float(f[m].attrs["backscatter_frac"]))
        ex_m = float(np.mean(bs_exact)) if bs_exact else np.nan
        ex_s = float(np.std(bs_exact)) if bs_exact else np.nan
        mdl = {m: (float(np.mean(v)) if v else np.nan,
                   float(np.std(v)) if v else np.nan) for m, v in bs_model.items()}
        return ex_m, ex_s, mdl

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.0), sharey=True)
    for panel_idx, ds in enumerate(["iso", "channel"]):
        ax = axes[panel_idx]
        x_base = np.arange(len(DELTAS))
        bar_w = 0.2

        for i, cat in enumerate(["exact"] + MODELS):
            means, stds = [], []
            for di in DELTAS:
                ex_m, ex_s, mdl = _collect(ds, di)
                if cat == "exact":
                    means.append(ex_m * 100.0)
                    stds.append(ex_s * 100.0)
                else:
                    m_m, m_s = mdl[cat]
                    means.append(m_m * 100.0)
                    stds.append(m_s * 100.0)
            pos = x_base + (i - 1.5) * bar_w
            ax.bar(pos, means, width=bar_w, yerr=stds, capsize=3,
                   color=MODEL_COLORS[cat], alpha=0.85,
                   edgecolor="white", linewidth=0.7,
                   label=MODEL_LABELS[cat])

        ax.axhspan(30, 40, color="#fff2cc", alpha=0.6, zorder=0,
                   label="Clark 1979: 30-40% (exact)")
        ax.axhline(50, color="gray", lw=0.8, ls=":", alpha=0.5)

        ax.set_xticks(x_base)
        ax.set_xticklabels([fr"$\Delta={r}\Delta x$" for r in FILTER_RATIOS])
        ax.set_xlabel(f"Filter width   ({'Isotropic' if ds == 'iso' else 'Channel'})")
        if panel_idx == 0:
            ax.set_ylabel("Backscatter fraction (%)")
        ax.set_ylim(0, 60)
        if panel_idx == 1:
            ax.legend(loc="upper right", fontsize=8, framealpha=0.92)

    fig.tight_layout()
    out = figures / "exact_vs_model_backscatter.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[exact_vs_model_backscatter] {out.name}")


# ---------------------------------------------------------------------------
# 13. Dynamic effective Cs² (iso: scalar vs Δ, channel: profile vs y⁺)
# ---------------------------------------------------------------------------

def plot_dynamic_cs2_effective(results: Path, data: Path, figures: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.0))

    # Iso panel: box-averaged scalar Cs² vs Δ
    ax = axes[0]
    cs2_all = [[] for _ in DELTAS]
    for i, di in enumerate(DELTAS):
        for t in SNAPSHOTS:
            fp = _rfile(results, "iso", t, di, "dynamic")
            if fp.exists():
                with h5py.File(fp, "r") as f:
                    if "Cs2_eff" in f:
                        cs2_all[i].append(float(f["Cs2_eff"][:].mean()))

    means = [np.mean(v) if v else np.nan for v in cs2_all]
    stds = [np.std(v) if v else np.nan for v in cs2_all]

    x = np.arange(len(DELTAS))
    ax.bar(x, means, yerr=stds, capsize=5,
           color="#1a9641", alpha=0.8, edgecolor="white", linewidth=0.7,
           label=r"$(C_s^2)_\mathrm{eff}$ box-averaged")
    for i, vals in enumerate(cs2_all):
        if vals:
            xp = np.full(len(vals), x[i]) + np.random.uniform(-0.06, 0.06, len(vals))
            ax.scatter(xp, vals, s=30, color="#0b5d28", zorder=5,
                       edgecolor="white", linewidth=0.7)
    ax.axhline(0.0289, color="#555555", lw=1.3, ls=":",
               label=r"Lilly 1967: $C_s^2 = 0.17^2 \approx 0.0289$")
    ax.axhline(0.0, color="k", lw=0.8, ls="--", alpha=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels([fr"$\Delta={r}\Delta x$" for r in FILTER_RATIOS])
    ax.set_xlabel(r"Filter width   (Isotropic: scalar $(C_s^2)_{\mathrm{eff}}$)")
    ax.set_ylabel(r"Effective $C_s^2$ (box average)")
    ax.legend(loc="upper right", fontsize=9)

    # Channel panel: Cs²(y⁺) profile at mid filter
    ax = axes[1]
    profiles = []
    y_plus = None
    for t in SNAPSHOTS:
        fp = _rfile(results, "channel", t, DELTA_IDX_MID, "dynamic")
        if fp.exists():
            with h5py.File(fp, "r") as f:
                if "Cs2_eff" in f:
                    profiles.append(f["Cs2_eff"][:][0, :, 0])
        if y_plus is None:
            vf = _vfile(data, "channel", t)
            if vf.exists():
                with h5py.File(vf, "r") as f:
                    y_plus = f["y_plus"][:]

    if profiles and y_plus is not None:
        prof = np.asarray(profiles)
        mean = prof.mean(axis=0)
        std = prof.std(axis=0)
        half = int(np.argmax(y_plus)) + 1
        yp = y_plus[:half]
        ax.semilogx(yp, mean[:half], "-", lw=2.0, color="#1a9641",
                    label=r"$(C_s^2)_\mathrm{eff}(y^+)$ (mean over 4 snapshots)")
        ax.fill_between(yp, (mean - std)[:half], (mean + std)[:half],
                        color="#1a9641", alpha=0.25, label=r"$\pm 1\sigma$")
    else:
        ax.text(0.5, 0.5, "Cs2_eff not found in channel results",
                ha="center", va="center", transform=ax.transAxes)

    ax.axhline(0.0289, color="#555555", lw=1.3, ls=":",
               label=r"Lilly 1967: $C_s^2 \approx 0.0289$")
    ax.axhline(0.0, color="k", lw=0.8, ls="--", alpha=0.6)
    ax.set_xlabel(r"$y^+$   (Channel: $(C_s^2)_{\mathrm{eff}}(y^+)$ at $\Delta=8\Delta x$)")
    ax.set_ylabel(r"Effective $C_s^2(y^+)$  (plane-averaged)")
    ax.legend(loc="lower right", fontsize=9)

    fig.tight_layout()
    out = figures / "dynamic_cs2_effective.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[dynamic_cs2_effective] {out.name}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def _default_paths() -> Tuple[Path, Path, Path]:
    """Defaults matching the prior_analysis/ GitHub layout."""
    here = Path(__file__).resolve().parent           # post_processing/
    root = here.parent                               # prior_analysis/
    return (root / "cpp_driver" / "results",
            root / "cpp_driver" / "data",
            root / "figures")


def main() -> int:
    default_results, default_data, default_figures = _default_paths()
    ap = argparse.ArgumentParser(
        description="Generate all single-snapshot / diagnostic figures."
    )
    ap.add_argument("--results", type=Path, default=default_results,
                    help="Directory of C++ HDF5 result files.")
    ap.add_argument("--data", type=Path, default=default_data,
                    help="Directory of DNS velocity HDF5 files.")
    ap.add_argument("--figures", type=Path, default=default_figures,
                    help="Output directory for PNG/PDF figures.")
    ap.add_argument("--snapshot", default="t1",
                    help="Snapshot tag used by single-snapshot figures "
                         "(default: t1).")
    args = ap.parse_args()

    args.figures.mkdir(parents=True, exist_ok=True)
    print(f"  results : {args.results}")
    print(f"  data    : {args.data}")
    print(f"  figures : {args.figures}")
    print(f"  snapshot: {args.snapshot}")
    print()

    # Order is deliberate: schematic first, then iso, then channel.
    plot_band_diagram(args.figures)
    plot_energy_spectrum(args.data, args.figures, args.snapshot)
    plot_scatter_cloud(args.results, args.figures, args.snapshot)
    plot_fat_tails(args.results, args.figures, args.snapshot)
    plot_backscatter_maps(args.results, args.figures, args.snapshot)
    plot_structure_vs_stress(args.results, args.data, args.figures, args.snapshot)
    plot_dynamic_cs2_pdf(args.results, args.figures, args.snapshot)

    plot_channel_law_of_wall(args.data, args.figures, args.snapshot)
    plot_channel_correlation_profile(args.results, args.data, args.figures,
                                     args.snapshot)
    plot_channel_dissipation_profile(args.results, args.data, args.figures,
                                     args.snapshot)
    plot_channel_reynolds_stresses(args.data, args.figures)

    plot_exact_vs_model_backscatter(args.results, args.figures)
    plot_dynamic_cs2_effective(args.results, args.data, args.figures)

    print()
    print(f"Done — figures written to {args.figures}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
