# Post-Processing & Visualization

The Python scripts in this folder read the HDF5 result files produced by
`../cpp_driver/` and produce every figure in `../figures/`. Nothing here
does any physics — all modeling and tensor arithmetic lives in C++.
Python only reads, aggregates, and plots.

> **Input assumption:** `cpp_driver/run_full_ensemble.sh` has finished and
> `../cpp_driver/results/` contains the 72 HDF5 files it writes. The
> isotropic/channel velocity HDF5s are expected in `../cpp_driver/data/`.

---

## 1. Install

```bash
cd post_processing
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Reusing the `data_download/` venv works too — the dependency list is
identical except you don't need the JHTDB client here.

---

## 2. Scripts in this folder

The post-processing code has been consolidated into three files:

| Script | Role |
|--------|------|
| [plot_figures.py](plot_figures.py) | All single-snapshot / diagnostic plots (13 figures). Consolidates what used to be 5 scripts: `advanced_plots.py`, `project_completeness_plots.py`, `validate_channel.py`, `make_band_diagram.py`, and the lone salvageable routine from `visualize.py`. |
| [plot_ensemble.py](plot_ensemble.py) | All ensemble plots (5 figures in `figures/ensemble/`) plus the summary CSV, AIAA LaTeX table, and snapshot-convergence stability report. Consolidates `ensemble_visualize.py`, `ensemble_stats.py`, and `advanced_plots_ensemble.py`. |
| [run_pipeline.py](run_pipeline.py) | Tiny orchestrator that invokes the two scripts above in sequence. |

Each script can also be run on its own — see the CLI flags at the top of
every `main()`.

---

## 3. End-to-end — C++ outputs → figures

```bash
cd post_processing
python run_pipeline.py \
    --results ../cpp_driver/results \
    --data    ../cpp_driver/data \
    --figures ../figures
```

This writes every PNG/PDF already checked into `../figures/` and
`../figures/ensemble/`, plus `ensemble_summary.csv`, `table_stats.tex`,
and `stability_report.csv` into `../cpp_driver/results/`.

Run the two scripts directly if you only want a subset:

```bash
python plot_figures.py  --results ../cpp_driver/results \
                        --data    ../cpp_driver/data \
                        --figures ../figures
python plot_ensemble.py --results ../cpp_driver/results \
                        --figures ../figures/ensemble
```

---

## 4. What gets produced

### `plot_figures.py` — single-snapshot figures in `../figures/`

| File | Source / content |
|------|-----------|
| `band_diagram.{png,pdf}` | Schematic E(k) with the three SGS filter bands shaded. |
| `iso_energy_spectrum.png` | 3-D radially-averaged E(k) from the DNS velocity field. |
| `iso_scatter_cloud_all_models.png` | Hexbin density, τ₁₂ model vs exact, one panel per model. |
| `iso_fat_tails_stress_pdf.png` | Semi-log PDF of τ₁₂ showing intermittent fat tails. |
| `iso_backscatter_{exact,dynamic,comparison}_zslice.png` | Π field on a z-slice; individual + side-by-side panels. |
| `iso_z_vorticity_slice.png` | ω₃ on the same slice. |
| `iso_sgs_stress_tau12_slice.png` | τ₁₂_exact on the same slice. |
| `iso_structure_vs_stress_zslice.png` | ω₃ / τ₁₂ side-by-side (co-location). |
| `iso_dynamic_cs2_pdf.png` | Pointwise PDF of the dynamic Smagorinsky Cs² field. |
| `channel_law_of_the_wall_validation.png` | U⁺ vs y⁺ with viscous and log-law references. |
| `channel_correlation_profile.png` | ρ(y⁺) for each model. |
| `channel_dissipation_profile.png` | ⟨Π⟩(y⁺) for each model + exact DNS. |
| `channel_reynolds_stresses.png` | ⟨uᵢ'uⱼ'⟩⁺(y⁺), with Moser 1999 benchmark overlay. |
| `exact_vs_model_backscatter.png` | Backscatter %: DNS vs each model, grouped by Δ. |
| `dynamic_cs2_effective.png` | (iso) scalar Cs² vs Δ + (channel) Cs²(y⁺). |

### `plot_ensemble.py` — ensemble figures in `../figures/ensemble/`

| File | Content |
|------|---------|
| `ensemble_correlation_{iso,channel}.png` | ρ vs Δ, mean ± std over 4 snapshots with snapshot scatter. |
| `ensemble_backscatter_{iso,channel}.png` | Backscatter % vs Δ, grouped bars with DNS reference band. |
| `ensemble_dissipation_pdf_iso.png` | Π PDF — min–max shaded band + mean line for each source. |
| `ensemble_alignment_pdf_iso.png` | cos θ PDF overlaid for all three Δ. |
| `ensemble_stability_summary.png` | ρ per model, both datasets, shared y-axis. |

Also written (tables, to `--results` or `--tables-out`):

| File | Content |
|------|---------|
| `ensemble_summary.csv` | Full mean/std table for every (dataset, model, Δ). |
| `table_stats.tex` | AIAA-style LaTeX table for mid-filter. |
| `stability_report.csv` | Per-snapshot ρ + STABLE/UNSTABLE flag. |

---

## 5. Troubleshooting

* **“File not found”** — check that `../cpp_driver/results/` contains the
  72 HDF5 files written by `run_full_ensemble.sh`. If a snapshot died
  mid-run, rerun only that snapshot.
* **Font warnings from Matplotlib** — harmless; the scripts fall back to
  a safe default font if the preferred one is missing.
* **NaN / empty PDFs** — usually a result file came from a failed
  run-time physics check. Rerun that C++ invocation with `--debug`.
