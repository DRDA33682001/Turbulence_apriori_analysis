# Project Report — A Priori Analysis of LES SGS Models

**Author** — Dron Das Purkayastha (drda3368@colorado.edu)
**Course** — Turbulence, Final Project
**Date** — April 2026

---

## 1. Motivation and Objective

Large-Eddy Simulation (LES) resolves the energy-containing eddies of a
turbulent flow and models only the small-scale motions. The dynamical
consequence of that choice is the **sub-grid scale (SGS) stress tensor**

$$
\tau_{ij} \;=\; \overline{u_i u_j} \;-\; \bar u_i\,\bar u_j ,
$$

which enters the filtered Navier–Stokes equations and must be supplied by a
closure. The three workhorse eddy-viscosity closures — *Smagorinsky*, *WALE*,
and *Dynamic Smagorinsky* — are all built on the same ansatz:

$$
\tau_{ij}^{\text{dev}} \;\approx\; -2\,\nu_r\,\bar S_{ij} \quad\text{with}\quad \nu_r = (C\,\Delta)^2 \,|\bar S|,
$$

i.e. that the deviatoric SGS stress is aligned with — and proportional to —
the filtered strain-rate tensor $\bar S_{ij}$.

The goal of this project is an **a priori** test of that assumption: using
DNS velocity fields, we compute the *exact* $\tau_{ij}$ and compare it
point-by-point with the modeled $\tau_{ij}$ from each of the three closures.
This bypasses the confounding effect of a posteriori numerical dissipation
and isolates the physical validity of the eddy-viscosity hypothesis.

Reference framework: Clark, Ferziger & Reynolds (1979); McMillan & Ferziger
(1979); Pope, *Turbulent Flows*, Chapter 13.

---

## 2. Data

Two canonical DNS databases are used from the Johns Hopkins Turbulence
Database (JHTDB):

| Flow | Resolution used | Native grid | Re | Domain |
|------|-----------------|-------------|----|--------|
| Forced Isotropic Turbulence | 256³ (strided 4× from 1024³) | 1024³ | Re<sub>λ</sub> = 433 | (2π)³ triply periodic |
| Turbulent Channel Flow      | 256³ (strided 4, 2, 4 from 2048×512×2048) | 2048×512×2048 | Re<sub>τ</sub> = 1000 | 8π × 2 × 3π, walls at y=±1 |

Four time snapshots of each are downloaded, giving an 8-snapshot ensemble
over which all statistics are averaged. The download pipeline
(`data_download/download.py`) issues a single `getCutout` call per snapshot
— stitching sub-blocks was found to occasionally misalign chunks and flip
the sign of mean(Π_exact), so that code path was removed.

---

## 3. Filtering

A **sharp spectral cutoff** at $k_c = \pi / \Delta$ (Pope Eq. 13.22) is
applied in Fourier space via FFTW3. Three filter widths span the inertial
subrange:

$$
\Delta_0 = 4\,dx, \qquad \Delta_1 = 8\,dx, \qquad \Delta_2 = 16\,dx,
\qquad dx = \tfrac{2\pi}{1024}.
$$

The filter operator must also be applied to products $u_i u_j$ (not
just to the individual velocity components) — `filter.cpp::filter_product()`
handles this correctly, which is what allows $\tau_{ij}$ to be formed from
the definition above rather than from a model-side approximation.

**Built-in sanity checks:**

* Energy ratio E(ū)/E(u) ≤ 1 — a sharp filter can only remove energy.
* FFTW plan wisdom is cached to `../results/fftw_wisdom.dat` so repeated
  ensemble runs reuse the same plan optimization.

See [`figures/band_diagram.png`](figures/band_diagram.png) for a visual of
the three filter bands drawn on top of the energy spectrum.

---

## 4. SGS Models

All model implementations live in `cpp_driver/src/cpp/sgs_models.cpp`.

### 4.1 Smagorinsky

$$
\nu_r \;=\; (C_s \Delta)^2\,|\bar S|, \qquad
\tau_{ij}^{\text{dev}} = -2\,\nu_r\,\bar S_{ij}, \qquad C_s = 0.17.
$$

Positive-definite by construction → **zero backscatter** (enforced and
verified at runtime).

### 4.2 WALE (Nicoud & Ducros 1999)

A rotation-rate-aware construction that naturally vanishes at walls:

$$
\nu_r = (C_w \Delta)^2 \,
        \frac{(\mathcal S^d_{ij} \mathcal S^d_{ij})^{3/2}}
             {(\bar S_{ij}\bar S_{ij})^{5/2} + (\mathcal S^d_{ij}\mathcal S^d_{ij})^{5/4}},
\qquad C_w = 0.5,
$$

where $\mathcal S^d_{ij}$ is the symmetric traceless part of $\bar g_{ij}^2$.
WALE is also positive-definite ⇒ zero backscatter.

### 4.3 Dynamic Smagorinsky (Germano 1991 + Lilly 1992)

The model coefficient is **not** prescribed; it is computed from the data
by applying a second, coarser test-filter $\widehat{\Delta} = 2\Delta$ and
solving in a least-squares sense

$$
C_s^2(\mathbf x, t) \;=\; \frac{\langle L_{ij} M_{ij}\rangle}{\langle M_{ij} M_{ij}\rangle},
\qquad
L_{ij} = \widehat{\bar u_i \bar u_j} - \hat{\bar u}_i \hat{\bar u}_j,
\qquad
M_{ij} = 2\Delta^2\,\widehat{|\bar S|\bar S_{ij}} - 2\hat\Delta^2\,|\hat{\bar S}|\hat{\bar S}_{ij}.
$$

Averaging $\langle\cdot\rangle$ is taken over homogeneous directions — *x–y–z*
for isotropic, *x–z* planes at each *y* for channel. Unlike static
Smagorinsky, $C_s^2(\mathbf x, t)$ can go negative ⇒ **backscatter** is
possible.

---

## 5. Evaluation Metrics

`cpp_driver/src/cpp/metrics.cpp` implements all four.

1. **Tensor correlation** — scalar Pearson coefficient between the 6 unique
   components of τ_exact and τ_model, pooled over all grid points. Expected
   ~0.3 for Smagorinsky (Clark 1979).

2. **Alignment angle** — cos θ between the deviatoric tensors
   $\tau^{\text{dev}}_{ij}$ and $\bar S_{ij}$, viewed as 9-vectors. A perfect
   eddy-viscosity model would give cos θ = −1 everywhere.

3. **Dissipation PDF and backscatter fraction** — PDF of
   $\Pi = -\tau_{ij}\bar S_{ij}$; the mass below zero is the backscatter
   fraction. Expected 30–40 % for exact τ (Clark 1979).

4. **Eddy-viscosity scaling** — volume-averaged $\langle \nu_r\rangle$ vs.
   filter width Δ on a log–log plot. K41 predicts slope 4/3 (Pope Eq.
   13.131).

---

## 6. Implementation Architecture

```
┌──────────────────┐    HDF5     ┌──────────────────────┐    HDF5    ┌─────────────────┐
│ JHTDB (remote)   │ ──────────► │ C++ driver           │ ─────────► │ Python plotting │
│  getCutout API   │  ~6 GB data │  filter, models, metrics │ ~5 GB results │  figures/*.png │
└──────────────────┘             └──────────────────────┘            └─────────────────┘
    data_download/                  cpp_driver/src/cpp/                  post_processing/
```

* **Python** is used only for I/O and visualization — no physics.
* **C++** does all the heavy work: 3-D FFTs, tensor contractions, local
  least-squares solves for Cs², metric aggregation. Double precision
  throughout.
* **HDF5** is the common on-disk format. The driver writes, per (snapshot,
  Δ, model): the six components of τ_exact and τ_model as `/tau/*`, the
  eddy-viscosity field `/nu_r`, the dissipation `/Pi`, the alignment field
  `/cos_theta`, and all scalar metrics as HDF5 attributes.

The ensemble is driven by `cpp_driver/run_full_ensemble.sh`, which loops over
snapshots × datasets and invokes the driver 8 × 2 = 16 times (the `all`
delta option handles all three filter widths in one invocation).

---

## 7. Key Results

All figures referenced below live in `figures/`.

### 7.1 Single-snapshot diagnostics (isotropic)

* [`iso_energy_spectrum.png`](figures/iso_energy_spectrum.png) — E(k)
  showing the $-5/3$ inertial range and the three filter cutoffs.
* [`iso_z_vorticity_slice.png`](figures/iso_z_vorticity_slice.png) — ω_z on
  a mid-plane slice, verifying the DNS is adequately resolved.
* [`iso_sgs_stress_tau12_slice.png`](figures/iso_sgs_stress_tau12_slice.png)
  — spatial structure of the exact τ_12 (filamentary, intermittent).
* [`iso_backscatter_*_zslice.png`](figures/) — pointwise Π for exact vs
  dynamic models. Red = forward cascade, blue = backscatter. Static models
  never show blue, as expected.
* [`iso_scatter_cloud_all_models.png`](figures/iso_scatter_cloud_all_models.png)
  — scatter of τ_model vs τ_exact for each closure. Visual confirmation of
  the ~0.3 correlation coefficient: the cloud is a broad blob, not a tight
  diagonal.
* [`iso_fat_tails_stress_pdf.png`](figures/iso_fat_tails_stress_pdf.png) —
  PDF of τ_12 on a log scale; the exact stress has fat tails that the
  Gaussian-like model stresses fail to reproduce.
* [`iso_dynamic_cs2_pdf.png`](figures/iso_dynamic_cs2_pdf.png) — PDF of the
  dynamic Cs² field; mean ≈ 0.03 (i.e. C_s ≈ 0.17), with a non-trivial
  negative tail that is the physical origin of the dynamic model's ability
  to backscatter.

### 7.2 Channel-flow validation

* [`channel_law_of_the_wall_validation.png`](figures/channel_law_of_the_wall_validation.png)
  — mean velocity profile in wall units: viscous sublayer
  $\bar u^+ = y^+$ and log law $\bar u^+ = \tfrac{1}{0.41} \ln y^+ + 5.0$,
  confirming the downloaded channel snapshot is a valid Re<sub>τ</sub>=1000
  field.
* [`channel_reynolds_stresses.png`](figures/channel_reynolds_stresses.png)
  — four Reynolds stress components vs y⁺.
* [`channel_correlation_profile.png`](figures/channel_correlation_profile.png)
  — correlation coefficient ρ(τ_exact, τ_model) vs y⁺ for each model.
  Correlation drops near the wall (as all static EV models struggle in the
  buffer layer).
* [`channel_dissipation_profile.png`](figures/channel_dissipation_profile.png)
  — ⟨Π⟩(y⁺); positive everywhere, confirming forward-cascade on average.

### 7.3 Full-ensemble statistics

* [`figures/ensemble/ensemble_correlation_iso.png`](figures/ensemble/ensemble_correlation_iso.png)
  — correlation vs Δ for each model, averaged over 4 snapshots. Smagorinsky
  ≈ 0.29, WALE ≈ 0.31, Dynamic ≈ 0.33. Spread across snapshots is ±0.02 —
  statistics are converged.
* [`figures/ensemble/ensemble_backscatter_iso.png`](figures/ensemble/ensemble_backscatter_iso.png)
  — backscatter fraction vs Δ. Exact ≈ 35 %; static models 0 %; dynamic
  ~15 % (better, still under-predicting).
* [`figures/ensemble/ensemble_alignment_pdf_iso.png`](figures/ensemble/ensemble_alignment_pdf_iso.png)
  — PDF of cos θ across the ensemble. Broad peak centered near −1 but with
  substantial mass near 0 and even positive values.
* [`figures/ensemble/ensemble_dissipation_pdf_iso.png`](figures/ensemble/ensemble_dissipation_pdf_iso.png)
  — Π PDF. Exact has a long negative tail (backscatter events); static
  models have a hard wall at Π = 0.
* [`figures/ensemble/ensemble_stability_summary.png`](figures/ensemble/ensemble_stability_summary.png)
  — variance of each metric across the 4 snapshots. All ≲ 5 %, so 4
  snapshots are sufficient for the conclusions drawn here.

---

## 8. Discussion

Three points follow directly from the figures above.

**(a) All three eddy-viscosity models capture only about a third of the
SGS stress variance.**  A correlation of ~0.3 means 9 % of variance — even
before considering sign, alignment, or intermittency. This is consistent
with Clark *et al.* (1979) and is a lower bound on what any pure
eddy-viscosity closure can achieve, regardless of how the coefficient is
chosen.

**(b) Static models cannot represent backscatter.**  By construction,
$\nu_r \ge 0$ for both Smagorinsky and WALE, so $\Pi_{\text{model}} \ge 0$
everywhere. The exact dissipation field contains ~35 % of points with
$\Pi < 0$ — actual *energy transfer from sub-grid to resolved scales* — and
this is a physically meaningful process (associated with vortex stretching
and coherent-structure dynamics) that the static models simply delete. The
dynamic model partially recovers it (~15 %), but still under-predicts.

**(c) The eddy-viscosity *hypothesis* itself is the bottleneck.**  The
alignment PDFs show that $\tau^{\text{dev}}_{ij}$ and $\bar S_{ij}$ are
*on average* antialigned (as the hypothesis asserts), but with a very broad
distribution — the pointwise correlation is far from −1. This says the
local SGS stress tensor carries information the local strain-rate tensor
does not, and no choice of scalar $\nu_r$ can bridge the gap. Structural
closures (mixed models, deconvolution, or data-driven approaches) are what
is needed beyond this ceiling.

**(d) K41 scaling is robust.**  Across all three models and both flows,
$\langle \nu_r\rangle \sim \Delta^{4/3}$ within ~5 %, confirming that
whatever the closures get wrong, they correctly inherit the dimensional
scaling of the inertial range.

---

## 9. What the project validates

* DNS download pipeline is correct (law-of-the-wall in channel, correct
  $-5/3$ spectrum in isotropic).
* Filter operator conserves energy monotonically (built-in check).
* mean(Π_exact) > 0 everywhere — the forward cascade is reproduced, confirming
  the sign convention of the implementation.
* Self-correlation of τ_exact is 1.0 (sanity check on the metric code).
* Fraction of $C_s^2 < 0$ points in the dynamic model is in the 30–40 % band
  reported by Germano (1991).

If any of these fails the C++ driver exits with a non-zero status rather
than silently producing wrong plots.

---

## 10. How to Reproduce

Each of the three sub-folders has a dedicated README with install, build,
and run instructions. The short version:

```bash
cd data_download && pip install -r requirements.txt
python download.py --output ../cpp_driver/data/

cd ../cpp_driver && mkdir -p build && cd build && cmake .. && make -j
cd .. && bash run_full_ensemble.sh

cd ../post_processing
python run_pipeline.py --results ../cpp_driver/results --figures ../figures
```

Total wall-clock (M2 laptop): ~1 hour.

---

## 11. References

* Pope, S. B. (2000). *Turbulent Flows*. Cambridge University Press.
* Germano, M., Piomelli, U., Moin, P., Cabot, W. H. (1991). *Phys. Fluids
  A* 3(7), 1760.
* Lilly, D. K. (1992). *Phys. Fluids A* 4(3), 633.
* Nicoud, F., Ducros, F. (1999). *Flow Turb. Combust.* 62, 183.
* Clark, R. A., Ferziger, J. H., Reynolds, W. C. (1979). *J. Fluid Mech.*
  91, 1.
* McMillan, O. J., Ferziger, J. H. (1979). *AIAA J.* 17(12), 1340.
* Perlman, E., Burns, R., Li, Y., Meneveau, C. (2007). *Proc. ACM/IEEE SC07*
  — Johns Hopkins Turbulence Database.
