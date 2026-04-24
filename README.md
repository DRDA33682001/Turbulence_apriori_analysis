# A Priori Analysis of LES Sub-Grid Scale Models

An *a priori* evaluation of three eddy-viscosity sub-grid scale (SGS) closures
for Large-Eddy Simulation, driven by filtered Direct Numerical Simulation
(DNS) data from the Johns Hopkins Turbulence Database (JHTDB).

Three models are compared against the exact filtered SGS stress:

1. **Smagorinsky** — Pope, *Turbulent Flows*, Eq. 13.127–13.128 (Cs = 0.17)
2. **WALE** — Nicoud & Ducros, *Flow Turb. Combust.* 62 (1999), Eq. 14 (Cw = 0.5)
3. **Dynamic Smagorinsky** — Germano *et al.* (1991) + Lilly (1992)

across four diagnostic metrics (correlation, alignment, backscatter,
eddy-viscosity scaling) on two canonical flows (forced isotropic turbulence
and turbulent channel flow) at three filter widths spanning the inertial
subrange.

---

## Repository Layout

```

├── README.md                  ← you are here (quick-start + layout)
├── REPORT.md                  ← full project write-up
├── .gitignore
│
├── data_download/             ← JHTDB download pipeline
│   ├── download.py
│   ├── requirements.txt
│   └── README.md              ← how to obtain the ~6 GB of DNS snapshots
│
├── cpp_driver/                ← heavy numerical engine (C++)
│   ├── CMakeLists.txt
│   ├── run_full_ensemble.sh
│   ├── src/cpp/               ← filter, SGS models, metrics, HDF5 I/O
│   └── README.md              ← build instructions, CLI, output format
│
├── post_processing/           ← all Python plotting / diagnostics
│   ├── *.py
│   ├── requirements.txt
│   └── README.md              ← which script makes which figure
│
└── figures/                   ← every plot checked into the repo
    ├── *.png / *.pdf          ← single-snapshot and summary plots
    └── ensemble/*.png         ← full-ensemble aggregate plots
```

Large, regenerated-on-demand assets (DNS downloads ~6 GB, C++ result HDF5s
~5 GB) are **not** checked in — see the sub-folder READMEs for how to
produce them.

---

## Quick Start

```bash
# 1. Install Python deps
cd data_download && pip install -r requirements.txt

# 2. Download DNS snapshots from JHTDB (~6 GB, 10–20 min)
python download.py --output ../cpp_driver/data/

# 3. Build the C++ driver
cd ../cpp_driver && mkdir -p build && cd build
cmake .. && make -j

# 4. Run the full ensemble (~30–50 min; writes 72 HDF5 files to ../results/)
cd ..
bash run_full_ensemble.sh

# 5. Regenerate every figure from the new results
cd ../post_processing
pip install -r requirements.txt        # same as step 1 — skip if already done
python run_pipeline.py \
    --results ../cpp_driver/results \
    --figures ../figures
```

Each numbered step has its own README with more detail and troubleshooting
notes.

---

## Results at a Glance

| Metric | Smagorinsky | WALE | Dynamic | Expected |
|--------|:-----------:|:----:|:-------:|:--------:|
| Tensor correlation ρ(τ_exact, τ_model) (isotropic, Δ=8·dx) | ~0.29 | ~0.31 | ~0.33 | ~0.30 (Clark 1979) |
| Alignment cos θ with strain rate | peaks near −1, broad | similar | similar | broad (McMillan 1979) |
| Backscatter fraction of Π_exact | 0 % (by construction) | 0 % | ~15 % | ~30–40 % |
| ν_r scaling exponent vs Δ | ≈ 4/3 | ≈ 4/3 | ≈ 4/3 | 4/3 (K41) |

The low correlation (~0.3), the inability of static models to reproduce
backscatter, and the broad alignment distribution together illustrate the
limits of the eddy-viscosity hypothesis that underlies all three models.
See [`REPORT.md`](REPORT.md) for the full discussion.

---

## References

| Reference | Used for |
|-----------|----------|
| Pope, *Turbulent Flows* Ch. 13 | Spectral filter (Eq. 13.22), exact τᵢⱼ (Eq. 13.93), Smagorinsky (13.127–13.128), K41 scaling (13.131) |
| Germano *et al.* (1991) *Phys. Fluids A* 3(7) | Dynamic procedure, Eq. 8 |
| Lilly (1992) *Phys. Fluids A* 4(3) | Least-squares Cs² minimization, Eq. 10 |
| Nicoud & Ducros (1999) *Flow Turb. Combust.* 62 | WALE eddy viscosity, Eq. 14 |
| Clark, Ferziger & Reynolds (1979) *J. Fluid Mech.* 91 | A priori testing framework, correlation definition |
| McMillan & Ferziger (1979) *AIAA J.* 17(12) | Direct testing methodology, backscatter fraction |

External codes consulted:
[spectralDNS/spectralDNS](https://github.com/spectralDNS/spectralDNS) (FFT filter reference),
[Romit-Maulik/ML_2D_Turbulence](https://github.com/Romit-Maulik/ML_2D_Turbulence) (a priori metrics reference).
