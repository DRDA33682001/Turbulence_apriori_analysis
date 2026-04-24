# C++ Driver — A Priori SGS Analysis

This folder contains the heavy numerical engine of the project. The executable
`les_apriori` reads a filtered DNS snapshot, applies three sub-grid scale
(SGS) models, and writes HDF5 result files that the Python scripts in
`../post_processing/` turn into figures.

> **Results are not checked into this repository.** A full ensemble produces
> **72 HDF5 files, ~5 GB total** (8 snapshots × 3 filter widths × 3 models).
> They are regenerated from the data in `data/` by running the driver as
> described below.

---

## 1. What each source file does

```
src/cpp/
├── config.h         Global compile-time constants: grid size N=256, dx=2π/1024,
│                    filter widths Δ ∈ {4dx, 8dx, 16dx}, Smagorinsky Cs=0.17,
│                    WALE Cw=0.5, helper functions idx(), check_scalar().
│
├── io.h / io.cpp    HDF5 I/O layer. Reads (u1, u2, u3) velocity from input
│                    files; writes 6-component symmetric tensors, scalar fields,
│                    and scalar metrics (as HDF5 attributes) for each run.
│
├── filter.h / .cpp  Spectral filtering via FFTW. Sharp cutoff at k_c = π/Δ
│                    (Pope Eq. 13.22). Caches FFTW plans to fftw_wisdom.dat so
│                    the second run re-uses prior plan optimization.
│                    Key API: spectral_filter_3d(), filter_product() — the
│                    latter correctly computes filter(uᵢuⱼ) needed for the
│                    exact SGS stress.
│
├── sgs_exact.h / .cpp
│                    Exact SGS stress from DNS:   τᵢⱼ = filter(uᵢuⱼ) − ūᵢ ūⱼ.
│                    Also computes the velocity gradient tensor and strain rate
│                    S̄ᵢⱼ of the filtered field.
│
├── sgs_models.h / .cpp
│                    Three model implementations:
│                       • Smagorinsky   (Pope Eq. 13.127–13.128)
│                       • WALE          (Nicoud & Ducros 1999, Eq. 14)
│                       • Dynamic Smag. (Germano 1991 + Lilly 1992
│                                        least-squares Cs²(x,t) field)
│
├── metrics.h / .cpp Four evaluation metrics per model:
│                      1. alignment angle cos θ between τᵢⱼ_dev and S̄ᵢⱼ
│                      2. tensor correlation coefficient ρ(τ_exact, τ_model)
│                      3. backscatter-fraction PDF of Π = −τᵢⱼ S̄ᵢⱼ
│                      4. eddy-viscosity scaling ν_r vs Δ^{4/3}   (K41)
│
└── main.cpp         CLI driver and top-level pipeline. See §3 below.
```

All numerics use double precision (FFTW3 double + HDF5 C API + Eigen for
small-matrix solves in the dynamic procedure).

---

## 2. Build

**Requirements** — CMake ≥ 3.15, a C++17 compiler, FFTW3 (double & single
precision), HDF5 (C API), Eigen3.

On macOS with Homebrew:

```bash
brew install cmake fftw hdf5 eigen
```

On Ubuntu:

```bash
sudo apt install cmake build-essential libfftw3-dev libhdf5-dev libeigen3-dev
```

Build:

```bash
cd cpp_driver
mkdir -p build && cd build
cmake ..
make -j$(nproc)        # on macOS: make -j$(sysctl -n hw.ncpu)
```

The executable `les_apriori` is placed in `build/`. CMake auto-detects all
`.cpp` files under `src/cpp/`.

---

## 3. Run

### 3.1 CLI usage

```
./les_apriori  <dataset>  <delta_index>  [options]

dataset       iso | channel
delta_index   0 | 1 | 2 | all          (filter widths 4dx, 8dx, 16dx)
options
  --input  <path>   override input HDF5 (default: ../data/<dataset>_256.h5)
  --tag    <label>  override output filename prefix (default: the dataset name)
  --debug           abort with exit code 1 if any physics check fails
```

### 3.2 Single run

From `build/`, filter isotropic data at Δ = 8·dx with all three models:

```bash
./les_apriori iso 1 --input ../data/isotropic_256_t1.h5 --tag iso_t1
```

This writes, under `../results/`:

```
iso_t1_delta1_smagorinsky.h5
iso_t1_delta1_wale.h5
iso_t1_delta1_dynamic.h5
```

Each file contains: filtered velocity ū, exact stress τ_exact, modeled stress
τ_model (6 symmetric components), eddy viscosity field ν_r, dissipation Π,
and scalar metrics (correlation ρ, mean(cos θ), backscatter fraction,
ν_r exponent) as HDF5 attributes.

### 3.3 Full ensemble

`run_full_ensemble.sh` loops over every combination:

```
8 snapshots  ×  3 filter widths (Δ = 4dx, 8dx, 16dx)  ×  3 SGS models
= 72 HDF5 files in ../results/
```

From the project root:

```bash
cd cpp_driver
bash run_full_ensemble.sh
```

A timestamped log is written to `../results/ensemble_run_<timestamp>.log`.

---

## 4. Physics sanity checks performed at run time

With `--debug`, the driver aborts on the first failure:

| Check | Expected | Reference |
|-------|----------|-----------|
| Energy ratio E(ū)/E(u) | < 1 | Pope Eq. 13.22 |
| Smagorinsky backscatter fraction | exactly 0 | Pope Eq. 13.127 (ν_r ≥ 0) |
| Dynamic Smagorinsky Cs² < 0 fraction | between 30 % and 40 % | Germano 1991 |
| Corr(τ_exact, τ_exact) | 1.0 | definition |
| mean(Π_exact) | > 0 | Pope p. 589 (forward cascade) |

A non-zero exit from these means the numerics or physics is broken and later
plots cannot be trusted.

---

## 5. Expected run time

On a M2-class laptop (single-threaded except for FFTW plans):

| Task | Wall time |
|------|-----------|
| One `./les_apriori iso 1` invocation | 40–70 s |
| Full 8-snapshot ensemble             | 30–50 min |

FFTW plan wisdom is cached, so the second and later runs are 20–30 % faster.
