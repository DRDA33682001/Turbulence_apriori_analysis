# Data Download

This folder contains the Python pipeline that fetches DNS velocity snapshots
from the **Johns Hopkins Turbulence Database (JHTDB)** and writes them to local
HDF5 files. These HDF5 files are the input to the C++ driver in
`../cpp_driver/`.

> **Why is the data not in this repository?**
> Each snapshot is ~768 MB (256³ grid × 3 velocity components × 4 time snapshots
> × 2 datasets ≈ 6 GB total). GitHub blocks files larger than 100 MB, and the
> JHTDB terms of use do not allow redistribution. You must download the data
> yourself with the instructions below.

---

## 1. What gets downloaded

| Dataset | File name | Size | Description |
|---------|-----------|------|-------------|
| Forced Isotropic (Re<sub>λ</sub>=433) | `isotropic_256_t{1..4}.h5` | ~768 MB each | 1024³ native, strided 4× to 256³, 4 time snapshots |
| Channel Flow (Re<sub>τ</sub>=1000)    | `channel_256_t{1..4}.h5`   | ~768 MB each | 2048×512×2048 native, strided to 256³, 4 snapshots |

Each HDF5 file contains three datasets — `u1`, `u2`, `u3` (x, y, z velocity
components) — on a uniform 256³ grid.

Total download: **~6 GB** (8 files).

---

## 2. Prerequisites

* Python 3.10 or newer
* A **JHTDB authorization token**. Register for free at
  <https://turbulence.pha.jhu.edu/authtoken.aspx>. Paste the token into
  `download.py` (replace the placeholder `"you need to get your api key and put it here"` 
  for the variable `JHTDB_TOKEN` at the top of the file) or export it
  as the environment variable `JHTDB_TOKEN` before running.

---

## 3. Install dependencies

```bash
cd data_download
python -m venv .venv
source .venv/bin/activate       # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

`requirements.txt` pins `numpy`, `scipy`, `h5py`, `matplotlib`, `seaborn`, and
`giverny` (the official JHTDB Python client).

---

## 4. Run the download

The script creates `../cpp_driver/data/` relative to the C++ driver (or
whichever `--output` path you pass) and drops the HDF5 files there.

```bash
# From the data_download/ folder
mkdir -p ../cpp_driver/data
python download.py --output ../cpp_driver/data/
```

The script will:

1. Issue a single `getCutout` call per snapshot (no chunk stitching — see the
   note in `download.py` about a historical bug from 8-chunk stitching).
2. Validate each returned array (shape, dtype, finite values, physical sanity
   such as RMS ratio across octants).
3. Atomically save to `<name>.h5.tmp` and rename on success, so a crash mid-run
   never overwrites a good file.

A full run takes about **10–20 minutes** on a typical home connection and
depends on the JHTDB queue.

---

## 5. Verifying the download

After the download finishes, you should see:

```
../cpp_driver/data/
    isotropic_256_t1.h5
    isotropic_256_t2.h5
    isotropic_256_t3.h5
    isotropic_256_t4.h5
    channel_256_t1.h5
    channel_256_t2.h5
    channel_256_t3.h5
    channel_256_t4.h5
```

A quick sanity check with `h5ls`:

```bash
h5ls ../cpp_driver/data/isotropic_256_t1.h5
# u1   Dataset {256, 256, 256}
# u2   Dataset {256, 256, 256}
# u3   Dataset {256, 256, 256}
```

---

## 6. Notes on physical setup

* **Isotropic dataset** — triply periodic cube, domain L = 2π, dx = 2π/1024.
  The 256³ sub-cube we extract preserves the first three octaves of the energy
  spectrum (k ≤ 128).
* **Channel dataset** — streamwise (x) and spanwise (z) are periodic with
  stride 4; wall-normal (y) uses stride 2 to keep enough resolution near the
  wall (y⁺ ≈ 1000 at channel centerline). Reynolds number Re<sub>τ</sub>=1000.

These conventions are required by the C++ driver's filter widths and model
constants. Changing the grid size or stride requires recompiling
`../cpp_driver/` with new `config.h` values.
