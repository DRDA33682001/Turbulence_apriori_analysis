#!/usr/bin/env python3
"""
run_pipeline.py — End-to-end orchestrator for the post-processing stage.

Assumes the C++ ensemble has already been run (so `cpp_driver/results/` is
populated with 72 HDF5 files). Calls, in order:
    1. plot_figures.py   — single-snapshot / diagnostic plots
    2. plot_ensemble.py  — ensemble figures + summary tables

Use this when you want every figure regenerated from a single command.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


def _run(cmd: list[str], label: str) -> bool:
    print(f"\n{'=' * 70}\n  {label}\n{'=' * 70}")
    print(f"> {' '.join(cmd)}")
    result = subprocess.run(cmd)
    if result.returncode != 0:
        print(f"\n[FAIL] {label}  (exit {result.returncode})")
        return False
    return True


def main() -> int:
    here = Path(__file__).resolve().parent
    root = here.parent

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results", type=Path,
                    default=root / "cpp_driver" / "results",
                    help="Directory of C++ HDF5 result files.")
    ap.add_argument("--data", type=Path,
                    default=root / "cpp_driver" / "data",
                    help="Directory of DNS velocity HDF5 files.")
    ap.add_argument("--figures", type=Path,
                    default=root / "figures",
                    help="Top-level figures directory.")
    ap.add_argument("--snapshot", default="t1",
                    help="Snapshot used by single-snapshot figures.")
    args = ap.parse_args()

    args.figures.mkdir(parents=True, exist_ok=True)
    (args.figures / "ensemble").mkdir(parents=True, exist_ok=True)

    py = sys.executable

    ok = _run(
        [py, str(here / "plot_figures.py"),
         "--results", str(args.results),
         "--data", str(args.data),
         "--figures", str(args.figures),
         "--snapshot", args.snapshot],
        "Stage 1/2 — single-snapshot plots",
    )
    if not ok:
        return 1

    ok = _run(
        [py, str(here / "plot_ensemble.py"),
         "--results", str(args.results),
         "--figures", str(args.figures / "ensemble")],
        "Stage 2/2 — ensemble plots + summary tables",
    )
    if not ok:
        return 1

    print(f"\nDone — figures in {args.figures}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
