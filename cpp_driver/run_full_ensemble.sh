#!/usr/bin/env bash
# run_full_ensemble.sh — Run the full 8-snapshot SGS ensemble sequentially.
#
# Each snapshot is processed at all 3 filter widths × 3 models = 9 result files.
# Total output: 8 × 9 = 72 HDF5 files in les_apriori/results/.
#
# Usage:
#   bash run_full_ensemble.sh            (run all 8 snapshots)
#   bash run_full_ensemble.sh 2>&1 | tee my.log   (custom log destination)

set -euo pipefail

# ---------------------------------------------------------------------------
# Paths (resolved relative to this script's location)
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="${SCRIPT_DIR}/les_apriori/build"
DATA_DIR="${SCRIPT_DIR}/les_apriori/data"
RESULTS_DIR="${SCRIPT_DIR}/les_apriori/results"
BINARY="${BUILD_DIR}/les_apriori"

# ---------------------------------------------------------------------------
# Log file
# ---------------------------------------------------------------------------
mkdir -p "${RESULTS_DIR}"
LOG_FILE="${RESULTS_DIR}/ensemble_run_$(date '+%Y%m%d_%H%M%S').log"
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "=============================================================="
echo "  LES SGS ENSEMBLE RUN"
echo "  Started : $(date)"
echo "  Binary  : ${BINARY}"
echo "  Results : ${RESULTS_DIR}"
echo "  Log     : ${LOG_FILE}"
echo "=============================================================="
echo ""

# ---------------------------------------------------------------------------
# Helper: run one snapshot and print a progress banner
# ---------------------------------------------------------------------------
run_snapshot() {
    local dataset="$1"   # iso | channel
    local tag="$2"       # iso_t1 | channel_t2 | ...
    local input="$3"     # full path to .h5 file

    echo ""
    echo "=============================================================="
    echo "  === Running ${tag} ==="
    echo "  dataset : ${dataset}"
    echo "  input   : ${input}"
    echo "  tag     : ${tag}"
    echo "  time    : $(date)"
    echo "=============================================================="
    echo ""

    (cd "${BUILD_DIR}" && "${BINARY}" "${dataset}" all \
        --input "${input}" \
        --tag   "${tag}")

    echo ""
    echo "  --- ${tag} COMPLETE at $(date) ---"
    echo ""
}

# ---------------------------------------------------------------------------
# Isotropic snapshots (t1-t4)
# ---------------------------------------------------------------------------
echo "##############################################################"
echo "  ISOTROPIC TURBULENCE ENSEMBLE (4 snapshots)"
echo "##############################################################"

for N in 1 2 3 4; do
    run_snapshot "iso" "iso_t${N}" "${DATA_DIR}/isotropic_256_t${N}.h5"
done

# ---------------------------------------------------------------------------
# Channel flow snapshots (t1-t4)
# ---------------------------------------------------------------------------
echo "##############################################################"
echo "  CHANNEL FLOW ENSEMBLE (4 snapshots)"
echo "##############################################################"

for N in 1 2 3 4; do
    run_snapshot "channel" "channel_t${N}" "${DATA_DIR}/channel_256_t${N}.h5"
done

# ---------------------------------------------------------------------------
# Completion summary
# ---------------------------------------------------------------------------
echo ""
echo "=============================================================="
echo "  ENSEMBLE RUN COMPLETE"
echo "  Finished : $(date)"
echo ""
echo "  Result files:"
ls -lh "${RESULTS_DIR}"/*.h5 2>/dev/null | awk '{print "    " $NF "  " $5}' || echo "    (none found)"
echo ""
echo "  Total HDF5 files: $(ls "${RESULTS_DIR}"/*.h5 2>/dev/null | wc -l | tr -d ' ')"
echo "  Log: ${LOG_FILE}"
echo "=============================================================="
