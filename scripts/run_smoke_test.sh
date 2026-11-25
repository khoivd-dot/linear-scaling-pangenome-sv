#!/usr/bin/env bash
# run_smoke_test.sh — Run the full smoke test pipeline.
# Each command is logged with a timestamp on start, end and failure; the first
# failure stops the run with that command's exit code.
#
# Chunk workflow overrides (unset = config/pipeline.yaml, now 400000 / 50000 / 8):
#   CHUNK_SIZE_BP=200000 CHUNK_OVERLAP_BP=0 PGGB_THREADS=16 bash scripts/run_smoke_test.sh
# PGGB_THREADS applies to the baseline build as well as the chunks.
set -euo pipefail

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

run() {
    local t0=$SECONDS rc=0
    log "START $*"
    "$@" || rc=$?
    if [ "$rc" -ne 0 ]; then
        log "FAILED (exit $rc) after $((SECONDS - t0))s: $*"
        exit "$rc"
    fi
    log "END $((SECONDS - t0))s: $*"
}

echo "=== Smoke Test Pipeline ==="
echo "Target: chr21:20000000-21000000"
echo ""

# Step 1: Check environment
echo "[1/8] Checking environment..."
run bash scripts/check_environment.sh

# Step 2: Fetch HPRC index
echo "[2/8] Fetching HPRC index..."
run python3 scripts/fetch_hprc_index.py

# Step 3: Download assemblies
echo "[3/8] Downloading assemblies..."
run python3 scripts/download_hprc.py

# Step 4: Prepare chr21
echo "[4/8] Preparing chr21..."
run python3 pipeline/prepare/prepare_sequences.py

# Step 5: Baseline PGGB graph
echo "[5/8] Building baseline graph..."
run bash pipeline/baseline/build_baseline.sh

# Step 6: Parallel chunks
echo "[6/8] Creating chunks..."
run python3 pipeline/parallel/make_chunks.py
run python3 pipeline/parallel/build_all_chunks.py

# Step 7: Merge
echo "[7/8] Merging graphs..."
run python3 pipeline/merge/merge_graphs.py

# Step 8: Benchmark
echo "[8/8] Running benchmark..."
run python3 pipeline/benchmark/graph_stats.py
run python3 pipeline/benchmark/compare_paths.py
run python3 pipeline/benchmark/build_report.py

echo ""
echo "=== Smoke Test Complete ==="
echo "Results in results/"
echo "Run: python3 scripts/sync_web_results.py"
