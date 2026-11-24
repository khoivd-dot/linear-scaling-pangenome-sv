#!/bin/bash
# Tier 1 driver as run on 2026-11-24 (paths are that machine's; see docs/TIER1_CONTROLS.md)
set -euo pipefail
export PIPELINE_CONFIG=/home/user/linear-scaling-pangenome-sv/config/pipeline.yaml
cd /home/user/work
S="--strategy windows=manifests/windows.tsv=fa/windows --strategy windows0=manifests/windows0.tsv=fa/windows0 --strategy hotspots300=manifests/hotspots300k.tsv=fa/hotspots300 --strategy hotspots50=manifests/hotspots.tsv=fa/hotspots50"
X="--strategy shift50=manifests/shift50.tsv=fa/shift50 --strategy shift100=manifests/shift100.tsv=fa/shift100"
for i in $(seq -w 1 20); do X="$X --strategy random_$i=manifests/random/random_$i.tsv=fa/random_$i"; done
echo "rep1 start $(date -u +%FT%TZ)"
python3 /home/user/linear-scaling-pangenome-sv/pipeline/benchmark/compare_chunking.py /home/user/work/runs/rep1 --mono fa/mono/_mono.fa $S $X --threads 4
for r in 2 3; do
  echo "rep$r start $(date -u +%FT%TZ)"
  python3 /home/user/linear-scaling-pangenome-sv/pipeline/benchmark/compare_chunking.py /home/user/work/runs/rep$r --mono fa/mono/_mono.fa $S --threads 4
done
echo "done $(date -u +%FT%TZ)"
