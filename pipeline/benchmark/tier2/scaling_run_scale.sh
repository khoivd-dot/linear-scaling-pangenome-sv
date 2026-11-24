#!/bin/bash
# Tier 2 scaling, as run on 2026-11-24 (paths are that machine's; see docs/TIER2.md)
# Tier 2 scaling driver: resumable (compare_chunking skips finished builds)
set -euo pipefail
export PIPELINE_CONFIG=/home/user/linear-scaling-pangenome-sv/config/pipeline.yaml
# replicate 1 ran before the -B fix; keep pggb's own default for both replicates (chr21 1 Mb is unaffected, see docs/TIER2.md)
export PGGB_TRANSCLOSE_BATCH=pggb-default
cd /home/user/work/scale
R=/home/user/linear-scaling-pangenome-sv; M=/home/user/work/manifests
for job in N10:1 N20:1 N40:1 N10:2 N20:2 N40:2; do
  N=${job%%:*}; r=${job#*:}
  [ -f runs/$N/rep$r/comparison.json ] && continue
  for d in runs/$N/rep$r/*/chunk_* runs/$N/rep$r/monolithic; do [ -d $d ] && [ ! -f $d/final.gfa ] && rm -rf $d; done
  echo "$N rep$r start $(date -u +%FT%TZ)"
  python3 $R/pipeline/benchmark/compare_chunking.py runs/$N/rep$r --mono fa/$N/mono/_mono.fa \
    --strategy windows0=$M/windows0.tsv=fa/$N/windows0 --strategy hotspots300=$M/hotspots300k.tsv=fa/$N/hotspots300 \
    --strategy hotspots50=$M/hotspots.tsv=fa/$N/hotspots50 --threads 4
  python3 $R/pipeline/benchmark/seam_metrics.py runs/$N/rep$r --manifests windows0=$M/windows0.tsv \
    hotspots300=$M/hotspots300k.tsv hotspots50=$M/hotspots.tsv --truvari /opt/truvari/bin/truvari > runs/$N/rep$r/seam.log 2>&1
done
echo "done $(date -u +%FT%TZ)"
