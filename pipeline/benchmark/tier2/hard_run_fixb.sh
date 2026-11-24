#!/bin/bash
# Tier 2 as run on 2026-11-24 (paths are that machine's; see docs/TIER2.md)
# Rerun with an explicit numeric seqwish batch (-B) so transitive closure runs in one batch
set -euo pipefail
export PIPELINE_CONFIG=/home/user/linear-scaling-pangenome-sv/config/pipeline.yaml
export PGGB_TRANSCLOSE_BATCH=100000000
R=/home/user/linear-scaling-pangenome-sv
cd /home/user/work/fixb
# 1-2) chr21 monolithic checks
for spec in N5:/home/user/work/fa/mono N40:/home/user/work/scale/fa/N40/mono; do
  n=${spec%%:*}; d=${spec#*:}
  [ -f chr21_$n/final.gfa ] || { echo "chr21 $n mono start $(date -u +%FT%TZ)"; python3 $R/scripts/run_pggb.py $d/_mono.fa $PWD/chr21_$n 4 > /dev/null; }
done
# 3-5) hard regions
cd /home/user/work/hard
for job in r17:1 r22:1 r17:2 r22:2; do
  reg=${job%%:*}; r=${job#*:}
  case $reg in r17) st=45400000; ch=chr17;; r22) st=18600000; ch=chr22;; esac
  W=runs_fixb/$reg/rep$r
  if [ ! -f $W/comparison.json ]; then
    for d in $W/*/chunk_* $W/monolithic; do [ -d $d ] && [ ! -f $d/final.gfa ] && rm -rf $d; done
    echo "$reg rep$r start $(date -u +%FT%TZ)"
    python3 $R/pipeline/benchmark/compare_chunking.py $W --mono fa/$reg/mono/_mono.fa --region-start $st --chrom $ch \
      --strategy windows0=$reg.windows0.tsv=fa/$reg/windows0 --strategy hotspots300=$reg.hotspots300.tsv=fa/$reg/hotspots300 --threads 4
  fi
  [ -f $W/seam_metrics.json ] || python3 $R/pipeline/benchmark/seam_metrics.py $W --region-start $st \
    --manifests windows0=$reg.windows0.tsv hotspots300=$reg.hotspots300.tsv --truvari /opt/truvari/bin/truvari > $W/seam.log 2>&1
done
echo "done $(date -u +%FT%TZ)"
