#!/bin/bash
# Tier 2 scaling, as run on 2026-11-24 (paths are that machine's; see docs/TIER2.md)
# project every chunking for N = 10, 20, 40 haplotypes (nested sample sets, order from haps.tsv)
set -euo pipefail
cd /home/user/work/scale
R=/home/user/linear-scaling-pangenome-sv; M=/home/user/work/manifests
for N in 10 20 40; do
  HAPS=""
  while IFS=$'\t' read -r i s h n c u; do
    [ $(( (i+1)*2 )) -le $N ] || continue
    HAPS="$HAPS --hap $s#$h=data/$n.paf=data/$n.contig.fa"
  done < haps.tsv
  for spec in mono=$M/mono.tsv windows0=$M/windows0.tsv hotspots300=$M/hotspots300k.tsv hotspots50=$M/hotspots.tsv; do
    name=${spec%%=*}; man=${spec#*=}
    [ -f fa/N$N/$name/projection.tsv ] && continue
    mkdir -p fa/N$N/$name
    python3 $R/pipeline/prepare/project_chunks.py $man /home/user/work/data/region.fa fa/N$N/$name \
      --region-start 20000000 --ref-name 'GRCh38#0#chr21' $HAPS | tail -1
  done
done
