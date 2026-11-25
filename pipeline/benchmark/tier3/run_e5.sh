#!/bin/bash
# Tier 3 E5: the same comparison with PGGB's recommended human settings (-p 98 -s 10k -k 47, pggb README)
# instead of the project's -p 90 -s 5k -k 29, on chr21, 17q21.31 and 22q11.2. Monolithic vs windows0.
# As run on 2026-11-25 (paths are that machine's; see docs/TIER3.md). Resumable.
set -euo pipefail
R=/home/user/linear-scaling-pangenome-sv; T=$R/pipeline/benchmark/tier3; W=/home/user/work/t3; cd $W
sed -e 's/minimum_identity: 90/minimum_identity: 98/' -e 's/segment_length: 5000/segment_length: 10000/' \
    -e 's/match_length: 29/match_length: 47/' $R/config/pipeline.yaml > pipeline_p98.yaml
grep -q "minimum_identity: 98" pipeline_p98.yaml && grep -q "segment_length: 10000" pipeline_p98.yaml && grep -q "match_length: 47" pipeline_p98.yaml
export PIPELINE_CONFIG=$W/pipeline_p98.yaml PATH=/home/user/work/bin:$PATH
PY=/home/user/work/venv/bin/python
for job in "r21 chr21 20000000" "r17 chr17 45400000" "r22 chr22 18600000"; do
  set -- $job; D=runs_p98/$1/rep1
  for d in $D/*/chunk_* $D/monolithic; do [ -d $d ] && [ ! -f $d/final.gfa ] && rm -rf $d; done
  [ -f $D/comparison.json ] || { echo "$1 p98 start $(date -u +%FT%TZ)"; python3 $R/pipeline/benchmark/compare_chunking.py $D --mono fa/$1/mono/_mono.fa \
    --region-start $3 --chrom $2 --strategy windows0=manifests/$1.windows0.tsv=fa/$1/windows0 --threads 4 | tail -2; }
  for g in monolithic/final.gfa windows0/merged.gfa; do
    out=$D/truth_$(dirname $g).json
    [ -f $out ] || $PY $T/score_vs_truth.py $D/$g $out --chrom $2 --region-start $3 --ref ref/$2.chr.fa \
      --truth truth/$1 --haps HG00673#1 HG00673#2 HG00733#1 HG00733#2 --segdups ref/segdups.chr.bed --truvari /home/user/work/venv/bin/truvari
  done
done
echo "e5 done $(date -u +%FT%TZ)"
