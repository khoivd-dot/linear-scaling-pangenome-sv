#!/bin/bash
# Sensitivity of the truth scores to vcfbub's bubble-size cut (-a 100000, the PGGB pipeline's, against -a 10000).
# Why: when a haplotype path, GRCh38 included, revisits nodes, vg deconstruct can report a short top-level
# bubble whose nested variants span hundreds of kb; vcfbub -l 0 then drops all of them (17q21.31 at -p 98).
# With -a 10000 such a parent is popped and its nested variants are kept. As run on 2026-11-25. Resumable.
set -euo pipefail
R=/home/user/linear-scaling-pangenome-sv; T=$R/pipeline/benchmark/tier3; W=/home/user/work/t3; cd $W
PY=/home/user/work/venv/bin/python
for job in "r17 chr17 45400000" "r22 chr22 18600000"; do
  set -- $job
  for D in runs/$1/rep1 runs_p98/$1/rep1; do
    for g in monolithic/final.gfa windows0/merged.gfa sdgroup/laced.gfa; do
      [ -f $D/$g ] || continue
      out=$D/truth_$(dirname $g)_a10k.json
      [ -f $out ] || $PY $T/score_vs_truth.py $D/$g $out --chrom $2 --region-start $3 --ref ref/$2.chr.fa \
        --truth truth/$1 --haps HG00673#1 HG00673#2 HG00733#1 HG00733#2 --segdups ref/segdups.chr.bed \
        --truvari /home/user/work/venv/bin/truvari --bub "-l 0 -a 10000"
    done
  done
done
echo "bub sensitivity done $(date -u +%FT%TZ)"
