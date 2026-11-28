#!/bin/bash
# Tier 3 E4: 10 Mb of chr21 (20-30 Mb), 4 haplotypes + GRCh38, monolithic vs 400 kb / 1 Mb / 2.5 Mb zero-overlap chunks.
# As run on 2026-11-24 (paths are that machine's; see docs/TIER3.md). Resumable.
set -euo pipefail
export PIPELINE_CONFIG=/home/user/linear-scaling-pangenome-sv/config/pipeline.yaml   # seqwish -B 100000000
export PATH=/home/user/work/bin:$PATH
R=/home/user/linear-scaling-pangenome-sv; T=$R/pipeline/benchmark/tier3; W=/home/user/work/t3; cd $W
PY=/home/user/work/venv/bin/python
HAPARGS=""
while IFS=$'\t' read -r reg s h asm contig strand; do
  [ $reg = r21x ] && HAPARGS="$HAPARGS --hap $s#$h=data/r21x.$asm.paf=data/$asm.$contig.fa"
done < contigs.tsv
for k in mono windows0 windows1m windows2500k; do
  [ -f fa/r21x/$k/projection.tsv ] && continue
  mkdir -p fa/r21x/$k
  python3 $R/pipeline/prepare/project_chunks.py manifests/r21x.$k.tsv data/r21x.fa fa/r21x/$k --region-start 20000000 --ref-name "GRCh38#0#chr21" $HAPARGS | tail -1
done
D=runs/r21x/rep1
for d in $D/*/chunk_* $D/monolithic; do [ -d $d ] && [ ! -f $d/final.gfa ] && rm -rf $d; done
[ -f $D/comparison.json ] || { echo "r21x start $(date -u +%FT%TZ)"; python3 $R/pipeline/benchmark/compare_chunking.py $D --mono fa/r21x/mono/_mono.fa \
  --region-start 20000000 --chrom chr21 --threads 4 --strategy windows0=manifests/r21x.windows0.tsv=fa/r21x/windows0 \
  --strategy windows1m=manifests/r21x.windows1m.tsv=fa/r21x/windows1m --strategy windows2500k=manifests/r21x.windows2500k.tsv=fa/r21x/windows2500k | tail -4; }
for g in monolithic/final.gfa windows0/merged.gfa windows1m/merged.gfa windows2500k/merged.gfa; do
  out=$D/truth_$(dirname $g).json
  [ -f $out ] || $PY $T/score_vs_truth.py $D/$g $out --chrom chr21 --region-start 20000000 --ref ref/chr21.chr.fa \
    --truth truth/r21x --haps HG00673#1 HG00673#2 HG00733#1 HG00733#2 --segdups ref/segdups.chr.bed --truvari /home/user/work/venv/bin/truvari
done
echo "e4 done $(date -u +%FT%TZ)"
