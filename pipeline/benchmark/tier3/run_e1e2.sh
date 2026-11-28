#!/bin/bash
# Tier 3 E1 (truth) + E2 (duplication-aware groups), as run on 2026-11-24 (paths are that machine's; see docs/TIER3.md)
# Resumable: finished builds are skipped, interrupted build dirs (no final.gfa) are removed.
set -euo pipefail
export PIPELINE_CONFIG=/home/user/linear-scaling-pangenome-sv/config/pipeline.yaml   # seqwish -B 100000000
export PATH=/home/user/work/bin:$PATH
R=/home/user/linear-scaling-pangenome-sv; T=$R/pipeline/benchmark/tier3; W=/home/user/work/t3; cd $W
TRUVARI=/home/user/work/venv/bin/truvari; PY=/home/user/work/venv/bin/python
HAPS="HG00673#1 HG00673#2 HG00733#1 HG00733#2"
REPS=${REPS:-1}
project() {  # region chrom start strategy
  [ -f fa/$1/$4/projection.tsv ] && return 0
  local HAPARGS="" rg s h asm contig strand
  while IFS=$'\t' read -r rg s h asm contig strand; do
    [ $rg = $1 ] && HAPARGS="$HAPARGS --hap $s#$h=data/$1.$asm.paf=data/$asm.$contig.fa"
  done < contigs.tsv
  mkdir -p fa/$1/$4
  python3 $R/pipeline/prepare/project_chunks.py manifests/$1.$4.tsv data/$1.fa fa/$1/$4 --region-start $3 --ref-name "GRCh38#0#$2" $HAPARGS | tail -1
}
for job in "r21 chr21 20000000 windows0,hotspots300" "r17 chr17 45400000 windows0" "r22 chr22 18600000 windows0"; do
  set -- $job; reg=$1; chrom=$2; st=$3
  for k in mono ${4//,/ }; do project $reg $chrom $st $k; done
  for r in $(seq 1 $REPS); do
    D=runs/$reg/rep$r
    for d in $D/*/chunk_* $D/monolithic $D/*/group_*/build; do [ -d $d ] && [ ! -f $d/final.gfa ] && rm -rf $d; done
    S=""; for k in ${4//,/ }; do S="$S --strategy $k=manifests/$reg.$k.tsv=fa/$reg/$k"; done
    [ -f $D/comparison.json ] || { echo "$reg rep$r chunking start $(date -u +%FT%TZ)"; python3 $R/pipeline/benchmark/compare_chunking.py $D --mono fa/$reg/mono/_mono.fa --region-start $st --chrom $chrom $S --threads 4 | tail -4; }
    # impg lace of the same windows0 chunks (all singletons), and duplication-aware groups
    [ -f manifests/$reg.single.tsv ] || { echo -e "group\tchunk_ids\tlinked_bp"; tail -n +2 manifests/$reg.windows0.tsv | awk '{printf "s_%s\t%s\t0\n",$1,$1}'; } > manifests/$reg.single.tsv
    [ -f $D/lace_windows0/comparison.json ] || $PY $T/group_build.py $D lace_windows0 manifests/$reg.single.tsv fa/$reg/windows0 --reuse windows0 --region-start $st --chrom $chrom
    if [ -f manifests/$reg.sdgroups.tsv ]; then
      [ -f $D/sdgroup/comparison.json ] || { echo "$reg rep$r sdgroup start $(date -u +%FT%TZ)"; $PY $T/group_build.py $D sdgroup manifests/$reg.sdgroups.tsv fa/$reg/windows0 --reuse windows0 --region-start $st --chrom $chrom; }
    fi
    # E1: every graph against the assembly-based truth calls
    for g in monolithic/final.gfa windows0/merged.gfa hotspots300/merged.gfa lace_windows0/laced.gfa sdgroup/laced.gfa; do
      [ -f $D/$g ] || continue
      out=$D/truth_$(dirname $g).json
      [ -f $out ] || $PY $T/score_vs_truth.py $D/$g $out --chrom $chrom --region-start $st --ref ref/$chrom.chr.fa \
        --truth truth/$reg --haps $HAPS --segdups ref/segdups.chr.bed --truvari $TRUVARI
    done
  done
done
echo "done $(date -u +%FT%TZ)"
