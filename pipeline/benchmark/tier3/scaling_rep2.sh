#!/bin/bash
# Tier 3 E3: second replicate of the Tier 2 scaling runs (N = 10, 20, 40), as run on 2026-11-24.
# Same 20 samples (select_samples.py, seed 1: 229 eligible, identical order), same manifests, pggb's own -B as in rep1.
# Usage: scaling_rep2.sh fetch | project | build
set -euo pipefail
export PATH=/home/user/work/bin:$PATH
R=/home/user/linear-scaling-pangenome-sv; W=/home/user/work/scale; cd $W
IMG=ghcr.io/pangenome/pggb@sha256:44f3563c5c3e2032c293fe309eed178b754021b07963ca4857e4ae5efda86474
REGION=/home/user/work/t3/data/r21.fa
case $1 in
fetch1)  # one haps.tsv row: i s h n c u
  n=$5; c=$6; u=$7
  [ -s data/$n.paf ] && exit 0
  t3=/home/user/work/t3/data/$n.$c.fa
  if [ -s $t3 ]; then ln -sf $t3 data/$n.contig.fa; ln -sf $t3.fai data/$n.contig.fa.fai
  elif [ ! -s data/$n.contig.fa ]; then
    name=$(curl -s $u.fai | cut -f1 | grep "#${c}$")
    docker run --rm --network host -e HTTPS_PROXY=$HTTPS_PROXY -e https_proxy=$HTTPS_PROXY -e CURL_CA_BUNDLE=/ca.crt \
      -v /root/.ccr/ca-bundle.crt:/ca.crt:ro -v $W/data:/d -w /d $IMG bash -c \
      "samtools faidx '$u' '$name' > $n.tmp && mv $n.tmp $n.contig.fa && samtools faidx $n.contig.fa"
  fi
  minimap2 -x asm5 -c --eqx $REGION data/$n.contig.fa > data/$n.paf.tmp && mv data/$n.paf.tmp data/$n.paf
  echo "$n aligned";;
fetch)
  awk '{print "fetch1",$0}' haps.tsv | xargs -P ${P:-2} -L 1 bash $0;;
project)
  for N in 10 20 40; do
    HAPS=""
    while IFS=$'\t' read -r i s h n c u; do
      [ $(( (i+1)*2 )) -le $N ] && HAPS="$HAPS --hap $s#$h=data/$n.paf=data/$n.contig.fa"
    done < haps.tsv
    for name in mono windows0 hotspots300k hotspots; do
      [ -f fa/N$N/$name/projection.tsv ] && continue
      mkdir -p fa/N$N/$name
      python3 $R/pipeline/prepare/project_chunks.py manifests/$name.tsv $REGION fa/N$N/$name \
        --region-start 20000000 --ref-name 'GRCh38#0#chr21' $HAPS | tail -1
    done
  done;;
build)
  export PIPELINE_CONFIG=$R/config/pipeline.yaml PGGB_TRANSCLOSE_BATCH=pggb-default
  for N in N10 N20 N40; do
    D=runs/$N/rep2
    [ -f $D/comparison.json ] && [ -f $D/seam_metrics.json ] && continue
    for d in $D/*/chunk_* $D/monolithic; do [ -d $d ] && [ ! -f $d/final.gfa ] && rm -rf $d; done
    echo "$N rep2 start $(date -u +%FT%TZ)"
    [ -f $D/comparison.json ] || python3 $R/pipeline/benchmark/compare_chunking.py $D --mono fa/$N/mono/_mono.fa \
      --strategy windows0=manifests/windows0.tsv=fa/$N/windows0 --strategy hotspots300=manifests/hotspots300k.tsv=fa/$N/hotspots300k \
      --strategy hotspots50=manifests/hotspots.tsv=fa/$N/hotspots --threads 4 | tail -3
    /home/user/work/venv/bin/python $R/pipeline/benchmark/seam_metrics.py $D --manifests windows0=manifests/windows0.tsv \
      hotspots300=manifests/hotspots300k.tsv hotspots50=manifests/hotspots.tsv --truvari /home/user/work/venv/bin/truvari > $D/seam.log 2>&1
  done
  echo "scaling rep2 done $(date -u +%FT%TZ)";;
esac
