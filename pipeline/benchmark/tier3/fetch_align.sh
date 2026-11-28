#!/bin/bash
# Tier 3 inputs, as run on 2026-11-24 (paths are that machine's; see docs/TIER3.md)
# Region FASTAs from Ensembl GRCh38, each haplotype's spanning contig fetched remotely, region aligned to it.
# minimap2 here is a mappy (minimap2 2.31 Python binding) stand-in with the same -x asm5 -c --eqx output.
set -euo pipefail
W=/home/user/work/t3; cd $W; mkdir -p data
IMG=ghcr.io/pangenome/pggb@sha256:44f3563c5c3e2032c293fe309eed178b754021b07963ca4857e4ae5efda86474
T=$(dirname $(readlink -f $0))
while read -r reg chrom start end; do
  [ -s data/$reg.fa ] || samtools faidx ref/$chrom.fa ${chrom#chr}:$((start+1))-$end | sed "1s/.*/>GRCh38#0#$chrom:$start-$end/" > data/$reg.fa
done < $T/regions.tsv
while IFS=$'\t' read -r reg s h asm contig strand; do
  case $s in HG00733) base=HPRC_PLUS;; *) base=HPRC;; esac
  url=https://s3-us-west-2.amazonaws.com/human-pangenomics/working/$base/$s/assemblies/release2/$asm.fa.gz
  if [ ! -s data/$asm.$contig.fa ]; then
    name=$(curl -s $url.fai | cut -f1 | grep "#${contig}$")
    docker run --rm --network host -e HTTPS_PROXY=$HTTPS_PROXY -e https_proxy=$HTTPS_PROXY -e CURL_CA_BUNDLE=/ca.crt \
      -v /root/.ccr/ca-bundle.crt:/ca.crt:ro -v $W/data:/d -w /d $IMG bash -c \
      "samtools faidx '$url' '$name' > $asm.$contig.tmp && mv $asm.$contig.tmp $asm.$contig.fa && samtools faidx $asm.$contig.fa"
  fi
  out=data/$reg.$asm.paf
  [ -s $out ] || { minimap2 -x asm5 -c --eqx -t 4 data/$reg.fa data/$asm.$contig.fa > $out.tmp && mv $out.tmp $out; }
  echo "$reg $asm $(wc -l < $out) alignments"
done < $W/contigs.tsv
