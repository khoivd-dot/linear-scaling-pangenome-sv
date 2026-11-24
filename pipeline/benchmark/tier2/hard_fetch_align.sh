#!/bin/bash
# Tier 2 hard regions, as run on 2026-11-24 (paths are that machine's; see docs/TIER2.md)
set -euo pipefail
cd /home/user/work/hard
IMG=ghcr.io/pangenome/pggb@sha256:44f3563c5c3e2032c293fe309eed178b754021b07963ca4857e4ae5efda86474
reg=$1; s=$2; n=$4; contig=$5
case $s in HG00733) base=HPRC_PLUS;; *) base=HPRC;; esac
url=https://s3-us-west-2.amazonaws.com/human-pangenomics/working/$base/$s/assemblies/release2/$n.fa.gz
out=$reg.$n
[ -s $out.paf ] && exit 0
name=$(curl -s $url.fai | cut -f1 | grep "#${contig}$")
docker run --rm --network host -e HTTPS_PROXY=$HTTPS_PROXY -e https_proxy=$HTTPS_PROXY -e CURL_CA_BUNDLE=/ca.crt \
  -v /root/.ccr/ca-bundle.crt:/ca.crt:ro -v $PWD:/d -w /d $IMG bash -c "samtools faidx '$url' '$name' > $out.contig.fa && samtools faidx $out.contig.fa"
minimap2 -x asm5 -c --eqx -t 2 $reg.fa $out.contig.fa > $out.paf.tmp 2>/dev/null && mv $out.paf.tmp $out.paf
