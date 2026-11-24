#!/bin/bash
# Tier 2 scaling, as run on 2026-11-24 (paths are that machine's; see docs/TIER2.md)
# fetch one haplotype's chr21 contig remotely and align the GRCh38 region to it
set -euo pipefail
cd /home/user/work/scale/data
IMG=ghcr.io/pangenome/pggb@sha256:44f3563c5c3e2032c293fe309eed178b754021b07963ca4857e4ae5efda86474
n=$1; contig=$2; url=$3
[ -s $n.paf ] && exit 0
name=$(curl -s $url.fai | cut -f1 | grep "#${contig}$")
docker run --rm --network host -e HTTPS_PROXY=$HTTPS_PROXY -e https_proxy=$HTTPS_PROXY -e CURL_CA_BUNDLE=/ca.crt \
  -v /root/.ccr/ca-bundle.crt:/ca.crt:ro -v $PWD:/d -w /d $IMG bash -c "samtools faidx '$url' '$name' > $n.contig.fa && samtools faidx $n.contig.fa"
minimap2 -x asm5 -c --eqx -t 2 ../../data/region.fa $n.contig.fa > $n.paf.tmp 2>/dev/null && mv $n.paf.tmp $n.paf
