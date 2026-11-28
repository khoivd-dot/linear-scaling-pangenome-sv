#!/bin/bash
# GRCh38 chromosomes from Ensembl release 112 (as Tier 3): $W/ref/<chrom>.fa named "chrN" and
# $W/ref/<chrom>.ens.fa named "N" (truth_calls.py reads Ensembl names), plus GIAB v3.3 segdups.
#   fetch_ref.sh chr21 [chr17 ...]
set -euo pipefail
PY=${PY:-python3}; W=${W:-$PWD/work}; mkdir -p "$W/ref"; cd "$W/ref"
for c in "$@"; do
  n=${c#chr}
  [ -s $c.ens.fa ] || curl -sS -f --retry 6 --retry-all-errors --retry-delay 10 https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.$n.fa.gz | gunzip > $c.ens.fa
  [ -s $c.fa ] || sed "1s/.*/>$c/" $c.ens.fa > $c.fa
  $PY -c "import pysam; pysam.faidx('$c.fa'); pysam.faidx('$c.ens.fa')"
done
[ -s segdups.bed ] || curl -sS -f --retry 6 --retry-all-errors --retry-delay 10 https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.3/GRCh38@all/SegmentalDuplications/GRCh38_segdups.bed.gz \
  | zcat | awk '$1=="chr17"||$1=="chr21"||$1=="chr22"' > segdups.bed
