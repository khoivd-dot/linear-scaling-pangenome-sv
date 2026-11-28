#!/bin/bash
# Tier 3 truth calls, as run on 2026-11-24: one haploid VCF + callable BED per region and haplotype
# from the HPRC release 2 winnowmap GRCh38_no_alt alignments (hap1 = paternal = PanSN #1).
set -euo pipefail
W=/home/user/work/t3; cd $W; mkdir -p truth
export CURL_CA_BUNDLE=/root/.ccr/ca-bundle.crt
T=$(dirname $(readlink -f $0))
while read -r reg chrom start end; do
  for s in HG00673:HPRC HG00733:HPRC_PLUS; do
    n=${s%%:*}; b=${s#*:}
    for h in 1 2; do
      out=truth/${reg}_${n}_$h
      [ -s $out.vcf ] && continue
      /home/user/work/venv/bin/python $T/truth_calls.py \
        https://s3-us-west-2.amazonaws.com/human-pangenomics/working/$b/$n/assemblies/release2/alignments/winnowmap/GRCh38_no_alt/$n.GRCh38_no_alt.hap$h.sorted.bam \
        ref/$chrom.fa $chrom $start $end $out --sample $n
    done
  done
done < $T/regions.tsv
