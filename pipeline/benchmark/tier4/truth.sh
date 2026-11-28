#!/bin/bash
# Tier 4 truth sets for the scored haplotypes, per region:
#   HG00673, HG00733  dipcall-style calls from the HPRC winnowmap GRCh38 alignments (tier3/truth_calls.py)
#   HG002             GIAB HG002 Q100 v5.0q benchmark (giab_truth.py); SV scoring uses its stvar regions
# Writes $W/truth/<region>_<sample>_<hap>.{vcf,bed} and $W/truthsv/<region>_<sample>_<hap>.bed.
#   truth.sh REGION [REGION ...]      (needs fetch_ref.sh for each chromosome first)
set -euo pipefail
T=$(dirname "$(readlink -f "$0")"); T3=$T/../tier3; PY=${PY:-python3}; W=${W:-$PWD/work}
mkdir -p "$W/truth" "$W/truthsv" "$W/giab"; cd "$W"
G=https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/v5.0q
for t in smvar stvar; do
  for f in vcf.gz vcf.gz.tbi benchmark.bed; do
    [ -s giab/$t.$f ] || curl -sS -f --retry 6 --retry-all-errors --retry-delay 10 -o giab/$t.$f $G/HG002_GRCh38_v5.0q_$t.$f
  done
done
S3=https://s3-us-west-2.amazonaws.com/human-pangenomics/working
for reg in "$@"; do
  read -r chrom start end < <(python3 "$T/sets.py" region "$reg")
  for s in HG00673:HPRC HG00733:HPRC_PLUS; do
    n=${s%%:*}; b=${s#*:}
    for h in 1 2; do
      out=truth/${reg}_${n}_$h
      [ -s $out.vcf ] || $PY "$T3/truth_calls.py" \
        $S3/$b/$n/assemblies/release2/alignments/winnowmap/GRCh38_no_alt/$n.GRCh38_no_alt.hap$h.sorted.bam \
        ref/$chrom.ens.fa $chrom $start $end $out --sample $n
      cp $out.bed truthsv/${reg}_${n}_$h.bed
    done
  done
  $PY "$T/giab_truth.py" giab/smvar.vcf.gz giab/stvar.vcf.gz giab/smvar.benchmark.bed giab/stvar.benchmark.bed \
    $chrom $start $end truth/$reg truthsv/$reg
  rm -f truthsv/${reg}_HG002_*.vcf
done
