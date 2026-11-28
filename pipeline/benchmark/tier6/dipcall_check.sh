#!/bin/bash
# Tier 6 E10: stand-alone dipcall on HG00673 and HG00733 against GRCh38 chr21, compared with the
# Tier 3 truth calls in one region. Inputs are each haplotype's chr21 contig from tier4/haps_arm.tsv,
# so dipcall sees the same sequence the graphs were built from.
#   dipcall_check.sh REGION OUTDIR      e.g. dipcall_check.sh mid out    (REGION from tier4/sets.py)
# Needs: run-dipcall on PATH (dipcall.kit), python with pysam ($PY), truvari ($TRUVARI), docker (bcftools),
# fetch_ref.sh chr21 and truth.sh REGION already run under $W.
set -euo pipefail
REG=$1; OUT=$2
T=$(dirname "$(readlink -f "$0")"); T4=$T/../tier4; PY=${PY:-python3}; W=${W:-$PWD/work}
read -r CHROM START END < <(python3 "$T4/sets.py" region "$REG")
mkdir -p "$OUT" "$W/dipcall"
for S in HG00673 HG00733; do
  for H in 1 2; do
    fa=$W/dipcall/${S}_$H.fa
    [ -s "$fa" ] && continue
    read -r contig url < <(awk -F'\t' -v s=$S -v h=$H 'NR > 1 && $1 == s && $2 == h {print $4, $6}' "$T4/haps_arm.tsv")
    $PY - "$url" "$contig" "$fa" <<'PY'
import sys, time, pysam
url, contig, out = sys.argv[1:]
for t in range(4):
    try:
        f = pysam.FastaFile(url); break
    except OSError:
        time.sleep(15 * (t + 1))
name = [n for n in f.references if n == contig or n.endswith('#' + contig)]
assert len(name) == 1, (contig, name)
s = f.fetch(name[0])
open(out, 'w').write(f'>{name[0]}\n' + '\n'.join(s[i:i + 80] for i in range(0, len(s), 80)) + '\n')
PY
  done
  ( cd "$W/dipcall" && run-dipcall $S "$W/ref/$CHROM.fa" ${S}_1.fa ${S}_2.fa > $S.mak && make -j2 -f $S.mak > $S.make.log 2>&1 )
  $PY "$T/dipcall_check.py" "$W/dipcall/$S.dip.vcf.gz" "$W/dipcall/$S.dip.bed" "$W/truth/$REG" "$CHROM" "$START" "$END" \
    "$W/ref/$CHROM.fa" "$OUT/dipcall_${REG}_$S.json" --sample $S --truvari "${TRUVARI:-truvari}" --segdups "$W/ref/segdups.bed"
done
