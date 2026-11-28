#!/bin/bash
# Tier 4 prep for one sample: fetch each haplotype's contig, align the GRCh38 region to it
# (minimap2 -x asm5 -c --eqx), and project every chunking's boundaries onto it (coordinates only).
#   prep_sample.sh SET N SAMPLE OUTDIR      e.g. prep_sample.sh arm 20 HG01106 out/prep
# Needs: python with pysam ($PY), minimap2 on PATH, $W/ref/<chrom>.fa (fetch_ref.sh).
set -euo pipefail
SET=$1; N=$2; S=$3; OUT=$4
T=$(dirname "$(readlink -f "$0")"); PY=${PY:-python3}; W=${W:-$PWD/work}
mkdir -p "$OUT" "$W/contigs" "$W/regions"
python3 "$T/sets.py" haps "$SET" "$N" "$S" | while IFS=$'\t' read -r reg prefix asm contig url; do
  read -r chrom start end < <(python3 "$T/sets.py" region "$reg")
  rfa=$W/regions/$reg.fa
  [ -s "$rfa" ] || $PY -c "
import pysam; f = pysam.FastaFile('$W/ref/$chrom.fa')
s = f.fetch(f.references[0], $start, $end).upper()
open('$rfa', 'w').write('>GRCh38#0#$chrom:$start-$end\n' + '\n'.join(s[i:i + 80] for i in range(0, len(s), 80)) + '\n')"
  cfa=$W/contigs/$asm.$contig.fa
  if [ ! -s "$cfa" ]; then
    $PY - "$url" "$contig" "$cfa" <<'PY'
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
with open(out + '.tmp', 'w') as fh:
    fh.write(f'>{name[0]}\n' + '\n'.join(s[i:i + 80] for i in range(0, len(s), 80)) + '\n')
import os; os.rename(out + '.tmp', out)
PY
  fi
  paf=$OUT/$reg.$prefix.paf; paf=${paf//#/_}
  [ -s "$paf" ] || { minimap2 -x asm5 -c --eqx -t "${THREADS:-4}" "$rfa" "$cfa" > "$paf.tmp" && mv "$paf.tmp" "$paf"; }
  margs=()
  for m in "$T"/manifests/$reg.*.tsv; do b=$(basename "$m" .tsv); margs+=(--manifest "${b#*.}=$m"); done
  $PY "$T/project_coords.py" "$paf" "$url" "$prefix" "$OUT/coords.$reg.$S.tsv" --region-start "$start" "${margs[@]}"
done
