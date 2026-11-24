#!/bin/bash
# Tier 2 hard regions, as run on 2026-11-24 (paths are that machine's; see docs/TIER2.md)
# 17q21.31 (chr17:45.4-46.6 Mb) and 22q11.2 (chr22:18.6-21.5 Mb), HG00673 + HG00733, both haplotypes.
set -euo pipefail
cd /home/user/work/hard
R=/home/user/linear-scaling-pangenome-sv
for c in 17 22; do curl -s -O https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.$c.fa.gz; done
for c in 17 22; do zcat Homo_sapiens.GRCh38.dna.chromosome.$c.fa.gz > chr$c.fa; done
samtools faidx chr17.fa 17:45400001-46600000 | sed '1s/.*/>GRCh38#0#chr17:45400000-46600000/' > r17.fa
samtools faidx chr22.fa 22:18600001-21500000 | sed '1s/.*/>GRCh38#0#chr22:18600000-21500000/' > r22.fa
# contig per haplotype: the highest-scoring minigraph-cactus chain that spans the region
python3 - <<'PY'
import gzip
regs = {'r17': ('chr17', 45400000, 46600000), 'r22': ('chr22', 18600000, 21500000)}
with open('contigs.tsv', 'w') as out:
    for s in ('HG00673', 'HG00733'):
        for h, hn in (('1', 'pat'), ('2', 'mat')):
            f = f'/home/user/work/chains/{s}_{hn}_hprc_r2_v1.0.1_vs_GRCh38.chain.gz'
            for k, (c, a, b) in regs.items():
                best = None
                for l in gzip.open(f, 'rt'):
                    if not l.startswith('chain'):
                        continue
                    p = l.split()
                    if p[7] != c:
                        continue
                    qs, qe = int(p[10]), int(p[11])
                    if p[9] == '-':
                        qs, qe = int(p[8]) - qe, int(p[8]) - qs
                    if qs <= a and qe >= b and (best is None or int(p[1]) > best[0]):
                        best = (int(p[1]), p[2], p[9])
                out.write(f"{k}\t{s}\t{h}\t{s}_{hn}_hprc_r2_v1.0.1\t{best[1]}\t{best[2]}\n")
PY
cut -f1-6 contigs.tsv | xargs -P 2 -n 6 ./fa.sh   # = hard_fetch_align.sh
# manifests: monolithic, fixed zero-overlap windows, haploblocks.org cuts with min 300 kb
# block_stats.tsv: haploblocks.org bidirectional_blast_samples/block_stats.tsv
for x in "r17 chr17 45400000 46600000 400000" "r22 chr22 18600000 21500000 414286"; do
  set -- $x
  awk -F'\t' -v c=$2 -v s=$3 -v e=$4 '$1==c && $4>s && $3<e {print $2}' block_stats.tsv > $1.blocks.txt
  python3 $R/pipeline/parallel/hotspot_chunks.py $1.blocks.txt $1.hotspots300.tsv --chrom $2 --start $3 --end $4 --min-size 300000
  python3 -c "
import sys; sys.path.insert(0,'$R')
from pipeline.parallel.make_chunks import create_chunks, write_manifest
write_manifest(create_chunks('$2',$3,$4,$5,0),'$1.windows0.tsv')
write_manifest(create_chunks('$2',$3,$4,$4-$3,0),'$1.mono.tsv')"
  sed -i 's/chunk_0001/_mono/' $1.mono.tsv
done
# projection (longest colinear chain, see project_chunks.py)
for x in "r17 chr17 45400000" "r22 chr22 18600000"; do
  set -- $x; HAPS=""
  while IFS=$'\t' read -r reg s h n c st; do [ $reg = $1 ] && HAPS="$HAPS --hap $s#$h=$1.$n.paf=$1.$n.contig.fa"; done < contigs.tsv
  for k in mono windows0 hotspots300; do
    mkdir -p fa/$1/$k
    python3 $R/pipeline/prepare/project_chunks.py $1.$k.tsv $1.fa fa/$1/$k --region-start $3 --ref-name "GRCh38#0#$2" $HAPS
  done
done
