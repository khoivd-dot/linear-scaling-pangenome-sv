#!/bin/bash
# Tier 3 shared inputs, as run on 2026-11-24 (paths are that machine's; see docs/TIER3.md):
# minigraph-cactus GRCh38 chains for the two pilot samples, Ensembl GRCh38 chr17/21/22, GIAB v3.3 segdups,
# contig choice per region, and chunk manifests (zero-overlap windows, duplication-aware groups).
set -euo pipefail
R=/home/user/linear-scaling-pangenome-sv; T=$R/pipeline/benchmark/tier3; W=/home/user/work/t3
PY=${PY:-/home/user/work/venv/bin/python}
mkdir -p $W/chains $W/ref $W/manifests; cd $W
S3=https://s3-us-west-2.amazonaws.com/human-pangenomics/working
for s in HG00673:HPRC HG00733:HPRC_PLUS; do n=${s%%:*}; b=${s#*:}
  for h in mat pat; do f=${n}_${h}_hprc_r2_v1.0.1_vs_GRCh38.chain.gz
    [ -s chains/$f ] || curl -s -f -o chains/$f $S3/$b/$n/assemblies/release2/annotation/chains/minigraph-cactus/$f; done; done
for c in 17 21 22; do
  [ -s ref/chr$c.fa ] || curl -s https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.$c.fa.gz | gunzip > ref/chr$c.fa
  [ -s ref/chr$c.chr.fa ] || sed "1s/.*/>chr$c/" ref/chr$c.fa > ref/chr$c.chr.fa
  $PY -c "import pysam; pysam.faidx('ref/chr$c.fa'); pysam.faidx('ref/chr$c.chr.fa')"
done
[ -s ref/segdups.chr.bed ] || curl -s https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.3/GRCh38@all/SegmentalDuplications/GRCh38_segdups.bed.gz \
  | zcat | awk '$1=="chr17"||$1=="chr21"||$1=="chr22"' > ref/segdups.chr.bed
python3 $T/pick_contigs.py $T/regions.tsv chains HG00673 HG00733 > contigs.tsv
mk() {  # name chrom start end size
  python3 -c "import sys; sys.path.insert(0,'$R')
from pipeline.parallel.make_chunks import create_chunks, write_manifest
write_manifest(create_chunks('$2',$3,$4,$5,0),'manifests/$1.tsv')"; }
for x in "r17 chr17 45400000 46600000 400000" "r22 chr22 18600000 21500000 414286" "r21 chr21 20000000 21000000 333334" "r21x chr21 20000000 30000000 400000"; do
  set -- $x; mk $1.windows0 $2 $3 $4 $5; mk $1.mono $2 $3 $4 $(($4-$3)); sed -i 's/chunk_0001/_mono/' manifests/$1.mono.tsv
done
mk r21x.windows1m chr21 20000000 30000000 1000000; mk r21x.windows2500k chr21 20000000 30000000 2500000
# haploblocks.org cuts with a 300 kb minimum on chr21 (Tier 1 manifest)
printf 'chunk_id\tchrom\treference_start\treference_end\tcore_start\tcore_end\toverlap_left\toverlap_right\tpairwise_overlap_bp\tstatus\nchunk_0001\tchr21\t20000000\t20347668\t20000000\t20347668\t0\t0\t0\tpending\nchunk_0002\tchr21\t20347668\t20650101\t20347668\t20650101\t0\t0\t0\tpending\nchunk_0003\tchr21\t20650101\t21000000\t20650101\t21000000\t0\t0\t0\tpending\n' > manifests/r21.hotspots300.tsv
# region FASTAs come from fetch_align.sh; duplication-aware groups need them
bash $T/fetch_align.sh
for r in r17:45400000 r22:18600000; do
  $PY $T/sd_groups.py data/${r%%:*}.fa manifests/${r%%:*}.windows0.tsv --region-start ${r#*:} \
    > manifests/${r%%:*}.sdgroups.tsv 2> manifests/${r%%:*}.sdpairs.tsv
done
