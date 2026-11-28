#!/bin/bash
# Tier 3 E3 inputs, as run on 2026-11-24: redraw the Tier 2 scaling samples (229 eligible, seed 1) and write haps.tsv.
set -euo pipefail
R=/home/user/linear-scaling-pangenome-sv; W=/home/user/work/scale
mkdir -p $W/chains $W/manifests $W/data; cd $W
[ -s index_raw.csv ] || curl -s -f -o index_raw.csv https://raw.githubusercontent.com/human-pangenomics/hprc_intermediate_assembly/refs/heads/main/data_tables/assemblies_release2_v1.0.index.csv
python3 - > chain_urls.txt <<'PY'
import csv
for r in csv.DictReader(open('index_raw.csv')):
    base = r['assembly'].rsplit('/', 1)[0].replace('s3://human-pangenomics/', 'https://s3-us-west-2.amazonaws.com/human-pangenomics/')
    print(f"{base}/annotation/chains/minigraph-cactus/{r['assembly_name']}_vs_GRCh38.chain.gz")
PY
(cd chains && xargs -P 4 -n 1 curl -s -f -O < ../chain_urls.txt || true)
python3 $R/pipeline/benchmark/tier2/select_samples.py chains index_raw.csv sel --chrom chr21 --start 20000000 --end 21000000 \
  --first HG00673 HG00733 --n 20 --seed 1
python3 - <<'PY'
import csv, json
el = json.load(open('sel/eligible.json')); order = json.load(open('sel/order.json'))
rows = {r['assembly_name']: r for r in csv.DictReader(open('index_raw.csv'))}
with open('haps.tsv', 'w') as out:
    for i, s in enumerate(order):
        for h in ('1', '2'):
            n, c, _ = el[s][h]
            u = rows[n]['assembly'].replace('s3://human-pangenomics/', 'https://s3-us-west-2.amazonaws.com/human-pangenomics/')
            out.write(f"{i}\t{s}\t{h}\t{n}\t{c}\t{u}\n")
PY
cp /home/user/work/t3/manifests/r21.windows0.tsv manifests/windows0.tsv
cp /home/user/work/t3/manifests/r21.hotspots300.tsv manifests/hotspots300k.tsv
# hotspots.tsv (haploblocks.org cuts, 50 kb minimum, 10 chunks) and mono.tsv as in Tier 1/2
python3 -c "
import sys; sys.path.insert(0,'$R')
from pipeline.parallel.cut_chunks import main
main(['manifests/hotspots.tsv','--chrom','chr21','--start','20000000','--end','21000000','--cuts','20271080,20347668,20405686,20508839,20583159,20650101,20701257,20795065,20899823'])"
cp /home/user/work/t3/manifests/r21.mono.tsv manifests/mono.tsv
