#!/usr/bin/env python3
"""For each region and haplotype, the contig whose minigraph-cactus GRCh38 chain spans the region (+pad).

    pick_contigs.py regions.tsv chains/ HG00673 HG00733 > contigs.tsv
Output: region, sample, hap (1=pat, 2=mat), assembly, contig, strand ('-' when no single chain spans it).
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "tier2"))
from select_samples import best_chain

regions_tsv, chains, *samples = sys.argv[1:]
regions = [l.split() for l in open(regions_tsv) if l.strip()]
for s in samples:
    for h, hn in (("1", "pat"), ("2", "mat")):
        asm = f"{s}_{hn}_hprc_r2_v1.0.1"
        path = os.path.join(chains, f"{asm}_vs_GRCh38.chain.gz")
        for reg, chrom, start, end in regions:
            b = best_chain(path, chrom, int(start) - 50000, int(end) + 50000)
            print("\t".join([reg, s, h, asm, b[1] if b else "NA", b[2] if b else "NA"]))
