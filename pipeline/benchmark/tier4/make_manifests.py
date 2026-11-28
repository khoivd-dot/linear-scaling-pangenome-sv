#!/usr/bin/env python3
"""Tier 4 chunk manifests (zero overlap), written into pipeline/benchmark/tier4/manifests/.

arm  chr21:15,000,000-46,000,000 (31 Mb of the long arm)
  mono        one chunk
  windows1m   fixed 1 Mb windows
  hotspots1m  haploblocks.org block boundaries, merged left to right until a chunk reaches 1 Mb;
              a last piece under 500 kb joins the chunk before it
mid  chr21:20,000,000-25,000,000 (Tier 6 E8): mono, windows1m
r17  chr17:45,400,000-46,600,000 and r22 chr22:18,600,000-21,500,000 (Tier 3 hard regions)
  mono, windows0 (the Tier 3 zero-overlap windows: 3 x 400 kb and 7 x 414,286 bp)

    make_manifests.py BLOCK_STATS.tsv   # data.haploblocks.org/bidirectional_blast_samples/block_stats.tsv
"""
import csv
import os
import sys

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..")
sys.path.insert(0, ROOT)
from pipeline.parallel.cut_chunks import main as cut  # noqa: E402
from pipeline.parallel.make_chunks import create_chunks, write_manifest  # noqa: E402

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "manifests")


def hotspot_cuts(block_stats, chrom, start, end, min_size, min_last):
    bounds = sorted({int(r[k]) for r in csv.DictReader(open(block_stats), delimiter="\t")
                     if r["chr"] == chrom for k in ("start", "end")})
    cuts, last = [], start
    for b in bounds:
        if start < b < end and b - last >= min_size:
            cuts.append(b)
            last = b
    if cuts and end - cuts[-1] < min_last:
        cuts.pop()
    return cuts


def window(name, chrom, start, end, size):
    write_manifest(create_chunks(chrom, start, end, size, 0), os.path.join(OUT, f"{name}.tsv"))


def main():
    os.makedirs(OUT, exist_ok=True)
    s, e = 15_000_000, 46_000_000
    cuts = hotspot_cuts(sys.argv[1], "chr21", s, e, 1_000_000, 500_000)
    cut([os.path.join(OUT, "arm.hotspots1m.tsv"), "--chrom", "chr21", "--start", str(s), "--end", str(e),
         "--cuts", ",".join(map(str, cuts))])
    window("arm.windows1m", "chr21", s, e, 1_000_000)
    window("mid.windows1m", "chr21", 20_000_000, 25_000_000, 1_000_000)
    window("mid.mono", "chr21", 20_000_000, 25_000_000, 5_000_000)
    for reg, chrom, rs, re_, size in (("arm", "chr21", s, e, None), ("r17", "chr17", 45_400_000, 46_600_000, 400_000),
                                      ("r22", "chr22", 18_600_000, 21_500_000, 414_286)):
        if size:
            window(f"{reg}.windows0", chrom, rs, re_, size)
        window(f"{reg}.mono", chrom, rs, re_, re_ - rs)
    print(f"arm hotspots1m: {len(cuts) + 1} chunks; cuts {cuts}")


if __name__ == "__main__":
    main()
