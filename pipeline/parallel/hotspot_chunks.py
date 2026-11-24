"""hotspot_chunks.py — chunk a region at recombination hotspots instead of fixed windows.

haploblocks.org splits each chromosome at recombination hotspots (1000 Genomes,
2,548 samples). Its block names, like `chr21_20023729-20271080`, give the cut
points. Consecutive blocks share their boundary coordinate, so chunks cut there
touch without overlapping, and the stitch joins each haplotype end to end
(seam_position accepts touching intervals).

Blocks shorter than --min-size are merged with their right neighbour, so no
chunk is shorter than PGGB's mapping block length (-l, 5 x -s by default).

Usage:
    python3 pipeline/parallel/hotspot_chunks.py BLOCKS.txt OUT.tsv \
        --chrom chr21 --start 20000000 --end 21000000 [--min-size 50000]
BLOCKS.txt holds one block name per line (aws/blocks/chr21_pilot_20-21Mb.txt).
"""
import argparse
import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from pipeline.parallel.make_chunks import write_manifest  # noqa: E402

BLOCK = re.compile(r"^(?P<chrom>[^_]+)_(?P<start>\d+)-(?P<end>\d+)$")


def read_cuts(names, chrom, start, end):
    """Sorted block boundaries strictly inside (start, end)."""
    cuts = set()
    for name in names:
        m = BLOCK.match(name.strip())
        if not m or m["chrom"] != chrom:
            continue
        for pos in (int(m["start"]), int(m["end"])):
            if start < pos < end:
                cuts.add(pos)
    return sorted(cuts)


def create_hotspot_chunks(chrom, start, end, cuts, min_size=50000):
    """Manifest rows (make_chunks schema) for chunks that touch at hotspot cuts."""
    edges = [start] + [c for c in cuts if start < c < end] + [end]
    spans, lo = [], start
    for hi in edges[1:]:
        if hi - lo >= min_size or hi == end:
            spans.append([lo, hi])
            lo = hi
    if len(spans) > 1 and spans[-1][1] - spans[-1][0] < min_size:
        tail = spans.pop()  # a short tail joins its left neighbour
        spans[-1][1] = tail[1]
    return [{
        "chunk_id": f"chunk_{i:04d}", "chrom": chrom,
        "reference_start": lo, "reference_end": hi,
        "core_start": lo, "core_end": hi,
        "overlap_left": 0, "overlap_right": 0,
        "pairwise_overlap_bp": 0, "status": "pending",
    } for i, (lo, hi) in enumerate(spans, 1)]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("blocks")
    ap.add_argument("out")
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--start", type=int, required=True)
    ap.add_argument("--end", type=int, required=True)
    ap.add_argument("--min-size", type=int, default=50000)
    a = ap.parse_args(argv)
    with open(a.blocks) as fh:
        cuts = read_cuts(fh, a.chrom, a.start, a.end)
    chunks = create_hotspot_chunks(a.chrom, a.start, a.end, cuts, a.min_size)
    write_manifest(chunks, a.out)
    sizes = [c["reference_end"] - c["reference_start"] for c in chunks]
    print(f"{len(cuts)} hotspot cuts -> {len(chunks)} chunks "
          f"({min(sizes):,}-{max(sizes):,} bp) -> {a.out}")


if __name__ == "__main__":
    main()
