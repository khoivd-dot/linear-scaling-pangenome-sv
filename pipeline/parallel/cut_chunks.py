"""cut_chunks.py — touching (zero-overlap) chunks from given or random cut points.

Controls for hotspot_chunks.py: the same end-to-end join, with cuts that do
not come from a recombination map.

    # explicit cuts (e.g. hotspot cuts shifted by 50 kb)
    python3 pipeline/parallel/cut_chunks.py OUT.tsv --chrom chr21 \
        --start 20000000 --end 21000000 --cuts 20397668,20700101
    # N random cut sets, K chunks each, every chunk >= --min-size
    python3 pipeline/parallel/cut_chunks.py OUTDIR --chrom chr21 \
        --start 20000000 --end 21000000 --random 20 --chunks 3 \
        --min-size 200000 --seed 1
"""
import argparse
import os
import random
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from pipeline.parallel.hotspot_chunks import create_hotspot_chunks  # noqa: E402
from pipeline.parallel.make_chunks import write_manifest  # noqa: E402


def random_cuts(rng, start, end, chunks, min_size):
    """Uniform cut set: sample the slack above min_size, split it at random."""
    slack = (end - start) - chunks * min_size
    if slack < 0:
        raise ValueError("region too short for that many chunks of --min-size")
    extra = sorted(rng.randint(0, slack) for _ in range(chunks - 1))
    cuts, pos, prev = [], start, 0
    for e in extra:
        pos += min_size + (e - prev)
        prev = e
        cuts.append(pos)
    return cuts


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("out")
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--start", type=int, required=True)
    ap.add_argument("--end", type=int, required=True)
    ap.add_argument("--cuts", help="comma-separated cut positions")
    ap.add_argument("--random", type=int, default=0, help="number of random cut sets")
    ap.add_argument("--chunks", type=int, default=3)
    ap.add_argument("--min-size", type=int, default=200000)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args(argv)
    if a.cuts:
        cuts = [int(c) for c in a.cuts.split(",")]
        write_manifest(create_hotspot_chunks(a.chrom, a.start, a.end, cuts, min_size=0), a.out)
        return
    os.makedirs(a.out, exist_ok=True)
    rng = random.Random(a.seed)
    for i in range(1, a.random + 1):
        cuts = random_cuts(rng, a.start, a.end, a.chunks, a.min_size)
        write_manifest(create_hotspot_chunks(a.chrom, a.start, a.end, cuts, min_size=0),
                       os.path.join(a.out, f"random_{i:02d}.tsv"))


if __name__ == "__main__":
    main()
