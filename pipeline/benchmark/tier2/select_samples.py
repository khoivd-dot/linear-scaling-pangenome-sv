#!/usr/bin/env python3
"""Pick HPRC release 2 samples whose region sits on one contig per haplotype.

Reads the minigraph-cactus GRCh38 chains (annotation/chains/minigraph-cactus/
<assembly>_vs_GRCh38.chain.gz) and the release index CSV. A haplotype is
eligible when a single chain covers the region plus `--pad` on both sides;
a sample is eligible when both haplotypes are. Writes eligible.json and an
ordered sample list (the fixed samples first, then a seeded random draw), so
that N = 10, 20, 40 haplotype sets are nested.

    select_samples.py chains/ index_raw.csv out/ --chrom chr21 \
        --start 20000000 --end 21000000 --first HG00673 HG00733 --n 20 --seed 1
"""
import argparse, collections, csv, glob, gzip, json, os, random


def best_chain(path, chrom, start, end):
    best = None
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if not line.startswith("chain"):
                continue
            p = line.split()
            if p[7] != chrom:
                continue
            qs, qe = int(p[10]), int(p[11])
            if p[9] == "-":
                qs, qe = int(p[8]) - qe, int(p[8]) - qs
            if qs <= start and qe >= end and (best is None or int(p[1]) > best[0]):
                best = (int(p[1]), p[2], p[9])
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("chains_dir"); ap.add_argument("index_csv"); ap.add_argument("outdir")
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--start", type=int, required=True); ap.add_argument("--end", type=int, required=True)
    ap.add_argument("--pad", type=int, default=50000)
    ap.add_argument("--first", nargs="*", default=[])
    ap.add_argument("--n", type=int, default=20); ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    rows = {r["assembly_name"]: r for r in csv.DictReader(open(a.index_csv))}
    by = collections.defaultdict(dict)
    for f in glob.glob(os.path.join(a.chains_dir, "*_vs_GRCh38.chain.gz")):
        name = os.path.basename(f).replace("_vs_GRCh38.chain.gz", "")
        b = best_chain(f, a.chrom, a.start - a.pad, a.end + a.pad)
        if b and name in rows:
            r = rows[name]
            by[r["sample_id"]][r["haplotype"]] = (name, b[1], b[2])
    eligible = {s: h for s, h in sorted(by.items()) if len(h) == 2}
    rest = sorted(s for s in eligible if s not in a.first)
    order = list(a.first) + random.Random(a.seed).sample(rest, a.n - len(a.first))
    os.makedirs(a.outdir, exist_ok=True)
    json.dump(eligible, open(os.path.join(a.outdir, "eligible.json"), "w"), indent=0)
    json.dump(order, open(os.path.join(a.outdir, "order.json"), "w"))
    print(f"{len(eligible)} samples with both haplotypes on one chain; order: {' '.join(order)}")


if __name__ == "__main__":
    main()
