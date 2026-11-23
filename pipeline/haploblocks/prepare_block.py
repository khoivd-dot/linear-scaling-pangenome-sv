"""prepare_block.py — turn a haploblocks.org block FASTA into PGGB input.

data.haploblocks.org splits 1000 Genomes haplotypes (2,548 samples x 2) at
recombination hotspots and serves one multi-FASTA per block, e.g.
    chr22_small_dataset/chr22_44303006-44307622_all_seqs.fasta
with headers like  >HG00096_chr22_region_44303006-44307622_hap0

This script:
  * drops empty records (the published files carry a few header-only lines),
  * renames to PanSN  SAMPLE#HAP#CHROM:START-END  (hap0 -> 1, hap1 -> 2),
  * optionally keeps a seeded random subset of N haplotypes (--n-haps),
  * optionally collapses identical sequences (--dedupe) and writes a TSV that
    maps every haplotype to the sequence it was collapsed into.

Usage:
    python3 pipeline/haploblocks/prepare_block.py IN.fasta OUT.fa [--n-haps N] [--seed S] [--dedupe]
"""
import argparse
import random
import re
import sys

HEADER = re.compile(r"^(?P<sample>[^_]+)_(?P<chrom>[^_]+)_region_(?P<start>\d+)-(?P<end>\d+)_hap(?P<hap>\d)$")


def read_fasta(path):
    name, buf = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        yield name, "".join(buf)


def pansn(name):
    m = HEADER.match(name)
    if not m:
        sys.exit(f"FATAL: unexpected header {name!r}")
    return f"{m['sample']}#{int(m['hap']) + 1}#{m['chrom']}:{m['start']}-{m['end']}"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("input")
    ap.add_argument("output")
    ap.add_argument("--n-haps", type=int, default=0, help="keep a random subset of N haplotypes (0 = all)")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--dedupe", action="store_true", help="collapse identical sequences")
    args = ap.parse_args(argv)

    records, seen, empty = [], set(), 0
    for name, seq in read_fasta(args.input):
        if not seq:
            empty += 1
            continue
        new = pansn(name)
        if new in seen:
            sys.exit(f"FATAL: duplicate non-empty record {name}")
        seen.add(new)
        records.append((new, seq))

    if args.n_haps and args.n_haps < len(records):
        records = sorted(random.Random(args.seed).sample(records, args.n_haps))

    mapping = None
    if args.dedupe:
        first, kept, mapping = {}, [], []
        for name, seq in records:
            if seq not in first:
                first[seq] = name
                kept.append((name, seq))
            mapping.append((name, first[seq]))
        records = kept

    with open(args.output, "w") as out:
        for name, seq in records:
            out.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                out.write(seq[i:i + 80] + "\n")
    if mapping is not None:
        with open(args.output + ".dedupe.tsv", "w") as out:
            out.write("haplotype\trepresentative\n")
            out.writelines(f"{a}\t{b}\n" for a, b in mapping)

    total = len(mapping) if mapping is not None else len(records)
    print(f"{args.input}: dropped {empty} empty records; {total} haplotypes -> {len(records)} sequences written to {args.output}")


if __name__ == "__main__":
    main()
