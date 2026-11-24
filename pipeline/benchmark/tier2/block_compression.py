#!/usr/bin/env python3
"""Do haplotypes inside a hotspot-bounded block repeat more often than inside random blocks?

Inputs are chunk FASTA directories from project_chunks.py: one for the hotspot
blocks and one per random layout (the same block sizes shuffled and rotated
around the region, see `random_layouts`). For each block and each nested
haplotype set (N = 10, 20, 40, in haps.tsv order) it counts distinct
haplotype sequences by exact identity (case-insensitive). The summary is the
duplicate fraction 1 - distinct/N, weighted by block length, over blocks of at
least `--min-bp`.

    block_compression.py haps.tsv comp/ --write-layouts 20   # comp/blocks.tsv -> comp/rand_XX.tsv
    (project every layout into comp/fa/<layout>/ with project_chunks.py, then)
    block_compression.py haps.tsv comp/ --min-bp 5000
"""
import argparse, csv, glob, hashlib, json, os, random, statistics as st


def random_layouts(sizes, start, n, seed=1):
    """Cut lists with the same block sizes in random order, rotated by a random offset."""
    rng, total, out = random.Random(seed), sum(sizes), []
    for _ in range(n):
        s = sizes[:]
        rng.shuffle(s)
        off = rng.randint(0, total - 1)
        pos = [start + (sum(s[:k]) + off) % total for k in range(1, len(s))] + ([start + off] if off else [])
        out.append(sorted({p for p in pos if start < p < start + total}))
    return out


def read_fasta(path):
    d, k = {}, None
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            k = line[1:].split(":")[0]
            d[k] = []
        else:
            d[k].append(line)
    return {k: "".join(v) for k, v in d.items()}


def count_distinct(comp_dir, sets):
    res = {}
    for d in sorted(glob.glob(os.path.join(comp_dir, "fa", "*"))):
        lay = os.path.basename(d)
        rows = {r["chunk_id"]: r for r in csv.DictReader(open(os.path.join(comp_dir, f"{lay}.tsv")), delimiter="\t")}
        res[lay] = {}
        for f in sorted(glob.glob(os.path.join(d, "chunk_*.fa"))):
            cid = os.path.basename(f)[:-3]
            seqs = read_fasta(f)
            e = {"bp": int(rows[cid]["reference_end"]) - int(rows[cid]["reference_start"])}
            for n, keep in sets.items():
                e[str(n)] = len({hashlib.md5(v.upper().encode()).hexdigest()
                                 for k, v in seqs.items() if "#".join(k.split("#")[:2]) in keep})
            res[lay][cid] = e
    return res


def dup_fraction(blocks, n, min_bp):
    xs = [(e["bp"], e[str(n)]) for e in blocks.values() if e["bp"] >= min_bp]
    tot = sum(b for b, _ in xs)
    return 1 - sum(b * c / n for b, c in xs) / tot


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("haps_tsv"); ap.add_argument("comp_dir")
    ap.add_argument("--min-bp", type=int, default=5000)
    ap.add_argument("--write-layouts", type=int, metavar="N", help="write N random layouts of comp/blocks.tsv and stop")
    a = ap.parse_args()
    if a.write_layouts:
        import sys
        sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))
        from pipeline.parallel.hotspot_chunks import create_hotspot_chunks
        from pipeline.parallel.make_chunks import write_manifest
        rows = list(csv.DictReader(open(os.path.join(a.comp_dir, "blocks.tsv")), delimiter="\t"))
        chrom, start, end = rows[0]["chrom"], int(rows[0]["reference_start"]), int(rows[-1]["reference_end"])
        sizes = [int(r["reference_end"]) - int(r["reference_start"]) for r in rows]
        for i, cuts in enumerate(random_layouts(sizes, start, a.write_layouts), 1):
            write_manifest(create_hotspot_chunks(chrom, start, end, cuts, min_size=0),
                           os.path.join(a.comp_dir, f"rand_{i:02d}.tsv"))
        return
    sets = {10: set(), 20: set(), 40: set()}
    for i, s, h, *_ in (l.rstrip("\n").split("\t") for l in open(a.haps_tsv)):
        for n in sets:
            if (int(i) + 1) * 2 <= n:
                sets[n].add(f"{s}#{h}")
    res = count_distinct(a.comp_dir, sets)
    summary = {}
    for n in sets:
        hot = dup_fraction(res["blocks"], n, a.min_bp)
        rnd = [dup_fraction(v, n, a.min_bp) for k, v in res.items() if k.startswith("rand_")]
        summary[n] = {"hotspot_blocks": round(hot, 4), "random_median": round(st.median(rnd), 4),
                      "random_min": round(min(rnd), 4), "random_max": round(max(rnd), 4),
                      "random_layouts_more_duplicated": sum(r > hot for r in rnd), "random_layouts": len(rnd)}
        print(n, summary[n])
    json.dump({"min_bp": a.min_bp, "summary": summary, "blocks": res},
              open(os.path.join(a.comp_dir, "compression.json"), "w"), indent=0)


if __name__ == "__main__":
    main()
