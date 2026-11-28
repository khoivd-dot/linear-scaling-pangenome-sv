#!/usr/bin/env python3
"""Group chunks that share a segmental duplication, so paralogous copies are built in one PGGB job.

The reference region is cut into 20 kb tiles (10 kb step) and each tile is mapped back to the
whole region with minimap2 (mappy, asm20). A hit away from the tile's own position, at least
--min-len long and --min-identity identical, links the chunk holding the tile to the chunk holding
the hit. Chunks joined by at least --min-link bp of such hits are merged (union-find); each group is
built as one PGGB job from its members' FASTAs, which need not be adjacent.

    sd_groups.py REGION.fa MANIFEST.tsv --region-start 18600000 > groups.tsv
Output columns: group, chunk_ids (comma-separated), linked_bp; pairs go to stderr.
"""
import argparse, collections, csv, sys
import mappy


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("region_fa"); ap.add_argument("manifest")
    ap.add_argument("--region-start", type=int, required=True)
    ap.add_argument("--tile", type=int, default=20000); ap.add_argument("--step", type=int, default=10000)
    ap.add_argument("--min-len", type=int, default=5000)
    ap.add_argument("--min-identity", type=float, default=0.90)
    ap.add_argument("--min-link", type=int, default=10000)
    a = ap.parse_args(argv)
    name, seq, _ = next(mappy.fastx_read(a.region_fa))
    rows = list(csv.DictReader(open(a.manifest), delimiter="\t"))
    bounds = [(int(r["reference_start"]) - a.region_start, int(r["reference_end"]) - a.region_start, r["chunk_id"]) for r in rows]

    def chunk_at(p):
        return next((c for s, e, c in bounds if s <= p < e), None)

    al = mappy.Aligner(seq=seq, preset="asm20")
    link = collections.Counter()
    for t0 in range(0, len(seq) - a.tile + 1, a.step):
        tile = seq[t0:t0 + a.tile]
        if tile.count("N") > a.tile // 2:
            continue
        for h in al.map(tile):
            if h.r_st < t0 + a.tile and h.r_en > t0:       # the tile itself (or overlapping it)
                continue
            if h.blen < a.min_len or h.mlen / h.blen < a.min_identity:
                continue
            c1 = chunk_at(t0 + (h.q_st + h.q_en) // 2); c2 = chunk_at((h.r_st + h.r_en) // 2)
            if c1 and c2 and c1 != c2:
                link[tuple(sorted((c1, c2)))] += h.blen * a.step // a.tile   # tiles overlap 2x: count half
    parent = {c: c for _, _, c in bounds}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x

    for (c1, c2), bp in sorted(link.items()):
        print(f"pair\t{c1}\t{c2}\t{bp}", file=sys.stderr)
        if bp >= a.min_link:
            parent[find(c1)] = find(c2)
    groups = collections.defaultdict(list)
    for _, _, c in bounds:
        groups[find(c)].append(c)
    print("group\tchunk_ids\tlinked_bp")
    for i, members in enumerate(sorted(groups.values(), key=lambda m: m[0]), 1):
        bp = sum(v for k, v in link.items() if k[0] in members and k[1] in members)
        print(f"group_{i:02d}\t{','.join(members)}\t{bp}")


if __name__ == "__main__":
    main()
