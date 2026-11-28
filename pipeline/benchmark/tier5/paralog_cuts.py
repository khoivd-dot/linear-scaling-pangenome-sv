#!/usr/bin/env python3
"""Add cuts that separate paralogous copies sharing a chunk, plus a control with the same number of cuts.

Tier 3 found that zero-overlap chunks merge fewer paralogs than a monolithic build because a chunk
cannot align two copies it does not hold. Copies that still share a chunk can be merged. This script
takes a base chunking (e.g. windows0), finds paralog pairs in the GRCh38 region by self-alignment
(20 kb tiles every 10 kb, minimap2 asm20 via mappy; hits >= --min-len bp at >= --min-identity, the
two copies >= --min-dist apart), and for each pair whose copies sit in the same chunk adds one cut
between them, in sequence outside every self-alignment hit (midpoint of the longest such stretch).
Pairs are taken in order of their right copy, and a pair already split by an earlier cut needs none.
A cut closer than --min-chunk to another is not added. The control ("ctrl") adds the same number of
cuts in chunks that hold no paralog pair, each at the same distance from the nearest original boundary
as its sdsplit counterpart (moved to the nearest unmasked base), so chunk count and sizes change the same
way without separating any copies.

    paralog_cuts.py REGION.fa BASE.tsv --region-start 45400000 --out-prefix manifests/r17 [--pairs pairs.tsv]
Writes OUT_PREFIX.sdsplit.tsv and OUT_PREFIX.ctrl.tsv (chunk manifests, zero overlap).
"""
import argparse, csv, sys
import mappy

FIELDS = ["chunk_id", "chrom", "reference_start", "reference_end", "core_start", "core_end",
          "overlap_left", "overlap_right", "pairwise_overlap_bp", "status"]


def self_hits(seq, tile, step, min_len, min_id):
    al = mappy.Aligner(seq=seq, preset="asm20")
    hits = []
    for t0 in range(0, len(seq) - tile + 1, step):
        s = seq[t0:t0 + tile]
        if s.count("N") > tile // 2:
            continue
        for h in al.map(s):
            if h.r_st < t0 + tile and h.r_en > t0:
                continue
            if h.blen < min_len or h.mlen / h.blen < min_id:
                continue
            a, b = (t0 + h.q_st, t0 + h.q_en), (h.r_st, h.r_en)
            hits.append(tuple(sorted((a, b))))
    return hits


def merge(iv):
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def free_runs(lo, hi, mask):
    """Stretches of [lo, hi) outside the merged mask."""
    runs, p = [], lo
    for s, e in mask:
        if e <= p or s >= hi:
            continue
        if s > p:
            runs.append((p, min(s, hi)))
        p = max(p, e)
    if p < hi:
        runs.append((p, hi))
    return runs


def write(path, chrom, start, cuts, end):
    b = [start] + sorted(cuts) + [end]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, FIELDS, delimiter="\t"); w.writeheader()
        for i, (s, e) in enumerate(zip(b, b[1:]), 1):
            w.writerow({"chunk_id": f"chunk_{i:04d}", "chrom": chrom, "reference_start": s, "reference_end": e,
                        "core_start": s, "core_end": e, "overlap_left": 0, "overlap_right": 0,
                        "pairwise_overlap_bp": 0, "status": "pending"})


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("region_fa"); ap.add_argument("base")
    ap.add_argument("--region-start", type=int, required=True); ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--tile", type=int, default=20000); ap.add_argument("--step", type=int, default=10000)
    ap.add_argument("--min-len", type=int, default=5000); ap.add_argument("--min-identity", type=float, default=0.90)
    ap.add_argument("--min-dist", type=int, default=50000); ap.add_argument("--min-chunk", type=int, default=100000)
    ap.add_argument("--pairs", help="write the paralog pairs and the decision for each")
    a = ap.parse_args(argv)
    name, seq, _ = next(mappy.fastx_read(a.region_fa))
    rs = a.region_start
    rows = list(csv.DictReader(open(a.base), delimiter="\t"))
    chrom, end = rows[0]["chrom"], int(rows[-1]["reference_end"])
    base = [int(r["reference_start"]) for r in rows[1:]]
    hits = self_hits(seq, a.tile, a.step, a.min_len, a.min_identity)
    mask = merge([iv for h in hits for iv in h])
    pairs = sorted({(x[0] + rs, x[1] + rs, y[0] + rs, y[1] + rs) for x, y in hits if y[0] - x[1] >= a.min_dist - (x[1] - x[0])
                    and (y[0] + y[1]) // 2 - (x[0] + x[1]) // 2 >= a.min_dist}, key=lambda p: (p[2], p[0]))
    mask_abs = [(s + rs, e + rs) for s, e in mask]
    cuts, log = list(base), []
    for p in pairs:
        if any(p[1] <= c <= p[2] for c in cuts):
            log.append((*p, "split")); continue
        runs = free_runs(p[1], p[2], mask_abs)
        if not runs:
            log.append((*p, "no unmasked gap")); continue
        s, e = max(runs, key=lambda r: r[1] - r[0])
        c = (s + e) // 2
        if min(abs(c - x) for x in cuts + [rs, end]) < a.min_chunk:
            log.append((*p, f"cut {c} too close")); continue
        cuts.append(c); log.append((*p, f"new cut {c}"))
    added = sorted(set(cuts) - set(base))
    # control: same number of extra cuts, in base chunks that hold no pair at all
    bounds = [rs] + sorted(base) + [end]
    def has_pair(lo, hi):
        return any(lo <= (p[0] + p[1]) // 2 < hi and lo <= (p[2] + p[3]) // 2 < hi for p in pairs)
    free = [(lo, hi) for lo, hi in zip(bounds, bounds[1:]) if not has_pair(lo, hi)]
    ctrl, used = list(base), set()
    for c in added:   # same distance from the nearest original boundary, in a chunk holding no pair
        lo, hi = next((lo, hi) for lo, hi in zip(bounds, bounds[1:]) if lo < c < hi)
        d, left = (c - lo, True) if c - lo <= hi - c else (hi - c, False)
        for flo, fhi in sorted(free, key=lambda b: abs((b[1] - b[0]) - (hi - lo))):
            if (flo, fhi) in used or fhi - flo <= 2 * min(d, a.min_chunk):
                continue
            want = flo + d if left else fhi - d
            runs = free_runs(flo + min(d, a.min_chunk) // 2, fhi - min(d, a.min_chunk) // 2, mask_abs)
            if runs:
                ctrl.append(min((max(s, min(want, e - 1)) for s, e in runs), key=lambda x: abs(x - want)))
                used.add((flo, fhi)); break
    if len(ctrl) - len(base) < len(added):
        sys.exit(f"control could place only {len(ctrl) - len(base)} of {len(added)} cuts")
    write(f"{a.out_prefix}.sdsplit.tsv", chrom, rs, cuts, end)
    write(f"{a.out_prefix}.ctrl.tsv", chrom, rs, ctrl, end)
    if a.pairs:
        with open(a.pairs, "w") as fh:
            fh.write("copy1_start\tcopy1_end\tcopy2_start\tcopy2_end\tdecision\n")
            for r in log:
                fh.write("\t".join(map(str, r)) + "\n")
    print(f"base cuts {sorted(base)}\nsdsplit adds {added}\nctrl adds {sorted(set(ctrl) - set(base))}\n"
          f"paralog pairs >= {a.min_dist} bp apart: {len(pairs)}", file=sys.stderr)


if __name__ == "__main__":
    main()
