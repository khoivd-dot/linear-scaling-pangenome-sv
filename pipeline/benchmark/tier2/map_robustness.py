#!/usr/bin/env python3
"""Do haploblock boundaries fall near recombination hotspots in independent maps?

For each map (stdpopsim HomSap genetic maps, HapMap-format text: chrom, pos,
rate cM/Mb, ...) a hotspot is a run of map intervals with rate >= FOLD times
the chromosome-wide length-weighted mean (runs within 5 kb merged). The test
statistic is the fraction of haploblock cuts within `--dist` of a hotspot; the
null is 1,000 sets of the same number of uniformly random positions (seed 1).

    map_robustness.py chr21_pilot_20-21Mb.txt out.json --start 20000000 --end 21000000 \
        --map HapMapII=HapMapII.chr21.txt --map PyrhoYRI=PyrhoYRI.chr21.txt ...
"""
import argparse, json, random, re, statistics as st


def load(path):
    pts = []
    for line in open(path):
        p = line.split()
        if p and p[1].isdigit():
            pts.append((int(p[1]), float(p[2])))
    pts.sort()
    return [(pts[i][0], pts[i + 1][0], pts[i][1]) for i in range(len(pts) - 1)]


def hotspots(iv, fold, r0, r1):
    total = sum(e - s for s, e, _ in iv)
    mean = sum((e - s) * r for s, e, r in iv) / total
    hs = []
    for s, e, r in iv:
        if r >= fold * mean and e > r0 - 50000 and s < r1 + 50000:
            if hs and s - hs[-1][1] <= 5000:
                hs[-1] = (hs[-1][0], e, max(hs[-1][2], r))
            else:
                hs.append((s, e, r))
    return [((s + e) // 2, s, e, r) for s, e, r in hs], mean


def frac_within(points, hs, d):
    return sum(1 for p in points if any(s - d <= p <= e + d for _, s, e, _ in hs)) / len(points)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("blocks"); ap.add_argument("out")
    ap.add_argument("--start", type=int, required=True); ap.add_argument("--end", type=int, required=True)
    ap.add_argument("--map", action="append", required=True, help="NAME=PATH")
    ap.add_argument("--fold", type=int, nargs="+", default=[5, 10])
    ap.add_argument("--dist", type=int, default=5000)
    a = ap.parse_args()
    r0, r1 = a.start, a.end
    blocks = [l.strip() for l in open(a.blocks) if l.strip()]
    cuts = sorted({int(x) for b in blocks for x in re.match(r"\w+_(\d+)-(\d+)", b).groups() if r0 < int(x) < r1})
    rng = random.Random(1)
    nulls = [sorted(rng.randint(r0, r1) for _ in cuts) for _ in range(1000)]
    res = {}
    for spec in a.map:
        name, path = spec.split("=", 1)
        iv = load(path)
        widths = [e - s for s, e, _ in iv if r0 <= s < r1]
        for fold in a.fold:
            hs, mean = hotspots(iv, fold, r0, r1)
            inreg = [h for h in hs if r0 <= h[0] < r1]
            obs = frac_within(cuts, hs, a.dist)
            null = [frac_within(n, hs, a.dist) for n in nulls]
            res[f"{name}_x{fold}"] = {
                "map_bins_in_region": len(widths), "median_bin_bp": st.median(widths) if widths else None,
                "chrom_mean_cM_per_Mb": round(mean, 3), "hotspots_in_region": len(inreg),
                "hotspot_bp_in_region": sum(min(e, r1) - max(s, r0) for _, s, e, _ in inreg),
                "haploblock_cuts_within_5kb": round(obs, 3), "random_cuts_within_5kb_mean": round(st.mean(null), 3),
                "p_value": sum(1 for x in null if x >= obs) / len(null)}
            print(name, fold, res[f"{name}_x{fold}"])
    json.dump({"cuts": cuts, "results": res}, open(a.out, "w"), indent=1)


if __name__ == "__main__":
    main()
