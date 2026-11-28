#!/usr/bin/env python3
"""Tier 4 experiment sets: regions, chunkings and haplotypes, and the GitHub job matrices for one run.

    sets.py plan SET N            JSON for the workflow: samples, build units, graphs, regions
    sets.py prefixes SET N        SAMPLE#HAP list in input order
    sets.py haps SET N SAMPLE     rows of the haplotype table for one sample (region, prefix, assembly, contig, url)
    sets.py region REGION         chrom start end

SET arm   chr21:15-46 Mb, the first N haplotypes of haps_arm.tsv (N = 6, 20 or 40);
          mono + hotspots1m, plus windows1m when N <= 6
SET hard  17q21.31 and 22q11.2 (Tier 3 regions), HG00673, HG00733, HG002 (N = 6); mono + windows0
"""
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REGIONS = {"arm": ("chr21", 15_000_000, 46_000_000), "r17": ("chr17", 45_400_000, 46_600_000),
           "r22": ("chr22", 18_600_000, 21_500_000)}


def haps(set_, n):
    if set_ == "arm":
        rows = list(csv.DictReader(open(os.path.join(HERE, "haps_arm.tsv")), delimiter="\t"))[:n]
        return [dict(r, region="arm") for r in rows]
    return list(csv.DictReader(open(os.path.join(HERE, "haps_hard.tsv")), delimiter="\t"))


def strategies(set_, n):
    if set_ == "arm":
        return [("arm", "mono"), ("arm", "hotspots1m")] + ([("arm", "windows1m")] if n <= 6 else [])
    return [(r, s) for r in ("r17", "r22") for s in ("mono", "windows0")]


def manifest(region, strategy):
    return os.path.join(HERE, "manifests", f"{region}.{strategy}.tsv")


def chunk_ids(region, strategy):
    return [r["chunk_id"] for r in csv.DictReader(open(manifest(region, strategy)), delimiter="\t")]


def prefixes(set_, n):
    out = []
    for r in haps(set_, n):
        p = f"{r['sample']}#{r['hap']}"
        if p not in out:
            out.append(p)
    return out


def main():
    cmd, *args = sys.argv[1:]
    if cmd == "region":
        print(*REGIONS[args[0]])
        return
    set_, n = args[0], int(args[1])
    if cmd == "prefixes":
        print(" ".join(prefixes(set_, n)))
    elif cmd == "haps":
        for r in haps(set_, n):
            if r["sample"] == args[2]:
                print("\t".join([r["region"], f"{r['sample']}#{r['hap']}", r["assembly"], r["contig"], r["url"]]))
    elif cmd == "plan":
        st = strategies(set_, n)
        units = [{"region": reg, "strategy": s, "chunk": c, "id": f"{reg}.{s}.{c}"}
                 for reg, s in st for c in chunk_ids(reg, s)]
        # the largest builds first, so they start before the runner queue fills
        units.sort(key=lambda u: u["strategy"] != "mono")
        samples = list(dict.fromkeys(p.split("#")[0] for p in prefixes(set_, n)))
        print(json.dumps({"samples": samples,
                          "units": units, "graphs": [{"region": reg, "strategy": s, "id": f"{reg}.{s}"} for reg, s in st],
                          "regions": sorted({reg for reg, _ in st})}))


if __name__ == "__main__":
    main()
