#!/usr/bin/env python3
"""Tier 6 E10: compare the Tier 3 truth calls with the stand-alone dipcall on the same haplotypes.

dipcall (run by pipeline/linear/run_dipcall.sh's method: whole contigs aligned to the reference) is the
base, the Tier 3 dipcall-style calls from HPRC's winnowmap alignments (tier3/truth_calls.py) are the
comparison. Each haplotype is scored inside the region, dipcall's confident BED and the truth's
callable BED, with the same truvari settings as score_vs_truth.py, on all sizes and on SVs (>= 50 bp).

    dipcall_check.py DIP.vcf.gz DIP.bed TRUTH_PREFIX CHROM START END REF.fa OUT.json \
        --sample HG00673 --truvari truvari [--segdups segdups.bed]
TRUTH_PREFIX_<sample>_<1|2>.{vcf,bed} are truth_calls.py outputs; dipcall's GT is hap1|hap2 in the
order its two FASTA files were given.
"""
import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "tier3"))
svt = None  # tier3/score_vs_truth.py, imported in main() because it needs pysam


def intervals(path, chrom, start, end):
    out = []
    for line in open(path):
        f = line.split()
        if f and f[0] == chrom:
            s, e = max(int(f[1]), start), min(int(f[2]), end)
            if s < e:
                out.append((s, e))
    out.sort()
    merged = []
    for s, e in out:
        if merged and s <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append((s, e))
    return merged


def intersect(a, b):
    out, i, j = [], 0, 0
    while i < len(a) and j < len(b):
        s, e = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if s < e:
            out.append((s, e))
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return out


def dip_haploid(vcf, chrom, start, end, hap, out):
    import pysam
    recs = set()
    with pysam.VariantFile(vcf) as f:
        for r in f.fetch(chrom, start, end) if f.index is not None else f:
            if r.chrom != chrom or not start <= r.pos - 1 < end:
                continue
            gt = r.samples[0]["GT"]
            if gt is None or len(gt) < hap or not gt[hap - 1]:
                continue
            alt = r.alleles[gt[hap - 1]]
            if alt == "*" or set(alt) - set("ACGTN"):
                continue
            recs.add((r.pos, r.ref, alt))
    svt.write_vcf(out, [], sorted(recs), chrom)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("dip_vcf"); ap.add_argument("dip_bed"); ap.add_argument("truth"); ap.add_argument("chrom")
    ap.add_argument("start", type=int); ap.add_argument("end", type=int); ap.add_argument("ref"); ap.add_argument("out")
    ap.add_argument("--sample", required=True); ap.add_argument("--truvari", required=True)
    ap.add_argument("--segdups")
    a = ap.parse_args(argv)
    global svt
    import score_vs_truth as svt
    svt.SD = svt.load_bed(a.segdups, a.chrom) if a.segdups else ([], [])
    work = a.out + "_work"
    os.makedirs(work, exist_ok=True)
    dbed = intervals(a.dip_bed, a.chrom, a.start, a.end)
    res = {"sample": a.sample, "region": f"{a.chrom}:{a.start}-{a.end}", "haps": {}}
    for h in (1, 2):
        tbed = intervals(f"{a.truth}_{a.sample}_{h}.bed", a.chrom, a.start, a.end)
        both = intersect(dbed, tbed)
        bed = os.path.join(work, f"bed_{h}.bed")
        open(bed, "w").write("".join(f"{a.chrom}\t{s}\t{e}\n" for s, e in both))
        base = svt.norm(dip_haploid(a.dip_vcf, a.chrom, a.start, a.end, h, os.path.join(work, f"dip_{h}.vcf")), a.ref, a.chrom)
        comp = svt.norm(svt.haploid(f"{a.truth}_{a.sample}_{h}.vcf", a.sample, 1, os.path.join(work, f"truth_{h}.vcf")),
                        a.ref, a.chrom)
        hv = {"bp": {"dipcall_bed": sum(e - s for s, e in dbed), "truth_bed": sum(e - s for s, e in tbed),
                     "both": sum(e - s for s, e in both)}}
        for m in ("all", "sv"):
            parts = svt.bench(a.truvari, base, comp, bed, os.path.join(work, f"truvari_{m}_{h}"), m == "sv")
            tot = {k: sum(v.values()) for k, v in parts.items()}
            hv[m] = {"total": svt.prf(tot["tp-base"], tot["fn"], tot["tp-comp"], tot["fp"]),
                     "sd": svt.prf(*(parts[k]["sd"] for k in ("tp-base", "fn", "tp-comp", "fp"))),
                     "other": svt.prf(*(parts[k]["other"] for k in ("tp-base", "fn", "tp-comp", "fp")))}
        res["haps"][f"{a.sample}#{h}"] = hv
        print(a.sample, h, hv["bp"], "F1 all", hv["all"]["total"]["f1"], "SV", hv["sv"]["total"]["f1"])
    json.dump(res, open(a.out, "w"), indent=1)


if __name__ == "__main__":
    main()
