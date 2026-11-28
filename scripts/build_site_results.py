#!/usr/bin/env python3
"""Collect the numbers the website shows into web/public/data/results.json.

Reads the result files in docs/ (tier*.json, recombination_chunks_comparison.json,
tier4_runs/*.json, tier3_path_revisits_*.tsv) and rounds each value to the precision
the page prints, so every number on the page can be traced to one entry here.
Numbers that exist only in a Markdown write-up are constants below, each naming
its source. Stdlib only.

    python3 scripts/build_site_results.py           # write
    python3 scripts/build_site_results.py --check   # exit 1 if results.json is stale
"""
import csv
import glob
import json
import os
import sys
from decimal import ROUND_HALF_UP, Decimal

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
DOCS = os.path.join(ROOT, "docs")
OUT = os.path.join(ROOT, "web", "public", "data", "results.json")

# docs/TIER4.md, intro: GitHub-hosted runners, 4 vCPU, 16 GB RAM, one runner per chunk.
RUNNER = {"vcpu": 4, "ram_gb": 16}
# docs/TIER4.md, E6: the arm is chr21:15-46 Mb (31 Mb); chunks start at 15 Mb.
ARM = {"chrom": "chr21", "start": 15_000_000, "end": 46_000_000, "mb": 31}
# docs/TIER4.md, Summary and E6 table: what happened to the monolithic arm builds
# that did not finish (the run JSONs only list the build as missing).
MONO_DNF = {
    20: {"status": "killed", "after_min": 63,
         "why": "killed by the kernel, out of memory (exit 137)"},
    40: {"status": "lost", "after_min": 300,
         "why": "runner stopped responding inside PGGB; GitHub kept no log"},
}
# docs/TIER4.md, E6: arm truth set has 1,700 truth SVs at 6 haplotypes; the SV
# differences between graphs are 25 to 40 calls; chunk builds ran up to 20 at a time.
ARM_TRUTH_SVS = 1700
ARM_SV_DIFF_CALLS = [25, 40]
MAX_CONCURRENT_JOBS = 20
# docs/TIER3.md, "Truth set" table and docs/TIER4.md, E7: region coordinates.
REGIONS = {
    "r21": {"name": "chr21 pilot", "locus": "chr21:20,000,000-21,000,000"},
    "r17": {"name": "17q21.31", "locus": "chr17:45,400,000-46,600,000"},
    "r22": {"name": "22q11.2", "locus": "chr22:18,600,000-21,500,000"},
}
# docs/TIER3.md, E1 table: 17q21.31 has 9 truth SVs ("0.900 (9)").
T3_R17_TRUTH_SVS = 9
# docs/TIER4.md, E7: GIAB v5.0q has benchmark regions over about 70% of 22q11.2
# and none in 17q21.31.
GIAB_22Q_COVERAGE_PCT = 70


def load(path):
    with open(os.path.join(DOCS, path)) as fh:
        return json.load(fh)


def r(x, nd):
    """Round half up, as the write-ups do (0.8975 -> 0.898, where round() gives 0.897)."""
    if x is None:
        return None
    q = Decimal(repr(x)).quantize(Decimal(1).scaleb(-nd), rounding=ROUND_HALF_UP)
    return int(q) if nd == 0 else float(q)


def gb(mb):
    """Megabytes to GB as the write-ups print them: 0.97, 1.50, 12.1."""
    g = mb / 1000
    return r(g, 1 if g >= 10 else 2)


def f1(block, key="truth"):
    """Overall F1 for all variant sizes from a scored graph block."""
    t = block.get(key)
    return None if t is None else t["all"]["total"]["f1"]


def tier4():
    runs = {}
    for path in sorted(glob.glob(os.path.join(DOCS, "tier4_runs", "*.json"))):
        d = load(os.path.relpath(path, DOCS))
        runs[(d["set"], d["n_haplotypes"], d["pggb"], d["rep"])] = d["graphs"]

    def cost(g, n):
        c = g["cost"]
        out = {
            "builds": g["chunks"],
            "cpu_h": r(c["cpu_s"] / 3600, 2),
            "longest_min": r(c["max_wall_s"] / 60, 1),
            "peak_gb": gb(c["max_peak_rss_mb"]),
            "cpu_s_per_mb_hap": r(c["cpu_s"] / ARM["mb"] / n, 1),
        }
        return out

    def acc(g):
        hap = g["truth_a100k"]["per_hap"]
        return {
            "f1_all": r(g["truth_a100k"]["all"]["total"]["f1"], 3),
            "f1_sv": r(g["truth_a100k"]["sv"]["total"]["f1"], 3),
            "hg002_f1": [r(hap["HG002#1"]["all"]["f1"], 3), r(hap["HG002#2"]["all"]["f1"], 3)],
        }

    scaling = []
    for n in (6, 20, 40):
        graphs = runs[("arm", n, "p90", "1")]
        row = {"n": n}
        mono = graphs["arm.mono"]
        if mono["cost"]["max_wall_s"] is None:
            row["mono"] = dict(MONO_DNF[n])
        else:
            row["mono"] = {"status": "done", **cost(mono, n), **acc(mono)}
        for key in ("hotspots1m", "windows1m"):
            g = graphs.get("arm." + key)
            if not g:
                continue
            row[key] = {**cost(g, n), **acc(g),
                        "stitch_s": r(g["stitch_s"], 0),
                        "seams_pass": g["boundaries_pass"], "seams": g["boundaries"],
                        "paths_identical": g["haplotypes_identical_to_monolithic"],
                        "paths": g["haplotypes"]}
            if row["mono"]["status"] == "done":
                row[key]["cpu_vs_mono"] = r(g["cost"]["cpu_s"] / mono["cost"]["cpu_s"], 2)
        scaling.append(row)

    # Per-chunk layout and cost at 6 haplotypes, for the hero tracks.
    g6 = runs[("arm", 6, "p90", "1")]
    chunks, pos = [], ARM["start"]
    for cid in sorted(g6["arm.hotspots1m"]["chunk"]):
        c = g6["arm.hotspots1m"]["chunk"][cid]
        chunks.append({"start": pos, "end": pos + c["bp"],
                       "wall_min": r(c["wall_s"] / 60, 2), "peak_gb": r(c["peak_rss_mb"] / 1000, 3)})
        pos += c["bp"]
    assert pos == ARM["end"], "hotspot chunks must tile the arm exactly"

    hard = []
    for reg in ("r22", "r17"):
        for p in ("p90", "p98"):
            row = {"region": REGIONS[reg]["name"], "pggb": "-p 90" if p == "p90" else "-p 98"}
            for key, name in (("mono", "mono"), ("windows0", "chunked")):
                reps = [runs[("hard", 6, p, rep)][f"{reg}.{key}"] for rep in ("1", "2")]
                row[name] = {
                    "f1_a10k": [r(g["truth_a10k"]["all"]["total"]["f1"], 3) for g in reps],
                    "f1_a100k": [r(g["truth_a100k"]["all"]["total"]["f1"], 3) for g in reps],
                    "f1_sd_a10k": [r(g["truth_a10k"]["all"]["sd"]["f1"], 3) for g in reps],
                    "peak_gb": gb(reps[0]["cost"]["max_peak_rss_mb"]),
                    "ref_revisit_pct": r(100 * next(v for k, v in reps[0]["path_revisit_share"].items()
                                                    if k.startswith("GRCh38")), 0),
                }
                if reg == "r22":
                    hap = reps[0]["truth_a10k"]["per_hap"]
                    row[name]["hg002_f1"] = [r(hap["HG002#1"]["all"]["f1"], 3), r(hap["HG002#2"]["all"]["f1"], 3)]
            hard.append(row)

    return {"runner": RUNNER, "arm": ARM, "arm_truth_svs": ARM_TRUTH_SVS, "arm_sv_diff_calls": ARM_SV_DIFF_CALLS,
            "max_concurrent_jobs": MAX_CONCURRENT_JOBS,
            "giab_22q_coverage_pct": GIAB_22Q_COVERAGE_PCT,
            "scaling": scaling, "chunks_n6": chunks, "hard": hard}


def revisit_shares():
    """Length-weighted share of each haplotype path that revisits a node, per graph."""
    out = {}
    for rep in ("rep1", "rep2"):
        with open(os.path.join(DOCS, f"tier3_path_revisits_{rep}.tsv")) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                _, reg, _, graph, _ = row["gfa"].split("/")
                acc = out.setdefault((reg, graph, rep), [0, 0])
                acc[0] += int(row["revisit_bp"])
                acc[1] += int(row["path_bp"])
    return {k: v[0] / v[1] for k, v in out.items()}


def tier3():
    d = load("tier3.json")["runs"]
    rv = revisit_shares()
    regions = []
    for reg in ("r21", "r17", "r22"):
        row = {"key": reg, **REGIONS[reg]}
        for graph, name in (("monolithic", "mono"), ("windows0", "chunked"), ("sdgroup", "sdgroup")):
            if graph not in d[f"{reg}/rep1"]:
                continue
            row[name] = {
                "f1": [r(f1(d[f"{reg}/{rep}"][graph]), 3) for rep in ("rep1", "rep2")],
                "revisit_pct": [r(100 * rv[(reg, graph, rep)], 1) for rep in ("rep1", "rep2")],
                "sv_f1": [r(d[f"{reg}/{rep}"][graph]["truth"]["sv"]["total"]["f1"], 2) for rep in ("rep1", "rep2")],
            }
            if d[f"{reg}/rep1"][graph].get("truth_a10k"):
                row[name]["f1_a10k"] = r(f1(d[f"{reg}/rep1"][graph], "truth_a10k"), 3)
            sd = d[f"{reg}/rep1"][graph]["truth"]["all"]["sd"]["f1"]
            if sd:
                row[name]["f1_sd"] = r(sd, 3)
        regions.append(row)

    # E5: PGGB's human settings (-p 98) scored with vcfbub -a 100000 and -a 10000.
    trap = []
    for reg in ("r17", "r22"):
        g = d[f"p98/{reg}/rep1"]
        trap.append({"region": REGIONS[reg]["name"],
                     "mono_a100k": r(f1(g["monolithic"]), 3), "mono_a10k": r(f1(g["monolithic"], "truth_a10k"), 3),
                     "chunked_a10k": r(f1(g["windows0"], "truth_a10k"), 3)})

    # E4: 10 Mb of chr21, monolithic against 25 x 400 kb chunks.
    x = d["r21x/rep1"]
    e4 = {"mono_wall_s": r(x["monolithic"]["wall_s"], 0), "mono_peak_gb": gb(x["monolithic"]["peak_rss_mb"]),
          "chunk_longest_s": r(x["windows0"]["max_wall_s"], 0), "chunk_peak_mb": r(x["windows0"]["peak_rss_mb"], 0),
          "chunks": x["windows0"]["builds"]}
    return {"regions": regions, "r17_truth_svs": T3_R17_TRUTH_SVS, "vcfbub_trap": trap, "e4_10mb": e4}


def where_to_cut():
    pilot = load("recombination_chunks_comparison.json")
    mono_bp = pilot["monolithic"]["total_bp"]
    out = {"pilot": {
        "windows_extra_bp": pilot["windows"]["total_bp"] - mono_bp,
        "hotspots_extra_bp": pilot["hotspots300k"]["total_bp"] - mono_bp,
        "overlap_kb": 50,  # docs/RECOMBINATION_CHUNKS.md, Question: 400 kb windows, 50 kb overlap
    }}

    t1 = load("tier1_controls.json")
    reps = t1["replicates"]
    keys = ("windows", "windows0", "hotspots300")
    mono = [reps[k]["monolithic"] for k in ("rep1", "rep2", "rep3")]
    out["tier1"] = {
        "f1": {k: [reps[rep][k]["truvari_f1"] for rep in ("rep1", "rep2", "rep3")] for k in keys},
        "near_cut_bp": {k: [reps[rep][k]["excess_bp_within_10kb_of_cut"] for rep in ("rep1", "rep2", "rep3")]
                        for k in keys},
        "random_f1": sorted(reps["rep1"][f"random_{i:02d}"]["truvari_f1"] for i in range(1, 21)),
        "random_f1_median": t1["random_cuts"]["truvari_f1_median"],
        "mono_vs_mono_f1": sorted(t1["monolithic_vs_monolithic_truvari_f1"].values()),
        "mono_drift_nodes": max(m["nodes"] for m in mono) - min(m["nodes"] for m in mono),
        "mono_drift_bp": max(m["node_bp"] for m in mono) - min(m["node_bp"] for m in mono),
    }

    t2 = load("tier2.json")
    bug = t2["seqwish_batch_bug"]["graph_length_bp"]["r22"]
    n40 = t2["scaling_rep1"]["N40/rep1"]
    out["tier2"] = {
        "seqwish_bug_bp": bug["10M"], "seqwish_fixed_bp": bug["30000000"],
        "n40": {"mono_cpu_s": r(n40["mono"]["cpu"], 0), "chunk_cpu_s": r(n40["windows0"]["cpu_sum"], 0),
                "longest_vs_mono_pct": r(100 * n40["windows0"]["max_wall"] / n40["mono"]["wall"], 0),
                "mono_peak_gb": gb(n40["mono"]["rss"]), "chunk_peak_mb": r(n40["windows0"]["max_rss"], 0)},
        "hard_mono_vs_mono_f1": [r(t2["hard_regions_mono_rep1_vs_rep2_f1"][k], 3) for k in ("r17", "r22")],
    }
    return out


def build():
    return {
        "generated_by": "scripts/build_site_results.py",
        "tier4": tier4(),
        "tier3": tier3(),
        "where_to_cut": where_to_cut(),
    }


def render(data):
    return json.dumps(data, indent=1) + "\n"


def main(argv):
    text = render(build())
    if "--check" in argv:
        with open(OUT) as fh:
            if fh.read() != text:
                sys.exit("web/public/data/results.json is stale: run scripts/build_site_results.py")
        return
    with open(OUT, "w") as fh:
        fh.write(text)
    print("wrote", os.path.relpath(OUT, ROOT))


if __name__ == "__main__":
    main(sys.argv[1:])
