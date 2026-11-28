#!/usr/bin/env python3
"""Collect one Tier 4 run (graph_eval.py and score_vs_truth.py outputs) into a compact JSON.

    summarize.py RESULTS_DIR OUT.json --set arm --n 20 --pggb p90 --rep 1
RESULTS_DIR/<region>.<strategy>/ holds graph.json and truth_a100k.json / truth_a10k.json.
"""
import argparse
import glob
import json
import os

STATS = ("nodes", "edges", "total_bp", "components", "paths")


def truth(path):
    if not os.path.exists(path):
        return None
    t = json.load(open(path))
    out = {m: {z: t[m][z] for z in ("total", "sd", "other")} for m in ("all", "sv")}
    out["per_hap"] = {}
    for h, v in t["haps"].items():
        o = {}
        for m in ("all", "sv"):
            s = {k: v[m]["sd"][k] + v[m]["other"][k] for k in ("tp_base", "fn", "tp_comp", "fp")}
            p = s["tp_comp"] / (s["tp_comp"] + s["fp"]) if s["tp_comp"] + s["fp"] else 0.0
            r = s["tp_base"] / (s["tp_base"] + s["fn"]) if s["tp_base"] + s["fn"] else 0.0
            o[m] = {**s, "f1": round(2 * p * r / (p + r), 4) if p + r else 0.0}
        out["per_hap"][h] = o
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("results"); ap.add_argument("out")
    ap.add_argument("--set", required=True); ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--pggb", default="p90"); ap.add_argument("--rep", default="1")
    a = ap.parse_args(argv)
    run = {"set": a.set, "n_haplotypes": a.n, "pggb": a.pggb, "rep": a.rep, "graphs": {}}
    for d in sorted(glob.glob(os.path.join(a.results, "*.*"))):
        gj = os.path.join(d, "graph.json")
        if not os.path.isdir(d) or not os.path.exists(gj):
            continue
        g = json.load(open(gj))
        st = g.get("stats", {})
        run["graphs"][os.path.basename(d)] = {
            "chunks": g["chunks"], "missing": g["missing"], "cost": g["cost"],
            "chunk": g["chunk"], "stitch_s": g.get("stitch_s"),
            "boundaries_pass": g.get("boundaries_pass"), "boundaries": g.get("boundaries"),
            "haplotypes_identical_to_monolithic": g.get("haplotypes_identical_to_monolithic"),
            "haplotypes": g.get("haplotypes"),
            "stats": {k: st.get(k) for k in STATS if k in st},
            "path_revisit_share": g.get("path_revisit_share"),
            "truth_a100k": truth(os.path.join(d, "truth_a100k.json")),
            "truth_a10k": truth(os.path.join(d, "truth_a10k.json")),
        }
    json.dump(run, open(a.out, "w"), indent=1)
    for k, v in run["graphs"].items():
        t1, t2 = v["truth_a100k"], v["truth_a10k"]
        f = lambda t, m: t[m]["total"]["f1"] if t else None  # noqa: E731
        print(f"{k}: cpu {v['cost']['cpu_s']} s, max wall {v['cost']['max_wall_s']} s, max rss {v['cost']['max_peak_rss_mb']} MB, "
              f"F1 all/sv a100k {f(t1, 'all')}/{f(t1, 'sv')} a10k {f(t2, 'all')}/{f(t2, 'sv')}, missing {v['missing']}")


if __name__ == "__main__":
    main()
