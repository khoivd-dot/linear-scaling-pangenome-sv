#!/usr/bin/env python3
"""Collect Tier 3 results into one JSON (docs/tier3.json) and print the tables used in docs/TIER3.md.

    summarize.py /home/user/work/t3/runs docs/tier3.json [--scale /home/user/work/scale/runs]
"""
import argparse, glob, json, os


def meta(d):
    p = os.path.join(d, "run_metadata.json")
    return json.load(open(p)) if os.path.exists(p) else {}


def builds(strategy_dir):
    """run_metadata of every PGGB build under a strategy (chunks, or groups' build dirs)."""
    ms = [meta(os.path.dirname(f)) for f in glob.glob(os.path.join(strategy_dir, "**", "run_metadata.json"), recursive=True)]
    return [m for m in ms if m]


def compute(ms):
    return {"cpu_s": round(sum(m.get("cpu_seconds", 0) for m in ms), 1),
            "max_wall_s": max((m.get("wall_seconds", 0) for m in ms), default=0),
            "peak_rss_mb": max((m.get("peak_rss_mb", 0) for m in ms), default=0), "builds": len(ms)}


def truth(path):
    if not os.path.exists(path):
        return None
    d = json.load(open(path))
    return {m: {z: {k: d[m][z][k] for k in ("precision", "recall", "f1", "tp_base", "fn", "fp")} for z in ("total", "sd", "other")}
            for m in ("all", "sv")}


def run_summary(rd):
    c = json.load(open(os.path.join(rd, "comparison.json")))
    mono = c["monolithic"]; mm = meta(os.path.join(rd, "monolithic"))
    out = {"monolithic": {"nodes": mono["nodes"], "bp": mono["total_bp"], "wall_s": mm.get("wall_seconds"),
                          "cpu_s": mm.get("cpu_seconds"), "peak_rss_mb": mm.get("peak_rss_mb"),
                          "truth": truth(os.path.join(rd, "truth_monolithic.json"))}}
    names = [k for k in c if isinstance(c[k], dict) and "nodes" in c[k] and k != "monolithic"]
    for extra in ("lace_windows0", "sdgroup"):
        if os.path.exists(os.path.join(rd, extra, "comparison.json")):
            names.append(extra)
    for n in names:
        e = c[n] if n in c else json.load(open(os.path.join(rd, n, "comparison.json")))
        if n == "sdgroup":
            ms = builds(os.path.join(rd, n)) + [meta(os.path.join(rd, "windows0", g.split(",")[0]))
                                                for g in e["groups"].values() if "," not in g]
        elif n == "lace_windows0":
            ms = builds(os.path.join(rd, "windows0"))
        else:
            ms = builds(os.path.join(rd, n))
        out[n] = {"nodes": e["nodes"], "bp": e["total_bp"], "excess_bp": e["total_bp"] - mono["total_bp"],
                  "haplotypes_identical": f'{e["haplotypes_identical_to_monolithic"]}/{e["haplotypes"]}',
                  "f1_vs_mono_exact": e["variants_vs_monolithic"]["f1"], **compute(ms),
                  "truth": truth(os.path.join(rd, f"truth_{n}.json"))}
        if "groups" in e:
            out[n]["groups"] = e["groups"]
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("runs"); ap.add_argument("out"); ap.add_argument("--scale")
    a = ap.parse_args()
    res = {}
    for rd in sorted(glob.glob(os.path.join(a.runs, "*", "rep*"))):
        if os.path.exists(os.path.join(rd, "comparison.json")):
            res[os.path.relpath(rd, a.runs)] = run_summary(rd)
    if a.scale:
        for rd in sorted(glob.glob(os.path.join(a.scale, "N*", "rep2"))):
            if os.path.exists(os.path.join(rd, "seam_metrics.json")):
                s = run_summary(rd)
                s["seam_metrics"] = json.load(open(os.path.join(rd, "seam_metrics.json")))
                res["scaling/" + os.path.relpath(rd, a.scale)] = s
    json.dump({"date": "2026-11-25", "runs": res}, open(a.out, "w"), indent=1)
    for k, v in res.items():
        print(f"## {k}")
        print("| strategy | nodes | bp | excess bp | haps | F1 vs mono | CPU s | max wall s | peak MB | truth F1 all (SD / other) | truth F1 SV |")
        print("|---|---|---|---|---|---|---|---|---|---|---|")
        for n, e in v.items():
            if n == "seam_metrics":
                continue
            t = e.get("truth")
            ts = f'{t["all"]["total"]["f1"]} ({t["all"]["sd"]["f1"]} / {t["all"]["other"]["f1"]})' if t else ""
            tv = str(t["sv"]["total"]["f1"]) if t else ""
            print(f'| {n} | {e["nodes"]} | {e["bp"]} | {e.get("excess_bp", "")} | {e.get("haplotypes_identical", "")} | '
                  f'{e.get("f1_vs_mono_exact", "")} | {e.get("cpu_s")} | {e.get("max_wall_s", e.get("wall_s"))} | {e.get("peak_rss_mb")} | {ts} | {tv} |')


if __name__ == "__main__":
    main()
