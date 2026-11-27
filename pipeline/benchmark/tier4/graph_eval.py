#!/usr/bin/env python3
"""Tier 4: assemble one graph from the GitHub build jobs and describe it.

For a chunked strategy the chunk graphs are stitched end to end (the project's overlap-aware stitch,
zero overlap here); for `mono` the single build is used as is. Writes GRAPH.gfa and a JSON with
build cost per chunk (wall, CPU, peak memory from GNU time), graph statistics, seam checks and,
when the monolithic graph of the same run is given, how many haplotypes spell the same sequence.

    graph_eval.py BUILDS_DIR REGION STRATEGY OUT_DIR [--mono-gfa MONO.gfa]
BUILDS_DIR/<REGION>.<STRATEGY>.<chunk_id>/ holds final.gfa and run_metadata.json for each chunk.
"""
import argparse
import csv
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..", "..")
sys.path.insert(0, ROOT)
from pipeline.benchmark.graph_stats import compute_stats  # noqa: E402
from pipeline.merge.gfa import GfaGraph  # noqa: E402
from pipeline.merge.merge_graphs import overlap_aware_stitch, walks_as_paths  # noqa: E402


def revisits(g):
    """Share of each path's bp spent re-entering nodes it already used (Tier 3 path_revisits.py)."""
    seglen = {k: len(v.sequence) for k, v in g.segments.items()}
    out = {}
    for name, p in g.paths.items():
        seen, tot, rep = set(), 0, 0
        for step in p.segment_names:
            n = step.rstrip("+-")
            tot += seglen.get(n, 0)
            if n in seen:
                rep += seglen.get(n, 0)
            seen.add(n)
        out[name] = round(rep / tot, 4) if tot else 0.0
    return out


def haplotype_sequences(g):
    return {tuple(n.split("#")[:2]): g.get_path_sequence(n) for n in g.paths}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("builds"); ap.add_argument("region"); ap.add_argument("strategy"); ap.add_argument("out")
    ap.add_argument("--mono-gfa")
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)
    manifest = os.path.join(HERE, "manifests", f"{a.region}.{a.strategy}.tsv")
    rows = {r["chunk_id"]: r for r in csv.DictReader(open(manifest), delimiter="\t")}
    res = {"region": a.region, "strategy": a.strategy, "chunks": len(rows), "chunk": {}, "missing": []}
    graphs = []
    for cid, r in rows.items():
        d = os.path.join(a.builds, f"{a.region}.{a.strategy}.{cid}")
        gfa, meta = os.path.join(d, "final.gfa"), os.path.join(d, "run_metadata.json")
        if not os.path.exists(gfa):
            res["missing"].append(cid)
            continue
        m = json.load(open(meta))
        res["chunk"][cid] = {"bp": int(r["reference_end"]) - int(r["reference_start"]), "paths": m.get("num_paths"),
                             "wall_s": round(m["wall_seconds"], 1), "cpu_s": m.get("cpu_seconds"),
                             "peak_rss_mb": m.get("peak_rss_mb"), "vcpus": m.get("vcpus")}
        graphs.append((cid, gfa))
    c = res["chunk"].values()
    res["cost"] = {"cpu_s": round(sum(x["cpu_s"] or 0 for x in c), 1), "sum_wall_s": round(sum(x["wall_s"] for x in c), 1),
                   "max_wall_s": max((x["wall_s"] for x in c), default=None),
                   "max_peak_rss_mb": max((x["peak_rss_mb"] or 0 for x in c), default=None)}
    out_gfa = os.path.join(a.out, "graph.gfa")
    if res["missing"]:
        print(f"{a.region}.{a.strategy}: missing chunks {res['missing']}; no graph")
    else:
        t0 = time.time()
        if a.strategy == "mono":
            g = walks_as_paths(GfaGraph.parse_file(graphs[0][1]))
            os.replace(graphs[0][1], out_gfa)
        else:
            chunks = [(cid, walks_as_paths(GfaGraph.parse_file(p))) for cid, p in graphs]
            g, bounds = overlap_aware_stitch(chunks, chunk_rows=rows)
            del chunks
            g.write_gfa(out_gfa)
            res["boundaries"] = len(bounds)
            res["boundaries_pass"] = sum(b["status"] == "PASS" for b in bounds)
        res["stitch_s"] = round(time.time() - t0, 1)
        res["stats"] = compute_stats(g, f"{a.region}.{a.strategy}")
        res["path_revisit_share"] = revisits(g)
        if a.mono_gfa and a.strategy != "mono" and os.path.exists(a.mono_gfa):
            mine = haplotype_sequences(g)
            del g
            mono = haplotype_sequences(walks_as_paths(GfaGraph.parse_file(a.mono_gfa)))
            res["haplotypes"] = len(mono)
            res["haplotypes_identical_to_monolithic"] = sum(mine.get(k) == v for k, v in mono.items())
    json.dump(res, open(os.path.join(a.out, "graph.json"), "w"), indent=1, default=str)
    print(json.dumps({k: res.get(k) for k in ("region", "strategy", "chunks", "missing", "cost", "boundaries_pass",
                                              "haplotypes_identical_to_monolithic")}, default=str))


if __name__ == "__main__":
    main()
