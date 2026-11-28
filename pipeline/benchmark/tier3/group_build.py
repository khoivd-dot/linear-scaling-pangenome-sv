#!/usr/bin/env python3
"""Build chunk groups (sd_groups.py) as single PGGB jobs, join everything with impg lace, compare to monolithic.

A group of non-adjacent chunks is built from the concatenation of its members' FASTAs, so
paralogous copies in different chunks are aligned to each other as in a monolithic build. Singleton
groups reuse the per-chunk graphs of a finished compare_chunking.py strategy when --reuse is given.
All graphs are joined with `impg lace` (subrange path names), the joined graph is compared with the
monolithic graph of the same run the way compare_chunking.py does it.

    group_build.py RUNDIR NAME GROUPS.tsv FASTA_DIR --reuse windows0 --region-start 18600000 --chrom chr22
Writes RUNDIR/NAME/{group_XX/, laced.gfa, comparison.json}.
"""
import argparse, csv, json, os, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(HERE, "..", "..", ".."))
from compare_chunking import build, variants, score, by_class, haplotype_sequences, load_config  # noqa: E402
from pipeline.benchmark.graph_stats import compute_stats  # noqa: E402
from pipeline.merge.gfa import GfaGraph  # noqa: E402
from pipeline.merge.merge_graphs import walks_as_paths  # noqa: E402


def lace(gfas, out):
    image = load_config()["image"]
    d = os.path.dirname(out)
    mounts = sorted({os.path.dirname(g) for g in gfas} | {d})
    args = ["docker", "run", "--rm"] + sum((["-v", f"{m}:{m}"] for m in mounts), []) + [image]
    subprocess.run(args + ["impg", "lace", "-f", *gfas, "-o", out, "--compress", "none", "-t", "4", "-v", "0"],
                   check=True, capture_output=True)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("rundir"); ap.add_argument("name"); ap.add_argument("groups"); ap.add_argument("fasta_dir")
    ap.add_argument("--reuse", help="strategy in RUNDIR whose chunk graphs serve the singleton groups")
    ap.add_argument("--region-start", type=int, required=True); ap.add_argument("--chrom", required=True)
    ap.add_argument("--threads", type=int, default=4)
    a = ap.parse_args(argv)
    run = os.path.abspath(a.rundir); out = os.path.join(run, a.name); os.makedirs(out, exist_ok=True)
    groups = list(csv.DictReader(open(a.groups), delimiter="\t"))
    gfas, walls, cpu = [], {}, {}
    for g in groups:
        members = g["chunk_ids"].split(",")
        if len(members) == 1 and a.reuse:
            bdir = os.path.join(run, a.reuse, members[0])
            gfa = os.path.join(bdir, "final.gfa")
        else:
            gdir = os.path.join(out, g["group"]); os.makedirs(gdir, exist_ok=True)
            fa = os.path.join(gdir, "input.fa")
            if not os.path.exists(fa):
                with open(fa + ".tmp", "w") as fh:
                    for m in members:
                        fh.write(open(os.path.join(a.fasta_dir, f"{m}.fa")).read())
                os.rename(fa + ".tmp", fa)
            bdir = os.path.join(gdir, "build")
            gfa, _ = build(fa, bdir, a.threads)
        meta = json.load(open(os.path.join(bdir, "run_metadata.json")))
        walls[g["group"]] = meta["wall_seconds"]
        cpu[g["group"]] = {k: meta.get(k) for k in ("cpu_seconds", "peak_rss_mb")}
        gfas.append(os.path.abspath(gfa))
    laced = lace(gfas, os.path.join(out, "laced.gfa"))

    mono_gfa = os.path.join(run, "monolithic", "final.gfa")
    mono = walks_as_paths(GfaGraph.parse_file(mono_gfa))
    merged = walks_as_paths(GfaGraph.parse_file(laced))
    ms, ss = haplotype_sequences(mono), haplotype_sequences(merged)
    truth = variants(mono_gfa, a.region_start, a.chrom)
    test = variants(laced, a.region_start, a.chrom)
    tc, sc = by_class(truth), by_class(test)
    res = {"groups": {g["group"]: g["chunk_ids"] for g in groups}, "group_wall_s": walls, "group_cpu_rss": cpu,
           "sum_wall_s": round(sum(walls.values()), 1), "max_group_wall_s": max(walls.values()),
           "haplotypes_identical_to_monolithic": sum(ss.get(k) == v for k, v in ms.items()), "haplotypes": len(ms),
           **compute_stats(merged, a.name),
           "variants_vs_monolithic": score(truth, test),
           "variants_by_class": {c: score(tc.get(c, set()), sc.get(c, set())) for c in sorted(set(tc) | set(sc))}}
    json.dump(res, open(os.path.join(out, "comparison.json"), "w"), indent=2, default=str)
    print(a.name, json.dumps({k: res[k] for k in ("haplotypes_identical_to_monolithic", "sum_wall_s", "max_group_wall_s")}),
          res["variants_vs_monolithic"])


if __name__ == "__main__":
    main()
