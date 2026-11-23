"""compare_chunking.py — monolithic vs chunked-and-stitched PGGB on the same input.

Builds one monolithic graph and, for each chunking strategy, one graph per
chunk; stitches the chunks; then compares every stitched graph with the
monolithic one on topology, haplotype paths and variants.

    python3 pipeline/benchmark/compare_chunking.py WORKDIR --mono MONO.fa \
        --strategy windows=MANIFEST.tsv=FASTA_DIR --strategy hotspots=... [--threads 4]

Each chunk FASTA is FASTA_DIR/<chunk_id>.fa (pipeline/prepare/project_chunks.py).
PGGB parameters come from config/pipeline.yaml (or PIPELINE_CONFIG), the same
for every build. Builds are skipped when their final.gfa already exists.

Variants: `vg deconstruct -P GRCh38` on each graph (pinned PGGB image), matched
exactly on (POS, REF, ALT) with the monolithic call set as truth. This is a
quick exact-match check, not a replacement for the truvari comparison.

Writes WORKDIR/comparison.json and WORKDIR/<strategy>/merged.gfa.
"""
import argparse
import csv
import json
import os
import subprocess
import sys

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from gen_pggb_config import load_config  # noqa: E402
from pipeline.benchmark.graph_stats import compute_stats  # noqa: E402
from pipeline.merge.gfa import GfaGraph  # noqa: E402
from pipeline.merge.merge_graphs import overlap_aware_stitch, walks_as_paths  # noqa: E402


def build(fasta, outdir, threads):
    """PGGB via scripts/run_pggb.py; returns (final.gfa, wall seconds)."""
    gfa, meta = os.path.join(outdir, "final.gfa"), os.path.join(outdir, "run_metadata.json")
    if not os.path.exists(gfa):
        subprocess.run([sys.executable, os.path.join(ROOT, "scripts", "run_pggb.py"),
                        fasta, outdir, str(threads)], check=True, cwd=ROOT,
                       stdout=subprocess.DEVNULL)
    return gfa, json.load(open(meta))["wall_seconds"]


def variants(gfa, region_start):
    """{(pos, ref, alt)} on the GRCh38 path, positions on the chromosome (1-based)."""
    image = load_config()["image"]
    d, name = os.path.split(os.path.abspath(gfa))
    out = subprocess.run(["docker", "run", "--rm", "-v", f"{d}:/g", image, "bash", "-c",
                          f"vg deconstruct -P GRCh38 -a /g/{name} 2>/dev/null"],
                         check=True, capture_output=True, text=True).stdout
    calls = set()
    for line in out.splitlines():
        if line.startswith("#"):
            continue
        f = line.split("\t")
        for alt in f[4].split(","):
            calls.add((int(f[1]) + region_start, f[3], alt))
    return calls


def score(truth, test):
    tp = len(truth & test)
    p = tp / len(test) if test else 0.0
    r = tp / len(truth) if truth else 0.0
    return {"truth": len(truth), "test": len(test), "tp": tp,
            "precision": round(p, 4), "recall": round(r, 4),
            "f1": round(2 * p * r / (p + r), 4) if p + r else 0.0}


def by_class(calls):
    def cls(ref, alt):
        d = abs(len(ref) - len(alt))
        return "snv" if len(ref) == len(alt) == 1 else ("sv_50bp+" if d >= 50 else "indel_or_mnp")
    out = {}
    for pos, ref, alt in calls:
        out.setdefault(cls(ref, alt), set()).add((pos, ref, alt))
    return out


def haplotype_sequences(g):
    """{(sample, hap): spelled sequence} with subranges and _rc suffixes ignored."""
    seqs = {}
    for name in g.paths:
        sample, hap = name.split("#")[:2]
        seqs[(sample, hap)] = g.get_path_sequence(name)
    return seqs


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("workdir")
    ap.add_argument("--mono", required=True)
    ap.add_argument("--strategy", action="append", required=True, help="NAME=MANIFEST=FASTA_DIR")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--region-start", type=int, default=20000000)
    ap.add_argument("--skip-variants", action="store_true")
    a = ap.parse_args(argv)
    os.makedirs(a.workdir, exist_ok=True)

    mono_gfa, mono_wall = build(os.path.abspath(a.mono), os.path.join(a.workdir, "monolithic"), a.threads)
    mono = walks_as_paths(GfaGraph.parse_file(mono_gfa))
    mono_seqs = haplotype_sequences(mono)
    result = {"pggb_params": {k: v for k, v in load_config().items() if k != "image"},
              "threads_per_build": a.threads,
              "monolithic": {"wall_s": mono_wall, **compute_stats(mono, "monolithic")}}
    truth = None if a.skip_variants else variants(mono_gfa, a.region_start)
    if truth is not None:
        result["monolithic"]["variants"] = {k: len(v) for k, v in by_class(truth).items()}

    for spec in a.strategy:
        name, manifest, fasta_dir = spec.split("=")
        with open(manifest) as fh:
            rows = {r["chunk_id"]: r for r in csv.DictReader(fh, delimiter="\t")}
        chunks, walls = [], {}
        for cid in rows:
            gfa, walls[cid] = build(os.path.abspath(os.path.join(fasta_dir, f"{cid}.fa")),
                                    os.path.join(a.workdir, name, cid), a.threads)
            chunks.append((cid, GfaGraph.parse_file(gfa)))
        merged, boundaries = overlap_aware_stitch(chunks, chunk_rows=rows)
        merged_path = os.path.join(a.workdir, name, "merged.gfa")
        merged.write_gfa(merged_path)

        seqs = haplotype_sequences(merged)
        same = {k: seqs.get(k) == v for k, v in mono_seqs.items()}
        entry = {
            "chunks": len(rows),
            "chunk_bp": [int(r["reference_end"]) - int(r["reference_start"]) for r in rows.values()],
            "chunk_wall_s": walls,
            "sum_wall_s": round(sum(walls.values()), 1),
            "max_chunk_wall_s": max(walls.values()),
            "boundaries_pass": sum(b["status"] == "PASS" for b in boundaries),
            "boundaries": len(boundaries),
            "haplotypes_identical_to_monolithic": sum(same.values()),
            "haplotypes": len(same),
            **compute_stats(merged, name),
        }
        if truth is not None:
            test = variants(merged_path, a.region_start)
            entry["variants_vs_monolithic"] = score(truth, test)
            tc, sc = by_class(truth), by_class(test)
            entry["variants_by_class"] = {c: score(tc.get(c, set()), sc.get(c, set()))
                                          for c in sorted(set(tc) | set(sc))}
        result[name] = entry
        print(f"{name}: {entry['chunks']} chunks, nodes {entry['nodes']}, "
              f"boundaries {entry['boundaries_pass']}/{entry['boundaries']} PASS, "
              f"haplotypes identical {entry['haplotypes_identical_to_monolithic']}/{entry['haplotypes']}")

    with open(os.path.join(a.workdir, "comparison.json"), "w") as fh:
        json.dump(result, fh, indent=2, default=str)
    print(f"-> {a.workdir}/comparison.json")


if __name__ == "__main__":
    main()
