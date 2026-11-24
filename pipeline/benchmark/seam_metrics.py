"""seam_metrics.py — where does a stitched graph differ from the monolithic one?

For each strategy in a compare_chunking.py workdir, places every node on the
reference (GRCh38) coordinate of the last reference node before it on any
haplotype path, then reports node count and node sequence (bp) in excess of the
monolithic graph: in total, within --window bp of each cut, and elsewhere.
With --truvari, also runs `truvari bench` on all sizes, SNVs included (monolithic VCF as base) and splits
false negatives and false positives by distance to the nearest cut.

    python3 pipeline/benchmark/seam_metrics.py WORKDIR --manifests NAME=MANIFEST.tsv ...
        [--window 10000] [--truvari /path/to/truvari] [--region-start 20000000]
Writes WORKDIR/seam_metrics.json.
"""
import argparse
import bisect
import csv
import gzip
import json
import os
import shutil
import subprocess
import sys

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
sys.path.insert(0, ROOT)
from pipeline.merge.gfa import GfaGraph  # noqa: E402
from pipeline.merge.merge_graphs import walks_as_paths  # noqa: E402


def node_positions(g, region_start):
    """{node: reference coordinate}; non-reference nodes take the last ref position seen."""
    ref = next(n for n in g.paths if n.startswith("GRCh38"))
    pos, off = {}, region_start
    for node, _o in g.path_steps(ref):
        pos.setdefault(node, off)
        off += g.segments[node].length
    for name in g.paths:
        if name == ref:
            continue
        last = region_start
        for node, _o in g.path_steps(name):
            if node in pos:
                last = pos[node]
            else:
                pos[node] = last
    return pos


def profile(g, region_start, cuts, window):
    """(nodes, bp) near cuts and elsewhere."""
    pos = node_positions(g, region_start)
    near = [0, 0]
    far = [0, 0]
    for node, p in pos.items():
        i = bisect.bisect_left(cuts, p)
        d = min([abs(p - cuts[j]) for j in (i - 1, i) if 0 <= j < len(cuts)] or [10**12])
        bucket = near if d <= window else far
        bucket[0] += 1
        bucket[1] += g.segments[node].length
    return near, far


def cuts_of(manifest):
    with open(manifest) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return sorted({int(r["core_start"]) for r in rows[1:]})


def split_multiallelic(vcf, out):
    """One ALT per record (truvari refuses multi-allelic lines); GT recoded per ALT."""
    with open(vcf) as fh, open(out, "w") as w:
        for line in fh:
            if line.startswith("#"):
                w.write(line)
                continue
            f = line.rstrip("\n").split("\t")
            alts = f[4].split(",")
            for i, alt in enumerate(alts, 1):
                g = f[:4] + [alt] + f[5:7] + ["."] + f[8:9]
                for gt in f[9:]:
                    sep = "|" if "|" in gt else "/"
                    g.append(sep.join("1" if a == str(i) else ("." if a == "." else "0")
                                      for a in gt.split(sep)))
                w.write("\t".join(g) + "\n")


def bgzip_index(vcf, truvari_bin):
    """Split, sort, bgzip and tabix a VCF with pysam from truvari's own environment."""
    split_multiallelic(vcf, vcf + ".split")
    vcf = vcf + ".split"
    py = os.path.join(os.path.dirname(truvari_bin), "python3")
    code = ("import pysam,sys\n"
            "src=sys.argv[1]; hdr=[]; body=[]\n"
            "for l in open(src):\n"
            "    (hdr if l.startswith('#') else body).append(l)\n"
            "body.sort(key=lambda l:int(l.split('\\t')[1]))\n"
            "h=[x for x in hdr if not x.startswith('#CHROM')]\n"
            "chrom=body[0].split('\\t')[0] if body else 'chr21'\n"
            "length={'chr17':83257441,'chr21':46709983,'chr22':50818468}.get(chrom,250000000)\n"
            "h.append(f'##contig=<ID={chrom},length={length}>\\n')\n"
            "h+= [x for x in hdr if x.startswith('#CHROM')]\n"
            "open(src+'.s','w').writelines(h+body)\n"
            "pysam.tabix_compress(src+'.s',src+'.gz',force=True)\n"
            "pysam.tabix_index(src+'.gz',preset='vcf',force=True)\n")
    subprocess.run([py, "-c", code, vcf], check=True)
    return vcf + ".gz"


def truvari_by_distance(base_vcf, comp_vcf, outdir, cuts, truvari_bin, window):
    """truvari bench (all sizes, sequence-aware); FN/FP split by distance to a cut."""
    if os.path.exists(outdir):
        shutil.rmtree(outdir)
    subprocess.run([truvari_bin, "bench", "-b", base_vcf, "-c", comp_vcf, "-o", outdir,
                    "--sizemin", "0", "--sizefilt", "0", "--pick", "multi",
                    "--refdist", "500", "--pctseq", "0.7", "--pctsize", "0.7"],
                   check=True, capture_output=True)
    summ = json.load(open(os.path.join(outdir, "summary.json")))
    out = {k: summ.get(k) for k in ("TP-base", "FN", "FP", "precision", "recall", "f1")}

    def dist_split(vcf_gz):
        near = far = 0
        sv_near = sv_far = 0
        with gzip.open(vcf_gz, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.split("\t")
                p = int(f[1])
                d = min(abs(p - c) for c in cuts) if cuts else 10**12
                sv = max(abs(len(f[3]) - len(a)) for a in f[4].split(",")) >= 50
                if d <= window:
                    near += 1
                    sv_near += sv
                else:
                    far += 1
                    sv_far += sv
        return {"near_cut": near, "elsewhere": far, "sv_near_cut": sv_near, "sv_elsewhere": sv_far}
    out["fn_by_distance"] = dist_split(os.path.join(outdir, "fn.vcf.gz"))
    out["fp_by_distance"] = dist_split(os.path.join(outdir, "fp.vcf.gz"))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("workdir")
    ap.add_argument("--manifests", nargs="+", required=True, help="NAME=MANIFEST.tsv")
    ap.add_argument("--window", type=int, default=10000)
    ap.add_argument("--region-start", type=int, default=20000000)
    ap.add_argument("--truvari")
    a = ap.parse_args(argv)

    mono = walks_as_paths(GfaGraph.parse_file(os.path.join(a.workdir, "monolithic", "final.gfa")))
    base_vcf = None
    if a.truvari:
        base_vcf = bgzip_index(os.path.join(a.workdir, "monolithic", "final.vcf"), a.truvari)
    result = {"window_bp": a.window}
    for spec in a.manifests:
        name, manifest = spec.split("=")
        merged = os.path.join(a.workdir, name, "merged.gfa")
        if not os.path.exists(merged):
            continue
        cuts = cuts_of(manifest)
        g = walks_as_paths(GfaGraph.parse_file(merged))
        (cn, cb), (fn_, fb) = profile(g, a.region_start, cuts, a.window)
        (mn, mb), (mfn, mfb) = profile(mono, a.region_start, cuts, a.window)
        entry = {"cuts": cuts,
                 "excess_nodes_near_cut": cn - mn, "excess_bp_near_cut": cb - mb,
                 "excess_nodes_elsewhere": fn_ - mfn, "excess_bp_elsewhere": fb - mfb,
                 "excess_nodes_total": g.node_count() - mono.node_count(),
                 "excess_bp_total": g.total_sequence_bp() - mono.total_sequence_bp()}
        if a.truvari:
            comp = bgzip_index(os.path.join(a.workdir, name, "merged.vcf"), a.truvari)
            entry["truvari"] = truvari_by_distance(base_vcf, comp, os.path.join(a.workdir, name, "truvari"),
                                                   cuts, a.truvari, a.window)
        result[name] = entry
        print(name, {k: v for k, v in entry.items() if k.startswith("excess")})
    with open(os.path.join(a.workdir, "seam_metrics.json"), "w") as fh:
        json.dump(result, fh, indent=2)


if __name__ == "__main__":
    main()
