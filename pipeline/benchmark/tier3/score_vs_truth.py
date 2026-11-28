#!/usr/bin/env python3
"""Score a graph's variant calls against the assembly-based truth calls, per haplotype.

Graph calls: vg deconstruct -a on the GRCh38 path, vcfbub -l 0 -a 100000 (top level, bubbles <= 100 kb; --bub to change) and
vcfwave (as in the PGGB pipeline) to break bubbles into primitive alleles, then one haploid VCF per
SAMPLE#HAP, left-aligned with bcftools norm. Truth: truth_calls.py output for the same haplotype,
normalised the same way. truvari bench runs inside each haplotype's callable BED, once on all sizes
(SNVs included) and once on SVs (>= 50 bp); calls are split by whether they sit in a GIAB
segmental duplication.

    score_vs_truth.py GFA OUT_JSON --chrom chr21 --region-start 20000000 --ref chr21.chr.fa \
        --truth truth/r21 --haps HG00673#1 HG00673#2 ... --segdups segdups.bed --truvari truvari
"""
import argparse, bisect, gzip, json, os, subprocess, tempfile
import pysam

IMG = "ghcr.io/pangenome/pggb@sha256:44f3563c5c3e2032c293fe309eed178b754021b07963ca4857e4ae5efda86474"


def docker(cmd, mounts):
    args = ["docker", "run", "--rm"]
    for m in mounts:
        args += ["-v", f"{m}:{m}"]
    return subprocess.run(args + [IMG, "bash", "-c", cmd], check=True, capture_output=True, text=True).stdout


def graph_vcf(gfa, out, chrom, region_start, ref, bub="-l 0 -a 100000"):
    d = os.path.dirname(os.path.abspath(gfa))
    raw = docker(f"vg deconstruct -P GRCh38 -a {os.path.abspath(gfa)} 2>/dev/null > /tmp/d.vcf && "
                 f"vcfbub {bub} --input /tmp/d.vcf > /tmp/b.vcf && vcfwave -I 1000 /tmp/b.vcf 2>/dev/null", [d])
    with open(out, "w") as fh:
        for line in raw.splitlines():
            if line.startswith("##contig"):
                continue
            if line.startswith("#"):
                fh.write(line + "\n"); continue
            f = line.split("\t")
            fh.write("\t".join([chrom, str(int(f[1]) + region_start)] + f[2:]) + "\n")
    return out


def haploid(vcf_in, sample, hap, out):
    recs, header = [], []
    for line in open(vcf_in):
        if line.startswith("##"):
            header.append(line); continue
        f = line.rstrip("\n").split("\t")
        if line.startswith("#"):
            idx = f.index(sample); continue
        gt = f[idx].split(":")[0].replace("/", "|").split("|")
        if len(gt) < hap or gt[hap - 1] in (".", "0"):
            continue
        alt = f[4].split(",")[int(gt[hap - 1]) - 1]
        if alt == "*" or set(alt) - set("ACGTN"):
            continue
        recs.append((int(f[1]), f[3], alt))
    write_vcf(out, header, sorted(set(recs)), f[0] if recs else "chr")
    return out


def write_vcf(out, header, recs, chrom):
    with open(out, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
        for p, r, a in recs:
            fh.write(f"{chrom}\t{p}\t.\t{r}\t{a}\t60\tPASS\t.\tGT\t1\n")


def norm(vcf, ref, chrom):
    """Left-align, split, dedup with bcftools; then sort, bgzip, tabix with pysam."""
    d = os.path.dirname(os.path.abspath(vcf)); rd = os.path.dirname(os.path.abspath(ref))
    fai = open(ref + ".fai").readline().split("\t")
    hdr = vcf + ".h.vcf"
    with open(vcf) as fi, open(hdr, "w") as fo:
        for line in fi:
            fo.write(line)
            if line.startswith("##fileformat"):
                fo.write(f"##contig=<ID={fai[0]},length={fai[1]}>\n")
    out = docker(f"bcftools norm -m- -f {os.path.abspath(ref)} -c s {os.path.abspath(hdr)} 2>/dev/null"
                 f" | bcftools norm -d exact 2>/dev/null", list({d, rd}))
    lines = out.splitlines()
    head = [l for l in lines if l.startswith("#")]
    body = sorted((l for l in lines if l and not l.startswith("#")), key=lambda l: int(l.split("\t")[1]))
    plain = vcf.replace(".vcf", ".norm.vcf")
    open(plain, "w").write("\n".join(head + body) + "\n")
    pysam.tabix_index(plain, preset="vcf", force=True)
    return plain + ".gz"


def load_bed(path, chrom):
    iv = sorted((int(l.split()[1]), int(l.split()[2])) for l in open(path) if l.split()[0] == chrom)
    return [s for s, _ in iv], [e for _, e in iv]


def in_bed(pos0, bed):
    i = bisect.bisect_right(bed[0], pos0) - 1
    return i >= 0 and pos0 < bed[1][i]


def count(vcf_gz, sd):
    c = {"sd": 0, "other": 0}
    if os.path.exists(vcf_gz):
        for line in gzip.open(vcf_gz, "rt"):
            if not line.startswith("#"):
                c["sd" if in_bed(int(line.split("\t")[1]) - 1, sd) else "other"] += 1
    return c


def bench(truvari, base, comp, bed, outdir, sv):
    args = [truvari, "bench", "-b", base, "-c", comp, "-o", outdir, "--includebed", bed,
            "--pick", "multi", "--refdist", "500", "--pctseq", "0.7", "--pctsize", "0.7"]
    args += ["--sizemin", "50", "--sizefilt", "30"] if sv else ["--sizemin", "0", "--sizefilt", "0"]
    subprocess.run(args, check=True, capture_output=True)
    parts = {k: count(os.path.join(outdir, f"{k}.vcf.gz"), SD) for k in ("tp-base", "fn", "tp-comp", "fp")}
    return parts


def prf(tpb, fn, tpc, fp):
    r = tpb / (tpb + fn) if tpb + fn else 0.0
    p = tpc / (tpc + fp) if tpc + fp else 0.0
    return {"tp_base": tpb, "fn": fn, "tp_comp": tpc, "fp": fp, "precision": round(p, 4),
            "recall": round(r, 4), "f1": round(2 * p * r / (p + r), 4) if p + r else 0.0}


def main(argv=None):
    global SD
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("gfa"); ap.add_argument("out_json")
    ap.add_argument("--chrom", required=True); ap.add_argument("--region-start", type=int, required=True)
    ap.add_argument("--ref", required=True); ap.add_argument("--truth", required=True, help="prefix: PREFIX_SAMPLE_HAP.vcf/.bed")
    ap.add_argument("--haps", nargs="+", required=True); ap.add_argument("--segdups", required=True)
    ap.add_argument("--truvari", required=True)
    ap.add_argument("--truth-sv", help="prefix for the SV (>= 50 bp) scoring; default --truth (Tier 4: GIAB stvar regions)")
    ap.add_argument("--bub", default="-l 0 -a 100000",
                    help="vcfbub options; the default is the PGGB pipeline's (sensitivity check: '-l 0 -a 10000')")
    a = ap.parse_args(argv)
    SD = load_bed(a.segdups, a.chrom)
    work = os.path.splitext(a.out_json)[0] + "_work"
    os.makedirs(work, exist_ok=True)
    gv = graph_vcf(a.gfa, os.path.join(work, "graph.vcf"), a.chrom, a.region_start, a.ref, a.bub)
    result = {"gfa": a.gfa, "vcfbub": a.bub, "haps": {}}
    totals = {m: {k: {"sd": 0, "other": 0} for k in ("tp-base", "fn", "tp-comp", "fp")} for m in ("all", "sv")}
    for h in a.haps:
        s, hp = h.split("#")
        t = f"{a.truth}_{s}_{hp}"
        tv = os.path.join(work, f"truth_{s}_{hp}.vcf")
        open(tv, "w").write(open(t + ".vcf").read())
        base = norm(tv, a.ref, a.chrom)
        comp = norm(haploid(gv, s, int(hp), os.path.join(work, f"graph_{s}_{hp}.vcf")), a.ref, a.chrom)
        result["haps"][h] = {}
        for m in ("all", "sv"):
            od = os.path.join(work, f"truvari_{m}_{s}_{hp}")
            subprocess.run(["rm", "-rf", od])
            bed = f"{a.truth_sv}_{s}_{hp}.bed" if m == "sv" and a.truth_sv else t + ".bed"
            parts = bench(a.truvari, base, comp, bed, od, m == "sv")
            result["haps"][h][m] = {z: prf(parts["tp-base"][z], parts["fn"][z], parts["tp-comp"][z], parts["fp"][z])
                                    for z in ("sd", "other")}
            for k in parts:
                for z in ("sd", "other"):
                    totals[m][k][z] += parts[k][z]
    for m in totals:
        T = totals[m]
        result[m] = {z: prf(T["tp-base"][z], T["fn"][z], T["tp-comp"][z], T["fp"][z]) for z in ("sd", "other")}
        result[m]["total"] = prf(*(T[k]["sd"] + T[k]["other"] for k in ("tp-base", "fn", "tp-comp", "fp")))
    json.dump(result, open(a.out_json, "w"), indent=1)
    print(a.out_json, json.dumps({m: result[m]["total"]["f1"] for m in ("all", "sv")}))


if __name__ == "__main__":
    main()
