"""project_chunks.py — cut each haplotype at the same reference coordinates.

For every haplotype, one minimap2 alignment of its assembly against the
reference region (PAF with a cg:Z CIGAR) projects every chunk boundary from
reference to haplotype coordinates. Adjacent chunks therefore get haplotype
intervals that meet exactly (or overlap exactly as much as the reference
windows do), which is what the stitch needs.

Records are named PanSN with a subrange, `SAMPLE#HAP#CONTIG:start-end`, in the
haplotype's own orientation along the reference. A haplotype aligned on the
minus strand is written reverse-complemented and named `CONTIG_rc`, with
coordinates counted on the reverse-complemented contig.

Usage:
    python3 pipeline/prepare/project_chunks.py MANIFEST.tsv REGION.fa OUTDIR \
        --region-start 20000000 --ref-name GRCh38#0#chr21 \
        --hap SAMPLE#HAP=ALIGN.paf=ASSEMBLY.fa.gz [--hap ...]
Writes OUTDIR/<chunk_id>.fa for each manifest row and OUTDIR/projection.tsv.
"""
import argparse
import csv
import re
import subprocess

CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


def _read_paf(paf_path, min_mapq):
    with open(paf_path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 12 or int(p[11]) < min_mapq or "tp:A:P" not in p[12:]:
                continue
            cg = next((t[5:] for t in p[12:] if t.startswith("cg:Z:")), None)
            if cg:
                yield {"qname": p[0], "qlen": int(p[1]), "qs": int(p[2]), "qe": int(p[3]),
                       "strand": p[4], "tlen": int(p[6]), "ts": int(p[7]), "te": int(p[8]),
                       "mapq": int(p[11]), "cigar": cg}


def chained_alignment(paf_path, min_mapq=20):
    """One colinear chain of primary alignments, merged into a single alignment.

    minimap2 often splits a 1 Mb region into several pieces at structural
    variants or repeats. Pieces on the same contig and strand that advance on
    both sequences are joined; the stretch between two pieces becomes a
    deletion (reference gap) plus an insertion (haplotype gap), so every
    reference position still projects to one haplotype position.
    """
    groups = {}
    for a in _read_paf(paf_path, min_mapq):
        groups.setdefault((a["qname"], a["strand"]), []).append(a)
    if not groups:
        raise SystemExit(f"FATAL: no primary alignment with MAPQ >= {min_mapq} in {paf_path}")
    best = None
    for (_q, strand), alns in groups.items():
        chain = []
        for a in sorted(alns, key=lambda x: (x["ts"], -(x["te"] - x["ts"]))):
            if chain:
                last = chain[-1]
                ahead = a["qs"] >= last["qe"] if strand == "+" else a["qe"] <= last["qs"]
                if a["ts"] < last["te"] or not ahead:
                    continue
            chain.append(a)
        span = sum(a["te"] - a["ts"] for a in chain)
        if best is None or span > best[0]:
            best = (span, chain)
    chain = best[1]
    first, last = chain[0], chain[-1]
    cigar = [first["cigar"]]
    for prev, cur in zip(chain, chain[1:]):
        tgap = cur["ts"] - prev["te"]
        qgap = cur["qs"] - prev["qe"] if prev["strand"] == "+" else prev["qs"] - cur["qe"]
        cigar += [f"{tgap}D" if tgap else "", f"{qgap}I" if qgap else "", cur["cigar"]]
    merged = dict(first, te=last["te"], cigar="".join(cigar), pieces=len(chain),
                  aligned_bp=best[0], mapq=min(a["mapq"] for a in chain))
    if first["strand"] == "+":
        merged["qe"] = last["qe"]
    else:
        merged["qs"] = last["qs"]
    return merged


def projector(aln):
    """f(ref_offset) -> haplotype coordinate in the reference-oriented frame.

    ref_offset is relative to the region start. Positions before or after the
    alignment clamp to its ends. Inside an insertion the cut goes before it, the
    same rule for both neighbouring chunks, so their intervals still meet.
    """
    q0 = aln["qs"] if aln["strand"] == "+" else aln["qlen"] - aln["qe"]
    tpts, qpts = [aln["ts"]], [q0]
    t, q = aln["ts"], q0
    for n, op in CIGAR.findall(aln["cigar"]):
        n = int(n)
        if op in "M=X":
            t, q = t + n, q + n
            tpts.append(t)
            qpts.append(q)
        elif op in "DN":
            tpts.append(t)
            qpts.append(q)
            t += n
            tpts.append(t)
            qpts.append(q)
        elif op == "I":
            q += n
            tpts.append(t)
            qpts.append(q)

    def f(r):
        r = min(max(r, aln["ts"]), aln["te"])
        lo, hi = 0, len(tpts) - 1
        while lo < hi:  # first index with tpts[i] >= r
            mid = (lo + hi) // 2
            if tpts[mid] < r:
                lo = mid + 1
            else:
                hi = mid
        i = lo
        if tpts[i] == r:
            return qpts[i]
        # inside a match run the haplotype advances with the reference;
        # inside a deletion it stays put
        return qpts[i - 1] + min(r - tpts[i - 1], qpts[i] - qpts[i - 1])
    return f


def fetch(assembly, contig, start, end, revcomp, qlen):
    """Sequence of [start, end) in the reference-oriented frame."""
    if revcomp:
        start, end = qlen - end, qlen - start
    args = ["samtools", "faidx", assembly, f"{contig}:{start + 1}-{end}"]
    if revcomp:
        args.insert(2, "-i")
    out = subprocess.run(args, check=True, capture_output=True, text=True).stdout
    return "".join(out.splitlines()[1:])


def read_region(path):
    with open(path) as fh:
        return "".join(l.strip() for l in fh if not l.startswith(">"))


def write_record(fh, name, seq):
    fh.write(f">{name}\n")
    for i in range(0, len(seq), 80):
        fh.write(seq[i:i + 80] + "\n")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("manifest")
    ap.add_argument("region_fa")
    ap.add_argument("outdir")
    ap.add_argument("--region-start", type=int, required=True)
    ap.add_argument("--ref-name", required=True, help="PanSN prefix, e.g. GRCh38#0#chr21")
    ap.add_argument("--hap", action="append", default=[], help="SAMPLE#HAP=ALIGN.paf=ASSEMBLY")
    a = ap.parse_args(argv)

    ref = read_region(a.region_fa)
    with open(a.manifest) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    haps = []
    for spec in a.hap:
        prefix, paf, asm = spec.split("=")
        aln = chained_alignment(paf)
        aln["cov"] = aln["aligned_bp"] / aln["tlen"]
        haps.append((prefix, asm, aln, projector(aln)))

    proj = []
    for row in rows:
        cid = row["chunk_id"]
        rs, re_ = int(row["reference_start"]), int(row["reference_end"])
        lo, hi = rs - a.region_start, re_ - a.region_start
        with open(f"{a.outdir}/{cid}.fa", "w") as out:
            write_record(out, f"{a.ref_name}:{rs}-{re_}", ref[lo:hi])
            for prefix, asm, aln, f in haps:
                qs, qe = f(lo), f(hi)
                rc = aln["strand"] == "-"
                contig = aln["qname"] + ("_rc" if rc else "")
                seq = fetch(asm, aln["qname"], qs, qe, rc, aln["qlen"])
                # HPRC release 2 contigs are already PanSN-named (SAMPLE#HAP#CONTIG)
                name = contig if contig.count("#") >= 2 else f"{prefix}#{contig}"
                write_record(out, f"{name}:{qs}-{qe}", seq)
                proj.append({"chunk_id": cid, "haplotype": prefix, "contig": contig,
                             "ref_start": rs, "ref_end": re_, "hap_start": qs, "hap_end": qe,
                             "strand": aln["strand"], "mapq": aln["mapq"], "pieces": aln["pieces"],
                             "region_coverage": round(aln["cov"], 4)})
    with open(f"{a.outdir}/projection.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(proj[0]), delimiter="\t")
        w.writeheader()
        w.writerows(proj)
    print(f"{len(rows)} chunks x {len(haps) + 1} sequences -> {a.outdir}")


if __name__ == "__main__":
    main()
