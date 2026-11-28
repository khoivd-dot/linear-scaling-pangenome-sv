#!/usr/bin/env python3
"""Write one chunk's PGGB input from projected coordinates: the GRCh38 slice plus each haplotype's
slice fetched remotely (bgzipped HPRC assembly over HTTPS, via pysam/htslib).

Record names and orientation match pipeline/prepare/project_chunks.py, so the stitch sees the
same input either way: `GRCh38#0#chr21:start-end`, then `SAMPLE#HAP#CONTIG[_rc]:start-end`.

    chunk_fasta.py COORDS_DIR STRATEGY CHUNK_ID REF_CHROM.fa OUT.fa --region arm --chrom chr21 --prefixes HG00673#1 ...
COORDS_DIR holds the project_coords.py tables (coords.<REGION>.<SAMPLE>.tsv from prep_sample.sh; strategy and
chunk names repeat across regions, so only the --region tables are read); REF_CHROM.fa is a whole-chromosome
FASTA (indexed with pysam if needed) whose single record is the chromosome.
"""
import argparse
import csv
import glob
import os
import time

import pysam

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def write_record(fh, name, seq):
    fh.write(f">{name}\n")
    for i in range(0, len(seq), 80):
        fh.write(seq[i:i + 80] + "\n")


def fetch(url, contig, start, end, revcomp, qlen, tries=4):
    if revcomp:
        start, end = qlen - end, qlen - start
    for t in range(tries):
        try:
            seq = pysam.FastaFile(url).fetch(contig, start, end).upper()
            break
        except (OSError, ValueError):
            if t == tries - 1:
                raise
            time.sleep(10 * (t + 1))
    if len(seq) != end - start:
        raise SystemExit(f"FATAL: {contig}:{start}-{end} returned {len(seq)} bp")
    return seq.translate(COMP)[::-1] if revcomp else seq


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("coords_dir"); ap.add_argument("strategy"); ap.add_argument("chunk_id")
    ap.add_argument("ref_fa"); ap.add_argument("out")
    ap.add_argument("--region", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--prefixes", nargs="+", required=True, help="SAMPLE#HAP in output order")
    a = ap.parse_args(argv)
    rows = {}
    for path in sorted(glob.glob(os.path.join(a.coords_dir, "**", f"coords.{a.region}.*.tsv"), recursive=True)):
        for r in csv.DictReader(open(path), delimiter="\t"):
            if r["strategy"] == a.strategy and r["chunk_id"] == a.chunk_id:
                rows[r["prefix"]] = r
    missing = [p for p in a.prefixes if p not in rows]
    if missing:
        raise SystemExit(f"FATAL: no coordinates for {missing}")
    ref = pysam.FastaFile(a.ref_fa)
    rs, re_ = int(rows[a.prefixes[0]]["ref_start"]), int(rows[a.prefixes[0]]["ref_end"])
    if any((int(rows[p]["ref_start"]), int(rows[p]["ref_end"])) != (rs, re_) for p in a.prefixes):
        raise SystemExit("FATAL: haplotypes disagree on the chunk's reference interval")
    with open(a.out, "w") as out:
        write_record(out, f"GRCh38#0#{a.chrom}:{rs}-{re_}", ref.fetch(ref.references[0], rs, re_).upper())
        for p in a.prefixes:
            r = rows[p]
            rc = r["strand"] == "-"
            qs, qe = int(r["hap_start"]), int(r["hap_end"])
            seq = fetch(r["url"], r["qname"], qs, qe, rc, int(r["qlen"]))
            contig = r["qname"] + ("_rc" if rc else "")
            name = contig if contig.count("#") >= 2 else f"{p}#{contig}"
            write_record(out, f"{name}:{qs}-{qe}", seq)
    print(f"{a.strategy}/{a.chunk_id}: {len(a.prefixes) + 1} sequences -> {a.out}")


if __name__ == "__main__":
    main()
