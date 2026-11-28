#!/usr/bin/env python3
"""Project chunk boundaries onto one haplotype and write coordinates only (no sequence).

The same projection as pipeline/prepare/project_chunks.py (longest colinear chain of the
minimap2 alignment of the haplotype's contig to the reference region), but the output is a
table, so each GitHub build job can fetch just its own chunk from the remote assembly
(chunk_fasta.py) instead of every job downloading every contig.

    project_coords.py PAF URL PREFIX OUT.tsv --region-start 15000000 --manifest NAME=MANIFEST.tsv [...]
Appends rows: strategy, chunk_id, ref_start, ref_end, prefix, qname, strand, qlen, hap_start, hap_end,
url, mapq, pieces, region_coverage.
"""
import argparse
import csv
import os
import sys

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..")
sys.path.insert(0, ROOT)
from pipeline.prepare.project_chunks import chained_alignment, projector  # noqa: E402

FIELDS = ["strategy", "chunk_id", "ref_start", "ref_end", "prefix", "qname", "strand", "qlen",
          "hap_start", "hap_end", "url", "mapq", "pieces", "region_coverage"]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("paf"); ap.add_argument("url"); ap.add_argument("prefix"); ap.add_argument("out")
    ap.add_argument("--region-start", type=int, required=True)
    ap.add_argument("--manifest", action="append", required=True, help="STRATEGY=MANIFEST.tsv")
    a = ap.parse_args(argv)
    aln = chained_alignment(a.paf)
    f = projector(aln)
    new = not os.path.exists(a.out)
    with open(a.out, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t")
        if new:
            w.writeheader()
        for spec in a.manifest:
            name, manifest = spec.split("=")
            for row in csv.DictReader(open(manifest), delimiter="\t"):
                rs, re_ = int(row["reference_start"]), int(row["reference_end"])
                qs, qe = f(rs - a.region_start), f(re_ - a.region_start)
                w.writerow({"strategy": name, "chunk_id": row["chunk_id"], "ref_start": rs, "ref_end": re_,
                            "prefix": a.prefix, "qname": aln["qname"], "strand": aln["strand"],
                            "qlen": aln["qlen"], "hap_start": qs, "hap_end": qe, "url": a.url,
                            "mapq": aln["mapq"], "pieces": aln["pieces"],
                            "region_coverage": round(aln["aligned_bp"] / aln["tlen"], 4)})
    print(f"{a.prefix}: {aln['qname']} {aln['strand']} coverage {aln['aligned_bp'] / aln['tlen']:.4f}")


if __name__ == "__main__":
    main()
