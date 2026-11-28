"""
extract_regions.py — Extract orthologous regions from assemblies.

Uses mapping from map_chromosome.py to extract correct intervals.
Never naively assumes GRCh38 coordinates apply to de novo assemblies.
prepare_sequences.py builds the full multi-haplotype FASTA (samtools only);
this module extracts one interval from a plain or bgzipped FASTA with
`samtools faidx` when it is on PATH, otherwise with pysam.

Usage (0-based half-open, as in PAF and sequence_mapping.tsv):
    python3 pipeline/prepare/extract_regions.py ASSEMBLY.fa[.gz] CONTIG START END [--strand -]
"""
import argparse
import os
import shutil
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from pipeline.prepare.faidx_utils import _revcomp
from pipeline.prepare.project_chunks import fetch, write_record


def extract_region(fasta: str, contig: str, start: int, end: int, strand: str = "+") -> str:
    """Sequence of [start, end) on contig, reverse-complemented when strand is '-'."""
    if shutil.which("samtools"):
        seq = fetch(fasta, contig, start, end, False, None)
    else:
        try:
            import pysam
        except ImportError:
            raise RuntimeError("extract_region needs samtools on PATH or the pysam module")
        with pysam.FastaFile(fasta) as fa:  # builds .fai (and .gzi for bgzip) if missing
            seq = fa.fetch(contig, start, end)
    return _revcomp(seq) if strand == "-" else seq


def main(argv=None):
    ap = argparse.ArgumentParser(description="Extract one 0-based half-open interval")
    ap.add_argument("fasta")
    ap.add_argument("contig")
    ap.add_argument("start", type=int)
    ap.add_argument("end", type=int)
    ap.add_argument("--strand", choices="+-", default="+")
    a = ap.parse_args(argv)
    write_record(sys.stdout, f"{a.contig}:{a.start}-{a.end}", extract_region(
        a.fasta, a.contig, a.start, a.end, a.strand))


if __name__ == "__main__":
    main()
