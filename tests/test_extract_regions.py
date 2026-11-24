"""Region extraction from plain and bgzipped FASTA (samtools or pysam)."""
import os
import shutil
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from pipeline.prepare import extract_regions
from pipeline.prepare.extract_regions import extract_region

try:
    import pysam
except ImportError:
    pysam = None

needs_backend = pytest.mark.skipif(not (shutil.which("samtools") or pysam),
                                   reason="needs samtools or pysam")
FASTA = ">ctgA\nAACCGGTTAC\nGTAC\n>ctgB\nTTTT\n"


@pytest.fixture(params=["plain", "bgzip"])
def fasta(request, tmp_path):
    p = tmp_path / "asm.fa"
    p.write_text(FASTA)
    if request.param == "plain":
        return str(p)
    if not pysam:
        pytest.skip("bgzip fixture needs pysam")
    pysam.tabix_compress(str(p), str(p) + ".gz")
    return str(p) + ".gz"


@needs_backend
def test_extracts_half_open_interval_across_lines(fasta):
    assert extract_region(fasta, "ctgA", 8, 12) == "ACGT"   # spans the line break
    assert extract_region(fasta, "ctgA", 0, 3, "-") == "GTT"  # revcomp of AAC
    assert extract_region(fasta, "ctgB", 0, 4) == "TTTT"


def test_clear_error_without_samtools_or_pysam(monkeypatch):
    monkeypatch.setattr(extract_regions.shutil, "which", lambda _: None)
    monkeypatch.setitem(sys.modules, "pysam", None)
    with pytest.raises(RuntimeError, match="samtools"):
        extract_region("missing.fa", "ctgA", 0, 1)
