"""prepare_block.py: haploblocks.org FASTA -> PanSN PGGB input."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from pipeline.haploblocks.prepare_block import main, pansn

FASTA = """>HG00096_chr22_region_100-110_hap0
>HG00096_chr22_region_100-110_hap0
ACGTACGTAC
>HG00096_chr22_region_100-110_hap1
ACGTACGTAC
>HG00097_chr22_region_100-110_hap0
ACGTTCGTAC
>HG00097_chr22_region_100-110_hap1
ACGTACGTAC
"""


def _run(tmp_path, *extra):
    src = tmp_path / "in.fasta"
    src.write_text(FASTA)
    out = tmp_path / "out.fa"
    main([str(src), str(out), *extra])
    return out, [l[1:].strip() for l in out.read_text().splitlines() if l.startswith(">")]


def test_pansn_name():
    assert pansn("HG00096_chr22_region_100-110_hap1") == "HG00096#2#chr22:100-110"


def test_drops_empty_and_renames(tmp_path):
    _, names = _run(tmp_path)
    assert names == ["HG00096#1#chr22:100-110", "HG00096#2#chr22:100-110",
                     "HG00097#1#chr22:100-110", "HG00097#2#chr22:100-110"]


def test_dedupe_maps_every_haplotype(tmp_path):
    out, names = _run(tmp_path, "--dedupe")
    assert names == ["HG00096#1#chr22:100-110", "HG00097#1#chr22:100-110"]
    rows = (tmp_path / "out.fa.dedupe.tsv").read_text().splitlines()[1:]
    assert len(rows) == 4
    assert "HG00097#2#chr22:100-110\tHG00096#1#chr22:100-110" in rows


def test_subset_is_seeded(tmp_path):
    _, a = _run(tmp_path, "--n-haps", "2", "--seed", "7")
    _, b = _run(tmp_path, "--n-haps", "2", "--seed", "7")
    assert a == b and len(a) == 2
