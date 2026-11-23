"""Reference -> haplotype projection used to cut every haplotype at the same places."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from pipeline.prepare.project_chunks import projector


def _aln(cigar, strand="+", qs=100, qe=None, qlen=1000, ts=0, te=None):
    return {"cigar": cigar, "strand": strand, "qs": qs, "qe": qe, "qlen": qlen, "ts": ts, "te": te}


def test_plus_strand_match_insertion_deletion():
    # 10= then 5 bp insertion, 10=, 4 bp deletion, 10=  (target 34 bp, query 35 bp)
    f = projector(_aln("10=5I10=4D10=", qs=100, qe=135, te=34))
    assert f(0) == 100
    assert f(10) == 100 + 10          # cut goes before the insertion
    assert f(11) == 100 + 10 + 5 + 1
    assert f(20) == 125               # end of second match run
    assert f(22) == 125               # inside the deletion the haplotype stays put
    assert f(24) == 125
    assert f(34) == 135


def test_minus_strand_counts_on_reverse_complement():
    f = projector(_aln("20=", strand="-", qs=300, qe=320, qlen=1000, te=20))
    assert f(0) == 1000 - 320
    assert f(20) == 1000 - 300


def test_clamps_outside_alignment():
    f = projector(_aln("10=", qs=50, qe=60, ts=5, te=15))
    assert f(0) == 50 and f(100) == 60


def test_adjacent_chunks_meet():
    f = projector(_aln("7=3I3=2D8=", qs=0, qe=21, te=20))
    cuts = [0, 7, 9, 12, 20]
    spans = [(f(a), f(b)) for a, b in zip(cuts, cuts[1:])]
    assert all(x[1] == y[0] for x, y in zip(spans, spans[1:]))
    assert spans[0][0] == 0 and spans[-1][1] == 21


def test_chain_joins_split_alignment(tmp_path):
    from pipeline.prepare.project_chunks import chained_alignment
    paf = tmp_path / "x.paf"
    rows = [  # two pieces with a 5 bp reference gap and a 7 bp haplotype gap
        ["h", 1000, 100, 120, "+", "r", 50, 0, 20, 20, 20, 60, "tp:A:P", "cg:Z:20="],
        ["h", 1000, 127, 152, "+", "r", 50, 25, 50, 25, 25, 60, "tp:A:P", "cg:Z:25="],
        ["h", 1000, 900, 910, "+", "r", 50, 30, 40, 10, 10, 60, "tp:A:P", "cg:Z:10="],  # out of order
    ]
    paf.write_text("".join("\t".join(map(str, r)) + "\n" for r in rows))
    a = chained_alignment(str(paf))
    assert a["pieces"] == 2 and a["cigar"] == "20=5D7I25="
    f = projector(a)
    assert f(0) == 100 and f(20) == 120 and f(25) == 120 and f(26) == 128 and f(50) == 152
