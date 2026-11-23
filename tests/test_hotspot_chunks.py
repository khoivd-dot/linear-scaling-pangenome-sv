"""Chunks cut at recombination hotspots: they tile the region and touch."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from pipeline.parallel.hotspot_chunks import create_hotspot_chunks, read_cuts

NAMES = ["chr21_19978426-20023729", "chr21_20023729-20271080", "chr21_20271080-20347668",
         "chr21_20347668-20373919", "chr21_20999271-21030086", "chr22_1-2"]


def _spans(chunks):
    return [(c["reference_start"], c["reference_end"]) for c in chunks]


def test_cuts_inside_region_only():
    assert read_cuts(NAMES, "chr21", 20000000, 21000000) == [
        20023729, 20271080, 20347668, 20373919, 20999271]


def test_chunks_tile_region_and_touch():
    cuts = read_cuts(NAMES, "chr21", 20000000, 21000000)
    spans = _spans(create_hotspot_chunks("chr21", 20000000, 21000000, cuts, 50000))
    assert spans[0][0] == 20000000 and spans[-1][1] == 21000000
    assert all(a[1] == b[0] for a, b in zip(spans, spans[1:]))
    assert all(hi - lo >= 50000 for lo, hi in spans)


def test_short_tail_joins_left_neighbour():
    spans = _spans(create_hotspot_chunks("c", 0, 1000, [400, 990], 100))
    assert spans == [(0, 400), (400, 1000)]


def test_no_min_size_keeps_every_block():
    spans = _spans(create_hotspot_chunks("c", 0, 1000, [400, 990], 0))
    assert spans == [(0, 400), (400, 990), (990, 1000)]
