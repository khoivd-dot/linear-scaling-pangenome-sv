"""Tier 6: interval handling for the dipcall check."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "pipeline", "benchmark", "tier6"))


def test_intervals_clip_merge_and_intersect(tmp_path):
    from dipcall_check import intersect, intervals
    bed = tmp_path / "a.bed"
    bed.write_text("chr21\t50\t150\nchr21\t140\t300\nchr22\t0\t999\nchr21\t400\t900\n")
    a = intervals(str(bed), "chr21", 100, 800)
    assert a == [(100, 300), (400, 800)]
    assert intersect(a, [(0, 120), (250, 450), (790, 1000)]) == [(100, 120), (250, 300), (400, 450), (790, 800)]
