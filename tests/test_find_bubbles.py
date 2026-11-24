"""Bubble detection and nesting on a tiny GFA."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from pipeline.merge.gfa import GfaGraph
from pipeline.export.find_bubbles import find_bubbles

# outer bubble 1 -> {2, 7} -> 10; inner bubble 2 -> {3, 4} -> 5 sits on the 2 allele
NESTED_GFA = "H\tVN:Z:1.0\n" + "".join(f"S\t{n}\tA\n" for n in [1, 2, 3, 4, 5, 7, 10]) + "".join(
    f"L\t{a}\t+\t{b}\t+\t0M\n" for a, b in
    [(1, 2), (1, 7), (2, 3), (2, 4), (3, 5), (4, 5), (5, 10), (7, 10)])


def test_nested_bubble_gets_parent_and_level():
    bubbles = {(b["source"], b["sink"]): b for b in find_bubbles(GfaGraph.parse(NESTED_GFA))}
    assert set(bubbles) == {("1", "10"), ("2", "5")}  # sink is the nearest join, not min("10", "5")
    outer, inner = bubbles[("1", "10")], bubbles[("2", "5")]
    assert outer["parent"] is None and outer["nesting_level"] == 0
    assert inner["nesting_level"] == 1
    assert find_bubbles(GfaGraph.parse(NESTED_GFA))[inner["parent"]]["source"] == "1"


def test_simple_bubble_is_top_level():
    gfa = "H\tVN:Z:1.0\n" + "".join(f"S\t{n}\tA\n" for n in "abcd") + "".join(
        f"L\t{x}\t+\t{y}\t+\t0M\n" for x, y in ["ab", "ac", "bd", "cd"])
    [b] = find_bubbles(GfaGraph.parse(gfa))
    assert (b["source"], b["sink"], b["parent"], b["nesting_level"]) == ("a", "d", None, 0)
