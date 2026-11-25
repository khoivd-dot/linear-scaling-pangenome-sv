"""
find_bubbles.py — First-pass bubble detection in GFA.

Each bubble records its `parent` (index of the smallest bubble whose interior
holds this bubble's source and sink, or None) and its `nesting_level` (0 = top).
Link orientation is ignored.
TODO (Ali): Integrate with web display.
"""
import argparse
import json
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from pipeline.merge.gfa import GfaGraph


def find_bubbles(graph: GfaGraph, max_depth: int = 5) -> list:
    out = {}
    for link in graph.links:
        out.setdefault(link.from_node, []).append(link.to_node)

    def distances(start: str, depth: int) -> dict:
        dist, frontier = {start: 0}, [start]
        for d in range(1, depth + 1):
            frontier = [n for f in frontier for n in out.get(f, []) if n not in dist]
            for n in frontier:
                dist.setdefault(n, d)
        return dist

    back = {}
    for link in graph.links:
        back.setdefault(link.to_node, []).append(link.from_node)

    def spread(start: str, stop: str, adj: dict) -> set:
        seen, todo = {start}, [start]
        while todo:
            for n in adj.get(todo.pop(), []):
                if n not in seen and n != stop:
                    seen.add(n)
                    todo.append(n)
        return seen

    bubbles = []
    for node, targets in out.items():
        if len(targets) >= 2:
            for i, t1 in enumerate(targets):
                for t2 in targets[i + 1:]:
                    d1, d2 = distances(t1, max_depth), distances(t2, max_depth)
                    common = d1.keys() & d2.keys()
                    if common:
                        # nearest shared node, not min() of string IDs ("10" < "5")
                        sink = min(common, key=lambda n: (max(d1[n], d2[n]), d1[n] + d2[n], n))
                        bubbles.append({
                            "source": node,
                            "sink": sink,
                            "alleles": [t1, t2],
                            "depth": max_depth,
                        })

    interiors = [(spread(b["source"], b["sink"], out) & spread(b["sink"], b["source"], back))
                 - {b["source"], b["sink"]} for b in bubbles]
    order = sorted(range(len(bubbles)), key=lambda i: -len(interiors[i]))
    for i in order:  # outer (larger) bubbles get their level first
        b = bubbles[i]
        holders = [j for j in range(len(bubbles))
                   if b["source"] in interiors[j] and b["sink"] in interiors[j]]
        parent = min(holders, key=lambda j: len(interiors[j]), default=None)
        b["parent"] = parent
        b["nesting_level"] = 0 if parent is None else bubbles[parent].get("nesting_level", 0) + 1
    return bubbles


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input")
    parser.add_argument("--output", "-o")
    args = parser.parse_args()
    graph = GfaGraph.parse_file(args.input)
    bubbles = find_bubbles(graph)
    output = args.output or args.input.rsplit(".", 1)[0] + "_bubbles.json"
    with open(output, "w") as f:
        json.dump(bubbles, f, indent=2)
    print(f"Found {len(bubbles)} bubbles -> {output}")


if __name__ == "__main__":
    main()