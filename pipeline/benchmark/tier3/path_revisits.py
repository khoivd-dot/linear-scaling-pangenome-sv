#!/usr/bin/env python3
"""How much of each haplotype path re-enters nodes it already used (collapsed repeats/paralogs).

A path that visits the same node twice spells two copies of a sequence through one node: the graph
has merged two loci of that haplotype. Reports, per path, the bp in repeated visits and their share.

    path_revisits.py GFA [GFA ...]  -> TSV: gfa, path, path_bp, revisit_bp, share
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))
from pipeline.merge.gfa import GfaGraph  # noqa: E402
from pipeline.merge.merge_graphs import walks_as_paths  # noqa: E402

print("gfa\tpath\tpath_bp\trevisit_bp\tshare")
for f in sys.argv[1:]:
    g = walks_as_paths(GfaGraph.parse_file(f))
    seglen = {k: len(s.sequence) for k, s in g.segments.items()}
    for name, p in g.paths.items():
        seen, tot, rep = set(), 0, 0
        for step in p.segment_names:
            n = step.rstrip("+-")
            L = seglen.get(n, 0); tot += L
            if n in seen:
                rep += L
            seen.add(n)
        print(f"{f}\t{name}\t{tot}\t{rep}\t{rep / tot if tot else 0:.4f}")
