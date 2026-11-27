"""Tier 4 kit: experiment sets, hotspot merging and coordinate-only projection."""
import csv
import json
import os
import subprocess
import sys

ROOT = os.path.join(os.path.dirname(__file__), "..")
T4 = os.path.join(ROOT, "pipeline", "benchmark", "tier4")
sys.path.insert(0, ROOT)
sys.path.insert(0, T4)


def _plan(set_, n):
    out = subprocess.run([sys.executable, os.path.join(T4, "sets.py"), "plan", set_, str(n)],
                         check=True, capture_output=True, text=True).stdout
    return json.loads(out)


def test_arm_plan_nested_haplotypes_and_units():
    p6, p40 = _plan("arm", 6), _plan("arm", 40)
    assert p6["samples"] == ["HG00673", "HG00733", "HG002"]
    assert p40["samples"][:3] == p6["samples"] and len(p40["samples"]) == 20
    # N = 6: mono + 29 hotspot chunks + 31 windows; N = 40: no windows
    assert len(p6["units"]) == 1 + 29 + 31 and len(p40["units"]) == 1 + 29
    assert p6["units"][0]["strategy"] == "mono"          # largest build queued first
    assert len({u["id"] for u in p6["units"]}) == len(p6["units"])


def test_hard_plan_regions():
    p = _plan("hard", 6)
    assert p["regions"] == ["r17", "r22"]
    assert {g["id"] for g in p["graphs"]} == {"r17.mono", "r17.windows0", "r22.mono", "r22.windows0"}


def test_manifests_tile_regions_without_overlap():
    from sets import REGIONS
    for name in os.listdir(os.path.join(T4, "manifests")):
        region = name.split(".")[0]
        rows = list(csv.DictReader(open(os.path.join(T4, "manifests", name)), delimiter="\t"))
        chrom, start, end = REGIONS[region]
        assert int(rows[0]["reference_start"]) == start and int(rows[-1]["reference_end"]) == end
        for a, b in zip(rows, rows[1:]):
            assert a["reference_end"] == b["reference_start"] and b["overlap_left"] == "0"


def test_hotspot_cuts_merge_to_min_size(tmp_path):
    from make_manifests import hotspot_cuts
    bs = tmp_path / "bs.tsv"
    rows = ["chr\tblock\tstart\tend"] + [f"chr21\tb\t{s}\t{s + 300}" for s in range(0, 3000, 300)]
    bs.write_text("\n".join(rows) + "\n")
    # boundaries every 300 bp; each chunk reaches 1000 bp; the 600 bp tail stays unless the minimum is higher
    assert hotspot_cuts(str(bs), "chr21", 0, 3000, 1000, 500) == [1200, 2400]
    assert hotspot_cuts(str(bs), "chr21", 0, 3000, 1000, 700) == [1200]


def test_project_coords_matches_project_chunks(tmp_path):
    from pipeline.prepare.project_chunks import chained_alignment, projector
    import project_coords
    paf = tmp_path / "x.paf"
    paf.write_text("\t".join(map(str, ["S#1#c", 1000, 100, 135, "+", "r", 34, 0, 34, 30, 35, 60,
                                       "tp:A:P", "cg:Z:10=5I10=4D10="])) + "\n")
    man = tmp_path / "m.tsv"
    man.write_text("chunk_id\treference_start\treference_end\n"
                   "chunk_0001\t500\t510\nchunk_0002\t510\t534\n")
    out = tmp_path / "c.tsv"
    project_coords.main([str(paf), "https://x/a.fa.gz", "S#1", str(out), "--region-start", "500",
                         "--manifest", f"w={man}"])
    rows = list(csv.DictReader(open(out), delimiter="\t"))
    f = projector(chained_alignment(str(paf)))
    assert [(int(r["hap_start"]), int(r["hap_end"])) for r in rows] == [(f(0), f(10)), (f(10), f(34))]
    assert rows[0]["hap_end"] == rows[1]["hap_start"] and rows[0]["qname"] == "S#1#c"
