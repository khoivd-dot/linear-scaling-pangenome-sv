"""web/public/data/results.json must match what scripts/build_site_results.py makes from docs/."""
import json
import os
import subprocess
import sys

ROOT = os.path.join(os.path.dirname(__file__), "..")
SCRIPT = os.path.join(ROOT, "scripts", "build_site_results.py")
RESULTS = os.path.join(ROOT, "web", "public", "data", "results.json")


def test_results_json_is_up_to_date():
    run = subprocess.run([sys.executable, SCRIPT, "--check"], capture_output=True, text=True)
    assert run.returncode == 0, run.stderr or run.stdout


def test_headline_numbers_match_the_write_ups():
    # Expected values copied from docs/TIER4.md (E6 tables) and docs/TIER3.md (E1), not recomputed.
    d = json.load(open(RESULTS))
    n6, n20, n40 = d["tier4"]["scaling"]
    hot = n6["hotspots1m"]
    assert (hot["builds"], hot["f1_all"], hot["cpu_vs_mono"], hot["longest_min"], hot["peak_gb"]) == (29, 0.976, 1.05, 3.4, 0.97)
    assert (n6["mono"]["f1_all"], n6["mono"]["longest_min"], n6["mono"]["peak_gb"]) == (0.976, 48.8, 12.1)
    assert [r["hotspots1m"]["cpu_s_per_mb_hap"] for r in (n6, n20, n40)] == [55.7, 48.3, 72.9]
    assert n20["mono"]["status"] == "killed" and n40["mono"]["status"] == "lost"
    assert (n20["hotspots1m"]["peak_gb"], n40["hotspots1m"]["peak_gb"]) == (1.5, 2.96)
    r22 = next(x for x in d["tier3"]["regions"] if x["key"] == "r22")
    assert (r22["mono"]["f1"][0], r22["chunked"]["f1"][0]) == (0.88, 0.962)
    assert (r22["mono"]["revisit_pct"][0], r22["chunked"]["revisit_pct"][0]) == (9.0, 2.4)
