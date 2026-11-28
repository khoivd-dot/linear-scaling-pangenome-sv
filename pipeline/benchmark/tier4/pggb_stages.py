#!/usr/bin/env python3
"""Time and memory per PGGB stage, from the log PGGB writes next to its graph.

PGGB runs each tool under GNU time and appends the report to its log. Both the short
form ("CMD", then "Us user Ss system P cpu Es total MKb max memory") and the verbose -v
form are read. Times of repeated tools (odgi, for instance) are summed; memory is the peak.

    pggb_stages.py PGGB.log [...] > stages.json
"""
import json
import re
import sys

TOOLS = ("wfmash", "seqwish", "smoothxg", "gfaffix", "odgi", "vg", "multiqc")
SHORT = re.compile(r"([\d.]+)s user\s+([\d.]+)s system.*?([\d.]+)s total\s+(\d+)\s*Kb max memory")


def tool_of(cmd):
    for word in cmd.replace('"', " ").split():
        name = word.rsplit("/", 1)[-1]
        if name in TOOLS:
            return name
    return None


def clock(s):
    parts = [float(x) for x in s.strip().split(":")]
    return sum(p * 60 ** i for i, p in enumerate(reversed(parts)))


def parse(text):
    out, last, verbose = {}, None, None

    def add(tool, wall, cpu, kb):
        if not tool:
            return
        s = out.setdefault(tool, {"wall_s": 0.0, "cpu_s": 0.0, "peak_rss_mb": 0.0, "runs": 0})
        s["wall_s"] = round(s["wall_s"] + wall, 1)
        s["cpu_s"] = round(s["cpu_s"] + cpu, 1)
        s["peak_rss_mb"] = round(max(s["peak_rss_mb"], kb / 1024), 1)
        s["runs"] += 1

    for line in text.splitlines():
        m = SHORT.search(line)
        if m:
            add(last, float(m.group(3)), float(m.group(1)) + float(m.group(2)), int(m.group(4)))
            continue
        t = line.strip()
        if t.startswith("Command being timed:"):
            verbose = {"tool": tool_of(t.split(":", 1)[1]), "user": 0.0, "sys": 0.0, "wall": 0.0}
        elif verbose is not None and t.startswith("User time (seconds):"):
            verbose["user"] = float(t.rsplit(":", 1)[1])
        elif verbose is not None and t.startswith("System time (seconds):"):
            verbose["sys"] = float(t.rsplit(":", 1)[1])
        elif verbose is not None and t.startswith("Elapsed (wall clock) time"):
            verbose["wall"] = clock(t.split("):", 1)[1])
        elif verbose is not None and t.startswith("Maximum resident set size (kbytes):"):
            add(verbose["tool"], verbose["wall"], verbose["user"] + verbose["sys"], int(t.rsplit(":", 1)[1]))
            verbose = None
        elif tool_of(t) and not t.startswith("["):
            last = tool_of(t)
    return out


def main(argv=None):
    paths = (argv if argv is not None else sys.argv[1:])
    text = "\n".join(open(p, errors="replace").read() for p in paths)
    json.dump(parse(text), sys.stdout, indent=1)
    print()


if __name__ == "__main__":
    main()
