"""
chunk_manifest.py — Utilities for chunk manifest management.

Tracks each chunk's `status` column (pending/running/done/failed) and prints
a progress summary.

Usage:
    python3 pipeline/parallel/chunk_manifest.py                      # list + summary
    python3 pipeline/parallel/chunk_manifest.py --set chunk_003 failed
    python3 pipeline/parallel/chunk_manifest.py --sync               # done where final.gfa exists
"""
import argparse
import csv
import os
from collections import Counter

MANIFEST = "work/chunks/chunk_manifest.tsv"
STATUSES = ("pending", "running", "done", "failed")


def load_manifest(path: str = MANIFEST) -> list:
    if not os.path.exists(path):
        return []
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def _write_manifest(chunks: list, path: str) -> None:
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=chunks[0].keys(), delimiter="\t")
        w.writeheader()
        w.writerows(chunks)


def update_chunk_status(chunk_id: str, status: str, path: str = MANIFEST) -> bool:
    """Set one chunk's status; returns False if the chunk is not in the manifest."""
    if status not in STATUSES:
        raise ValueError(f"status must be one of {STATUSES}, got {status!r}")
    chunks = load_manifest(path)
    for c in chunks:
        if c["chunk_id"] == chunk_id:
            c["status"] = status
            _write_manifest(chunks, path)
            return True
    return False


def mark_built(path: str = MANIFEST, chunks_dir: str = "work/chunks") -> int:
    """Mark every chunk whose <chunks_dir>/<chunk_id>/final.gfa exists as done."""
    chunks = load_manifest(path)
    n = 0
    for c in chunks:
        if c["status"] != "done" and os.path.exists(os.path.join(chunks_dir, c["chunk_id"], "final.gfa")):
            c["status"] = "done"
            n += 1
    if n:
        _write_manifest(chunks, path)
    return n


def progress(chunks: list) -> str:
    n = Counter(c["status"] for c in chunks)
    rest = ", ".join(f"{s} {n[s]}" for s in STATUSES if s != "done")
    return f"{n['done']}/{len(chunks)} done ({rest})"


def main(argv=None):
    ap = argparse.ArgumentParser(description="Chunk status tracking")
    ap.add_argument("--manifest", default=MANIFEST)
    ap.add_argument("--set", nargs=2, metavar=("CHUNK_ID", "STATUS"))
    ap.add_argument("--sync", action="store_true", help="mark chunks with final.gfa as done")
    a = ap.parse_args(argv)
    if a.set and not update_chunk_status(*a.set, path=a.manifest):
        ap.error(f"chunk {a.set[0]} not in {a.manifest}")
    if a.sync:
        mark_built(a.manifest, os.path.dirname(a.manifest) or ".")
    chunks = load_manifest(a.manifest)
    print(f"Chunks in manifest: {len(chunks)}")
    for c in chunks:
        print(f"  {c['chunk_id']}: {c['status']}")
    print(progress(chunks))


if __name__ == "__main__":
    main()
