"""Chunk status tracking and progress summary."""
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from pipeline.parallel.chunk_manifest import load_manifest, mark_built, progress, update_chunk_status


@pytest.fixture
def manifest(tmp_path):
    p = tmp_path / "chunk_manifest.tsv"
    p.write_text("chunk_id\treference_start\treference_end\tstatus\n"
                 "chunk_000\t0\t100\tpending\nchunk_001\t100\t200\tpending\nchunk_002\t200\t300\tpending\n")
    return str(p)


def test_set_status_and_progress(manifest):
    assert update_chunk_status("chunk_001", "failed", manifest)
    assert update_chunk_status("chunk_002", "done", manifest)
    assert not update_chunk_status("chunk_999", "done", manifest)
    with pytest.raises(ValueError):
        update_chunk_status("chunk_000", "finished", manifest)
    chunks = load_manifest(manifest)
    assert [c["status"] for c in chunks] == ["pending", "failed", "done"]
    assert chunks[0]["reference_end"] == "100"  # other columns survive the rewrite
    assert progress(chunks) == "1/3 done (pending 1, running 0, failed 1)"


def test_mark_built_uses_final_gfa(manifest, tmp_path):
    (tmp_path / "chunk_001").mkdir()
    (tmp_path / "chunk_001" / "final.gfa").write_text("H\tVN:Z:1.0\n")
    assert mark_built(manifest, str(tmp_path)) == 1
    assert [c["status"] for c in load_manifest(manifest)] == ["pending", "done", "pending"]
