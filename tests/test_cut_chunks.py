import random

from pipeline.parallel.cut_chunks import random_cuts
from pipeline.parallel.hotspot_chunks import create_hotspot_chunks
from pipeline.benchmark.seam_metrics import split_multiallelic


def test_random_cuts_respect_min_size_and_touch():
    rng = random.Random(7)
    for _ in range(200):
        cuts = random_cuts(rng, 0, 1_000_000, 3, 200_000)
        rows = create_hotspot_chunks("chr21", 0, 1_000_000, cuts, min_size=0)
        assert len(rows) == 3
        assert all(r["reference_end"] - r["reference_start"] >= 200_000 for r in rows)
        assert all(a["reference_end"] == b["reference_start"] for a, b in zip(rows, rows[1:]))


def test_split_multiallelic(tmp_path):
    src = tmp_path / "in.vcf"
    src.write_text("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n"
                   "chr21\t10\t>1>4\tG\tGCA,GC\t60\t.\tAC=1,1\tGT\t1|0\t0|2\n")
    out = tmp_path / "out.vcf"
    split_multiallelic(str(src), str(out))
    body = [l.split("\t") for l in out.read_text().splitlines() if not l.startswith("#")]
    assert [b[4] for b in body] == ["GCA", "GC"]
    assert [b[9:] for b in body] == [["1|0", "0|0"], ["0|0", "0|1"]]
