#!/usr/bin/env python3
"""Assembly-based truth calls per haplotype from the HPRC winnowmap GRCh38 alignments (dipcall-style).

For one region and one haplotype BAM (HPRC release 2, assembly vs GRCh38_no_alt), writes
  OUT.vcf  haploid calls (SNVs, indels, SVs) read off the alignment CIGAR against the reference
  OUT.bed  callable region: covered by exactly one alignment >= --min-cov-len, and that alignment
           is >= --min-call-len with MAPQ >= --min-mapq (paftools call / dipcall defaults: 10 kb, 50 kb, 5)

    truth_calls.py BAM_URL REF.fa CHROM START END OUT [--sample HG00673]
REF.fa uses Ensembl names ("21"); output uses CHROM ("chr21").
"""
import argparse
import pysam


def callable_intervals(alns, start, end, min_cov_len, min_call_len, min_mapq):
    events = []
    for a in alns:
        if a.query_alignment_length >= min_cov_len:
            events += [(max(a.reference_start, start), 1, a), (min(a.reference_end, end), -1, a)]
    events.sort(key=lambda e: (e[0], e[1]))
    out, active, prev = [], set(), None
    for pos, kind, a in events:
        if prev is not None and pos > prev and len(active) == 1:
            (only,) = active
            if only.query_alignment_length >= min_call_len and only.mapping_quality >= min_mapq:
                if out and out[-1][1] == prev and out[-1][2] is only:
                    out[-1][1] = pos
                else:
                    out.append([prev, pos, only])
        (active.add if kind == 1 else active.discard)(a)
        prev = pos
    return [(s, e, a) for s, e, a in out if e > s]


def calls(aln, ref, lo, hi):
    """(pos0, REF, ALT) for one alignment, restricted to reference [lo, hi)."""
    q = aln.query_sequence.upper()
    rpos, qpos = aln.reference_start, 0
    # query_sequence includes soft clips (op 4 advances qpos) and excludes hard clips (op 5 is skipped)
    for op, n in aln.cigartuples:
        if op in (0, 7, 8):
            for i in range(n):
                r = rpos + i
                if lo <= r < hi:
                    rb, qb = ref[r], q[qpos + i]
                    if rb != qb and rb != "N" and qb != "N":
                        yield r, rb, qb
            rpos += n; qpos += n
        elif op == 1:
            if lo <= rpos - 1 < hi and lo <= rpos < hi:
                ins = q[qpos:qpos + n]
                if "N" not in ins:
                    yield rpos - 1, ref[rpos - 1], ref[rpos - 1] + ins
            qpos += n
        elif op == 2:
            if lo <= rpos - 1 and rpos + n < hi:
                dl = ref[rpos - 1:rpos + n]
                if "N" not in dl:
                    yield rpos - 1, dl, dl[0]
            rpos += n
        elif op == 3:
            rpos += n
        elif op == 4:
            qpos += n


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("bam"); ap.add_argument("ref"); ap.add_argument("chrom")
    ap.add_argument("start", type=int); ap.add_argument("end", type=int); ap.add_argument("out")
    ap.add_argument("--sample", default="TRUTH")
    ap.add_argument("--min-cov-len", type=int, default=10000)
    ap.add_argument("--min-call-len", type=int, default=50000)
    ap.add_argument("--min-mapq", type=int, default=5)
    a = ap.parse_args(argv)
    fa = pysam.FastaFile(a.ref)
    rname = a.chrom if a.chrom in fa.references else a.chrom.removeprefix("chr")
    ref = fa.fetch(rname).upper()
    bam = pysam.AlignmentFile(a.bam)
    alns = [x for x in bam.fetch(a.chrom, a.start, a.end) if not x.is_unmapped and not x.is_secondary]
    regions = callable_intervals(alns, a.start, a.end, a.min_cov_len, a.min_call_len, a.min_mapq)
    recs = set()
    for s, e, aln in regions:
        recs.update(calls(aln, ref, s, e))
    with open(a.out + ".bed", "w") as fh:
        for s, e, _ in regions:
            fh.write(f"{a.chrom}\t{s}\t{e}\n")
    with open(a.out + ".vcf", "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write(f"##contig=<ID={a.chrom},length={len(ref)}>\n")
        fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        fh.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{a.sample}\n")
        for p, r, alt in sorted(recs):
            fh.write(f"{a.chrom}\t{p + 1}\t.\t{r}\t{alt}\t60\tPASS\t.\tGT\t1\n")
    cov = sum(e - s for s, e, _ in regions)
    print(f"{a.out}: {len(recs)} calls, callable {cov} bp of {a.end - a.start} ({cov / (a.end - a.start):.3f})")


if __name__ == "__main__":
    main()
