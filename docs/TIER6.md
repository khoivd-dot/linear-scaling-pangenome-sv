# Tier 6: haplotype count against chunking, and cost per PGGB stage

Tier 4 ([TIER4.md](TIER4.md)) left three questions. (Tier 5, on paralog-aware cuts, is a separate experiment.)

1. **Is the accuracy loss at more haplotypes caused by chunking?** On the chr21
   arm, the same three samples scored F1 0.976, 0.972 and 0.969 at 6, 20 and
   40 haplotypes. The monolithic arm did not finish at 20 or 40 haplotypes, so
   chunking and haplotype count could not be told apart.
2. **Why does cost per haplotype rise at 40?** Chunked builds took 56, 48 and
   73 CPU-seconds per Mb per haplotype at 6, 20 and 40 haplotypes. PGGB's
   all-to-all alignment is the usual suspect, but nothing measured the stages
   one by one.
3. **How far can the truth set be trusted?** The truth calls for HG00673 and
   HG00733 come from our own dipcall-style caller on HPRC's alignments. The
   project's stand-alone dipcall wrapper (`pipeline/linear/run_dipcall.sh`)
   has never been run against them.

<!-- RESULTS -->

## E8: chunking against haplotype count, with time per stage (running)

The region is chr21:20-25 Mb. It is small enough that a monolithic build
should fit on a 16 GB runner at 40 haplotypes; Tier 3 measured 2.2 GB for
1 Mb at 40 haplotypes and 2.8 GB for 10 Mb at 6. It lies inside the Tier 4
arm, and uses the same nested haplotype sets: N = 6, 20 and 40 from
[`haps_arm.tsv`](../pipeline/benchmark/tier4/haps_arm.tsv), with HG00673,
HG00733 and HG002 in every set. Each N gets a monolithic graph and five
zero-overlap 1 Mb windows, scored against the same truth as Tier 4 (Tier 3
calls for HG00673 and HG00733, GIAB v5.0q for HG002).

Each build now keeps PGGB's own log. `pggb_stages.py` reads GNU time's report
for every tool (wfmash, seqwish, smoothxg, gfaffix, odgi) from it, so the
summary JSON holds wall time, CPU and peak memory per stage for every chunk.

What would settle the questions:

- If monolithic and chunked graphs score the same at each N, but both fall
  as N grows, the Tier 4 drop is a property of larger graphs (or of calling
  variants from them), not of chunking.
- If the chunked graph falls further than the monolithic one, chunking
  costs accuracy at high N, and the size of the gap says how much.
- The stage table should show which stage's share of CPU grows from 6 to 40
  haplotypes.

Run with stage `tier6-mid` and `n` = 6, 20 and 40 in `experiment.yml`. The
summaries land in `docs/tier4_runs/tier4_mid_n*_p90_rep1.json`.

## Proposed next

None of these has started. Each runs on free GitHub runners unless noted.

| ID | Question | Design | Cost |
|----|----------|--------|------|
| E9 | Do scale and segmental duplications interact? | The chr22 long arm, which holds 22q11.2, at 6 haplotypes, with HG002 scored against GIAB. Needs a scan for assemblies that span the LCR22s on one contig. | Moderate; the scan is the main work |
| E10 | Does the stand-alone dipcall agree with the Tier 3 truth calls? | Run `run_dipcall.sh` on HG00673 and HG00733 in chr21:20-21 Mb and 22q11.2, then compare with truvari. | Small |
| E11 | What chunk size is best at 40 haplotypes? | chr21:20-25 Mb with 500 kb, 1 Mb and 2.5 Mb windows at N = 40. | Small once E8 is in |
| E12 | Monolithic arm at 20 and 40 haplotypes | The Tier 4 arm on a 64 GB machine, to compare with the chunked graphs directly. | Paid compute; needs approval |
