# Tier 4: a chromosome arm, 40 haplotypes, and a curated truth set

Tiers 1 to 3 ([TIER1_CONTROLS.md](TIER1_CONTROLS.md), [TIER2.md](TIER2.md),
[TIER3.md](TIER3.md)) tested chunking on regions of 1 to 10 Mb with at most 40
haplotypes, and scored accuracy against truth calls made from HPRC's own
assembly-to-reference alignments. Tier 4 closes the two items the plan still
had open:

- **E6, chromosome scale.** 31 Mb of the chr21 long arm (chr21:15-46 Mb) at
  6, 20 and 40 haplotypes, chunked at haploblocks.org boundaries into 29
  chunks of about 1 Mb, against a monolithic build of the same input. This
  measures the cost per Mb per haplotype that the repository's name refers to.
- **E7, a curated truth set.** HG002 joins HG00673 and HG00733 in 22q11.2 and
  17q21.31, and its calls are scored against the GIAB HG002 Q100 v5.0q
  benchmark, at PGGB `-p 90` and `-p 98`, with two replicates.

All runs used GitHub-hosted runners (4 vCPU, 16 GB RAM, 14 GB SSD for public
repositories) through `.github/workflows/experiment.yml`, one runner per
chunk.

## Summary

- **On the chr21 long arm, chunking costs nothing in accuracy.** With 6
  haplotypes the monolithic graph of 31 Mb and both chunked graphs score the
  same against the truth: F1 0.976 for all variants, and 0.855 to 0.875 for
  SVs of 50 bp or more (1,700 truth SVs). All 28 or 30 seams pass, and every
  haplotype path in the chunked graphs is identical to the monolithic one.
  HG002 scores 0.980 against GIAB in all three graphs.
- **The monolithic arm does not fit on a standard runner past 6 haplotypes.**
  At 6 haplotypes it took 49 minutes and 12.1 GB. At 20 the kernel killed it
  (exit 137, out of memory) after 63 minutes on the 16 GB runner. At 40 the
  runner stopped responding after 5 hours inside PGGB and GitHub kept no log.
  The 29 hotspot chunks needed at most 8.2 minutes and 1.5 GB each at 20
  haplotypes, and 25 minutes and 3.0 GB at 40.
- **Cost per Mb per haplotype is near 50 CPU-seconds up to 20 haplotypes and
  rises after.** Chunked builds took 56 CPU-seconds per Mb per haplotype at 6
  haplotypes, 48 at 20 and 73 at 40 (53 for the monolithic build at 6). Total
  CPU grew 2.9 times from 6 to 20 haplotypes and 3.0 times from 20 to 40, as
  expected once PGGB's all-to-all alignment, which grows with the square of
  the number of sequences, takes over. At 6 haplotypes the chunked graphs used
  1.05 times (hotspots) and 1.08 times (fixed windows) the CPU of the
  monolithic build.
- **The same three samples score a little lower as haplotypes are added.**
  F1 for all variants is 0.976, 0.972 and 0.969 at 6, 20 and 40 haplotypes,
  and SV F1 0.855, 0.844 and 0.813. Without a monolithic graph at 20 or 40 we
  cannot say how much of this is chunking.
- **With HG002 added, the hard regions repeat Tier 3.** At PGGB `-p 90` the
  chunked 22q11.2 graph scores 0.964 to 0.965 against 0.898 monolithic, and the chunked
  17q21.31 graph 0.802 against 0.69 to 0.70. At `-p 98 -s 10000 -k 47` the
  monolithic graphs improve (0.94 and 0.85 once the vcfbub trap is avoided)
  but the chunked graphs still score higher (0.97 and 0.94). Both replicates
  agree to within 0.01.
- **Against the curated GIAB benchmark the gap is smaller.** In 22q11.2,
  HG002's F1 is 0.982 to 0.983 for the chunked graph at both settings. The
  monolithic graph reaches 0.943 and 0.976 for the two haplotypes at `-p 90`
  and 0.982 at `-p 98`. GIAB's benchmark regions leave out much of the LCR22
  duplications, where the two graphs differ most, so the curated truth can
  judge only part of the difference.
- **Adding HG002 made the monolithic 17q21.31 graph merge more.** Its GRCh38
  path revisits nodes over 17% of its length (11% in Tier 3 with 4
  haplotypes). The chunked graph's figure stays at 8%. On the chr21 arm
  every path's figure is 0.1% or less, so there is nothing to merge there.

## E6: chromosome scale

Cost of each graph of chr21:15-46 Mb. CPU is the sum over all builds; the
longest build and peak memory are those of the largest single job, which is
what a runner has to hold. "Per Mb per haplotype" divides CPU by 31 Mb and by
N (the GRCh38 path is not counted).

| N | Graph | Builds | CPU-h | vs monolithic | Longest build | Peak memory | CPU-s per Mb per haplotype |
|---|---|---|---|---|---|---|---|
| 6 | monolithic | 1 | 2.75 | 1.00 | 48.8 min | 12.1 GB | 53.2 |
| 6 | hotspots1m | 29 | 2.88 | 1.05 | 3.4 min | 0.97 GB | 55.7 |
| 6 | windows1m | 31 | 2.96 | 1.08 | 3.0 min | 0.97 GB | 57.4 |
| 20 | monolithic | 1 | killed | | > 63 min | > 16 GB | |
| 20 | hotspots1m | 29 | 8.31 | | 8.2 min | 1.50 GB | 48.3 |
| 40 | monolithic | 1 | runner lost | | > 5 h | | |
| 40 | hotspots1m | 29 | 25.11 | | 24.6 min | 2.96 GB | 72.9 |

Stitching the chunks took 6 s at 6 haplotypes, 16 s at 20 and 34 s at 40. The chunk
builds ran up to 20 at a time, the free plan's limit on concurrent jobs.

Accuracy against the truth, `vcfbub -a 100000` (`-a 10000` gives the same
figures to within 0.002):

| N | Graph | F1 all | F1 SV | HG002 F1 (GIAB), hap 1 / hap 2 | Seams passed | Paths identical to monolithic |
|---|---|---|---|---|---|---|
| 6 | monolithic | 0.976 | 0.862 | 0.980 / 0.980 | | |
| 6 | hotspots1m | 0.976 | 0.855 | 0.980 / 0.980 | 28 / 28 | 7 / 7 |
| 6 | windows1m | 0.976 | 0.875 | 0.980 / 0.980 | 30 / 30 | 7 / 7 |
| 20 | hotspots1m | 0.972 | 0.844 | 0.975 / 0.976 | 28 / 28 | no monolithic graph |
| 40 | hotspots1m | 0.969 | 0.813 | 0.972 / 0.973 | 28 / 28 | no monolithic graph |

The SV differences at 6 haplotypes are 25 to 40 calls out of 1,700 and run
in no consistent direction (fixed windows score highest, hotspots lowest).
The same three samples score 0.004 lower in the 20-haplotype graph than in
the 6-haplotype one, and 0.007 lower at 40; SV F1 falls from 0.855 to 0.813.
With no monolithic graph at 20 or 40 haplotypes we cannot tell whether this
comes from chunking or from calling variants in a graph with more
haplotypes. Tier 3's scaling series (10 to 40 haplotypes on 1 Mb)
found that agreement between chunked and monolithic graphs falls slowly with
N, from about 0.996 to 0.986-0.993, so some of the drop may be chunking.

## E7: HG002 and the GIAB benchmark

HG00673, HG00733 and HG002 (6 haplotypes) in 22q11.2 (chr22:18.6-21.5 Mb,
7 chunks of 414 kb) and 17q21.31 (chr17:45.4-46.6 Mb, 3 chunks of 400 kb).
F1 for all variants, both replicates (rep 1 / rep 2). HG00673 and HG00733 are
scored against the assembly-based truth calls, HG002 against GIAB v5.0q.

| Region | Settings | Graph | F1, `-a 100000` | F1, `-a 10000` | Inside duplications (`-a 10000`) | HG002 hap 1 / hap 2 (`-a 10000`, rep 1) |
|---|---|---|---|---|---|---|
| 22q11.2 | `-p 90` | monolithic | 0.898 / 0.898 | 0.914 / 0.914 | 0.789 / 0.787 | 0.943 / 0.976 |
| 22q11.2 | `-p 90` | windows0 | 0.964 / 0.965 | 0.964 / 0.965 | 0.935 / 0.936 | 0.982 / 0.982 |
| 22q11.2 | `-p 98` | monolithic | 0.189 / 0.189 | 0.944 / 0.943 | 0.869 / 0.868 | 0.982 / 0.982 |
| 22q11.2 | `-p 98` | windows0 | 0.966 / 0.966 | 0.968 / 0.968 | 0.947 / 0.946 | 0.982 / 0.982 |
| 17q21.31 | `-p 90` | monolithic | 0.694 / 0.704 | 0.694 / 0.704 | 0.565 / 0.579 | not scored |
| 17q21.31 | `-p 90` | windows0 | 0.802 / 0.802 | 0.802 / 0.802 | 0.726 / 0.725 | not scored |
| 17q21.31 | `-p 98` | monolithic | 0.856 / 0.853 | 0.857 / 0.854 | 0.789 / 0.785 | not scored |
| 17q21.31 | `-p 98` | windows0 | 0.939 / 0.938 | 0.939 / 0.938 | 0.919 / 0.917 | not scored |

- **The vcfbub trap is back at `-p 98`.** With PGGB's `-a 100000` the
  monolithic 22q11.2 VCF loses almost every call (F1 0.19, and no HG002 call
  survives), as Tier 3 found. The `-a 10000` column is the fair comparison
  for monolithic graphs.
- **SVs.** 22q11.2 has 159 truth SVs over the six haplotypes. SV F1 is 0.78 to
  0.84 for the monolithic and 0.78 to 0.86 for the chunked graphs (leaving
  out the trapped `-p 98` monolithic VCF); the
  differences are within a few calls. 17q21.31 has 9 truth SVs, too few to
  rank the graphs (chunked 0.95 to 1.00, monolithic 0.45 to 0.80).
- **Collapsed paralogs.** In 17q21.31 the monolithic graph holds 1.35 Mb of
  node sequence at `-p 90` (1.48 Mb at `-p 98`) for a 1.2 Mb region with six
  haplotypes, against 2.08 Mb (2.27 Mb) chunked. Up to 32% of a haplotype's
  length revisits nodes in the monolithic graph, against at most 10% chunked.
  In 22q11.2 the figures are 9 to 11% against 1.5 to 5.5%, as in Tier 3.
- **Cost.** The hard regions are small enough that both kinds of graph build
  in a few minutes; chunked builds took 0.67 to 1.55 times the CPU of the
  monolithic ones and at most 1.04 GB, against 1.9 GB monolithic in 22q11.2.
  The two replicates of the same build differed by up to 1.7 times in CPU
  time, so shared GitHub runners are too noisy for CPU ratios this close
  to 1 in small regions.

## Limits

- One run per haplotype count on the chr21 arm. The hard-region replicates
  agreed to within 0.01 F1, but the arm's small SV differences could move
  with a second build.
- The monolithic arm at 20 and 40 haplotypes did not finish on a 16 GB runner,
  so at those sizes chunked graphs are compared to the truth only, not to a
  monolithic graph. A larger machine would give that comparison; the plan
  limited this round to free runners.
- GIAB v5.0q has no benchmark regions in 17q21.31 and leaves out much of the
  LCR22 duplications in 22q11.2, which is where monolithic and chunked graphs
  differ. Its HG002 calls come from the same Q100 assembly that went into the
  graph, so they test whether the graph reproduces the assembly, not whether
  the assembly is right.
- Timings come from shared GitHub runners (4 vCPU, 16 GB). Replicate builds
  differed by up to 1.7 times in CPU time.
- The 17q21.31 SV counts are small (9).

## Reproduce

Every run is a dispatch of `.github/workflows/experiment.yml` on this branch:
stage `tier4-arm` with `n` = 6, 20 or 40, and stage `tier4-hard` with
`pggb` = `p90` or `p98` and `rep` = 1 or 2. Each run commits its summary to
[`docs/tier4_runs/`](tier4_runs/) (one JSON per run, the data behind every
table here).
The `from_run` input re-summarizes an earlier run's results without
rebuilding, which is how the three hard-region summaries lost to a commit
conflict were recovered. Scripts are in `pipeline/benchmark/tier4/`; the
chunk manifests are in `pipeline/benchmark/tier4/manifests/`.


## Inputs

**Haplotypes.** For the arm, every HPRC release 2 sample whose two haplotypes
each cover chr21:14.95-46.05 Mb with one minigraph-cactus GRCh38 chain is
eligible (212 samples). The order is HG00673, HG00733, HG002, then a random
draw (seed 1, as in Tier 2) from the rest
([`haps_arm.tsv`](../pipeline/benchmark/tier4/haps_arm.tsv)). N = 6, 20 and 40
take the first 3, 10 and 20 samples, so the sets are nested and the three
truth samples are in every graph. HG002 is the Q100 v1.1 assembly in PanSN form
from the HPRC bucket (`working/HPRC_PLUS/HG002/assemblies/Q100/pansn/`,
`HG002#1` paternal). The hard regions use HG00673, HG00733 and HG002 (6
haplotypes) with the Tier 3 coordinates
([`haps_hard.tsv`](../pipeline/benchmark/tier4/haps_hard.tsv)).

**Chunks.** [`make_manifests.py`](../pipeline/benchmark/tier4/make_manifests.py)
writes the manifests. On the arm, `hotspots1m` walks the haploblocks.org block
boundaries (`block_stats.tsv`, which covers chr21:14.22-46.16 Mb without gaps)
and cuts at the first boundary at least 1 Mb past the previous cut: 29 chunks
of 1.01 to 1.27 Mb. `windows1m` is 31 fixed 1 Mb windows (N = 6 only). The
hard regions reuse the Tier 3 `windows0` cuts. Every chunking has zero overlap.

**Projection and inputs.** One job per sample fetches each haplotype's contig
from the HPRC bucket, aligns the GRCh38 region to it with minimap2 2.28
(`-x asm5 -c --eqx`) and projects every chunk boundary onto it with the same
longest-colinear-chain projection as `project_chunks.py`, writing coordinates
only. Each build job then fetches just its own slice of every haplotype over
HTTPS ([`chunk_fasta.py`](../pipeline/benchmark/tier4/chunk_fasta.py)); its
output is byte-identical to `project_chunks.py` on a tested 17q21.31 chunk.

**PGGB.** The pinned image (pggb commit 4225c6c) with the project's settings (`-p 90 -s 5000
-K 19 -k 29`, numeric seqwish `-B 100000000`), or `-p 98 -s 10000 -k 47` where
marked, 4 threads, one build per runner. CPU time and peak memory come from GNU
`time -v` around `pggb`, as in Tiers 1 to 3.

**Truth.** HG00673 and HG00733: the Tier 3 dipcall-style calls from HPRC's
winnowmap GRCh38 alignments ([`truth_calls.py`](../pipeline/benchmark/tier3/truth_calls.py)).
HG002: the GIAB v5.0q GRCh38 benchmark (small variants and structural
variants, PASS records only), split into haplotypes by the phased genotype.
GIAB built v5.0q with dipcall from the same Q100 assembly; its `hap1` is the
paternal assembly (read names `chr*_PATERNAL` in GIAB's dipcall BAM), so
`HG002#1` takes the first allele
([`giab_truth.py`](../pipeline/benchmark/tier4/giab_truth.py)). All-size scores
use the small-variant benchmark regions and SV scores (50 bp or more) use the
SV benchmark regions (`score_vs_truth.py --truth-sv`). Because the benchmark
comes from the same assembly that went into the graph, an HG002 error is a
place where the graph, or the call made from it, does not reproduce the
assembly; GIAB's curation decides where that can be judged. v5.0q has
benchmark regions over about 70% of 22q11.2 and none in chr17:45.4-46.6 Mb, so
HG002 is scored on chr21 and in 22q11.2 only.

Graph calls and scoring are as in Tier 3: `vg deconstruct`, `vcfbub`, `vcfwave`,
per-haplotype VCFs, `bcftools norm`, truvari 5.5.0. Every graph is scored with
PGGB's `vcfbub -l 0 -a 100000` and with `-a 10000`, because Tier 3 found that
the first drops nested sites when the reference path loops.
