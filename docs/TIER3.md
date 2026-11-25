# Tier 3: accuracy against a truth set, duplication-aware chunks, 10 Mb

Tiers 1 and 2 ([TIER1_CONTROLS.md](TIER1_CONTROLS.md), [TIER2.md](TIER2.md))
compared every chunked graph with a monolithic build of the same input. That
shows whether chunking changes the graph, not whether the change is for better
or worse. In segmental duplications the two disagreed (truvari F1 0.78 to 0.88),
and Tier 2 left open which one was right. Tier 3 scores both against variant
calls made straight from the assemblies, tests whether keeping paralogous
copies in one PGGB job brings chunked graphs back in line, and runs 10 Mb to see
how cost and accuracy behave at a larger scale.

## Summary

- **Outside duplications, chunking costs no accuracy.** On chr21 (1 Mb and
  10 Mb) the monolithic and every chunked graph score the same against the
  truth calls: F1 0.979 to 0.981 for all variants and 0.928 to 0.943 for SVs of
  50 bp or more. Chunk size (400 kb, 1 Mb or 2.5 Mb) and cut placement
  (fixed or haploblocks.org hotspots) do not change this.
- **In 22q11.2 the chunked graph is closer to the truth than the monolithic
  one.** F1 is 0.962 against 0.880 overall, and 0.937 against 0.719 inside the
  LCR22 duplications. The monolithic graph loses recall there (0.66 against
  0.93) because it routes about 9% of each haplotype through nodes the same
  haplotype has already used, which merges paralogous copies; the chunked
  graph does this for 2.4%.
- **17q21.31 is mixed.** Chunking raises small-variant precision inside the
  duplicated flanks (0.63 against 0.50) but adds SV false positives there
  (7 against 2, of 5 true SVs), so SV F1 drops from 0.90 to 0.72. The SV counts
  are small.
- **Duplication-aware groups undo the gain.** Building the chunks that share a
  duplication in one job (found by self-alignment of the reference: 17q21.31
  chunks 1 and 3; 22q11.2 chunks 1, 5 and 7, which hold LCR22A, B and D) makes
  the graph more like the monolithic one (F1 against monolithic 0.65 to 0.79
  in 17q21.31, 0.66 to 0.72 in 22q11.2) and less like the truth (22q11.2 F1
  0.889). The grouped job collapses paralogs even more than the monolithic
  build (11 to 14% of haplotype bp revisit nodes) and, in 22q11.2, takes as
  long as the whole monolithic build.
- **At 10 Mb chunking pays off on wall time and memory.** The monolithic build
  of chr21:20-30 Mb took 582 s and 2.8 GB; 25 chunks of 400 kb took at most
  75 s and 970 MB each, for 1.5 times the CPU. With 1 Mb chunks the CPU
  overhead is 1.22 times, with 2.5 Mb chunks 1.06 times.
- **impg lace and our stitch give identical graphs** on all 4 regions (same
  nodes, sequence and truth scores).

PENDING: the recommended-settings control (E5). Both replicates of E1, E2 and
the scaling series (E3) are in, and they agree.

## Truth set

For each haplotype of HG00673 and HG00733 we read the HPRC release 2 winnowmap
alignments of the assembly to GRCh38 (no alt) from the HPRC S3 bucket
(`working/<project>/<sample>/assemblies/release2/alignments/winnowmap/GRCh38_no_alt/`,
hap1 is paternal, PanSN `#1`). These alignments were made by HPRC, not by us,
and with a different aligner from the minimap2 projection and the wfmash
alignments inside PGGB. `truth_calls.py` follows dipcall and paftools call:

- a reference base is callable when exactly one alignment of at least 10 kb
  covers it, and that alignment is at least 50 kb long with MAPQ 5 or more;
- inside callable regions, SNVs, insertions and deletions are read off the
  CIGAR against the reference.

| Region | Haplotype calls | Callable share |
| --- | --- | --- |
| chr21:20-21 Mb | 1,939 to 2,127 | 0.998 to 0.999 |
| chr17:45.4-46.6 Mb (17q21.31) | 721 to 894 | 0.75 to 0.99 |
| chr22:18.6-21.5 Mb (22q11.2) | 3,311 to 4,336 | 0.83 to 0.90 |
| chr21:20-30 Mb | 13,681 to 15,395 | 0.992 to 0.993 |

Graph calls: `vg deconstruct -a` on the GRCh38 path, then `vcfbub -l 0 -a
100000` and `vcfwave` as in the PGGB pipeline, split into one haploid VCF per
haplotype, and both sides left-aligned with `bcftools norm`. truvari 5.5.0
compares them inside each haplotype's callable BED, once on all sizes
(SNVs included; `--pick multi --refdist 500 --pctseq 0.7 --pctsize 0.7`, as in
Tiers 1 and 2) and once on SVs of 50 bp or more. Calls are split by whether
they fall in a GIAB v3.3 segmental duplication (`GRCh38_segdups.bed`).

The truth set has limits worth stating up front. Inside duplications the
assembly alignment itself can be wrong, and the callable filter drops 10 to
25% of the duplicated regions. The comparison says which graph agrees with
HPRC's assembly-to-reference alignments, which is the standard dipcall-style
benchmark, not which graph is right base by base.

## E1: chunked against monolithic, scored against the truth

PGGB settings as in Tier 2 (`-p 90 -s 5000 -K 19 -k 29`, numeric seqwish
`-B`), 4 haplotypes plus GRCh38. Rep1 numbers; rep2 is below. Chunks are zero-overlap
windows (`windows0`: 3 x 333 kb on chr21, 3 x 400 kb in 17q21.31,
7 x 414 kb in 22q11.2).

| Region, graph | Precision | Recall | F1 | F1 in SDs | F1 elsewhere | SV F1 (n) |
| --- | --- | --- | --- | --- | --- | --- |
| chr21 1 Mb, monolithic | 0.977 | 0.985 | 0.981 | (no SDs) | 0.981 | 0.943 (27) |
| chr21 1 Mb, windows0 | 0.977 | 0.985 | 0.981 | | 0.981 | 0.943 |
| chr21 1 Mb, hotspots300 | 0.977 | 0.985 | 0.981 | | 0.981 | 0.943 |
| 17q21.31, monolithic | 0.632 | 0.903 | 0.743 | 0.627 | 0.973 | 0.900 (9) |
| 17q21.31, windows0 | 0.726 | 0.902 | 0.804 | 0.720 | 0.973 | 0.720 |
| 22q11.2, monolithic | 0.909 | 0.853 | 0.880 | 0.719 | 0.973 | 0.761 (103) |
| 22q11.2, windows0 | 0.964 | 0.960 | 0.962 | 0.937 | 0.978 | 0.818 |

Outside duplications every graph scores 0.973 to 0.981, so the differences
are all in the duplicated blocks. There the two regions behave differently:

- **22q11.2.** Inside the LCR22s the monolithic graph misses a third of the
  truth calls (recall 0.660) and the chunked graph misses 7% (0.933), with
  higher precision as well (0.941 against 0.790). SV recall rises from 0.68 to
  0.89; SV precision falls from 0.86 to 0.75.
- **17q21.31.** Both graphs have low precision in the duplicated flanks
  (0.50 and 0.63) at the same recall (0.85). The SV numbers rest on 5 true
  SVs in duplications, where the chunked graph adds 5 more false positives.

### Why: collapsed paralogs

`path_revisits.py` counts, for every haplotype path, the bases it spells
through a node it has already passed. A path can only revisit a node if the
graph merged two copies of a sequence from the same haplotype.

| Region | Monolithic | windows0 | sdgroup |
| --- | --- | --- | --- |
| chr21 1 Mb | 0.000 | 0.000 | |
| 17q21.31 | 0.105 | 0.044 | 0.138 |
| 22q11.2 | 0.090 | 0.024 | 0.115 |

With `-p 90`, wfmash maps segments between paralogous copies that are 90% or
more identical, and seqwish then merges them, so a haplotype's variants in one
copy are expressed relative to a merged node rather than to its own GRCh38
locus. Chunks that separate the copies cannot do this. E5 tests whether the
stricter settings PGGB recommends for human data avoid the problem without
chunking.

## E2: duplication-aware chunk groups

`sd_groups.py` maps 20 kb tiles of the reference region back to the region
(minimap2 asm20 via mappy) and links two chunks when hits of 5 kb or more at
90% identity or more join them with at least 10 kb in total. Linked chunks are
built as one PGGB job from the concatenation of their FASTAs; the other
chunks reuse the windows0 graphs; `impg lace` joins everything
(`group_build.py`).

| | 17q21.31 windows0 | 17q21.31 sdgroup | 22q11.2 windows0 | 22q11.2 sdgroup |
| --- | --- | --- | --- | --- |
| Groups | 3 singletons | chunks 1+3 (78 kb linked), 2 | 7 singletons | chunks 1+5+7 (466 kb linked), 2, 3, 4, 6 |
| Node sequence vs monolithic (bp) | +90,473 | -53,620 | +284,710 | -106,654 |
| Exact-match F1 vs monolithic | 0.646 | 0.794 | 0.660 | 0.719 |
| Truth F1, all / in SDs | 0.804 / 0.720 | 0.750 / 0.637 | 0.962 / 0.937 | 0.889 / 0.760 |
| Truth F1, SVs | 0.720 | 0.706 | 0.818 | 0.789 |
| CPU, all builds (s) | 381 | 452 | 1,284 | 1,860 |
| Longest build (s) | 94 | 118 | 159 | 398 |
| Monolithic CPU / wall (s) | 320 / 120 | | 1,130 / 383 | |

Grouping does what it was designed to do, bringing the chunked graph closer to
the monolithic one, and that makes it worse against the truth. The grouped job
holds the paralogs and little of the unique sequence between them, so it
collapses them even harder than the monolithic build (node sequence below the
monolithic graph in both regions). In 22q11.2 the group of LCR22A, B and D
takes as long as the whole monolithic build. On this evidence the duplication-
aware cut rule proposed at the end of Tier 2 should not be adopted.

## Second replicate of E1 and E2

Every E1 and E2 graph was rebuilt from the same inputs and scored the same
way. Truth F1, rep1 / rep2:

| Region, graph | F1 all | F1 in SDs | SV F1 | Path revisits |
| --- | --- | --- | --- | --- |
| chr21 1 Mb, monolithic | 0.981 / 0.981 | | 0.943 / 0.943 | 0.000 / 0.000 |
| chr21 1 Mb, windows0 | 0.981 / 0.981 | | 0.943 / 0.943 | 0.000 / 0.000 |
| chr21 1 Mb, hotspots300 | 0.981 / 0.981 | | 0.943 / 0.943 | |
| 17q21.31, monolithic | 0.743 / 0.749 | 0.627 / 0.637 | 0.900 / 0.900 | 0.105 / 0.105 |
| 17q21.31, windows0 | 0.804 / 0.802 | 0.720 / 0.717 | 0.720 / 0.720 | 0.044 / 0.044 |
| 17q21.31, sdgroup | 0.750 / 0.750 | 0.637 / 0.636 | 0.706 / 0.706 | 0.138 / 0.138 |
| 22q11.2, monolithic | 0.880 / 0.880 | 0.719 / 0.713 | 0.761 / 0.755 | 0.090 / 0.090 |
| 22q11.2, windows0 | 0.962 / 0.962 | 0.937 / 0.936 | 0.818 / 0.814 | 0.024 / 0.024 |
| 22q11.2, sdgroup | 0.889 / 0.889 | 0.760 / 0.760 | 0.789 / 0.774 | 0.115 / 0.115 |

All differences between the replicates are 0.015 or less, and most are under
0.006, against gaps of 0.05 to 0.22 between graphs in the duplication regions.
The two monolithic builds of 22q11.2 differ by 183 nodes and 167 bp, but their
inside-SD recall is the same (0.660 and 0.651), as is the revisit share. So the
E1 and E2 conclusions do not rest on one lucky build. Rep2 ran faster on the
same machine (22q11.2 monolithic 291 s against 383 s wall), and the grouped
22q11.2 job again took longer than the whole monolithic build (307 s against
291 s). Per-path revisit counts are in
[`tier3_path_revisits_rep1.tsv`](tier3_path_revisits_rep1.tsv) and
[`tier3_path_revisits_rep2.tsv`](tier3_path_revisits_rep2.tsv).

## E4: 10 Mb of chr21

chr21:20,000,000-30,000,000, the same 4 haplotypes plus GRCh38, one replicate.
All four haplotypes sit on one contig across the 10 Mb.

| | Monolithic | 25 x 400 kb | 10 x 1 Mb | 4 x 2.5 Mb |
| --- | --- | --- | --- | --- |
| Node sequence vs monolithic (bp) | 10,066,965 | +278 | -180 | -303 |
| Seams PASS, haplotypes identical | | 24/24, 5/5 | 9/9, 5/5 | 3/3, 5/5 |
| Exact-match F1 vs monolithic | | 0.995 | 0.997 | 0.998 |
| Truth F1, all / SVs | 0.979 / 0.928 | 0.979 / 0.936 | 0.979 / 0.928 | 0.979 / 0.928 |
| CPU, all builds (s) | 1,681 | 2,544 | 2,057 | 1,785 |
| CPU vs monolithic | 1.00 | 1.51 | 1.22 | 1.06 |
| Longest build (s) | 582 | 75 | 109 | 189 |
| Peak memory per build (MB) | 2,840 | 969 | 971 | 979 |

Monolithic CPU rose 6.2 times from 1 Mb to 10 Mb (271 to 1,681 s), less than
the 10 times in length, because the 1 Mb build is dominated by fixed costs.
Each chunk pays its own fixed cost too, which is why 400 kb chunks cost half as
much CPU again; at 1 Mb and above the overhead falls to 6 to 22%. With one
worker per chunk the 10 Mb region finishes 3 to 8 times sooner and never needs
more than about 1 GB.

## E3: second scaling replicate

Tier 2's haplotype-scaling series (chr21:20-21 Mb, GRCh38 plus 10, 20 and 40
HPRC haplotypes, the same samples in the same order, seed 1) ran once. E3
repeats it with the same inputs and the same seqwish `-B` default as rep1
(`PGGB_TRANSCLOSE_BATCH=pggb-default`; at 1 Mb the bug does not apply). Each
chunked graph is scored against the monolithic build of its own replicate with
`seam_metrics.py`, as in Tiers 1 and 2.

| N | Graph | Extra bp near cuts, rep1 / rep2 | truvari F1 vs monolithic, rep1 / rep2 |
| --- | --- | --- | --- |
| 10 | windows0 | 0 / 0 | 0.9982 / 0.9956 |
| 10 | hotspots300 | 0 / 0 | 0.9977 / 0.9962 |
| 10 | hotspots50 | 101 / 100 | 0.9951 / 0.9963 |
| 20 | windows0 | 0 / 4 | 0.9978 / 0.9982 |
| 20 | hotspots300 | -2 / -2 | 0.9930 / 0.9931 |
| 20 | hotspots50 | 50 / 7 | 0.9929 / 0.9923 |
| 40 | windows0 | 0 / 0 | 0.9961 / 0.9934 |
| 40 | hotspots300 | 1 / 0 | 0.9930 / 0.9861 |
| 40 | hotspots50 | 30 / -35 | 0.9905 / 0.9888 |

Every seam passed and every haplotype path was identical to the monolithic
one in both replicates (2/2 or 9/9 seams, 11/11, 21/21 and 41/41 paths).

- **The ranking holds.** Fixed zero-overlap windows scored higher than both
  hotspot layouts in 5 of the 6 replicate-by-N comparisons; the exception
  (N = 10, rep2) is 0.9956 against 0.9962 and 0.9963. At N = 20 and 40
  windows0 is ahead in both replicates, by 0.003 to 0.007.
- **Agreement with monolithic falls slowly with N** for every chunking, from
  about 0.996 at N = 10 to 0.986 to 0.993 at N = 40. None of the extra
  disagreement sits at the cuts for the three-chunk layouts; with ten chunks
  (hotspots50) up to 46 of 99 false positives at N = 40 are within 10 kb of a
  cut. SVs of 50 bp or more missed against monolithic rose with N in rep2
  as in rep1: none for windows0 at N = 10 and 20, 4 of the monolithic
  graph's SVs at N = 40 (rep1: 0), and 0 to 7 for the hotspot layouts.
- **Compute repeats closely.** At N = 40 the monolithic build took 525 s wall,
  1,745 s CPU and 2.2 GB (rep1: 536 s, 1,775 s, 2.3 GB). Three windows0 chunks
  used 1,779 s CPU (1.02 times), the longest finished in 211 s (40% of
  monolithic) and none used more than 971 MB.
- **Monolithic builds drift between replicates** by 4 to 222 nodes (13,362
  against 13,442 at N = 10; 22,898 against 23,120 at N = 40), the same
  run-to-run noise seen in Tier 1, now growing with N.

## E5: PGGB's recommended human settings (running)

## What this changes

- The accuracy question from Tier 2 has an answer for the regions tested:
  where chunked and monolithic graphs disagree, the chunked graph is at least
  as close to assembly-based truth, and in 22q11.2 much closer. Chunk
  boundaries that separate paralogous copies act as a guard against collapsing
  them.
- Duplication-aware grouping is not worth its cost.
- For runs longer than a few Mb, chunks of 1 Mb or more keep the CPU overhead
  near 1.06 to 1.22 times while capping memory near 1 GB.

## Limits

- Two samples (four haplotypes) for the truth comparison; two replicates in
  the duplication regions, one for the 10 Mb runs.
- The truth set covers 75 to 99% of each region and is weakest inside
  duplications.
- SV counts in 17q21.31 duplications are small (5).
- Timings come from a 4 vCPU cloud machine; the container restarted once
  during the runs, and interrupted builds were rerun from scratch.
- The minimap2 used for projection here is mappy 2.31 (minimap2's own Python
  binding) behind a small command-line shim, not the 2.28 binary used in
  Tier 2.

## Reproduce

Scripts as run are in `pipeline/benchmark/tier3/`; `setup_inputs.sh`,
`run_truth.sh`, `run_e1e2.sh`, `run_e4.sh`, `scaling_select.sh`,
`scaling_rep2.sh`, `run_e5.sh`, and `run_queue.sh` for the order. The same
stages run on a GitHub runner from `.github/workflows/experiment.yml`.
`summarize.py` writes [`tier3.json`](tier3.json).
