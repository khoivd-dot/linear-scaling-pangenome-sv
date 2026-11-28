# Tier 2: more haplotypes, harder regions, other maps

Tier 1 ([TIER1_CONTROLS.md](TIER1_CONTROLS.md)) found that on 1 Mb of chr21
with five haplotypes, any zero-overlap cut works as well as a hotspot cut. Tier 2
asked whether that still holds with more haplotypes and in segmental
duplications, whether the hotspot positions are real, whether hotspot blocks
compress better, and how the stitch compares with impg. It also turned up a
bug in the pinned PGGB image that affects any build with a sequence longer
than about 1 Mb.

## Summary

- **A seqwish bug under-collapses large inputs.** PGGB v0.6.0 hands its
  default `-B 10M` to the bundled seqwish (v0.7.11-35-gaab5bee) without
  converting the suffix. seqwish reads it as 1,000,000, and its batched
  transitive closure then loses merges across batch edges. A 2.9 Mb
  monolithic build of 22q11.2 came out with 12.3 Mb of node sequence. With a
  numeric `-B` the same alignments give 4.3 Mb. The 1 Mb chr21 builds are not
  affected (checked at 5 and 40 haplotypes). `config/pipeline.yaml` now sets
  `transclose_batch: 100000000` and every runner passes it.
- **Seams stay clean up to 40 haplotypes.** Three zero-overlap chunks add 0 to
  1 bp of extra node sequence within 10 kb of a cut at N = 10, 20 and 40, and
  every haplotype path matches the monolithic graph base for base. Ten small
  chunks add up to about 100 bp near cuts. Chunked graphs differ from
  monolithic about as much as two monolithic builds differ from each other.
- **Hotspot placement still shows no advantage.** Fixed zero-overlap windows
  scored higher than hotspot chunks in four of six runs, and every gap was
  smaller than run-to-run noise (two replicates per N).
- **Chunking starts to pay at 40 haplotypes.** Monolithic CPU time grows faster
  than linearly (364, 708, 1,775 s at N = 10, 20, 40 in rep1) and peak memory
  rises from about 970 MB to 2.3 GB. Three chunks use about the same total CPU
  as monolithic at N = 40 (0.98 to 1.04 times), stay at about 970 MB each, and
  the longest chunk finishes in about 40% of the monolithic wall time.
- **In segmental duplications chunking changes the graph, not the seams.** In
  17q21.31 and 22q11.2 the seams are clean, but chunked graphs carry about
  6.5% more sequence than monolithic and truvari F1 against monolithic is 0.78
  to 0.79 and 0.87 to 0.88 (two monolithic builds agree at 0.98 and 0.97).
  At least 92% of the disagreeing calls sit in the
  duplicated blocks, where the monolithic graph aligns paralogous copies that
  fall into different chunks. Hotspot and fixed cuts give the same numbers.
- **The cut positions are real hotspots,** but blocks between them do not hold
  more identical haplotypes than random blocks of the same sizes, and impg
  already does the stitch.

## 1. The seqwish batch bug

The first pass on the hard regions showed monolithic graphs far larger than
the chunked ones: 2.26 Mb against 1.49 Mb of node sequence for 17q21.31, 12.3
Mb against 6.4 Mb for 22q11.2. The wfmash alignments were complete in both
cases (every haplotype aligned end to end to the reference at 99.8% identity
in 22q11.2), so the loss happened in seqwish. Rerunning seqwish from the same
image on the saved alignments isolated it:

| seqwish `-B` | 17q21.31 graph (bp) | 22q11.2 graph (bp) |
| --- | ---: | ---: |
| `10M` (what pggb passes) | 2,281,318 | 12,455,138 |
| `1000000` | 2,281,318 | |
| `2000000` | 1,849,589 | |
| `10000000` | 1,379,791 | |
| `30000000` | 1,379,791 | 4,366,871 |

`10M` behaves exactly like `1000000`, and the log shows 1 Mb batches. A single
2.8 Mb reference-to-haplotype alignment, induced alone, merged only reference
positions 109,564 to 1,000,000; everything after the first batch edge stayed
unmerged. In `/usr/local/bin/pggb` the user-supplied `-B` goes through
`parse_numeric`, but the default `TRANSCLOSE_BATCH=10M` does not. The DNAnexus
chunk applet in this repository also called seqwish with `-B 10M` directly. No
existing report of this was found.

Effect on earlier results: rebuilding the chr21 monolithic graphs with
`-B 100000000` gave 10,385 nodes / 1,006,527 bp at N = 5 (Tier 1 replicates:
10,386 to 10,395 / 1,006,528 to 1,006,575) and 22,830 / 1,013,665 at N = 40
(rep1: 22,898 / 1,013,701), both within run-to-run noise. The pilot, Tier 1
and the scaling runs below stand. The hard-region numbers in section 6 come
from builds with the fix; the first pass is discarded.

## 2. Do the cuts sit on hotspots in other recombination maps?

The haploblocks.org boundaries are the only source for the cuts, and its
method is not published. Five maps from the stdpopsim catalogue (GRCh38) were
checked: HapMap II, pyrho YRI, CEU and CHB (Spence and Song 2019), and the
deCODE sex-averaged map, which in that catalogue is the 2010 pedigree map
lifted to GRCh38 in 10 kb bins, not the 2019 map. A hotspot is a stretch with
at least 5 times the chromosome-wide mean rate. The statistic is the share of
the 19 haploblocks.org cuts in chr21:20-21 Mb within 5 kb of a hotspot, against
1,000 sets of 19 random positions.

| Map | Hotspots in region | Cuts within 5 kb | Random cuts | p |
| --- | ---: | ---: | ---: | ---: |
| HapMap II | 23 | 0.74 | 0.31 | 0.001 |
| pyrho YRI | 19 | 0.68 | 0.23 | < 0.001 |
| pyrho CEU | 14 | 0.58 | 0.18 | < 0.001 |
| pyrho CHB | 16 | 0.63 | 0.21 | < 0.001 |
| deCODE 2010 (10 kb bins) | 5 | 0.32 | 0.11 | 0.016 |

The four LD-based maps agree most, which is expected if haploblocks.org also
works from population LD. The pedigree map is the independent check; at 10 kb
resolution it finds only 5 hotspots here but still puts the cuts near them
more often than chance. At a 10-fold threshold it has one hotspot in the region
and no power (p = 0.32). The positions are real hotspots; the rest of this
document asks whether that matters.

## 3. Our stitch against impg

impg 0.3.3 (bioconda) does the same two jobs as this pipeline. Both were run on
the Tier 1 data.

- **Cut coordinates.** `impg partition -w 333334 -d 0` on PGGB's own wfmash
  alignments (not our minimap2 alignments) put the first cut at the same
  haplotype coordinate as `project_chunks.py` on all four haplotypes. It merged
  the last two windows into one 667 kb partition, so only one cut could be
  compared.
- **Joining.** `impg lace` on the same chunk graphs gave the same node, edge
  and base counts as our stitch, the same multiset of segment sequences, and
  haplotypes identical to the monolithic graph, for all three zero-overlap
  chunkings in all three replicates (9 of 9).

The stitch adds nothing impg does not already do. What this repository adds
is the cut choice and the scoring against a monolithic build.

## 4. Do hotspot blocks compress better?

Tier 1 left one argument for hotspot cuts untested: haplotypes inside a
recombination-bounded block should repeat more often, so a builder could
collapse identical copies first. For each block and N = 10, 20 and 40
haplotypes, the number of distinct haplotype sequences (exact identity) gives a
length-weighted duplicate fraction, 1 - distinct / N, over blocks of at least
5 kb. The 19 haploblocks.org cuts make 20 blocks; the null is 20 layouts with
the same block sizes shuffled and rotated around the region.

| N | Hotspot blocks | Random layouts, median (range) | Random layouts with more duplication |
| ---: | ---: | ---: | ---: |
| 10 | 0.011 | 0.013 (0.003-0.030) | 14 of 20 |
| 20 | 0.010 | 0.018 (0.011-0.037) | 20 of 20 |
| 40 | 0.029 | 0.038 (0.026-0.054) | 18 of 20 |

Exact duplicates are rare in assemblies: at N = 40 about 3% of block sequence
could be dropped, and hotspot blocks do worse than most random layouts.
haploblocks.org reports much stronger grouping for the same blocks (5,096
haplotypes in 3 to 1,002 MMseqs2 clusters per block, median about 120), but
those are clusters of sequences built from phased variant calls, not exact
copies of assembled sequence. Near-identity clustering was not tested.

## 5. Scaling to 40 haplotypes

Samples were drawn from HPRC release 2: HG00673 and HG00733 (the pilot pair),
then 18 more chosen at random (seed 1) among the 229 samples whose
chr21:19.95-21.05 Mb sits on a single minigraph-cactus chain in both
haplotypes. Haplotype sets are nested (the first 5, 10 or 20 samples). Each
contig was fetched remotely, aligned with minimap2 2.28, and cut with the same
manifests as Tier 1. Two replicates per N, each scored against its own
monolithic build; PGGB with pggb's default `-B` throughout, so the runs match
Tier 1 (1 Mb inputs are not affected by the bug, section 1). Numbers are
rep1 / rep2.

Monolithic builds:

| N haplotypes (+ GRCh38) | Nodes | Node sequence (bp) | Wall (s) | CPU (s) | Peak memory (MB) | truvari F1, rep1 vs rep2 |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 4 (Tier 1, 3 reps) | 10,386-10,395 | 1,006,528-1,006,575 | 112-135 | 221-266 | about 970 | 0.9995-1.0 |
| 10 | 13,362 / 13,449 | 1,008,347 / 1,008,401 | 163 / 167 | 364 / 401 | 972 / 973 | 0.9960 |
| 20 | 18,214 / 18,196 | 1,011,283 / 1,011,251 | 250 / 330 | 708 / 940 | 1,197 / 1,203 | 0.9914 |
| 40 | 22,898 / 22,904 | 1,013,701 / 1,013,633 | 536 / 671 | 1,775 / 2,201 | 2,335 / 2,359 | 0.9975 |

Chunked builds against the monolithic build of the same N and replicate
(windows0, hotspots300, hotspots50):

| | | N = 10 | N = 20 | N = 40 |
| --- | --- | --- | --- | --- |
| Extra node sequence within 10 kb of a cut (bp) | windows0 | 0 / 0 | 0 / 1 | 0 / 0 |
| | hotspots300 | 0 / 0 | -2 / 5 | 1 / 0 |
| | hotspots50 | 101 / 100 | 50 / 60 | 30 / -36 |
| Extra node sequence, total (bp) | windows0 | 1 / -57 | 7 / 36 | -134 / 61 |
| | hotspots300 | 20 / -163 | -24 / 39 | 71 / -271 |
| | hotspots50 | -44 / -127 | -12 / -181 | -47 / 172 |
| truvari F1 vs monolithic | windows0 | 0.9982 / 0.9964 | 0.9978 / 0.9933 | 0.9961 / 0.9924 |
| | hotspots300 | 0.9977 / 0.9966 | 0.9930 / 0.9936 | 0.9930 / 0.9895 |
| | hotspots50 | 0.9951 / 0.9961 | 0.9929 / 0.9923 | 0.9905 / 0.9921 |
| truvari FN within 10 kb of a cut | windows0 | 0 / 0 | 7 / 3 | 0 / 0 |
| | hotspots300 | 0 / 0 | 6 / 9 | 0 / 0 |
| | hotspots50 | 3 / 2 | 21 / 26 | 32 / 27 |
| SVs (50 bp+) missed, of 40 / 67 / 112 | windows0 | 0 / 0 | 0 / 2 | 0 / 5 |
| | hotspots300 | 0 / 0 | 1 / 2 | 2 / 4 |
| | hotspots50 | 3 / 0 | 1 / 0 | 5 / 5 |
| CPU, all chunks (s) | windows0 | 456 / 602 | 798 / 1,039 | 1,790 / 2,282 |
| | hotspots300 | 469 / 567 | 780 / 1,006 | 1,742 / 2,297 |
| | hotspots50 | 663 / 876 | 1,016 / 1,280 | 1,771 / 2,320 |
| Longest chunk, wall (s) | windows0 | 93 / 129 | 123 / 156 | 215 / 278 |
| | hotspots300 | 96 / 115 | 123 / 155 | 210 / 284 |
| | hotspots50 | 79 / 101 | 104 / 135 | 163 / 215 |

Every chunked build kept all haplotypes identical to monolithic and passed
every seam check. Peak memory per chunk was 968 to 970 MB in every build
except the smallest chunk at N = 40 rep1 (853 MB).

What the numbers say:

- **Chunked graphs differ from monolithic about as much as two monolithic
  builds differ from each other.** Monolithic rep1 against rep2 scores F1
  0.996, 0.991 and 0.9975 at N = 10, 20 and 40. Chunked against monolithic
  scores 0.990 to 0.998 across all 18 comparisons. Only hotspots300 at N = 40
  rep2 (0.9895) falls below the replicate range, and only just.
- **The seams hold.** With three chunks the near-cut excess is 0 to 5 bp in
  every run. The truvari errors counted as "within 10 kb of a cut" are not at
  the cut, and they repeat across replicates: at N = 20 the near-cut errors of
  windows0 and hotspots300 are one locus at about 20.657 Mb, 9 kb before one
  cut and 7 kb after the other. With ten chunks, 88 of the 111 near-cut errors
  across both replicates sit 3 to 5 kb after the cut at 20,583,159. This
  matches the Tier 1 observation that a few loci move with chunk context; the
  scale (a few kb, near the 5 kb wfmash segment) suggests the first or last
  segment of a chunk maps differently, but this was not tested.
- **Hotspot placement shows no advantage.** windows0 scored higher F1 than
  hotspots300 in four of six runs (both at N = 40, one each at N = 10 and 20),
  and every gap is smaller than the gap between two monolithic builds.
  Ten small chunks lose up to 5 SVs.
- **F1 drifts down with N, and so does replicate agreement.** For three
  chunks F1 falls from 0.996-0.998 at N = 10 to 0.990-0.996 at N = 40, but
  monolithic replicates also agree less at N = 20 and 40 than at N = 5. The
  drift is at least partly PGGB's own run-to-run variation growing with N.
- **Compute.** At N = 40 the three-chunk builds use 0.98 to 1.04 times the
  monolithic CPU (against 1.3 to 1.5 at N = 5 to 10), the longest chunk takes
  39 to 42% of the monolithic wall time, and peak memory per build stays near
  970 MB while monolithic reaches 2.3 to 2.4 GB. Replicate 2 ran about 25%
  slower throughout (shared machine), so compare ratios within a replicate.

## 6. Hard regions: 17q21.31 and 22q11.2

Regions: GRCh38 chr17:45.4-46.6 Mb (the 17q21.31 inversion region) and
chr22:18.6-21.5 Mb (22q11.2 with LCR22A to D; includes 50 kb of N), with
HG00673 and HG00733, both haplotypes. Contigs come from the release chains as
before. Haplotype lengths show the copy-number variation: 1.20 to 1.49 Mb for
the 1.2 Mb 17q21.31 region, and 2.78 to 4.70 Mb for the 2.9 Mb 22q11.2 region,
almost all of the extra in the first 414 kb window, which holds LCR22A (0.34
to 2.24 Mb per haplotype).

The projection had to change. The greedy chain in `project_chunks.py` locked
onto pieces that mapped to the wrong copy of an LCR and covered 0.2 to 2% of
22q11.2. It now takes the colinear chain with the most aligned reference bases
(dynamic programming); coverage is 0.89 to 0.96, no cut falls in a chain gap,
every chunking concatenates back to the monolithic haplotypes, and all 40
chr21 projections are byte-identical to the old ones.

Chunkings: monolithic, fixed zero-overlap windows (3 x 400 kb and 7 x 414 kb),
and haploblocks.org cuts with a 300 kb minimum (3 and 8 chunks). Two
replicates, all with the seqwish fix. Numbers are rep1 / rep2.

| | 17q21.31 windows0 | 17q21.31 hotspots300 | 22q11.2 windows0 | 22q11.2 hotspots300 |
| --- | --- | --- | --- | --- |
| Monolithic node sequence (bp) | 1,366,429 / 1,366,418 | | 4,347,015 / 4,346,516 | |
| Extra node sequence, total (bp) | 90,448 / 90,443 | 92,126 / 92,124 | 281,553 / 284,137 | 284,659 / 285,919 |
| Extra within 10 kb of a cut (bp) | 0 / 0 | 0 / 0 | 1 / 1 | -543 / -543 |
| Seams PASS, haplotypes identical | 2/2, 5/5 | 2/2, 5/5 | 6/6, 5/5 | 7/7, 5/5 |
| truvari F1 vs monolithic | 0.794 / 0.781 | 0.795 / 0.786 | 0.872 / 0.879 | 0.874 / 0.875 |
| truvari F1, monolithic rep1 vs rep2 | 0.983 | | 0.973 | |
| Share of truvari errors in duplicated blocks | 0.95 / 0.93 | 0.94 / 0.92 | 0.93 / 0.92 | 0.92 / 0.92 |
| SVs (50 bp+) missed | 14 / 14 | 14 / 14 | 45 / 42 | 39 / 37 |
| CPU, all chunks vs monolithic (s) | 336 / 324 vs 267 / 281 | 323 / 327 | 1,147 / 1,128 vs 1,014 / 1,028 | 1,168 / 1,202 |
| Peak memory, chunk vs monolithic (MB) | 967 vs 446 / 431 | 967 | 969 vs 1,185 / 1,123 | 970 |

"Duplicated blocks" are the 100 kb bins where the errors cluster:
chr17:45.5-45.6 and 46.5-46.6 Mb (the duplicated segments at both ends of the
inversion region), and chr22:18.6-19.0, 20.3-20.4 and 21.1-21.4 Mb, which
line up with LCR22A, B and D. Together they are 17% of the 17q21.31 region and
28% of 22q11.2.

What this shows:

- **The seams are fine.** Extra sequence near cuts is 0 to 1 bp except one
  hotspot cut in 22q11.2 (-543 bp, the same in both replicates, where the
  scoring places a large LCR insertion next to the cut).
- **Chunks cannot align paralogs that sit in different chunks.** The
  monolithic graph aligns copies of a duplicated block to each other; chunks
  of 300 to 400 kb separate LCR22A, B and D, and the two 17q21.31 flanking
  copies, so those copies stay apart. That is where the 6.5% extra sequence and
  at least 92% of the truvari disagreement sit. Whether the collapsed or the
  separated representation is closer to the truth is not settled by this
  comparison; truvari here measures agreement with monolithic, not accuracy.
- **This is beyond run-to-run noise, which is also larger here.** Two
  monolithic builds of the same input agree at F1 0.983 (17q21.31) and 0.973
  (22q11.2), against 0.9995 on chr21. Chunked against monolithic is 0.78 to
  0.88.
- **Cut placement does not matter here either.** Fixed and hotspot cuts give
  F1 within 0.01 of each other, and the errors are at the same loci.
- **Small inputs do not save memory at this size.** The 17q21.31 monolithic
  build peaks at about 440 MB while every chunk sits at the 970 MB floor.

## What this changes

- **Fix in place.** Any build with a sequence over about 1 Mb, including a
  chunk, needs the numeric `-B`. It is now the default in every runner.
- **Placement.** Across two regions, three haplotype counts and the Tier 1
  controls, recombination hotspots have not beaten fixed zero-overlap cuts
  once. The useful rules are: no overlap, chunks of a few hundred kb, and
  replicates for every comparison.
- **Where chunking helps.** On chr21 the payoff appears between 20 and 40
  haplotypes: equal CPU, lower wall time when chunks run in parallel, and a
  memory ceiling that does not grow with N.
- **Where chunking needs care.** In segmental duplications the cut should
  keep paralogous copies in one chunk, or the graph will represent them
  differently from a monolithic build. A duplication-aware cut rule (from a
  segmental duplication track) would be the next thing to test, ahead of
  anything recombination-based.

## Limits

- Scaling: two replicates per N, one region of chr21.
- Hard regions: two samples (four haplotypes) per region.
- truvari compares against the monolithic graph, not a truth set.
- Compression used exact identity only.
- The near-cut loci and the paralog explanation in section 6 come from where
  the errors fall, not from inspecting the graphs base by base.
- Timings are from a shared 4 vCPU cloud machine; the container restarted
  twice during the runs and interrupted builds were rerun from scratch.

## Reproduce

Scripts as run are in `pipeline/benchmark/tier2/`:

- `select_samples.py`: eligible samples and order from the release chains.
- `scaling_fetch_align.sh`, `scaling_project_all.sh`, `scaling_run_scale.sh`:
  scaling inputs, projection and builds.
- `map_robustness.py`: section 2.
- `block_compression.py`: section 4 (`--write-layouts` makes the random
  layouts).
- `hard_prepare.sh`, `hard_fetch_align.sh`, `hard_run_fixb.sh`: section 6.

Raw numbers: [`tier2.json`](tier2.json).
