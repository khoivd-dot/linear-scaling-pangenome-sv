# Building pangenome graphs in chunks: zero-overlap seams are clean, cut placement does not matter, and chunk boundaries can limit paralog collapse

Authors: Nhat Quang Tran, Dang Khoi Vu, Michael Olufemi, Alexander W. Lehr, Ali Rizvi

Draft of 2026-11-25. Every number below is copied from the result documents and JSON files in `docs/` of the code repository. All planned experiments are complete; figures are specified in captions but not yet rendered.

## Abstract

**Motivation.** PGGB graph construction gets more expensive as haplotypes are added and regions grow. Cutting a region into chunks, building each separately and joining the graphs is an obvious workaround, but in the sources we searched we found no measurement of how such graphs differ from a monolithic build.

**Results.** With HPRC release 2 assemblies and GRCh38, we built PGGB v0.6.0 graphs of chr21 (1 and 10 Mb), 17q21.31 and 22q11.2, monolithically and in chunks cut at fixed, shifted, random or recombination-hotspot positions. A pilot advantage of hotspot cuts came from the 50 kb overlap of the fixed windows. Without overlap, three-chunk builds added at most 2 bp of node sequence within 10 kb of a cut at 4 to 40 haplotypes, haplotype paths matched the monolithic graph, and hotspot cuts never gave higher variant agreement than fixed cuts. PGGB v0.6.0 passes its default seqwish `-B 10M` unparsed, leaving a 2.9 Mb build with 12.3 Mb of node sequence instead of 4.3 Mb. At 40 haplotypes, three chunks used the same CPU time as a monolithic build and peaked near 970 MB of memory, against 2.3 GB. In segmental duplications, chunked and monolithic graphs disagreed (truvari F1 0.78 to 0.88). Against assembly-based truth calls, the chunked graph scored higher in 22q11.2 (F1 0.962 against 0.880, reproduced in a second replicate) and mixed in 17q21.31, and monolithic haplotype paths revisited nodes for 9 to 11% of their length, against 2 to 4% when chunked, which points to paralog collapse. PGGB's recommended human settings (`-p 98 -s 10000 -k 47`) reduced but did not remove the collapse. They also exposed a calling trap: the reference path loops in monolithic graphs of duplicated regions, and the PGGB pipeline's `vcfbub -l 0 -a 100000` then drops hundreds of kilobases of nested variants (F1 0.24 to 0.47). With a 10 kb bubble cut the monolithic graphs scored 0.86 and 0.93, still below chunked graphs at the same settings (0.94 and 0.97).

**Availability.** https://github.com/khoivd-dot/linear-scaling-pangenome-sv

## 1. Introduction

A pangenome graph stores many haplotypes of a region as one structure, with each haplotype a path through shared and divergent sequence. PGGB builds such graphs without a reference backbone in three stages: all-to-all alignment with wfmash, graph induction by transitive closure of aligned bases with seqwish, and smoothing with smoothxg, followed by normalisation and graph statistics [1-4]. Its all-to-all alignment grows roughly with the square of the number of haplotypes, and induction and smoothing run over the full length of the input. Both limits matter for human cohorts. HPRC release 1 reported 47 phased diploid assemblies [5], and the release 2 preprint reports 460 haplotypes from 232 individuals [6].

Large graph resources already partition the genome coarsely. Minigraph-Cactus and the HPRC release 1 graphs are built per chromosome [5,7], and Progressive Cactus decomposes multiple-genome alignment into subproblems along a guide tree [8]. Finer partitioning has two published forms that we know of. nf-core/pangenome splits the wfmash alignment step into jobs across a cluster but induces one graph from the combined alignments [9], so it chunks the alignment, not the graph. impg partitions a cohort into windows (`partition`), builds each window's graph and merges the per-window graphs (`lace`) [10]. Its documentation states that it makes no accuracy comparison against a whole-graph build and asks users to measure on their own data [10]. A companion manuscript may contain such a validation; we could not read it [10].

Chunking raises three practical questions. Do the joins (seams) between chunk graphs leave artefacts? Does it matter where the cuts fall? Population genetics suggests cutting at recombination hotspots, where haplotype blocks end, so that the haplotypes inside a chunk share more sequence. And what happens in segmental duplications (SDs), where paralogous copies of a sequence lie hundreds of kilobases apart and a cut can place them in different chunks? SDs carry high SNV density and interlocus gene conversion [11,12], and the two regions we test, 17q21.31 and 22q11.2, are known for structural diversity and repeat-mediated rearrangement [13-17].

We answer these questions with controlled comparisons on real assemblies. Every chunked graph is compared with a monolithic graph built from the same input with the same parameters, and, for four regions, both are scored against variant calls taken directly from HPRC's assembly-to-reference alignments. The tests ran in tiers. Each tier was designed after the previous one corrected or narrowed its claims, and the Results follow that order. Along the way we found a parameter-parsing bug in the pinned PGGB image that affects any input sequence longer than about 1 Mb.

In the sources searched on 2026-11-24 (Section 6 lists them), we did not find a controlled accuracy comparison of chunked against monolithic PGGB-style graphs, a test of hotspot against fixed or random cut positions, a measured CPU and memory break-even for chunking, a description of the SD failure mode, or a report of the seqwish `-B` problem. These are statements about the sources searched, not claims that no such work exists.

## 2. Methods

### 2.1 Data

Assemblies came from HPRC release 2 [6]. The pilot and Tier 1 used both haplotypes of HG00673 and HG00733, which were also the two samples in the duplication-region and truth-set experiments. For each haplotype, the contig holding the region was identified from the release's Minigraph-Cactus GRCh38 chains and extracted with remote `samtools faidx` [26], so no whole assembly was downloaded. GRCh38 sequence came from Ensembl release 112. Records follow the PanSN naming convention (sample#haplotype#contig) [18].

For the scaling experiment, 18 further samples were drawn at random (seed 1) from the 229 release 2 samples whose chr21:19.95-21.05 Mb lies on a single Minigraph-Cactus chain in both haplotypes. The sets were nested: the first 5, 10 or 20 samples gave N = 10, 20 or 40 haplotypes. GRCh38 was added to every build, so a build at N haplotypes has N + 1 paths. Tier 1 used four haplotypes plus GRCh38 (five paths).

### 2.2 Regions

| Region | GRCh38 coordinates | Used in | Notes |
| --- | --- | --- | --- |
| chr21 pilot | chr21:20,000,000-21,000,000 | pilot, Tiers 1-3 | single-copy region; haplotypes 1,001,586 to 1,001,742 bp |
| chr21 10 Mb | chr21:20,000,000-30,000,000 | Tier 3 (E4) | all four haplotypes on one contig across the region |
| 17q21.31 | chr17:45.4-46.6 Mb | Tiers 2-3 | inversion region with duplicated flanks; haplotypes 1.20 to 1.49 Mb |
| 22q11.2 | chr22:18.6-21.5 Mb | Tiers 2-3 | LCR22A to D; includes 50 kb of N; haplotypes 2.78 to 4.70 Mb, almost all of the excess in the first 414 kb window, which holds LCR22A (0.34 to 2.24 Mb per haplotype) |

Table 1. Regions. Haplotype lengths are the projected per-haplotype sequences for the four HG00673 and HG00733 haplotypes.

### 2.3 Projecting cuts onto haplotypes

Each haplotype contig was aligned to the reference region with minimap2 2.28 (`-x asm5 -c --eqx`) [27]. In Tier 3 the same alignments were made with mappy 2.31, minimap2's Python binding, behind a small command-line shim. `project_chunks.py` keeps primary alignments with MAPQ of at least 20, groups them by contig and strand, and takes the colinear chain with the most aligned reference bases, found by an O(n²) dynamic programme over compatible predecessors. The pieces of the chain are merged into one alignment, with the stretch between two pieces recorded as a reference gap plus a haplotype gap, so that every reference coordinate projects to exactly one haplotype coordinate. Each chunk boundary is projected through this alignment, and adjacent chunks therefore meet exactly on every haplotype.

The dynamic programme replaced a greedy left-to-right chain during Tier 2. In 22q11.2 the greedy chain locked onto pieces that mapped to the wrong LCR copy and covered 0.2 to 2% of the region. The chained projection covers 0.89 to 0.96, no cut falls in a chain gap, and all 40 chr21 projections are byte-identical to those made with the greedy chain. Every zero-overlap chunking was checked to concatenate back to the monolithic haplotype sequences base for base.

### 2.4 Chunkings

| Name | Region | Cuts | Overlap | Purpose |
| --- | --- | --- | ---: | --- |
| windows | chr21 1 Mb | 20.4, 20.8 Mb (chunks 425, 450, 225 kb) | 50 kb | the original pipeline's fixed windows |
| windows0 | chr21 1 Mb | 20.333, 20.667 Mb | 0 | fixed windows without overlap |
| hotspots300 | chr21 1 Mb | 20.348, 20.650 Mb | 0 | haploblocks.org cuts, blocks merged to at least 300 kb |
| hotspots50 | chr21 1 Mb | 9 cuts (10 chunks, 51 to 271 kb) | 0 | haploblocks.org cuts, minimum 50 kb |
| shift50, shift100 | chr21 1 Mb | hotspot cuts + 50 kb, + 100 kb | 0 | same chunk sizes, other places |
| random_01..20 | chr21 1 Mb | 2 random cuts each, chunks at least 200 kb, seed 1 | 0 | null distribution |
| windows0 | 17q21.31; 22q11.2 | 3 x 400 kb; 7 x 414 kb | 0 | fixed windows |
| hotspots300 | 17q21.31; 22q11.2 | haploblocks.org cuts, 300 kb minimum (3; 8 chunks) | 0 | hotspot cuts |
| sdgroup | 17q21.31; 22q11.2 | windows0 chunks, those sharing a duplication built as one job | 0 | duplication-aware grouping (Section 2.12) |
| 400 kb, 1 Mb, 2.5 Mb | chr21 10 Mb | 25, 10 and 4 equal windows | 0 | chunk size at 10 Mb |

Table 2. Chunkings. Coordinates are GRCh38.

The hotspot cuts are the 19 haploblocks.org boundaries inside chr21:20-21 Mb, which split 1000 Genomes haplotypes into recombination-bounded blocks [19]. `hotspot_chunks.py` merges adjacent blocks shorter than a minimum size. The method behind haploblocks.org is not published, which is why we checked its boundaries against independent recombination maps (Section 2.10).

### 2.5 Graph construction

All graphs were built with PGGB v0.6.0 [1] from the pinned image `ghcr.io/pangenome/pggb@sha256:44f3563c5c3e2032c293fe309eed178b754021b07963ca4857e4ae5efda86474`, which bundles seqwish v0.7.11-35-gaab5bee. Parameters were `-p 90 -s 5000 -K 19 -k 29` with 4 threads per build. After the bug described in Section 3.3 was found, a numeric seqwish batch size, `-B 100000000`, became the default in every runner; all duplication-region builds reported here and all Tier 3 builds used it. The pilot, Tier 1 and the scaling runs used PGGB's default `-B`; we checked that the 1 Mb chr21 builds are unaffected (Section 3.3). A control with PGGB's recommended settings for human data, `-p 98 -s 10000 -k 47` (E5), was running at the time of writing.

The monolithic baseline is one PGGB run over the full multi-FASTA of the region. Each chunk is a separate PGGB run on the projected sub-sequences of all haplotypes, with identical parameters.

### 2.6 Joining chunk graphs

Chunk graphs were joined with the repository's stitch, which walks each haplotype across the boundary between adjacent chunks. For overlapping chunks it cuts inside the overlap; for chunks that touch exactly, `seam_position` joins them end to end at the shared coordinate, and a gap between chunks fails the boundary rather than being bridged. A boundary passes when every haplotype crosses it, and the joined graph passes when every haplotype path spells the same sequence as in the monolithic graph. We also joined the same chunk graphs with `impg lace` (impg 0.3.3, bioconda) [10] and compared the results. In the Tier 3 duplication-aware experiment, `impg lace` was the joiner.

### 2.7 Comparison with the monolithic graph

Graph size was measured as node count and total node sequence (bp). `seam_metrics.py` places every node at the reference coordinate of the last GRCh38 node before it on any haplotype path, and splits the node sequence in excess of the monolithic graph into a part within 10 kb of a cut and the rest. Variants were called with `vg deconstruct -P GRCh38 -a` [20] on every graph and compared with `truvari bench` [21] with the monolithic calls as the base, on all sizes including SNVs, after splitting multi-allelic records, with `--refdist 500 --pctseq 0.7 --pctsize 0.7`. truvari errors were also split by distance to the nearest cut. The pilot instead matched calls exactly on (POS, REF, ALT), which is stricter than truvari, and Tier 3 also reports this exact-match F1.

This comparison measures agreement with the monolithic graph. It does not measure accuracy: where the two disagree, it cannot say which one is right.

### 2.8 Truth set

To score accuracy we built per-haplotype truth calls for HG00673 and HG00733 from HPRC's release 2 winnowmap alignments of each assembly to GRCh38 without alternate contigs. HPRC made these alignments with a different aligner from both our minimap2 projection and the wfmash alignments inside PGGB. `truth_calls.py` follows the dipcall and paftools approach [22]. A reference base is callable when exactly one alignment of at least 10 kb covers it and that alignment is at least 50 kb long with MAPQ of 5 or more. Inside callable regions, SNVs, insertions and deletions are read from the CIGAR string.

| Region | Truth calls per haplotype | Callable share |
| --- | --- | --- |
| chr21:20-21 Mb | 1,939 to 2,127 | 0.998 to 0.999 |
| 17q21.31 | 721 to 894 | 0.75 to 0.99 |
| 22q11.2 | 3,311 to 4,336 | 0.83 to 0.90 |
| chr21:20-30 Mb | 13,681 to 15,395 | 0.992 to 0.993 |

Table 3. Truth set coverage.

Graph calls were made with `vg deconstruct -a` on the GRCh38 path, then `vcfbub -l 0 -a 100000` and `vcfwave` as in the PGGB pipeline [28, 29], split into one haploid VCF per haplotype, and left-aligned together with the truth calls by `bcftools norm` [26]. truvari 5.5.0 [21] compared them inside each haplotype's callable regions, once on all sizes (`--pick multi --refdist 500 --pctseq 0.7 --pctsize 0.7`) and once on SVs of 50 bp or more. Calls were stratified by overlap with the GIAB v3.3 segmental duplication track (`GRCh38_segdups.bed`) [23]. Because the PGGB pipeline's `vcfbub -l 0 -a 100000` turned out to drop nested sites under short top-level bubbles (Section 3.7), every graph in the duplication regions was also scored with `vcfbub -l 0 -a 10000` (`run_bub_sensitivity.sh`). E5 repeated the monolithic and windows0 builds of the three 1 Mb-scale regions with PGGB's recommended human settings, `-p 98 -s 10000 -k 47` (`run_e5.sh`).

The truth set is weakest where the question is hardest. Inside duplications the assembly alignment itself can be wrong, and the callable filter drops 10 to 25% of the duplicated sequence. The comparison therefore says which graph agrees with HPRC's assembly-to-reference alignments, the usual dipcall-style benchmark. It does not establish which graph is correct base by base.

### 2.9 Path revisits

`path_revisits.py` counts, for each haplotype path, the bases it spells through a node it has already visited, and reports their share of the path length. A path can only revisit a node if the graph has merged two copies of a sequence from the same haplotype, so the share measures how much of each haplotype has been collapsed onto its own paralogs or repeats.

### 2.10 Recombination maps and block compression

To check whether the haploblocks.org cuts lie on hotspots, we used five GRCh38 maps from the stdpopsim catalogue: HapMap phase II [24], pyrho maps for YRI, CEU and CHB [25], and the deCODE sex-averaged map. In that catalogue the deCODE map is the 2010 pedigree map lifted from GRCh36 to GRCh38 by stdpopsim [30], not the 2019 sequence-level map. A hotspot was a stretch with at least 5 times the chromosome-wide mean rate (10 times in a sensitivity check). The statistic was the share of the 19 cuts within 5 kb of a hotspot, compared with 1,000 sets of 19 random positions in the region.

To test whether hotspot-bounded blocks hold more identical haplotypes, we counted, for every block of at least 5 kb and N = 10, 20 and 40, the distinct haplotype sequences by exact identity, and computed a length-weighted duplicate fraction, 1 - distinct / N. The 19 cuts make 20 blocks. The null distribution was 20 layouts with the same block sizes shuffled and rotated around the region.

### 2.11 Replicates

PGGB was not deterministic between runs of the pinned image (Section 3.1), so replicates were planned for the main comparisons.

| Experiment | Replicates |
| --- | --- |
| Pilot | 1 |
| Tier 1: monolithic, windows, windows0, hotspots300, hotspots50 | 3 |
| Tier 1: shift50, shift100, 20 random cut sets | 1 |
| Tier 2: scaling at N = 10, 20, 40 | 2 (rep1 in Tier 2, rep2 in Tier 3 as E3) |
| Tier 2: duplication regions against monolithic | 2 |
| Tier 3: truth comparison and duplication-aware groups (E1, E2) | 2 |
| Tier 3: 10 Mb (E4) | 1 |
| Tier 3: recommended settings (E5) | 1 |

Table 4. Replicate design.

Before running Tier 1 we fixed a decision rule: a placement effect would be claimed only if the hotspot chunking beat at least 19 of the 20 random cut sets. No other formal tests were run apart from the permutation tests in Section 2.10.

### 2.12 Duplication-aware grouping

`sd_groups.py` maps 20 kb tiles of the reference region back to the region (minimap2 asm20 preset through mappy) and links two chunks when hits of at least 5 kb and at least 90% identity join them with at least 10 kb in total. Linked chunks are built as one PGGB job from the concatenation of their FASTA files, the other chunks reuse the windows0 graphs, and `impg lace` joins everything (`group_build.py`). In 17q21.31 this linked chunks 1 and 3 (78 kb of linking hits); in 22q11.2 it linked chunks 1, 5 and 7 (466 kb), which hold LCR22A, B and D.

### 2.13 Compute measurement

All builds ran one at a time on a shared cloud machine with 4 vCPU and 15 GB of memory, 4 threads per build. From Tier 1 on, CPU time and peak resident memory came from GNU `time -v` inside the container; wall time for chunked builds is reported as the sum over chunks and as the longest chunk, the latter being the wall time with one worker per chunk. The container restarted during Tiers 2 and 3, and interrupted builds were rerun from scratch. Timings vary between replicates on this machine (Tier 1 monolithic wall time was 135 s in replicate 1 and 112 to 122 s in replicates 2 and 3).

## 3. Results

### 3.1 The pilot's hotspot advantage came from the window overlap; zero-overlap seams are clean

The pilot compared three chunkings of the 1 Mb chr21 region with a monolithic build, one run each (Table 5). Three chunks cut at hotspots came closer to the monolithic graph than three fixed 400 kb windows: 5 extra nodes instead of 14, and 14 bp of extra node sequence instead of 4,723 bp. But the windows overlapped by 50 kb and the hotspot chunks did not overlap at all.

| | Monolithic | Windows, 400 kb, 50 kb overlap | Hotspots, min 300 kb | Hotspots, min 50 kb |
|---|---:|---:|---:|---:|
| Chunks | 1 | 3 | 3 | 10 |
| Nodes | 10,382 | 10,396 | 10,387 | 10,406 |
| Node sequence (bp) | 1,006,549 | 1,011,272 | 1,006,563 | 1,006,553 |
| Boundaries PASS | | 2/2 | 2/2 | 9/9 |
| Haplotypes identical to monolithic | | 5/5 | 5/5 | 5/5 |
| Exact-match variant F1 vs monolithic (3,754 calls) | | 0.9975 | 0.9979 | 0.9981 |
| SVs (50 bp+) matched exactly | | 24/24 | 24/24 | 21/24 |
| PGGB wall, longest chunk (s) | 141.6 | 100.0 | 94.9 | 87.3 |
| PGGB wall, all chunks summed (s) | 141.6 | 276.9 | 274.6 | 712.8 |

Table 5. Pilot, chr21:20-21 Mb, four haplotypes plus GRCh38, one run per build. The three SVs missed by the 10-chunk build are the same 130 to 156 bp insertions reported 5 bp apart, 41 kb from the nearest cut.

Tier 1 separated overlap from placement with zero-overlap fixed windows (windows0), shifted hotspot cuts and 20 random cut sets, and repeated the main chunkings three times (Table 6). The overlap explained the pilot. The 50 kb overlap added 4,723 bp of node sequence within 10 kb of the seams in every replicate, and all 7 of its truvari false negatives sat there. Without overlap, the near-cut excess was 0 bp in every replicate of windows0, hotspots300 and hotspots50, and in 18 of the 20 random sets. The four exceptions (two random sets and both shifted chunkings) were +21 or -22 bp, within the 47 bp by which two monolithic builds differ. Every zero-overlap seam passed and every haplotype path was identical to the monolithic one.

| | windows (50 kb overlap) | windows0 | hotspots300 | hotspots50 |
| --- | --- | --- | --- | --- |
| Extra node sequence (bp) | 4,744 / 4,697 / 4,697 | -22 / -48 / -69 | 8 / -39 / -12 | 5 / 5 / -48 |
| ...of which within 10 kb of a cut | 4,723 / 4,723 / 4,723 | 0 / 0 / 0 | 0 / 0 / 0 | 0 / 0 / 0 |
| truvari F1 vs monolithic | 0.9980 / 0.9975 / 0.9975 | 0.9993 / 0.9987 / 0.9988 | 0.9985 / 0.9980 / 0.9976 | 0.9992 / 0.9992 / 0.9989 |
| truvari errors within 10 kb of a cut (FN) | 7 / 7 / 7 | 0 / 0 / 0 | 0 / 0 / 0 | 0 / 0 / 0 |
| SVs (50 bp+) missed by truvari | 0 / 0 / 0 | 0 / 0 / 0 | 0 / 0 / 0 | 0 / 0 / 0 |
| Seams PASS, haplotypes identical | 2/2, 5/5 | 2/2, 5/5 | 2/2, 5/5 | 9/9, 5/5 |
| CPU, all chunks (s) | 393 / 326 / 317 | 359 / 308 / 314 | 361 / 309 / 312 | 651 / 586 / 584 |
| Longest chunk, wall (s) | 99 / 87 / 79 | 91 / 78 / 77 | 90 / 79 / 78 | 85 / 71 / 73 |

Table 6. Tier 1, chr21:20-21 Mb, four haplotypes plus GRCh38. Values are replicates 1 / 2 / 3, each scored against its own monolithic build. Monolithic builds had 10,386 / 10,395 / 10,395 nodes, 1,006,528 / 1,006,575 / 1,006,575 bp of node sequence, 135 / 122 / 112 s wall time and 266 / 229 / 221 s CPU time.

The monolithic builds themselves varied. Two of the three were identical; the third differed by 9 nodes and 47 bp, including 3 SVs reported at shifted positions, and truvari F1 between monolithic builds was 0.9995 (replicate 1 against 2 and 3) and 1.0 (replicate 2 against 3). Differences of a few nodes or a few tens of base pairs between single runs are therefore within run-to-run noise, and the remaining differences between chunked and monolithic graphs (tens of base pairs) sit at that scale.

What remained after removing the overlap was not at the seams. It clustered at a few loci 25 to 75 kb from the nearest cut. Random cut sets with a cut at 20.20 to 20.25 Mb (sets 08, 10, 14 and 17) shifted the graph at 20.275 to 20.300 Mb by -51 to -74 nodes and +102 to +155 bp relative to the replicate 2 monolithic graph, and all five chunkings inspected this way also differed from monolithic at 20.725 to 20.750 Mb. A plausible cause is that the graph built around a complex locus depends on where its chunk starts, because wfmash maps in 5 kb segments counted from the start of each sequence; we did not test this.

`impg lace` produced the same graphs as our stitch. On the Tier 1 chunk graphs it gave the same node, edge and base counts, the same multiset of segment sequences and haplotypes identical to the monolithic graph for all three zero-overlap chunkings in all three replicates (9 of 9). In Tier 3 the two joiners gave identical graphs, with the same truth scores, on all four regions. `impg partition -w 333334 -d 0`, run on PGGB's own wfmash alignments rather than our minimap2 alignments, placed the first cut at the same haplotype coordinate as `project_chunks.py` on all four haplotypes; it merged the last two windows into one 667 kb partition, so only one cut could be compared.

### 3.2 Cut placement: hotspots are real but give no advantage

On the Tier 1 controls, the hotspot chunking beat 14 of 20 random cut sets on extra node sequence and sat at the random median on truvari F1 (0.9985 against 0.9984). The random sets ranged from -93 to +163 bp of extra node sequence (median absolute value 26 bp) and from 0.9960 to 0.9999 in F1; the shifted chunkings scored 2 bp and F1 0.9965 (shift50) and 2 bp and F1 0.9964 (shift100). The fixed zero-overlap windows had higher F1 than the hotspot chunks in all three replicates, while the hotspot chunks had slightly less extra node sequence in all three (8, -39 and -12 bp against -22, -48 and -69 bp), a difference within the 47 bp that separated monolithic replicates. By the rule fixed before the run, there was no placement effect (Figure 1).

The same held in every later comparison. With 10, 20 and 40 haplotypes, windows0 had higher truvari F1 against the monolithic graph than hotspots300, and total extra node sequence for both stayed within tens of base pairs at N = 10 and 20 and within 134 bp at N = 40 (Section 3.4, one replicate per N). In both duplication regions, fixed and hotspot cuts gave truvari F1 within 0.01 of each other, with errors at the same loci (Section 3.5). Against the truth set on chr21, monolithic, windows0 and hotspots300 graphs all scored F1 0.981 (Table 14).

The null result is not because the cuts miss hotspots. Four LD-based maps and the pedigree-based deCODE map all placed more of the 19 haploblocks.org cuts near a hotspot than random positions (Table 7). The pedigree map is the most independent check; at 10 kb resolution it finds only 5 hotspots in the region but still puts the cuts near them more often than chance. At the 10-fold threshold it has one hotspot in the region and no power (p = 0.32).

| Map | Hotspots in region | Cuts within 5 kb | Random cuts | p |
| --- | ---: | ---: | ---: | ---: |
| HapMap II | 23 | 0.74 | 0.31 | 0.001 |
| pyrho YRI | 19 | 0.68 | 0.23 | < 0.001 |
| pyrho CEU | 14 | 0.58 | 0.18 | < 0.001 |
| pyrho CHB | 16 | 0.63 | 0.21 | < 0.001 |
| deCODE 2010 (10 kb bins) | 5 | 0.32 | 0.11 | 0.016 |

Table 7. Share of the 19 haploblocks.org cuts in chr21:20-21 Mb within 5 kb of a hotspot (at least 5 times the chromosome mean rate), against the mean of 1,000 random sets of 19 positions.

The remaining argument for hotspot cuts was compression: inside a recombination-bounded block, more haplotypes should be identical, and a builder could collapse identical copies first. Exact duplicates were rare in these assemblies, and hotspot blocks held fewer of them than most random layouts of the same block sizes (Table 8). At N = 40 about 3% of block sequence could be dropped. haploblocks.org reports much stronger grouping for the same blocks (5,096 haplotypes in 3 to 1,002 MMseqs2 clusters per block, median about 120) [19], but those are clusters of sequences built from phased variant calls, not exact copies of assembled sequence. We did not test near-identity clustering.

| N | Hotspot blocks | Random layouts, median (range) | Random layouts with more duplication |
| ---: | ---: | ---: | ---: |
| 10 | 0.011 | 0.013 (0.003-0.030) | 14 of 20 |
| 20 | 0.010 | 0.018 (0.011-0.037) | 20 of 20 |
| 40 | 0.029 | 0.038 (0.026-0.054) | 18 of 20 |

Table 8. Length-weighted duplicate fraction (1 - distinct / N, exact identity) over blocks of at least 5 kb.

### 3.3 PGGB v0.6.0 passes an unparsed seqwish batch size

In the first pass on the duplication regions, monolithic graphs came out far larger than chunked ones: 2.26 Mb against 1.49 Mb of node sequence for 17q21.31, and 12.3 Mb against 6.4 Mb for 22q11.2. The wfmash alignments were complete (in 22q11.2 every haplotype aligned end to end to the reference at 99.8% identity), so the loss happened during induction. Rerunning seqwish from the same image on the saved alignments isolated it (Table 9).

| seqwish `-B` | 17q21.31 graph (bp) | 22q11.2 graph (bp) |
| --- | ---: | ---: |
| `10M` (what pggb passes) | 2,281,318 | 12,455,138 |
| `1000000` | 2,281,318 | |
| `2000000` | 1,849,589 | |
| `10000000` | 1,379,791 | |
| `30000000` | 1,379,791 | 4,366,871 |

Table 9. Graph length from seqwish v0.7.11-35-gaab5bee on the same alignments, by batch size.

`10M` behaves exactly like `1000000`, and the log shows 1 Mb batches. A single 2.8 Mb reference-to-haplotype alignment, induced alone, merged only reference positions 109,564 to 1,000,000; everything after the first batch edge stayed unmerged. In the image's `pggb` script, a user-supplied `-B` passes through `parse_numeric`, but the default `TRANSCLOSE_BATCH=10M` does not, and seqwish reads the string as 1,000,000. Its batched transitive closure [4] then loses merges across batch edges. The DNAnexus chunk applet in our own repository also called seqwish with `-B 10M` directly. We did not find an existing report of this problem in the sources searched.

The 1 Mb chr21 results are not affected. Rebuilt with `-B 100000000`, the chr21 monolithic graph had 10,385 nodes and 1,006,527 bp at four haplotypes plus GRCh38 (Tier 1 replicates: 10,386 to 10,395 nodes and 1,006,528 to 1,006,575 bp) and 22,830 nodes and 1,013,665 bp at N = 40 (default `-B`: 22,898 and 1,013,701), both within run-to-run noise. All duplication-region and 10 Mb numbers below come from builds with the numeric `-B`; the first pass was discarded. The repository configuration now sets `transclose_batch: 100000000` and every runner passes it.

### 3.4 Scaling with haplotype count and region length

Monolithic CPU time grew faster than linearly with haplotype count, from 364 s at N = 10 to 708 s at N = 20 and 1,775 s at N = 40, and peak memory rose from about 970 MB to 2.3 GB (Table 10, Figure 3). Below N = 20, every build, monolithic or chunk, sat at a memory floor near 970 MB, so at that scale chunking cannot save memory.

| N haplotypes (+ GRCh38) | Nodes | Node sequence (bp) | Wall (s) | CPU (s) | Peak memory (MB) |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 4 (Tier 1, 3 reps) | 10,386-10,395 | 1,006,528-1,006,575 | 112-135 | 221-266 | about 970 |
| 10 | 13,362 | 1,008,347 | 163 | 364 | 972 |
| 20 | 18,214 | 1,011,283 | 250 | 708 | 1,197 |
| 40 | 22,898 | 1,013,701 | 536 | 1,775 | 2,335 |
| 10, rep2 | 13,442 | 1,008,293 | 153 | 365 | 972 |
| 20, rep2 | 18,218 | 1,011,246 | 249 | 712 | 1,181 |
| 40, rep2 | 23,120 | 1,013,737 | 525 | 1,745 | 2,210 |

Table 10. Monolithic builds of chr21:20-21 Mb. Rows without a replicate label are rep1 (Tier 2); rep2 is E3, with the same inputs and settings.

With three zero-overlap chunks the seams held at every N: extra node sequence within 10 kb of a cut was between -2 and +1 bp, and every haplotype path matched the monolithic graph (Table 11). The truvari errors counted as near a cut were not at the cut. At N = 20 the 6 to 7 near-cut errors of windows0 and hotspots300 are one locus at about 20.657 Mb, 9 kb before one cut and 7 kb after the other. Ten small chunks were worse: they added 30 to 101 bp near cuts, lost 1 to 5 SVs, and 46 of their 56 near-cut errors across the three values of N sat 3 to 5 kb after the cut at 20,583,159. That distance is close to the 5 kb wfmash segment length, which suggests that the first or last segment of a chunk maps differently; we did not test this.

| | N = 10 | N = 20 | N = 40 |
| --- | --- | --- | --- |
| Extra node sequence within 10 kb of a cut (bp) | 0 / 0 / 101 | 0 / -2 / 50 | 0 / 1 / 30 |
| Extra node sequence, total (bp) | 1 / 20 / -44 | 7 / -24 / -12 | -134 / 71 / -47 |
| Haplotypes identical to monolithic | all | all | all |
| truvari F1 vs monolithic | 0.9982 / 0.9977 / 0.9951 | 0.9978 / 0.9930 / 0.9929 | 0.9961 / 0.9930 / 0.9905 |
| truvari FN within 10 kb of a cut | 0 / 0 / 3 | 7 / 6 / 21 | 0 / 0 / 32 |
| SVs (50 bp+) missed, of 40 / 67 / 112 | 0 / 0 / 3 | 0 / 1 / 1 | 0 / 2 / 5 |
| CPU, all chunks (s) | 456 / 469 / 663 | 798 / 780 / 1,016 | 1,790 / 1,742 / 1,771 |
| Longest chunk, wall (s) | 93 / 96 / 79 | 123 / 123 / 104 | 215 / 210 / 163 |
| Peak memory, largest chunk (MB) | 969 / 968 / 968 | 970 / 969 / 968 | 970 / 970 / 853 |

Table 11. Chunked builds against the monolithic build of the same N (windows0 / hotspots300 / hotspots50), rep1.

| | N = 10 | N = 20 | N = 40 |
| --- | --- | --- | --- |
| Extra node sequence within 10 kb of a cut (bp) | 0 / 0 / 100 | 4 / -2 / 7 | 0 / 0 / -35 |
| Extra node sequence, total (bp) | 54 / 21 / 132 | 7 / 23 / -79 | -15 / -81 / 144 |
| Haplotypes identical to monolithic | all | all | all |
| truvari F1 vs monolithic | 0.9956 / 0.9962 / 0.9963 | 0.9982 / 0.9931 / 0.9923 | 0.9934 / 0.9861 / 0.9888 |
| truvari FN within 10 kb of a cut | 0 / 0 / 2 | 3 / 6 / 22 | 0 / 0 / 36 |
| SVs (50 bp+) missed | 0 / 2 / 0 | 0 / 1 / 0 | 4 / 7 / 3 |
| CPU, all chunks (s) | 431 / 435 / 660 | 793 / 810 / 976 | 1,779 / 1,733 / 1,772 |
| Longest chunk, wall (s) | 86 / 86 / 79 | 122 / 123 / 103 | 211 / 222 / 155 |
| Peak memory, largest chunk (MB) | 969 / 968 / 968 | 970 / 969 / 968 | 971 / 970 / 777 |

Table 11b. The same comparison for rep2 (E3), each chunked graph against the rep2 monolithic build.

Agreement with the monolithic graph drifted down with N. For three chunks, F1 fell from 0.998 at N = 10 to 0.993-0.996 at N = 40. Part of this drift is run-to-run noise: at four haplotypes, monolithic builds already agreed with each other at F1 0.9995, and the second replicate shows that the noise grows with N too. Rep2 gave 0.993-0.996 at N = 10 and 0.986-0.993 at N = 40, and the two monolithic builds at each N differed by 80, 4 and 222 nodes at N = 10, 20 and 40. The decline with N is therefore present in both replicates, but its size varies between them by up to 0.007, which is as large as the differences between chunkings at a given N. What does repeat is the order: windows0 scored higher than both hotspot layouts in 5 of the 6 replicate-by-N comparisons (the exception, N = 10 rep2, is 0.9956 against 0.9962 and 0.9963), and the three-chunk layouts kept the extra sequence within 10 kb of a cut between -2 and +4 bp in both replicates.

The compute balance shifted with N. At four haplotypes, three chunks used 1.3 to 1.4 times the monolithic CPU time and ten chunks 2.4 to 2.6 times, because each PGGB run carries a fixed cost of roughly a minute at this size. At N = 40, three chunks used the same CPU as monolithic (ratio 0.98 to 1.01; windows0 1,790 s against 1,775 s), the longest chunk finished in 40% of the monolithic wall time, and peak memory per build stayed near 970 MB while the monolithic build reached 2.3 GB.

At 10 Mb (chr21:20-30 Mb, four haplotypes plus GRCh38, one replicate), monolithic CPU time was 6.2 times that of the 1 Mb build (1,681 s against 271 s), less than the 10-fold increase in length because fixed costs dominate the 1 Mb build (Table 12). Each chunk also pays its own fixed cost, so 25 chunks of 400 kb used 1.51 times the monolithic CPU, while chunks of 1 Mb and 2.5 Mb used 1.22 and 1.06 times. With one worker per chunk, the longest build took 75 to 189 s against 582 s, and no build needed more than 979 MB, against 2,840 MB for the monolithic build. Seams passed and haplotypes were identical in all three chunkings, and truth F1 was 0.979 for every graph.

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

Table 12. chr21:20-30 Mb, one replicate, numeric seqwish `-B`.

### 3.5 Segmental duplications

**Chunked and monolithic graphs disagree.** In 17q21.31 and 22q11.2 the seams were clean: extra node sequence near cuts was 0 to 1 bp, except at one hotspot cut in 22q11.2 (-543 bp in both replicates, where the scoring places a large LCR insertion next to the cut). But the chunked graphs carried about 6.5% more node sequence than the monolithic graphs, and truvari F1 against the monolithic calls was 0.781 to 0.795 in 17q21.31 and 0.872 to 0.879 in 22q11.2 (Table 13). This is beyond run-to-run noise, which is also larger in these regions: two monolithic builds of the same input agreed at F1 0.983 (17q21.31) and 0.973 (22q11.2), against 0.9995 on chr21. At least 92% of the disagreeing calls fell in a few 100 kb bins that cover 17% of the 17q21.31 region and 28% of 22q11.2: chr17:45.5-45.6 and 46.5-46.6 Mb (the duplicated segments at both ends of the inversion region) and chr22:18.6-19.0, 20.3-20.4 and 21.1-21.4 Mb, which line up with LCR22A, B and D. Chunks of 300 to 400 kb separate these copies, so a chunked build cannot align them to each other, while the monolithic build can. This explanation comes from where the errors fall, not from base-level inspection of the graphs.

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

Table 13. Duplication regions against the monolithic graph (Tier 2), HG00673 and HG00733 plus GRCh38, numeric seqwish `-B`, replicates 1 / 2.

At 1.2 Mb and four haplotypes, chunking saved no memory: the 17q21.31 monolithic build peaked at about 440 MB, while every chunk sat at the 970 MB floor.

**Against the truth set, the chunked graph is closer in 22q11.2 and mixed in 17q21.31.** The comparison with the monolithic graph cannot say which representation is right, so Tier 3 scored both against the assembly-based truth calls (Table 14, Figure 4). Table 14 gives replicate 1; a full second replicate agreed within 0.015 on every F1 value and within 0.006 on most (Table 14b). Outside duplications, every graph scored F1 0.973 to 0.981, so the differences were all inside the duplicated blocks.

| Region, graph | Precision | Recall | F1 | F1 in SDs | F1 elsewhere | SV F1 (n) |
| --- | --- | --- | --- | --- | --- | --- |
| chr21 1 Mb, monolithic | 0.977 | 0.985 | 0.981 | (no SDs) | 0.981 | 0.943 (27) |
| chr21 1 Mb, windows0 | 0.977 | 0.985 | 0.981 | | 0.981 | 0.943 |
| chr21 1 Mb, hotspots300 | 0.977 | 0.985 | 0.981 | | 0.981 | 0.943 |
| 17q21.31, monolithic | 0.632 | 0.903 | 0.743 | 0.627 | 0.973 | 0.900 (9) |
| 17q21.31, windows0 | 0.726 | 0.902 | 0.804 | 0.720 | 0.973 | 0.720 |
| 22q11.2, monolithic | 0.909 | 0.853 | 0.880 | 0.719 | 0.973 | 0.761 (103) |
| 22q11.2, windows0 | 0.964 | 0.960 | 0.962 | 0.937 | 0.978 | 0.818 |

Table 14. Accuracy against truth calls from HPRC winnowmap assembly alignments, HG00673 and HG00733, `-p 90 -s 5000 -K 19 -k 29`, replicate 1. n is the number of truth SVs (50 bp or more).

| Region, graph | F1 all | F1 in SDs | SV F1 |
| --- | --- | --- | --- |
| chr21 1 Mb, monolithic / windows0 / hotspots300 | 0.981 / 0.981 / 0.981 | | 0.943 / 0.943 / 0.943 |
| 17q21.31, monolithic | 0.749 | 0.637 | 0.900 |
| 17q21.31, windows0 | 0.802 | 0.717 | 0.720 |
| 17q21.31, sdgroup | 0.750 | 0.636 | 0.706 |
| 22q11.2, monolithic | 0.880 | 0.713 | 0.755 |
| 22q11.2, windows0 | 0.962 | 0.936 | 0.814 |
| 22q11.2, sdgroup | 0.889 | 0.760 | 0.774 |

Table 14b. Replicate 2 of the truth comparison: every graph rebuilt from the same inputs. The two monolithic 22q11.2 builds differ by 183 nodes and 167 bp of node sequence, yet their recall inside SDs is 0.660 and 0.651.

In 22q11.2, inside the LCR22s, the monolithic graph missed a third of the truth calls (recall 0.660) and the chunked graph 7% (0.933), and the chunked graph also had higher precision there (0.941 against 0.790). About half of that missing recall is a calling artefact described in Section 3.7: with a 10 kb vcfbub bubble cut, the monolithic graph's recall inside the LCR22s is 0.821 and its F1 0.903, against 0.960 for the chunked graph. Across the whole region, SV recall rose from 0.68 to 0.89 and SV precision fell from 0.86 to 0.75. In 17q21.31 both graphs had low precision in the duplicated flanks (0.50 and 0.63) at the same recall (0.85). Chunking raised small-variant precision there but added SV false positives: 7 against 2, where the truth set holds only 5 SVs in duplications. SV F1 for the region fell from 0.90 to 0.72. With 9 truth SVs in the whole region, the 17q21.31 SV result rests on very few calls.

So under the settings used throughout (`-p 90`), where chunked and monolithic graphs disagreed, the chunked graph was at least as close to the assembly-based truth calls overall, and in 22q11.2 much closer. Section 3.7 tests whether the advantage remains under PGGB's recommended human settings (`-p 98 -s 10000 -k 47`); it does.

**Mechanism: the monolithic graph collapses paralogs.** Haplotype paths in the monolithic graphs revisited nodes they had already used for 10.5% of their length in 17q21.31 and 9.0% in 22q11.2, against 4.4% and 2.4% in the chunked graphs (Table 15). On chr21 neither graph had any revisits. A revisit means the graph merged two copies of a sequence from the same haplotype. With `-p 90`, wfmash maps segments between paralogous copies that are 90% or more identical, seqwish merges them, and a haplotype's variants in one copy are then expressed against a merged node rather than its own GRCh38 locus, which matches the loss of recall in the 22q11.2 monolithic graph. Chunks that separate the copies cannot merge them. The revisit shares were the same to three decimals in the second replicate. At `-p 98` they fall but stay well above the chunked graphs (Section 3.7).

| Region | Monolithic | windows0 | sdgroup |
| --- | --- | --- | --- |
| chr21 1 Mb | 0.000 | 0.000 | |
| 17q21.31 | 0.105 | 0.044 | 0.138 |
| 22q11.2 | 0.090 | 0.024 | 0.115 |

Table 15. Share of path length (all paths, GRCh38 included) that revisits a node already used by the same path. Replicates 1 and 2 agree to three decimals.

### 3.6 Duplication-aware grouping does not help

Tier 2 had proposed keeping paralogous copies in one chunk. Tier 3 tested this by building the chunks linked by a duplication as one PGGB job (Table 16). Grouping did what it was designed to do: it brought the chunked graph closer to the monolithic one (exact-match F1 against monolithic from 0.646 to 0.794 in 17q21.31 and from 0.660 to 0.719 in 22q11.2). That made it worse against the truth set. In 22q11.2, truth F1 fell from 0.962 to 0.889 overall and from 0.937 to 0.760 inside SDs; in 17q21.31 it fell from 0.804 to 0.750. The grouped job holds the paralogs and little of the unique sequence between them, and it collapsed them more than the monolithic build did: 13.8% and 11.5% of haplotype length revisited nodes (Table 15), and node sequence fell below that of the monolithic graph in both regions. In 22q11.2 the group of LCR22A, B and D took as long as the whole monolithic build (398 s against 383 s wall time).

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

Table 16. Duplication-aware groups (Tier 3, E2), replicate 1. In replicate 2 the grouped 22q11.2 job again took longer than the whole monolithic build (307 s against 291 s wall), and truth F1 was 0.889 (22q11.2) and 0.750 (17q21.31). These builds are separate from the Tier 2 builds in Table 13, so node-sequence differences vary slightly between the two tables.

### 3.7 Recommended human settings, and a calling artefact

E5 rebuilt the monolithic and windows0 graphs of chr21:20-21 Mb, 17q21.31 and 22q11.2 with PGGB's recommended human settings, `-p 98 -s 10000 -K 19 -k 47` (one replicate), and scored them against the same truth calls. On chr21 nothing changed: both graphs scored F1 0.981 and SV F1 0.943. In the duplication regions the `-p 98` monolithic graphs first scored F1 0.468 (17q21.31) and 0.239 (22q11.2), with recall outside duplications of 0.023 and 0.003. The graphs were not broken; wfmash had aligned every haplotype end to end. The loss happened while making the VCF. In both monolithic graphs the GRCh38 path itself revisits nodes (11.2% and 9.1% of its length), and `vg deconstruct` reports a top-level bubble with a short reference allele (20 kb at chr17:46.49 Mb) whose 990 nested sites span about 700 kb. `vcfbub -l 0 -a 100000`, as in the PGGB pipeline, keeps top-level sites and pops a bubble only when an allele exceeds 100 kb, so the nested sites were dropped; the 17q21.31 monolithic VCF had 880 records against 1,800 for windows0.

We rescored every graph in the two duplication regions with `vcfbub -l 0 -a 10000`, which pops such parents (Table 17).

| Graph | 17q21.31 F1 (in SDs) | 17q21.31 SV F1 | 22q11.2 F1 (in SDs) | 22q11.2 SV F1 |
| --- | --- | --- | --- | --- |
| `-p 90` monolithic | 0.743 / 0.743 (0.627 / 0.628) | 0.90 / 0.90 | 0.880 / 0.903 (0.719 / 0.796) | 0.76 / 0.79 |
| `-p 90` windows0 | 0.804 / 0.804 (0.720 / 0.720) | 0.72 / 0.72 | 0.962 / 0.960 (0.937 / 0.932) | 0.82 / 0.85 |
| `-p 90` sdgroup | 0.750 / 0.750 (0.637 / 0.637) | 0.71 / 0.71 | 0.889 / 0.891 (0.760 / 0.766) | 0.79 / 0.78 |
| `-p 98` monolithic | 0.468 / 0.862 (0.664 / 0.794) | 0.36 / 1.00 | 0.239 / 0.929 (0.505 / 0.855) | 0.07 / 0.78 |
| `-p 98` windows0 | 0.941 / 0.941 (0.917 / 0.917) | 1.00 / 1.00 | 0.962 / 0.965 (0.937 / 0.943) | 0.80 / 0.85 |

Table 17. Truth F1 with `vcfbub -a 100000` / `-a 10000`, replicate 1. Path revisit shares at `-p 98`: monolithic 0.064 (17q21.31) and 0.077 (22q11.2), windows0 0.026 and 0.023.

Three things follow. Chunked graphs hardly depend on the bubble cut (overall F1 within 0.003), while monolithic graphs of duplicated regions can depend on it entirely. With the artefact removed, `-p 98` improves the monolithic graph over `-p 90` (0.862 against 0.743 in 17q21.31; 0.929 against 0.903 in 22q11.2) but leaves it below windows0 at `-p 98` (0.941 and 0.965). And the chunked `-p 98` graph of 17q21.31, which found all 9 truth SVs with no false positives, is the most accurate graph of that region we built. At `-p 98` the 22q11.2 chunks used 1.15 times the monolithic CPU time (928 s against 806 s), and the longest chunk took 120 s against 276 s.

### Figures

**Figure 1: Seams and cut placement on chr21:20-21 Mb.** (a) Extra node sequence relative to the monolithic graph for windows, windows0, hotspots300 and hotspots50, one bar per replicate, split into the part within 10 kb of a cut and the rest. The 4,723 bp near-cut bar of the overlapping windows sets the scale; the zero-overlap bars are near zero. Fields: `tier1_controls.json`, fields `replicates.rep{1,2,3}.{windows,windows0,hotspots300,hotspots50}.excess_bp_within_10kb_of_cut` and `.excess_bp_elsewhere`. (b) truvari F1 against monolithic for the 20 random cut sets as a strip or histogram, with hotspots300 (rep1), windows0 (rep1), shift50 and shift100 marked and the monolithic-against-monolithic F1 drawn as a reference line. Fields: `replicates.rep1.random_01..random_20.truvari_f1`, `replicates.rep1.{hotspots300,windows0,shift50,shift100}.truvari_f1`, `random_cuts.truvari_f1_median`, `monolithic_vs_monolithic_truvari_f1`. (c) Same layout for total extra node sequence: `replicates.rep1.random_XX.excess_bp`. Caption should state n = 20 random sets, one run each, and three replicates for the four main chunkings.

**Figure 2: Are the cuts on hotspots, and do hotspot blocks compress?** (a) For each of the five maps at the 5-fold threshold, the share of haploblocks.org cuts within 5 kb of a hotspot against the random-set mean, with p values annotated. Fields: `tier2.json`, fields `map_robustness.results.{HapMapII,PyrhoYRI,PyrhoCEU,PyrhoCHB,DeCodeSexAveraged}_x5.haploblock_cuts_within_5kb`, `.random_cuts_within_5kb_mean`, `.p_value`; the `_x10` entries as a supplementary panel. (b) Duplicate fraction against N for hotspot blocks, with the random-layout median and range as a band. Fields: `block_compression.summary.{10,20,40}.hotspot_blocks`, `.random_median`, `.random_min`, `.random_max`.

**Figure 3: Cost of monolithic and chunked builds.** (a) CPU time against N haplotypes (10, 20, 40, plus the Tier 1 point at 4) for monolithic, windows0, hotspots300 and hotspots50. (b) Wall time: monolithic against longest chunk. (c) Peak memory: monolithic against largest chunk, with the 970 MB floor visible. Fields: `tier2.json`, fields `scaling_rep1.N{10,20,40}/rep1.mono.{cpu,wall,rss}` and `.{windows0,hotspots300,hotspots50}.{cpu_sum,max_wall,max_rss}`; Tier 1 point from `tier1_controls.json`, fields `replicates.rep*.monolithic.{cpu_s,wall_s,peak_rss_mb}`. Second replicate points: `tier3.json`, fields `runs.scaling/N{10,20,40}/rep2.{monolithic,windows0,hotspots300,hotspots50}.{cpu_s,max_wall_s,peak_rss_mb}` (monolithic uses `wall_s`). (d) Region length: CPU ratio, longest build and peak memory for the 10 Mb chunkings. Fields: `tier3.json`, fields `runs.r21x/rep1.{monolithic,windows0,windows1m,windows2500k}.{cpu_s,max_wall_s,peak_rss_mb}` (monolithic uses `wall_s`). Note that the JSON key `windows0` in `r21x` is the 25 x 400 kb chunking.

**Figure 4: Segmental duplications against the truth set.** (a) Truth F1 inside SDs and elsewhere for monolithic, windows0 and sdgroup graphs in 17q21.31 and 22q11.2, with chr21 as a control. Fields: `tier3.json`, fields `runs.{r17,r22,r21}/rep1.{monolithic,windows0,sdgroup}.truth.all.sd.f1`, `.truth.all.other.f1`; precision and recall from the same nodes for a supplementary panel. (b) SV F1 with the number of truth SVs annotated: `.truth.sv.total.f1`, `.truth.sv.total.tp_base + .truth.sv.total.fn`. (c) Path revisit share per graph from Table 15. These values are not in `tier3.json`; they come from the `path_revisits.py` output behind TIER3.md and need to be exported before plotting. Plot replicate 1 with replicate 2 as a second point (`runs.{r17,r22,r21}/rep2`, same fields), and add panel (d): truth F1 at `-p 90` and `-p 98` for monolithic and windows0 under both vcfbub cuts (`runs.{r17,r22}/rep1.*.truth`, `.truth_a10k`; `runs.p98/{r17,r22}/rep1.*.truth`, `.truth_a10k`).

## 4. Discussion

For builders who want to chunk a PGGB run, the results support a short set of rules. Cut without overlap, and project each cut onto every haplotype through one colinear chain so that neighbouring chunks meet exactly. Overlap was the only source of seam artefacts we found; without it, seams added at most a few base pairs at every haplotype count and region tested, and haplotype paths were identical to the monolithic ones. Use chunks of 1 Mb or more when the region allows it. Below a few hundred kilobases each PGGB run is mostly fixed cost, and at 10 Mb the CPU overhead fell from 1.51 times with 400 kb chunks to 1.22 and 1.06 times with 1 Mb and 2.5 Mb chunks. Replicate every comparison, because the pinned image changed by up to 9 nodes and 47 bp between runs on chr21 and agreed with itself at only F1 0.973 to 0.983 in the duplication regions. Pass seqwish a numeric `-B` larger than the total input length whenever any input sequence, including a chunk, exceeds about 1 Mb. And do not spend effort on hotspot cuts: the haploblocks.org cuts do lie on hotspots in five maps, but across two regions, three haplotype counts, random and shifted controls and the truth comparison, they never gave higher variant agreement or accuracy than fixed zero-overlap cuts, and their differences in node sequence were within run-to-run noise.

When chunking pays depends on scale. At four haplotypes on 1 Mb it cost 1.3 to 1.4 times the CPU and saved no memory, because every build sat at a floor near 970 MB. At 40 haplotypes on 1 Mb, CPU was equal, wall time with parallel chunks fell to 40% and memory stayed at the floor while the monolithic build reached 2.3 GB; at 10 Mb with four haplotypes the monolithic build reached 2,840 MB. impg's documentation says that partitioning does not cap peak memory [10]. In our runs the per-chunk peak did not grow with N or region length, but the floor itself was higher than the monolithic peak for the 1.2 Mb 17q21.31 region (about 440 MB), so a memory saving appears only once the monolithic build outgrows the floor.

The stitch itself is not new. `impg lace` produced the same graphs as our stitch in every comparison, and `impg partition` placed a comparable cut at the same haplotype coordinate. What this work adds is the evaluation: which cuts to make, and how a chunked graph compares with a monolithic one and with assembly-based truth. nf-core/pangenome takes a different route [9]. It distributes the alignment and keeps one induction step, so in principle it can still merge paralogs that a chunked build separates. Its accuracy treatment could not be read for this draft.

The segmental duplication result changes how chunk boundaries should be viewed. Tier 2 read the 6.5% extra sequence of chunked graphs as a defect and proposed keeping paralogous copies in one chunk. Scored against truth calls, the chunked graph was at least as accurate in both regions overall and much more accurate in 22q11.2, where the monolithic graph routed about 9% of each haplotype through nodes it had already used. Grouping paralogs into one job reproduced and worsened that collapse. A chunk boundary that separates paralogous copies acts as a guard against merging them. This held at PGGB's recommended human settings too: `-p 98` lowered the monolithic revisit share (to 6.4% and 7.7%) but did not bring it near the chunked graphs (2.6% and 2.3%), and chunked graphs at `-p 98` were the most accurate graphs we built (F1 0.941 in 17q21.31 with all 9 truth SVs found, 0.965 in 22q11.2). So the practical advice is both: use stricter identity in duplications, and chunk. The second result is a warning for anyone who uses the VCF from a monolithic PGGB run of a duplicated region. When the reference path itself revisits nodes, `vg deconstruct` can report a short top-level bubble whose nested sites cover hundreds of kilobases, and the pipeline's `vcfbub -l 0 -a 100000` silently drops them. In our `-p 98` builds this removed most calls; at `-p 90` it removed about half of the missing recall in 22q11.2. A path-revisit count, or a second VCF with a smaller bubble cut, catches it. 17q21.31 does not fit either reading cleanly: chunking improved small-variant precision in the duplicated flanks but added SV false positives, on 5 truth SVs.

### Limitations

- Few samples. The truth comparison and duplication regions used two samples (four haplotypes). Scaling reached 40 haplotypes on one 1 Mb region of chr21. Hotspot placement was tested on chr21 and in the two duplication regions only.
- Replicates. Scaling at N = 10, 20 and 40 has two replicates; the truth comparison, revisit shares and grouping in the duplication regions have two; random and shifted controls and the 10 Mb runs ran once each.
- A weak truth set where it matters. The truth calls come from HPRC's winnowmap assembly alignments with a dipcall-style callable filter. Inside duplications those alignments can be wrong, and 10 to 25% of the duplicated sequence was not callable. The SV results in 17q21.31 rest on 5 truth SVs in duplications.
- One aligner, two parameter sets. Most graphs used PGGB at `-p 90 -s 5000 -K 19 -k 29`; the recommended human settings (`-p 98 -s 10000 -k 47`) were tested once, on monolithic and windows0 graphs of three regions. The 10 kb vcfbub cut is one alternative chosen to expose the artefact, not a validated calling pipeline, and it moves SV F1 by up to 0.05 even for chunked graphs. We did not test Minigraph-Cactus or other builders [7].
- Timings. All timings come from a shared 4 vCPU, 15 GB cloud machine, one build at a time, with variation between replicates and container restarts. Absolute times and the memory floor will differ on other hardware; the ratios are the more transferable result.
- PGGB non-determinism. Monolithic replicates differed by up to 9 nodes and 47 bp on chr21 and agreed with each other at truvari F1 0.983 and 0.973 in the duplication regions, so single-run differences of that size cannot be read as effects.
- Untested explanations. The locus-specific differences 25 to 75 kb from cuts and the errors 3 to 5 kb after a cut may come from wfmash's 5 kb segmentation, and the paralog explanation comes from where errors fall and from revisit counts, not from base-level inspection of the graphs. Compression was tested with exact identity only.
- Recombination map provenance. The cuts come from haploblocks.org, whose method is unpublished [19]; the maps used to check them were derived from 1000 Genomes, HapMap and deCODE populations, not from the HPRC samples.

## 5. Data and code availability

Code, configuration, result documents and raw numbers are at https://github.com/khoivd-dot/linear-scaling-pangenome-sv. Input data are public: HPRC release 2 assemblies and winnowmap alignments from the HPRC S3 bucket, GRCh38 from Ensembl release 112, recombination maps from the stdpopsim catalogue, and the GIAB v3.3 segmental duplication track. Scripts by experiment:

- Pilot: `pipeline/parallel/hotspot_chunks.py`, `pipeline/prepare/project_chunks.py`, `pipeline/benchmark/compare_chunking.py`; results in `docs/RECOMBINATION_CHUNKS.md` and `docs/recombination_chunks_comparison.json`.
- Tier 1: `pipeline/parallel/cut_chunks.py`, `pipeline/benchmark/run_tier1_controls.sh`, `pipeline/benchmark/seam_metrics.py`; results in `docs/TIER1_CONTROLS.md` and `docs/tier1_controls.json`.
- Tier 2: `pipeline/benchmark/tier2/` (`select_samples.py`; `scaling_fetch_align.sh`, `scaling_project_all.sh`, `scaling_run_scale.sh`; `map_robustness.py`; `block_compression.py`; `hard_prepare.sh`, `hard_fetch_align.sh`, `hard_run_fixb.sh`); results in `docs/TIER2.md` and `docs/tier2.json`.
- Tier 3: `pipeline/benchmark/tier3/` (`setup_inputs.sh`, `truth_calls.py`, `score_vs_truth.py`, `run_truth.sh`, `run_e1e2.sh`, `sd_groups.py`, `group_build.py`, `path_revisits.py`, `run_e4.sh`, `scaling_select.sh`, `scaling_rep2.sh`, `run_e5.sh`, `run_queue.sh`, `summarize.py`); results in `docs/TIER3.md` and `docs/tier3.json`. The same stages run on a GitHub runner from `.github/workflows/experiment.yml`.

Funding, competing interests, author contributions and the statement on AI assistance in writing and analysis: [to be supplied by the authors].

## 6. References

Numbering is by first citation in this draft. All entries are taken from the project's citation base (`docs/paper/related_work.md`, compiled 2026-11-24); entries that were not verified there are marked [citation to verify]. References 26-30 were added on 2026-11-25 after checking each record online (PMC, PubMed, the vcfbub repository and the stdpopsim catalogue source) and are numbered after the rest; renumber by first citation before submission.

1. Garrison E, Guarracino A, et al. 2024. Building pangenome graphs. Nature Methods 21(11):2008-2012. DOI 10.1038/s41592-024-02430-3. [full author list to verify]
2. wfmash (Guarracino A, Garrison E, et al.). Software, no dedicated paper. https://github.com/waveygang/wfmash ; Zenodo DOI 10.5281/zenodo.6949373.
3. Guarracino A, Heumos S, Nahnsen S, Prins P, Garrison E. 2022. ODGI: understanding pangenome graphs. Bioinformatics 38(13):3319-3326. DOI 10.1093/bioinformatics/btac308.
4. Garrison E, Guarracino A. 2023. Unbiased pangenome graphs. Bioinformatics 39(1):btac743. DOI 10.1093/bioinformatics/btac743.
5. Liao W-W, Asri M, Ebler J, Doerr D, Haukness M, Hickey G, et al.; Human Pangenome Reference Consortium. 2023. A draft human pangenome reference. Nature 617(7960):312-324. DOI 10.1038/s41586-023-05896-x.
6. Human Pangenome Reference Consortium. 2026. HPRC2: A human pangenome reference with near-complete coverage of common genetic variation. Preprint, DOI 10.64898/2026.07.21.739710. [citation to verify]
7. Hickey G, et al. 2023. Pangenome graph construction from genome alignments with Minigraph-Cactus. Nature Biotechnology 42(4):663-673. DOI 10.1038/s41587-023-01793-w. [author list to verify]
8. Armstrong J, Hickey G, Diekhans M, Fiddes IT, Novak AM, et al. 2020. Progressive Cactus is a multiple-genome aligner for the thousand-genome era. Nature 587(7833):246-251. DOI 10.1038/s41586-020-2871-y.
9. Heumos S, Heuer ML, Hanssen F, Heumos L, Guarracino A, Heringer P, Ehmele P, Prins P, Garrison E, Nahnsen S. 2024. Cluster-efficient pangenome graph construction with nf-core/pangenome. Bioinformatics 40(11):btae609. DOI 10.1093/bioinformatics/btae609.
10. Guarracino A, Kille B, Garrison E. impg: implicit pangenome graph. Software. https://github.com/pangenome/impg . Companion manuscript "Implicit pangenome graphs" [citation to verify].
11. Vollger MR, Guitart X, Dishuck PC, et al. 2022. Segmental duplications and their variation in a complete human genome. Science 376:eabj6965. [DOI to verify]
12. Vollger MR, Dishuck PC, Harvey WT, et al. 2023. Increased mutation and gene conversion within human segmental duplications. Nature 617:325-334. DOI 10.1038/s41586-023-05895-y. [author list to verify]
13. Zody MC, Jiang Z, Fung H-C, et al. 2008. Evolutionary toggling of the MAPT 17q21.31 inversion region. Nature Genetics 40(9):1076-1083. DOI 10.1038/ng.193.
14. Steinberg KM, Antonacci F, Sudmant PH, et al. 2012. Structural diversity and African origin of the 17q21.31 inversion polymorphism. Nature Genetics 44(8). DOI 10.1038/ng.2335.
15. Boettger LM, Handsaker RE, Zody MC, McCarroll SA. 2012. Structural haplotypes and recent evolution of the human 17q21.31 region. Nature Genetics. DOI 10.1038/ng.2334. [citation to verify]
16. Pastor S, Tran O, Jin A, Carrado D, et al. 2020. Optical mapping of the 22q11.2DS region reveals complex repeat structures and preferred locations for non-allelic homologous recombination (NAHR). Scientific Reports 10:12235. DOI 10.1038/s41598-020-69134-4.
17. Vervoort L, et al. 2025. Multiple paralogs and recombination mechanisms contribute to the high incidence of 22q11.2 deletion syndrome. Genome Research 35(4):786. DOI 10.1101/gr.279331.124.
18. PanSN-spec: sample#haplotype#contig naming. https://github.com/pangenome/PanSN-spec . [citation to verify]
19. haploblocks.org (https://data.haploblocks.org). No associated paper located. [citation to verify]
20. Garrison E, Sirén J, Novak AM, Hickey G, Eizenga JM, et al. 2018. Variation graph toolkit improves read mapping by representing genetic variation in the reference. Nature Biotechnology 36(9):875-879. DOI 10.1038/nbt.4227. [author list to verify]
21. English AC, Menon VK, Gibbs RA, Metcalf GA, Sedlazeck FJ. 2022. Truvari: refined structural variant comparison preserves allelic diversity. Genome Biology 23:271. DOI 10.1186/s13059-022-02840-6.
22. Li H, Bloom JM, Farjoun Y, Fleharty M, Gauthier L, Neale B, MacArthur D. 2018. A synthetic-diploid benchmark for accurate variant-calling evaluation. Nature Methods 15:595-597. DOI 10.1038/s41592-018-0054-7. [citation to verify: that dipcall is introduced here]
23. Dwarshuis N, Kalra D, McDaniel J, et al. 2024. The GIAB genomic stratifications resource for human reference genomes. Nature Communications 15:9029. DOI 10.1038/s41467-024-53260-y. [author list to verify]
24. International HapMap Consortium. 2007. A second generation human haplotype map of over 3.1 million SNPs. Nature 449:851-861. DOI 10.1038/nature06258.
25. Spence JP, Song YS. 2019. Inference and analysis of population-specific fine-scale recombination maps across 26 diverse human populations. Science Advances 5(10):eaaw9206.
26. Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, Whitwham A, Keane T, McCarthy SA, Davies RM, Li H. 2021. Twelve years of SAMtools and BCFtools. GigaScience 10(2):giab008. DOI 10.1093/gigascience/giab008.
27. Li H. 2018. Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics 34(18):3094-3100. DOI 10.1093/bioinformatics/bty191.
28. vcfbub: popping bubbles in vg deconstruct VCFs. Software, no dedicated paper. https://github.com/pangenome/vcfbub .
29. Garrison E, Kronenberg ZN, Dawson ET, Pedersen BS, Prins P. 2022. A spectrum of free software tools for processing the VCF variant call format: vcflib, bio-vcf, cyvcf2, hts-nim and slivar. PLOS Computational Biology 18(5):e1009123. DOI 10.1371/journal.pcbi.1009123. (vcfwave is part of vcflib.)
30. Kong A, Thorleifsson G, Gudbjartsson DF, Masson G, Sigurdsson A, et al. 2010. Fine-scale recombination rate differences between sexes, populations and individuals. Nature 467:1099-1103. DOI 10.1038/nature09525.
