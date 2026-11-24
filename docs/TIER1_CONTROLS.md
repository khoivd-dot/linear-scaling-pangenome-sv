# Hotspot chunks: controls and replicates

**Question.** The chr21 pilot ([RECOMBINATION_CHUNKS.md](RECOMBINATION_CHUNKS.md))
found that three chunks cut at recombination hotspots matched the monolithic
graph better than three fixed 400 kb windows. But the windows overlapped by
50 kb and the hotspot chunks did not overlap at all. Was the gain from cutting
at hotspots, or from dropping the overlap?

**Answer.** From dropping the overlap. With zero overlap, chunks cut at fixed
positions, at shifted hotspot positions or at 20 random positions all leave no
extra sequence at the seams, and the hotspot cuts do no better than the others.
Every bit of the extra sequence in the overlapping windows (4,723 bp, the same
in all three replicates) sits within 10 kb of the two seams. What remains after
removing the overlap is small and comes from a few loci away from the cuts. The
other finding is that PGGB is not deterministic run to run, even with the pinned
image, so differences of a few nodes or a few tens of bp between single runs
are noise.

## Setup

- **Input:** the pilot region again, GRCh38 chr21:20,000,000-21,000,000 plus
  HG00673 and HG00733, both haplotypes, from HPRC release 2. Rebuilt from
  scratch on 2026-11-24: the chr21 contig of each haplotype was found with the
  release's minigraph-cactus GRCh38 chains and pulled with remote `samtools faidx`
  (no whole-assembly download), then aligned with minimap2 2.28 (`-x asm5 -c --eqx`).
  GRCh38 comes from Ensembl release 112. The contigs and coordinates match the pilot.
- **PGGB:** pinned v0.6.0, `config/pipeline.yaml` parameters, 4 threads, one
  build at a time on a 4 vCPU, 15 GB machine. CPU time and peak memory now come
  from GNU `time -v` inside the container.
- **Chunkings** (all projected onto every haplotype with `project_chunks.py`;
  every zero-overlap chunking was checked to concatenate back to the monolithic
  haplotypes base for base):

  | Name | Cuts (Mb) | Overlap | What it tests |
  | --- | --- | ---: | --- |
  | windows | 20.4, 20.8 | 50 kb | the pilot's fixed windows |
  | windows0 | 20.333, 20.667 | 0 | fixed windows without overlap |
  | hotspots300 | 20.348, 20.650 | 0 | the pilot's hotspot chunks (min 300 kb) |
  | hotspots50 | 9 hotspot cuts | 0 | the pilot's 10 small chunks (min 50 kb) |
  | shift50, shift100 | hotspot cuts + 50 kb, + 100 kb | 0 | same sizes, wrong places |
  | random_01..20 | 2 random cuts each, chunks >= 200 kb, seed 1 | 0 | null distribution |

- **Scoring** (`pipeline/benchmark/seam_metrics.py`): every node is placed at the
  reference coordinate of the last GRCh38 node before it on any haplotype path;
  extra nodes and extra node sequence relative to the monolithic graph are then
  split into "within 10 kb of a cut" and "elsewhere". Variants from
  `vg deconstruct` are compared with `truvari bench` (all sizes, SNVs included,
  multi-allelic records split, `--refdist 500 --pctseq 0.7 --pctsize 0.7`),
  monolithic as the base, and the errors are split by distance to a cut.
- **Replicates:** monolithic plus the four main chunkings three times; the shifted
  and random chunkings once.

## Results

Three replicates, each scored against its own monolithic build. Numbers are
rep1 / rep2 / rep3.

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

Monolithic builds: 10,386 / 10,395 / 10,395 nodes, 1,006,528 / 1,006,575 /
1,006,575 bp of node sequence, 135 / 122 / 112 s wall, 266 / 229 / 221 s CPU.
truvari F1 between monolithic builds: 0.9995 (rep1 vs rep2 and vs rep3) and
1.0 (rep2 vs rep3).

Controls, one run each:

| | Extra node sequence (bp) | truvari F1 |
| --- | ---: | ---: |
| shift50 | 2 | 0.9965 |
| shift100 | 2 | 0.9964 |
| 20 random cut sets: range | -93 to +163 | 0.9960 to 0.9999 |
| 20 random cut sets: median | 26 (absolute) | 0.9984 |
| hotspots300, rep1 | 8 | 0.9985 |

What the numbers say:

- **The overlap explains the pilot.** The 50 kb overlap adds 4,723 bp at the
  seams in every replicate, and all 7 of its truvari false negatives sit there.
  Remove the overlap and the seams are clean wherever the cut falls: the
  near-cut excess is 0 bp in 18 of 20 random sets and in every fixed, shifted
  and hotspot chunking (the two exceptions are +21 and -22 bp). The end-to-end join
  added in the pilot is the useful part.
- **Hotspot placement shows no advantage at N = 5.** On extra sequence, the
  hotspot run beats 14 of 20 random sets; on truvari F1 it sits at the random
  median (0.9985 against 0.9984). The fixed zero-overlap windows scored higher
  F1 than the hotspot chunks in all three replicates. By the rule set before the
  run (beat at least 19 of 20 random sets), there is no placement effect.
- **The remaining differences are not at the seams.** They cluster at a few
  loci. For example, cuts at 20.20 to 20.25 Mb (random sets 08, 10, 14 and 17)
  shift the graph at 20.275 to 20.300 Mb by -51 to -74 nodes and +102 to
  +155 bp (against the rep2 monolithic graph), 25 to 75 kb from the nearest
  cut. All five chunkings inspected this way also differ from monolithic at
  20.725 to 20.750 Mb. A plausible cause is that the graph
  induced around a complex locus depends on where its chunk starts (wfmash maps
  in 5 kb segments counted from the sequence start), but this was not tested.
- **PGGB is not deterministic run to run.** Two of the three monolithic builds
  were identical and the third differed by 9 nodes and 47 bp, including 3 SVs
  reported at shifted positions. The pilot's statement that topology is
  deterministic for the pinned image was wrong, and single-run differences of
  that size should not be read as effects. Under truvari every SV matched in
  every run; the exact-match "21/24" counts are the same left-alignment shifts.
- **Compute.** Three zero-overlap chunks used 1.3 to 1.4 times the monolithic
  CPU time; ten chunks used 2.4 to 2.6 times. Peak memory was 966 to 971 MB in 123 of
  126 builds, whatever the chunk size; only the smallest chunk (67 kb) came in
  lower, at 796 to 873 MB. At 1 Mb and 5 haplotypes memory is a fixed floor and the memory argument for chunking cannot be tested here.

## What this changes

The case for recombination-aware cuts now rests on what was not tested here:
whether haplotypes inside a block compress (deduplicate) better as the number of
haplotypes grows. At N = 5, any zero-overlap cut works as well as a hotspot cut.
Two things follow for the next tier. Every comparison needs replicates, since
PGGB's run-to-run noise is about the size of the effects measured so far. And
the loci that move with chunk context (20.275 to 20.300 Mb and 20.725 to
20.750 Mb here) are the ones to inspect in the hard-region experiment.

## Limits

- One region, five haplotypes. Random and shifted controls ran once each.
- truvari compares graphs against the monolithic graph, not against a truth
  set.
- Timings come from a shared cloud machine; replicate 1 ran slower than 2 and 3
  (for example, monolithic 135 s against 112 to 122 s).

## Reproduce

```bash
# manifests
python3 pipeline/parallel/cut_chunks.py windows0.tsv --chrom chr21 --start 20000000 --end 21000000 \
    --cuts 20333334,20666668
python3 pipeline/parallel/cut_chunks.py random/ --chrom chr21 --start 20000000 --end 21000000 \
    --random 20 --chunks 3 --min-size 200000 --seed 1
# per haplotype: contig from the HPRC chain, remote faidx, then align
samtools faidx https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC/HG00673/assemblies/release2/HG00673_mat_hprc_r2_v1.0.1.fa.gz \
    'HG00673#2#JAHBBY020000010.1' > HG00673_mat.contig.fa
minimap2 -x asm5 -c --eqx -t 4 region.fa HG00673_mat.contig.fa > HG00673_mat.paf
python3 pipeline/prepare/project_chunks.py windows0.tsv region.fa fa/windows0 --region-start 20000000 \
    --ref-name 'GRCh38#0#chr21' --hap 'HG00673#2=HG00673_mat.paf=HG00673_mat.contig.fa'  # one --hap per haplotype
# build, compare, score (PIPELINE_CONFIG must be an absolute path if run outside the repo)
python3 pipeline/benchmark/compare_chunking.py runs/rep1 --mono fa/mono/_mono.fa \
    --strategy windows0=windows0.tsv=fa/windows0 --threads 4
python3 pipeline/benchmark/seam_metrics.py runs/rep1 --manifests windows0=windows0.tsv --truvari $(which truvari)
```

Raw numbers for every build: [`tier1_controls.json`](tier1_controls.json).
