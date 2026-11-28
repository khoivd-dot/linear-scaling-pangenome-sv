
# Parallel Pangenome Graph Construction (Linear Scaling SV)

[![Tests](https://img.shields.io/badge/tests-206%20passed-brightgreen)](tests/)
[![Website](https://img.shields.io/badge/website-live-blue)](https://web-lemon-zeta-40.vercel.app)
[![GitHub Pages](https://img.shields.io/badge/GitHub%20Pages-live-2ea44f)](https://collaborativebioinformatics.github.io/linear-scaling-pangenome-sv/)
[![DNAnexus](https://img.shields.io/badge/DNAnexus-Group11__2026-orange)]()

> **Research question:** can a pangenome graph be built in independent chunks and stitched back together without losing haplotype paths, topology or variant calls, compared with building it in one piece?

**Short answer so far:** yes. Chunks that touch end to end (no overlap) stitch cleanly. Where the cuts go matters little. Chunked graphs are as accurate as monolithic ones outside segmental duplications and more accurate inside them. And chunking is what lets a 31 Mb chromosome arm with 20 or 40 haplotypes build on a 16 GB machine at all.

**Live site:** [collaborativebioinformatics.github.io/linear-scaling-pangenome-sv](https://collaborativebioinformatics.github.io/linear-scaling-pangenome-sv/) ([Vercel mirror](https://web-lemon-zeta-40.vercel.app))

## Contents

1. [Approach](#1-approach)
2. [Results](#2-results)
3. [Quick start](#3-quick-start)
4. [Running the experiments](#4-running-the-experiments)
5. [Pipeline](#5-pipeline)
6. [Repository layout](#6-repository-layout)
7. [Configuration](#7-configuration)
8. [Web explorer](#8-web-explorer)
9. [Documentation index](#9-documentation-index)
10. [Project status](#10-project-status)
11. [Team](#11-team)
12. [Data sources and license](#12-data-sources-and-license)

---

## 1. Approach

A pangenome graph stores many haplotypes in one structure. Shared sequence becomes common nodes, and insertions, deletions and rearrangements become alternative paths. PGGB, the graph builder this project uses, works in three stages: all-against-all alignment of the input haplotypes (wfmash), graph induction from those alignments (seqwish), and smoothing (smoothxg). The alignment stage grows roughly with the square of the number of haplotypes, and every stage grows with the length of the region. Past a few haplotypes over a few megabases, a monolithic build stops fitting on ordinary machines.

This project splits the problem:

1. **Cut** the GRCh38 region into chunks and project every cut onto each haplotype's own assembly (minimap2 alignment, longest colinear chain), so neighbouring chunks meet at exactly the same base on every haplotype.
2. **Build** one PGGB graph per chunk, in parallel, with identical settings.
3. **Stitch** the chunk graphs back into one graph with an in-house GFA library that joins each haplotype's path across every seam and reports per-haplotype PASS/FAIL at each boundary.
4. **Validate** the stitched graph against a monolithic build of the same input (topology, path identity, variant calls compared with truvari) and against truth sets that use no graph at all: dipcall-style calls from HPRC's assembly alignments, and the curated GIAB HG002 Q100 benchmark.

---

## 2. Results

All experiments use real HPRC release 2 assemblies, GRCh38 and the pinned PGGB v0.6.0 image. Each tier answers the question the previous one raised.

| Tier | Question | Finding | Write-up |
|------|----------|---------|----------|
| Pilot | Do cuts at recombination hotspots (haploblocks.org) beat fixed 400 kb windows? | Hotspot chunks looked closer to monolithic (+14 bp of node sequence vs +4,723 bp). | [RECOMBINATION_CHUNKS.md](docs/RECOMBINATION_CHUNKS.md) |
| 1 | Was that the hotspots, or the 50 kb overlap the windows had? | The overlap. Zero-overlap cuts anywhere (fixed, shifted or 20 random sets) leave no extra sequence at the seams, and hotspot cuts score at the random median (F1 0.9985 vs 0.9984). | [TIER1_CONTROLS.md](docs/TIER1_CONTROLS.md) |
| 2 | Does that hold at more haplotypes and in segmental duplications? | Seams stay clean up to 40 haplotypes, and fixed windows still match or beat hotspots. In 17q21.31 and 22q11.2, chunked graphs carry about 6.5% more sequence because paralogous copies land in different chunks. | [TIER2.md](docs/TIER2.md) |
| 3 | Which graph is right where they disagree? Do duplication-aware groups help? What happens at 10 Mb? | Against assembly-based truth calls, chunked graphs match monolithic on chr21 (F1 0.98) and beat it in 22q11.2 (0.96 vs 0.88), because monolithic PGGB merges paralogous LCR22 copies. Grouping chunks that share a duplication undoes the gain. PGGB's human settings (`-p 98 -s 10k -k 47`) reduce the merging, but chunked graphs still score higher (0.94-0.97 vs 0.86-0.93). At 10 Mb, 1 Mb chunks cost 1.22x the CPU and stay under 1 GB each, against 2.8 GB monolithic. | [TIER3.md](docs/TIER3.md) |
| 4 | Does chunking hold at chromosome scale, and against a curated truth set? | On 31 Mb of chr21 with 6 haplotypes, 29 hotspot chunks score the same as the monolithic graph (F1 0.976; HG002 0.980 against GIAB Q100) for 1.05x the CPU, at most 3.4 min and 0.97 GB per chunk against 49 min and 12.1 GB. At 20 and 40 haplotypes the monolithic build did not finish on a 16 GB runner, while chunks needed at most 1.5 and 3.0 GB. Cost is 48 to 56 CPU-seconds per Mb per haplotype up to 20 haplotypes and 73 at 40. With HG002 added, chunked graphs still beat monolithic in 22q11.2 (0.965 vs 0.898) and 17q21.31 (0.80 vs 0.70) at `-p 90`, and at `-p 98` (0.97 vs 0.94, 0.94 vs 0.85). Against GIAB alone, the 22q11.2 gap closes at `-p 98`. | [TIER4.md](docs/TIER4.md) |

### What to take from it

- **Cut without overlap.** Chunks that touch end to end stitch cleanly. `project_chunks.py` projects each cut onto every haplotype so neighbours meet exactly. A few hundred kb to about 1 Mb per chunk is a useful size; at 50-270 kb each PGGB run is mostly fixed cost.
- **Where you cut matters less than expected.** Recombination hotspots have not beaten fixed zero-overlap cuts in any comparison so far.
- **In segmental duplications, chunk boundaries help.** At `-p 90`, monolithic PGGB routes about 9% of each 22q11.2 haplotype through nodes it has already used (collapsed paralogs). Zero-overlap chunks bring that down to 2.4% and score closer to the truth. Keeping paralogous chunks together in one job makes it worse.
- **Chunking is what makes large builds possible.** A monolithic 31 Mb build needs 12 GB at 6 haplotypes and does not finish on a 16 GB machine at 20 or 40. The largest chunk at 40 haplotypes needs 3 GB.
- **Check monolithic VCFs in duplicated regions.** When a reference path loops through collapsed paralogs, `vg deconstruct` can put hundreds of kb of variants under one short top-level bubble, and PGGB's `vcfbub -l 0 -a 100000` drops them all. At `-p 98` the monolithic 22q11.2 VCF lost 86% of the truth calls. `pipeline/benchmark/tier3/path_revisits.py` shows the loops; scoring with `--bub "-l 0 -a 10000"` recovers the calls.
- **Replicate everything.** The pinned PGGB image is not deterministic: two monolithic builds of the same input differ by up to 9 nodes and 47 bp, which is the size of most effects at 1 Mb. Shared GitHub runners also vary by up to 1.7x in CPU time between replicate builds.
- **seqwish bug, fixed here.** PGGB v0.6.0 passes its default `-B 10M` to seqwish unparsed, so any input sequence over about 1 Mb comes out under-collapsed (22q11.2: 12.3 Mb of node sequence instead of 4.3 Mb). `config/pipeline.yaml` sets `transclose_batch: 100000000` and every runner passes it; set `PGGB_TRANSCLOSE_BATCH=pggb-default` for the old behaviour.

### The first benchmark (DNAnexus, August 2026)

| Metric | Value |
|--------|-------|
| Region | GRCh38 chr21:20,000,000-21,000,000 (1,000,000 bp) |
| Haplotypes | 5 (GRCh38 + HG00673/HG00733 paternal and maternal) |
| Monolithic graph | 10,389 nodes, 14,185 edges, 327 s on `mem3_ssd1_v2_x16` |
| Parallel build | 272 s (1.20x faster) with 3 overlapping 400 kb chunks |
| Stitch | Rebuilt in November 2026 with zero-overlap chunks: every haplotype path identical to monolithic, truvari F1 about 0.998 ([Tier 1](docs/TIER1_CONTROLS.md)) |

Figures from that benchmark: [results 1](slides/results_1.png), [results 2](slides/results_2.png), [speed comparison](slides/speed_comparison.html), [speed vs accuracy](slides/speed_vs_accuracy.html), [scaling by sample count](slides/scaling_by_sample_count.html) and [by genome length](slides/scaling_by_genome_length.html) (the last two are projections from a fitted model, not measurements). Slide decks: [pilot_pipeline.pptx](pilot_pipeline.pptx), [methods_results_slide.pptx](slides/methods_results_slide.pptx).

### Paper

The paper draft built from Tiers 1-4 is [docs/paper/manuscript.md](docs/paper/manuscript.md), with related work in [docs/paper/related_work.md](docs/paper/related_work.md).

---

## 3. Quick start

### Local demo (no cloud account, no HPRC data)

```bash
git clone https://github.com/collaborativebioinformatics/linear-scaling-pangenome-sv.git
cd linear-scaling-pangenome-sv
pip install pyyaml pytest
make test            # 206 tests
make demo            # synthetic graphs -> stitch -> web JSON; prints "Verdict: EQUIVALENT", "Boundaries: 2/2 PASS"
make web             # http://localhost:3000
```

`make demo` overwrites the real-data JSON in `web/public/data/` with synthetic data. Run `git checkout -- web/public/data` afterwards to get the real data back before committing.

`make check` lists the tools the real-data pipeline needs (samtools, minimap2, bcftools, truvari, and dx or aws for the cloud kits). The demo and the tests need only Python with `pyyaml` and `pytest`.

### Real data

There are three ways to run on real HPRC data:

| Where | Use it for | How |
|-------|-----------|-----|
| GitHub Actions | Tiers 3-4 and anything that takes hours | [Section 4](#4-running-the-experiments) |
| DNAnexus | The original 1 Mb benchmark | `bash dnanexus/setup_workstation.sh` then `bash dnanexus/run_pipeline.sh --upload` ([DNANEXUS.md](docs/DNANEXUS.md)) |
| AWS Batch (Spot) | haploblocks.org blocks at scale | [aws/README.md](aws/README.md) |

---

## 4. Running the experiments

Long runs go on GitHub-hosted runners (4 vCPU, 16 GB RAM, free for public repositories) through [`.github/workflows/experiment.yml`](.github/workflows/experiment.yml). Open the Actions tab, pick **Experiment**, choose a branch and a stage:

| Stage | What it runs | Inputs | Output |
|-------|-------------|--------|--------|
| `tier3-e1e2` | Tier 3 E1/E2: monolithic, chunked and SD-grouped graphs scored against the truth | `reps` | Artifact with JSON, TSV and logs |
| `tier3-e4` | Tier 3 E4: 10 Mb of chr21 | | Artifact |
| `tier3-scaling` | Tier 3 E3: 10/20/40-haplotype scaling replicate | | Artifact |
| `tier4-arm` | Tier 4 E6: chr21:15-46 Mb, monolithic vs 29 hotspot chunks (plus 31 fixed windows at N = 6) | `n` = 6, 20 or 40; `pggb`; `rep` | `docs/tier4_runs/*.json` committed to the branch |
| `tier4-hard` | Tier 4 E7: 22q11.2 and 17q21.31 with HG002, scored against GIAB v5.0q | `pggb` = `p90` or `p98`; `rep` | `docs/tier4_runs/*.json` |
| `tier5-mid` | Tier 5 E8: chr21:20-25 Mb, monolithic vs five 1 Mb windows at the same haplotype counts as Tier 4, with time and memory per PGGB stage | `n` = 6, 20 or 40; `pggb`; `rep` | `docs/tier4_runs/tier4_mid_*.json` |
| `verify-image` | Prints the PGGB and tool versions inside the pinned image | | Job log |

Tier 4 and 5 stages fan out: one job per sample (align and project the cuts), one per chunk and one for the monolithic graph (build), one per graph and `vcfbub` setting (stitch and score), then a summary job. On a `claude/` branch the summary is committed to `docs/tier4_runs/`. The `from_run` input re-summarizes an earlier run's results without rebuilding.

Scripts for each tier live in `pipeline/benchmark/tier2/` to `tier4/`, and each write-up ends with a Reproduce section.

---

## 5. Pipeline

```
HPRC release 2 assemblies + GRCh38
  -> minimap2: project each cut onto every haplotype   (pipeline/prepare/project_chunks.py)
  -> PGGB per chunk, in parallel                        (scripts/run_pggb.py)
     PGGB once over the whole region (baseline)         (pipeline/baseline/build_baseline.sh)
  -> stitch the chunk graphs                            (pipeline/merge/merge_graphs.py)
  -> validate: paths, topology, variants vs monolithic and vs truth
  -> compact JSON -> web explorer
```

| Component | Status | Where |
|-----------|--------|-------|
| GFA model and parser (S, L, P and W lines, round trip) | ✅ | `pipeline/merge/gfa.py` |
| Chunking: overlapping windows, zero-overlap windows, hotspot cuts | ✅ | `pipeline/parallel/` |
| Cut projection onto each haplotype | ✅ | `pipeline/prepare/project_chunks.py`, `pipeline/benchmark/tier4/project_coords.py` |
| Overlap-aware and touching-chunk stitch with boundary report | ✅ | `pipeline/merge/merge_graphs.py`, `pipeline/merge/validate_merge.py` |
| Monolithic vs chunked comparison (topology, paths, truvari) | ✅ | `pipeline/benchmark/compare_chunking.py`, `seam_metrics.py` |
| Truth sets: dipcall-style calls from HPRC alignments; GIAB HG002 Q100 | ✅ | `pipeline/benchmark/tier3/truth_calls.py`, `tier4/giab_truth.py` |
| Truth scoring (vg deconstruct, vcfbub, vcfwave, truvari, segdup split) | ✅ | `pipeline/benchmark/tier3/score_vs_truth.py` |
| Collapsed-paralog metric (path revisits) | ✅ | `pipeline/benchmark/tier3/path_revisits.py` |
| GitHub Actions fan-out for long runs | ✅ | `.github/workflows/experiment.yml` |
| DNAnexus applets (`pggb_chunk`, `pggb_baseline`) | ✅ | `dnanexus/` |
| haploblocks.org blocks and AWS Batch kit | ✅ | `pipeline/haploblocks/`, `aws/` |
| Web dataset export with genomic-data guard | ✅ | `pipeline/export/build_web_dataset.py`, `scripts/sync_web_results.py` |
| Stand-alone dipcall wrapper | 🟡 written, not yet run on real data | `pipeline/linear/run_dipcall.sh` |

---

## 6. Repository layout

```
config/                     pipeline.yaml (PGGB image digest and parameters), haploblocks.yaml
pipeline/
  prepare/                  minimap2 mapping, cut projection, interval extraction, multi-FASTA builder
  parallel/                 chunk manifests (windows, hotspots), chunk FASTA builder
  baseline/                 monolithic PGGB via Docker
  merge/                    GFA model, disjoint union, overlap-aware stitch, boundary report
  linear/                   dipcall and SVIM-asm wrappers (graph-free calls)
  benchmark/                graph stats, path and variant comparison, seam metrics
    tier2/ tier3/ tier4/    scripts behind each tier's write-up
  haploblocks/              haploblocks.org block -> PanSN FASTA, optional dedupe
  export/                   GFA -> bounded web JSON
scripts/                    run_pggb.py (canonical PGGB runner), HPRC index and download, demo, environment check
dnanexus/                   applets and orchestration for DNAnexus
aws/                        AWS Batch (Spot) stack for haploblocks.org blocks
web/                        Next.js + React + TypeScript explorer (Cytoscape.js)
tests/                      pytest suite (206 tests)
docs/                       write-ups, run summaries (tier*.json, tier4_runs/), paper draft
slides/                     benchmark figures and decks
.github/workflows/          ci.yml (tests, typecheck), deploy-pages.yml, experiment.yml (long runs)
```

---

## 7. Configuration

| Parameter | File | Default | Description |
|-----------|------|---------|-------------|
| `pggb.image` | config/pipeline.yaml | `ghcr.io/pangenome/pggb@sha256:44f3563c…` | Pinned PGGB v0.6.0 image; never `:latest` |
| `pggb.params` | config/pipeline.yaml | `-p 90 -s 5000 -K 19 -k 29` | Tier 3 and 4 also run `-p 98 -s 10000 -k 47` |
| `pggb.params.transclose_batch` | config/pipeline.yaml | 100000000 | seqwish `-B`; must exceed the total input bp ([TIER2.md](docs/TIER2.md)) |
| `pggb.threads` | config/pipeline.yaml | 8 | PGGB threads (the experiments use 4) |
| `target.chromosome`, `target.start/end` | config/pipeline.yaml | chr21, 20000000-21000000 | Default interval (0-based, half-open) |
| `mapping.min_query_coverage`, `mapping.min_mapq` | config/pipeline.yaml | 0.90, 20 | Interval mapping filters |
| `parallel.chunk_size_bp`, `parallel.overlap_bp` | config/pipeline.yaml | 400000, 50000 | Overlapping-window chunker (the tiers use zero overlap) |
| `pggb.params.segment_length` | config/haploblocks.yaml | 1000 | `-s` for haploblocks.org blocks, which can be under 5 kb |

---

## 8. Web explorer

A Next.js + React + TypeScript app with Cytoscape.js. Its tabs are Explore (interactive graph), Results, Chunks, Compare (two haplotypes side by side) and Speed.

```bash
cd web
npm ci
npm run dev       # http://localhost:3000
npm run build
npx vercel --prod # deploy
```

Data flow: GFA results are turned into bounded compact JSON offline by `pipeline/export/build_web_dataset.py` (node ids, lengths, degree and path membership, **no nucleotide sequences**). Then `scripts/sync_web_results.py` copies them into `web/public/data/`, refusing GFA, FASTA, VCF, BAM and CRAM files and anything over the size limit. No raw genomic files or credentials reach the browser. See [STATIC_RESULTS_ARCHITECTURE.md](docs/STATIC_RESULTS_ARCHITECTURE.md).

---

## 9. Documentation index

| Document | Contents |
|----------|----------|
| [ARCHITECTURE.md](docs/ARCHITECTURE.md) | Components and data flow |
| [METHODS.md](docs/METHODS.md) | Methods summary |
| [RECOMBINATION_CHUNKS.md](docs/RECOMBINATION_CHUNKS.md) | Pilot: hotspot vs fixed windows |
| [TIER1_CONTROLS.md](docs/TIER1_CONTROLS.md) | Overlap vs cut placement, random-cut controls |
| [TIER2.md](docs/TIER2.md) | 40 haplotypes, segmental duplications, recombination maps |
| [TIER3.md](docs/TIER3.md) | Truth set, duplication-aware groups, 10 Mb, `-p 98`, the vcfbub trap |
| [TIER4.md](docs/TIER4.md) | chr21 long arm at 6/20/40 haplotypes, HG002 vs GIAB Q100 |
| [TIER5.md](docs/TIER5.md) | Plan and results for the next round: chunking vs haplotype count, cost per PGGB stage |
| [paper/manuscript.md](docs/paper/manuscript.md) | Paper draft |
| [DNANEXUS.md](docs/DNANEXUS.md) | Running on DNAnexus |
| [DEMO.md](docs/DEMO.md) | The synthetic demo |
| [TASK_BREAKDOWN.md](docs/TASK_BREAKDOWN.md) | Hackathon task breakdown (historical) |
| [FINAL_INTEGRATION_AUDIT.md](docs/FINAL_INTEGRATION_AUDIT.md), [FINAL_LOGIC_AUDIT.md](docs/FINAL_LOGIC_AUDIT.md), [ZERO_EDGE_ROOT_CAUSE.md](docs/ZERO_EDGE_ROOT_CAUSE.md) | Audits and bug write-ups |

Machine-readable results: `docs/tier1_controls.json`, `docs/tier2.json`, `docs/tier3.json`, `docs/tier4_runs/*.json`.

---

## 10. Project status

The hackathon checklist, updated:

| Item | Status | Evidence |
|------|--------|----------|
| Q1-Q3: stitch, stitch tests, boundary report | ✅ | `make demo`: EQUIVALENT, 2/2 boundaries PASS |
| M1: PGGB image verified | ✅ | The pinned digest is ghcr.io tag `202607241950514225c6`, which `latest` also points to; ghcr.io carries no `v0.6.0` tag, so the version comes from `pggb --version` inside the image |
| K1: graph-free calls | ✅ via Tier 3 | Dipcall-style calls from HPRC's assembly-to-GRCh38 alignments (`tier3/truth_calls.py`); the stand-alone `run_dipcall.sh` has not been run |
| K2: variant comparison | ✅ | vg deconstruct + truvari in every tier (`compare_chunking.py`, `score_vs_truth.py`) |
| A1-A6: DNAnexus applets, web export and guard | ✅ | |
| A7-A10: benchmark display, chunk view, compare mode, Vercel deploy | ✅ | Results, Chunks and Compare tabs; site live on Vercel and GitHub Pages |
| D1-D4: documentation | ✅ | ARCHITECTURE.md, METHODS.md, decks in `slides/`, this README |
| S1: integration test | ✅ | Synthetic (`make demo`) and real data (Tiers 1-4) |
| S2: tests on CI | ✅ | `ci.yml`: 206 passed |
| S3: tagged release | ☐ | Not yet tagged |

---

## 11. Team

| Person | Area |
|--------|------|
| Quang | Pipeline architecture, graph stitching |
| Michael | Pangenome graph construction (PGGB) |
| Khoi | Linear and assembly-based workflow, validation |
| Ali | DNAnexus, integration, web |
| Alexander | Documentation and presentation |
| Lex | PGGB benchmarking |

Where to work: graph construction in `pipeline/baseline/` and `pipeline/parallel/`; stitching in `pipeline/merge/`; cloud runs in `dnanexus/`, `aws/` and `.github/workflows/`; benchmarking in `pipeline/benchmark/`; the site in `web/`; writing in `docs/`.

---

## 12. Data sources and license

All inputs are public: HPRC release 2 assemblies and alignments (HPRC S3 bucket), the HG002 Q100 assembly and GIAB v5.0q benchmark (NIST), GRCh38 (Ensembl release 112), GIAB v3.3 genome stratifications, haploblocks.org recombination blocks, and stdpopsim recombination maps.

MIT, see [LICENSE](LICENSE). Repository: https://github.com/collaborativebioinformatics/linear-scaling-pangenome-sv
