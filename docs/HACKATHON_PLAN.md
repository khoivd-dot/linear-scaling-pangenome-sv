# Hackathon task plan

The original hackathon task plan, kept for the record. For the current state of the project, see the [README](../README.md).

## Who Is Doing What?

**Full task breakdown with commands, file paths, and acceptance criteria:**
👉 **[docs/TASK_BREAKDOWN.md](TASK_BREAKDOWN.md)**

| Person | Primary Area | What They Are Working On |
|--------|-------------|--------------------------|
| **Quang** | Pipeline architecture + graph merging | ✅ **Q1–Q3:** overlap-aware stitch, `tests/test_merge.py`, boundary report. `make demo` → Strategy: overlap_aware, EQUIVALENT, 5 components. |
| **Michael** | Pangenome graph construction | M1: Run real PGGB applet on DNAnexus (micro-test). M2/M3: Build baseline + chunk graphs on DNAnexus. M4: Full 1 Mb smoke test. |
| **Khoi** | Pipeline + linear/assembly workflow | K1: Dipcall wrapper validation. K2: SVIM-asm wrapper. K3: Truvari variant comparison. K4: Method documentation. |
| **Ali** | Pipeline + DNAnexus + integration + web | ✅ A1-A5 complete. ✅ A6 (web dataset export) complete. 🔵 **A7-A10 (current):** benchmark display polish, chunk visualization, compare mode, Vercel deploy. |
| **Alexander** | Documentation + presentation | D1: Update architecture doc. D2: DNAnexus operations guide. D3: Web visualization guide. D4: Final presentation. |
| **Lex** | PGGB benchmarking | L1: Local PGGB cost benchmarks. L2: Graph cost model. L3: Document findings. L4: Merge into pipeline. |

> **Important:** The pipeline is a shared technical responsibility. Everyone collaboratively implements and troubleshoots.

## Where Should I Work?

```
Want to work on PGGB / graph construction?
  pipeline/baseline/ and pipeline/parallel/  (Michael + Khoi)

Want to work on graph stitching?
  pipeline/merge/  (Quang + Michael + Khoi + Ali)

Want to work on DNAnexus / cloud execution?
  dnanexus/  (Ali + Quang + Khoi)

Want to work on benchmarking / validation?
  pipeline/benchmark/  (Ali + Khoi + Michael)

Want to work on the website / visualization?
  web/  (Ali + anyone)

Want to work on writing / presentation?
  docs/  (Alexander)
```

## What To Do Next (Pre-Submission)

All 19 P0 audit items are fixed. Below is the exact pre-submission task list broken down per person. Stitch is implemented. Independent per-chunk PGGB on DNAnexus (`/data/prepared/chunks/*.fa`) is unblocked.

### Quang — Merge Algorithm (P0, highest priority)

Stitch is the single critical unimplemented component.

| # | Task | Area | Instructions |
|---|------|------|-------------|
| Q1 | **Overlap-aware stitch** | `pipeline/merge/merge_graphs.py` | **DONE.** `overlap_aware_stitch()`; `make demo` writes `results/merge/merged.gfa`. |
| Q2 | **Stitch tests** | `tests/test_merge.py` | **DONE.** Path continuity, sequence identity, no invented gap edges, demo EQUIVALENT. |
| Q3 | **Boundary report** | `pipeline/merge/validate_merge.py` | **DONE.** Per-haplotype PASS/FAIL + `results/merge/boundary_report.tsv`. |

### Michael — PGGB + Graph Validation (P0)

| # | Task | Area | Instructions |
|---|------|------|-------------|
| M1 | **PGGB version validation** | `config/pipeline.yaml` | Verify the pinned image digest matches actual `v0.6.0` tag on ghcr.io. Update if mismatch. |
| M2 | **Chunk FASTA review** | `pipeline/parallel/build_all_chunks.py` | Review parent-locus constraint (lines 180-202): contig FATAL, strand warning, margin calc. |
| M3 | **Baseline graph review** | `pipeline/baseline/build_baseline.sh` | Confirm `-w 50000` used by both baseline and chunks. Verify command matches manual pggb invocation. |
| M4 | **Test data creation** | `tests/` | Create a small (5 Kb) test FASTA with 2-3 haplotypes for stitch testing without real HPRC data. |

### Khoi — Linear Pipeline + Validation (P0/P1)

| # | Task | Area | Instructions |
|---|------|------|-------------|
| K1 | **Linear dipcall wrapper** | `pipeline/linear/run_dipcall.sh` | Download/verify dipcall, run on 4 HPRC haplotypes vs GRCh38 chr21, produce VCF. Needed for benchmark. |
| K2 | **Variant comparison** | `pipeline/benchmark/benchmark_variants.sh` | `vg deconstruct` on baseline + merged GFAs, `truvari` compare baseline vs merged VCF. Record precision/recall/F1. |
| K3 | **Benchmark audit** | `pipeline/benchmark/benchmark_variants.sh` | Fix existing stub. Wire results into `results/benchmark/`. |
| K4 | **Real data smoke test** | — | Once stitch ready (Q1), run full pipeline on chr21:20000000-21000000. Report baseline/chunks/stitch/equivalence status. |

### Ali — DNAnexus + Integration + Web (P0/P1)

| # | Task | Area | Instructions |
|---|------|------|-------------|
| A1 | **Baseline applet test** | `dnanexus/applets/pggb_baseline/` | Build: `dx build --brief`. Verify same instance type as chunk (`mem3_ssd1_v2_x16`). |
| A2 | **Orchestration dry run** | `dnanexus/run_parallel_chunks.sh` | Verify expected=local=uploaded=submitted=downloaded. Simulate missing FASTA, verify FATAL. |
| A3 | **Timing dashboard** | `dnanexus/run_parallel_chunks.sh` | Review `graph_parallel_wall_seconds` (max stop - min start) and `sum_worker_seconds`. |
| A4 | **Web stitch status** | `web/app/page.tsx` | Reads `latest.json` stitch/equivalence badges (`PASS` / `EQUIVALENT` after `make demo`). |
| A5 | **Web JSON guard** | `scripts/sync_web_results.py` | Prevent accidental GFA copy into `web/public/data/`. Only compact JSON allowed. |

### Alexander — Documentation + Presentation (P0)

| # | Task | Area | Instructions |
|---|------|------|-------------|
| D1 | **Architecture diagram** | `docs/ARCHITECTURE.md` | Draw pipeline architecture (ASCII/Mermaid). Mark merge as SYNTHETIC_VALIDATED, equivalence NOT_RUN. |
| D2 | **Methods summary** | `docs/` | 1-2 pages: parallelization strategy, parent-locus constraint, provenance hashing, benchmark design. |
| D3 | **Presentation slides** | `docs/presentation/` | Research question, architecture, status, timeline, next steps. |
| D4 | **README review** | `README.md` | Verify "Who Is Doing What" table. Read for clarity. |

### Shared — Critical Path

| # | Task | Who | Instructions |
|---|------|-----|-------------|
| S1 | **Integration test** | Quang + Michael + Ali | **Quang part done:** `make demo` → Strategy: overlap_aware, 5 components, EQUIVALENT. |
| S2 | **Final validation** | Everyone | `make check && make test` from clean checkout. All 55+ tests pass. CI on the merge commit. |
| S3 | **Tag release** | Ali | Squash-merge working branch to main. Tag `v0.1.0-pre`. Do NOT launch real PGGB. |

### Timeline

```
Now ___ Q1 (stitch: 1-2d) ___ S1 (integration: 1d) ___ S2 (validation: 1d) ___ Submission
        M1-M4, A1-A5 (parallel)   D1-D4 (parallel)        K1-K4 (parallel)
```

### Stitch Algorithm Reference (for Q1)

1. **Input:** Two adjacent chunk GFAs with known overlap region from `chunk_manifest.tsv`.
2. **Find common subpaths:** Extract each haplotype's path segment covering the overlap.
3. **Create weld mapping:** Map segment IDs from chunk_B's overlap to chunk_A's by sequence identity.
4. **Merge nodes:** Replace chunk_B IDs with chunk_A IDs in links/paths. Drop duplicate segments.
5. **Concatenate remaining:** Keep both graphs' content outside the overlap.
6. **Validate:** For every haplotype, verify path continuity and sequence identity to monolithic baseline.

Start with exact-match welding before attempting alignment-based welding for near-matches.

### Submission Checklist

```
☑ Q1: Stitch implemented and tested
☑ Q2: Stitch tests pass
☑ Q3: Boundary report
☐ M1: PGGB version verified
☐ K1: Dipcall wrapper complete
☐ K2: Variant comparison working
☑ A1-A5: DNAnexus + web updates done
☑ A6: Web dataset export (build_web_dataset.py)
☐ A7-A10: benchmark display, chunk viz, compare mode, Vercel deploy
☐ D1-D4: Documentation complete
☐ S1: Full integration test passes
☐ S2: All 55+ tests pass on CI
☐ S3: Tagged v0.1.0-pre
```
