# Cutting chunks at recombination hotspots

**Question.** The parallel build cuts the 1 Mb pilot region into fixed 400 kb
windows with 50 kb of overlap. Does cutting at recombination hotspots instead
change how close the stitched graph comes to the monolithic one?

**Short answer.** On the chr21 pilot with the same five haplotypes, three chunks
cut at hotspots (no overlap) gave a graph closer to the monolithic build than
three fixed windows: 5 more nodes than monolithic instead of 14, 14 more bp of
node sequence instead of 4,723, the same variant agreement, and all 24 SVs
matched. Ten small hotspot chunks also stitched cleanly, but each PGGB run
carries roughly a minute of fixed cost at this size, so total compute was 2.6
times that of three chunks and 5 times the monolithic run. These are single runs on one region; see the limits below.

## Setup

- **Region and input:** GRCh38 chr21:20,000,000-21,000,000 plus HG00673 and
  HG00733, both haplotypes, from HPRC release 2 (`config/samples.yaml`). Each
  haplotype's assembly holds the region on one contig, aligned in 3 to 5
  colinear minimap2 pieces (asm5, MAPQ 60) that `project_chunks.py` chains.
  Haplotype lengths are 1,001,586 to 1,001,742 bp.
- **Hotspot cuts:** the 19 haploblocks.org boundaries inside the region
  (`aws/blocks/chr21_pilot_20-21Mb.txt`; recombination-split blocks from 1000
  Genomes). `hotspot_chunks.py` merges blocks shorter than `--min-size`.
- **Projection:** every chunk boundary is projected onto every haplotype through
  the same chained alignment, so neighbouring chunks meet exactly on each
  haplotype. The hotspot chunks concatenate back to each full haplotype base
  for base (checked).
- **Stitch:** the existing overlap-aware stitch. `seam_position` now joins
  intervals that touch exactly (zero overlap) end to end; any gap still fails.
- **PGGB:** pinned v0.6.0, `config/pipeline.yaml` parameters (-p 90 -s 5000
  -K 19 -k 29), 4 threads per build, one build at a time on a 4 vCPU, 15 GB
  machine, 2026-11-23.
- **Variants:** `vg deconstruct -P GRCh38 -a` on every graph, matched exactly
  on (POS, REF, ALT) against the monolithic calls. This is stricter than
  truvari, which allows position and sequence tolerance.

## Results

| | Monolithic | Windows 400 kb, 50 kb overlap | Hotspots, min 300 kb | Hotspots, min 50 kb |
|---|---:|---:|---:|---:|
| Chunks | 1 | 3 | 3 | 10 |
| Chunk sizes (kb) | 1,000 | 425, 450, 225 | 348, 302, 350 | 51 to 271 |
| Nodes | 10,382 | 10,396 | 10,387 | 10,406 |
| Edges | 14,176 | 14,195 | 14,182 | 14,207 |
| Node sequence (bp) | 1,006,549 | 1,011,272 | 1,006,563 | 1,006,553 |
| Components | 1 | 1 | 1 | 1 |
| Boundaries PASS | | 2/2 | 2/2 | 9/9 |
| Haplotypes spelled identical to monolithic | | 5/5 | 5/5 | 5/5 |
| Variant F1 vs monolithic (3,754 calls) | | 0.9975 | 0.9979 | 0.9981 |
| SVs (50 bp+) matched exactly | | 24/24 | 24/24 | 21/24 |
| PGGB wall, longest chunk (s) | 141.6 | 100.0 | 94.9 | 87.3 |
| PGGB wall, all chunks summed (s) | 141.6 | 276.9 | 274.6 | 712.8 |

The monolithic graph here (10,382 nodes, 14,176 edges) is within 0.1% of the
DNAnexus baseline in the README (10,389 and 14,185), on 4 threads instead of 16.

What the numbers say:

- **Fixed windows leave about 4.7 kb of extra node sequence** that the hotspot
  chunks do not. This is consistent with each 50 kb overlap being aligned twice,
  once in each window, and the two copies not collapsing perfectly around the seam.
  Hotspot chunks share no sequence, so there is nothing to collapse.
- **The three SVs the 10-chunk run "missed" are not seam artifacts.** They
  are the same 130 to 156 bp insertions reported 5 bp away (20,742,472 vs
  20,742,477), 41 kb from the nearest cut, which is a left-alignment
  difference inside a repeat. All of that run's mismatches sit at two loci,
  both more than 10 kb from any cut. A tolerant comparison such as truvari should match
  them; exact matching does not.
- **Small chunks cost more in total.** Every 10-chunk run took 60 to 87 s
  regardless of size (51 to 271 kb), so most of it is fixed per-run cost at this
  scale, not alignment. Fully parallel, the 10-chunk build would finish in
  about 87 s against 142 s monolithic, but it uses 5 times the compute. With
  5 haplotypes and 1 Mb, the useful unit is a few hundred kb. Hotspots should
  decide where to cut, and a minimum size should decide how many cuts to make.

## Limits

- One region, five haplotypes, one run per build. Timings vary by a few
  seconds between runs, and so does the graph: replicates in
  [TIER1_CONTROLS.md](TIER1_CONTROLS.md) show the monolithic build itself
  changing by up to 9 nodes and 47 bp between runs of the pinned image. That
  follow-up also shows the hotspot advantage here came from dropping the 50 kb
  overlap, not from where the cuts fall.
- The hotspot map comes from 1000 Genomes short-read data, while the
  haplotypes here are HPRC long-read assemblies. The cut positions are
  reference coordinates, so they transfer, but they were not derived from these
  samples.
- The variant check is exact matching from `vg deconstruct`; the truvari
  comparison (K3) is still the one to quote.
- Scaling with many haplotypes (where the all-pairs cost dominates) is not
  tested here; that is what the AWS kit in `aws/` is for.

## Reproduce

```bash
# chunk manifests
python3 pipeline/parallel/hotspot_chunks.py aws/blocks/chr21_pilot_20-21Mb.txt hotspots.tsv \
    --chrom chr21 --start 20000000 --end 21000000 --min-size 300000
# windows: pipeline/parallel/make_chunks.py create_chunks(..., 400000, 50000)

# per haplotype: locate the region, pull the contig, align with CIGAR
minimap2 -x asm5 -t 4 region.fa HG00673_mat_hprc_r2_v1.0.1.fa.gz > approx.paf
samtools faidx HG00673_mat_hprc_r2_v1.0.1.fa.gz 'HG00673#2#JAHBBY020000010.1' > mat.contig.fa
minimap2 -x asm5 -c --eqx -t 4 region.fa mat.contig.fa > HG00673_mat.paf

python3 pipeline/prepare/project_chunks.py hotspots.tsv region.fa hotspots/ \
    --region-start 20000000 --ref-name 'GRCh38#0#chr21' \
    --hap x=HG00673_mat.paf=mat.contig.fa   # one --hap per haplotype

python3 pipeline/benchmark/compare_chunking.py run/ --mono mono/_mono.fa \
    --strategy windows=windows.tsv=windows --strategy hotspots=hotspots.tsv=hotspots --threads 4
```

Raw output: [`recombination_chunks_comparison.json`](recombination_chunks_comparison.json).
