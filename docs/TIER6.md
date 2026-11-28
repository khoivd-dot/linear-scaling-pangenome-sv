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

## Which PGGB the pinned image holds (checklist item M1)

Every tier so far used `ghcr.io/pangenome/pggb@sha256:44f3563c…` and called it
PGGB v0.6.0. The new `verify-image` stage ran each tool's version command
inside that image (run 46159869291):

| Tool | Version reported |
|------|------------------|
| pggb | `4225c6c` (a commit, not a release tag) |
| wfmash | `d47b7e3` |
| seqwish | v0.7.11-35-gaab5bee |
| smoothxg | v0.8.2-9-g0ea0470 |
| vg | v1.71.0 "Cera" |

On ghcr.io the digest is the image tagged `202607241950514225c6` (built
2026-07-24 from pggb commit 4225c6c), which `latest` also points to today.
There is no `v0.6.0` tag. The results stand, since they all come from this
one pinned image, but "PGGB v0.6.0" in the write-ups should read "PGGB at
commit 4225c6c (ghcr.io build of 2026-07-24)". The seqwish `-B` finding in
Tier 2 also belongs to this build; whether the v0.6.0 release behaves the
same was not tested.

## Results so far (N = 40 still running)

**Haplotype count, not chunking, lowers accuracy.** In chr21:20-25 Mb the
monolithic graph loses accuracy as haplotypes are added, just as the chunked
one does. At each N the two score within 0.001 of each other on all variants,
and the chunked graph scores higher on SVs. So the Tier 4 drop from 6 to 20
haplotypes is a property of larger graphs (or of calling from them), not a
cost of cutting the chromosome.

**wfmash is what grows.** Going from 7 to 21 paths multiplied wfmash's CPU by
about 10, close to the square of the path count, while smoothxg grew by 2.8.
wfmash took about 10% of chunked CPU at N = 6 and a third at N = 20, which
explains why cost per haplotype rises at N = 40 in Tier 4.

**The truth set agrees with dipcall.** The stand-alone dipcall and the Tier 3
calls agree at F1 0.998 on all variants for all four haplotypes. Of dipcall's
133 SVs, the truth calls find 131 and add 2 of their own (E10 below).

### E8: accuracy and cost

HG00673, HG00733 and HG002 scored together; vcfbub `-a 100000`. The region has
no segmental duplications, so only the totals mean anything.

| N | Graph | F1 all | F1 SV | HG002 F1 (hap 1 / 2) | CPU s | Max wall s | Peak MB | Seams pass | Paths identical to monolithic |
|---|-------|--------|-------|----------------------|-------|------------|---------|------------|-------------------------------|
| 6 | monolithic | 0.9793 | 0.879 | 0.9836 / 0.9817 | 1219 | 385 | 1978 | | |
| 6 | 1 Mb windows | 0.9789 | 0.900 | 0.9837 / 0.9815 | 1435 | 152 | 972 | 4/4 | 7/7 |
| 20 | monolithic | 0.9733 | 0.855 | 0.9781 / 0.9782 | 5115 | 1471 | 5639 | | |
| 20 | 1 Mb windows | 0.9742 | 0.887 | 0.9787 / 0.9790 | 4909 | 411 | 1214 | 4/4 | 21/21 |
| 40 | both | running | | | | | | | |

With vcfbub `-a 10000` every F1 changes by less than 0.003. Chunking cost
18% more CPU at N = 6 and 4% less at N = 20, cut the longest wall time by a
factor of 2.5 to 3.6, and kept peak memory near 1 GB while the monolithic
peak rose from 2.0 to 5.6 GB.

### E8: CPU seconds per stage

Chunked rows are the sum over the five windows; peak memory in brackets.

| N | Graph | wfmash | seqwish | smoothxg | odgi |
|---|-------|--------|---------|----------|------|
| 6 | monolithic | 169 | 28 (1978 MB) | 873 | 148 (994 MB) |
| 6 | 1 Mb windows | 149 | 27 (464 MB) | 962 | 295 (972 MB) |
| 20 | monolithic | 1754 | 172 (5639 MB) | 2458 | 731 (1008 MB) |
| 20 | 1 Mb windows | 1596 | 159 (1214 MB) | 2367 | 784 (974 MB) |

seqwish sets the monolithic memory peak. In the chunked builds the floor is
odgi at about 1 GB per window, whatever the window holds. gfaffix takes under
a second everywhere.

### E10: stand-alone dipcall against the Tier 3 truth calls

Stage `dipcall-check` ran dipcall 0.3 on the chr21 contigs of HG00673 and
HG00733 and scored the Tier 3 calls against it in chr21:20-25 Mb, inside both
callable BEDs (about 4.998 of the 5 Mb). Same truvari settings as everywhere
else. Results: [`docs/tier6_runs/`](tier6_runs/).

| Haplotype | F1 all | F1 SV | SV calls |
|-----------|--------|-------|----------|
| HG00673#1 | 0.9981 | 1.000 | 34 of 34 found |
| HG00673#2 | 0.9979 | 1.000 | 33 of 33 found |
| HG00733#1 | 0.9977 | 0.952 | 30 of 32 found, 1 extra |
| HG00733#2 | 0.9982 | 0.986 | 34 of 34 found, 1 extra |

The two callers align the same assemblies in different ways (dipcall maps
whole contigs with minimap2, Tier 3 reads HPRC's winnowmap alignments), so
agreement at this level says the truth calls are not an artefact of one
aligner. This region has no segmental duplications, so it does not test the
truth in 17q21.31 or 22q11.2.


## E8 design: chunking against haplotype count, with time per stage

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

E10 is done (above). The rest have not started. Each runs on free GitHub runners unless noted.

| ID | Question | Design | Cost |
|----|----------|--------|------|
| E9 | Do scale and segmental duplications interact? | The chr22 long arm, which holds 22q11.2, at 6 haplotypes, with HG002 scored against GIAB. Needs a scan for assemblies that span the LCR22s on one contig. | Moderate; the scan is the main work |
| E10 | Does the stand-alone dipcall agree with the Tier 3 truth calls? | Done in chr21:20-25 Mb (F1 0.998). Still open: 22q11.2 and 17q21.31, where segmental duplications could make the two callers disagree. | Small |
| E11 | What chunk size is best at 40 haplotypes? | chr21:20-25 Mb with 500 kb, 1 Mb and 2.5 Mb windows at N = 40. | Small once E8 is in |
| E12 | Monolithic arm at 20 and 40 haplotypes | The Tier 4 arm on a 64 GB machine, to compare with the chunked graphs directly. | Paid compute; needs approval |
