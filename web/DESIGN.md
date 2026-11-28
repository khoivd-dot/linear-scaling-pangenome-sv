# Direction contract: partitioned pangenome graphs, results site

## Step 0
- **Mode:** Read first (reviewers, collaborators and hackathon judges follow the results in order), Operate for the graph explorer.
- **Level:** 2. One idea that belongs to the subject, motion only as feedback.
- **Kind:** a research project (BCM SV Hackathon 2026) with real measurements and an interactive pilot graph.
- **Job:** understand the results, then open the write-ups or explore the graph.

## The rut
The category default is a near-black hero with neon green and blue, a big-number metric row and bento cards. Its predictable opposite is a cream editorial page with a serif display face. Neither is used.

## Direction kept, and why
The genome-browser session from the previous pass still fits best: the audience reads UCSC/IGV views every day, and a browser track is the honest way to show where chunks fall on a chromosome. Candidates reconsidered: karyotype ideogram (too decorative, no place for measurements), journal figure panel (static, no first-viewport image), lab notebook (costume), transit map (the explorer already uses that idea). This pass moves the browser from the 1 Mb pilot to the result that matters.

## Thesis
The hero is a browser view of the 31 Mb chr21 arm built from real data: a coordinate ruler, the 29 hotspot chunks at their true positions, and two peak-memory signal tracks on one 0-16 GB scale (the runner's RAM). The monolithic build is one orange slab at 12.1 GB; the chunks are a low petrol comb under 1 GB. It is the image-first hero and the argument in one picture.

## World
- **Ground:** browser white, light mode only. Readers sit at lab desks or watch a projector in a lit room.
- **Colour:** cool neutrals. Colour carries data meaning only. Charts: petrol `--chunked` oklch(45% 0.12 230) for chunked graphs, orange `--mono` oklch(56% 0.15 45) for monolithic graphs, grey `--noise` for random controls and run-to-run noise. The pair passed the dataviz validator in light mode (CVD dE 17.5, normal dE 27.0, contrast >= 3:1). Graph explorer: green = shared, blue = unique, amber = the clicked node. Every colour is an oklch token on `:root`.
- **Type:** Atkinson Hyperlegible Next for text (slashed zero, so HG00673 cannot be misread) and JetBrains Mono for coordinates, IDs, parameters and table figures. Both OFL 1.1, self-hosted with next/font/local, licences in `app/fonts/`. No CDN, no Google Fonts.
- **Components:** locus boxes naming the region each section's evidence comes from, hairline panels with 4-6 px radii, tables for comparisons, ruler ticks for ordered steps.

## Story
1. What it is: reference-projected partitioning (positional GRCh38 windows, zero overlap, exact stitch), and what the project adds (measurement, not the mechanism; impg does the same partition and lace).
2. Tier 4: the chr21 arm, cost and accuracy by haplotype count, with builds that did not finish marked as such.
3. Hard regions: the fair -p 98 / vcfbub -a 10000 comparison first, 17q21.31 shown as mixed, -p 90 as context, path revisits, the vcfbub trap.
4. Where to cut: pilot, Tier 1 controls against run-to-run noise, Tier 2 and the seqwish bug.
5. Explore the pilot graph, methods and reproduction, team, write-ups.

## Charts
HTML and CSS marks positioned in percent, so labels keep their size on a phone and the page reads with JavaScript off. Lines drawn with a non-scaling SVG stroke. Every chart has a legend (or labelled rows), direct labels on the points that matter, a native title tooltip per mark and a "Show the numbers" table. Axes that do not start at zero say so.

## Data
Every number comes from `public/data/results.json`, built from `docs/` by `scripts/build_site_results.py`; `tests/test_site_results.py` fails when it drifts. Numbers that exist only in a write-up are constants in that script, each naming its doc and section.

## First viewport
At 1440x900: locus box, two-line headline, lede, both actions and the top three browser tracks. On a phone the headline, lede, actions and the ruler and chunk tracks.

## Motion
160 ms colour feedback on controls only. Everything rests under prefers-reduced-motion.
