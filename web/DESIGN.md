# Direction contract: Parallel Pangenome Graph Explorer

## Step 0
- **Mode:** Operate, with Read for the explanatory bands. Researchers explore a graph; hackathon judges read the result.
- **Level:** 2. The task comes first, with one idea that makes the site belong to its subject. Motion is feedback only.
- **Kind:** research software, a scientific project from the BCM SV Hackathon 2026.
- **Job:** explore the graph, then judge the result (equivalence, chunks, speed).

## The rut
The category default is a near-black hero with neon green and blue, a big-number metric row and bento cards. That was the previous pass. The predictable opposite is a cream editorial page with a serif display face. Neither one is used here.

## Candidate worlds
A genome browser session (UCSC/IGV ruler and tracks); a Sanger chromatogram printout; a karyotype ideogram; a lab notebook on graph paper; a journal figure panel; a transit map (rejected because the graph explorer already uses the subway idea literally).

## Thesis
The page reads like a genome browser looking at the benchmark region. A locus box (GRCh38 chr21:20,000,000–21,000,000) sits above the headline. Under it is a track panel: a coordinate ruler, the real chunk layout from overview.json (cores dark, overlaps grey), and a sketch of the graph with a shared backbone and a unique branch. It replaces the generic dark dashboard hero.

## World
- **Ground:** browser white, light mode. Readers sit at lab desks, and judges watch laptops or projectors in lit rooms.
- **Colour:** restrained, with cool-tinted neutrals. Two hues carry data meaning only: green means shared, blue means unique/selected. Orange appears only for the baseline in the speed chart, and that line is labelled directly so it never relies on colour alone. Amber marks the clicked node. Every colour is an oklch token on :root, and canvas and Cytoscape read the same tokens.
- **Type:** Atkinson Hyperlegible Next for text. It is built for legibility on projectors and has a slashed zero, so sample IDs such as HG00673 can't be misread. JetBrains Mono is used only for coordinates, IDs and counts. Both are OFL and self-hosted via next/font/local.
- **Components:** hairline-ruled panels, 4–6 px radii, labelled tracks, tables for comparisons (no big-number tiles), and tick marks on rules for ordered steps.

## Story
What is being tested (the hero and its track panel) → the concepts → five steps on a ruler → explore the graph → the results table, chunks, compare and speed tabs.

## First viewport
At 1280×800 the locus box, the two-line headline, the lede, "Explore the graph" and the top of the track panel are visible.

## Signature and motion
The signature is the browser track panel built from real chunk coordinates. Motion is limited to 160 ms colour and transform feedback. Everything rests under prefers-reduced-motion, including the graph search zoom.
