import data from "../public/data/results.json";
import Explorer from "../components/Explorer";
import { Bars, DotRows, Figure, Legend, ScalingPlot } from "../components/Charts";

// Every number on this page comes from public/data/results.json, built from docs/ by
// scripts/build_site_results.py (tests/test_site_results.py fails if it drifts).
interface Build {
  status?: string; after_min?: number; why?: string;
  builds?: number; cpu_h?: number; longest_min?: number; peak_gb?: number; cpu_s_per_mb_hap?: number;
  f1_all?: number; f1_sv?: number; hg002_f1?: number[]; cpu_vs_mono?: number; stitch_s?: number;
  seams_pass?: number; seams?: number; paths_identical?: number | null; paths?: number | null;
}
interface Scale { n: number; mono: Build; hotspots1m: Build; windows1m?: Build }
interface HardGraph { f1_a10k: number[]; f1_a100k: number[]; f1_sd_a10k: number[]; peak_gb: number; ref_revisit_pct: number; hg002_f1?: number[] }
interface T3Graph { f1: number[]; revisit_pct: number[]; sv_f1: number[]; f1_sd?: number; f1_a10k?: number }
interface Region { key: string; name: string; locus: string; mono: T3Graph; chunked: T3Graph; sdgroup?: T3Graph }
const R = data as unknown as {
  tier4: { runner: { vcpu: number; ram_gb: number }; arm: { start: number; end: number; mb: number };
    arm_truth_svs: number; arm_sv_diff_calls: number[]; max_concurrent_jobs: number; giab_22q_coverage_pct: number; scaling: Scale[];
    chunks_n6: { start: number; end: number; wall_min: number; peak_gb: number }[];
    hard: { region: string; pggb: string; mono: HardGraph; chunked: HardGraph }[] };
  tier3: { regions: Region[]; r17_truth_svs: number; vcfbub_trap: { region: string; mono_a100k: number; mono_a10k: number; chunked_a10k: number }[];
    e4_10mb: { mono_wall_s: number; mono_peak_gb: number; chunk_longest_s: number; chunk_peak_mb: number; chunks: number } };
  where_to_cut: {
    pilot: { windows_extra_bp: number; hotspots_extra_bp: number; overlap_kb: number };
    tier1: { f1: Record<"windows" | "windows0" | "hotspots300", number[]>; near_cut_bp: Record<"windows" | "windows0" | "hotspots300", number[]>;
      random_f1: number[]; random_f1_median: number; mono_vs_mono_f1: number[]; mono_drift_nodes: number; mono_drift_bp: number };
    tier2: { seqwish_bug_bp: number; seqwish_fixed_bp: number; hard_mono_vs_mono_f1: number[];
      n40: { mono_cpu_s: number; chunk_cpu_s: number; longest_vs_mono_pct: number; mono_peak_gb: number; chunk_peak_mb: number } };
  };
};

const GH = "https://github.com/khoivd-dot/linear-scaling-pangenome-sv";
const doc = (name: string) => `${GH}/blob/main/docs/${name}`;
const n = (x: number) => x.toLocaleString("en-US");
const f3 = (x: number) => x.toFixed(3);
const f2 = (x: number) => x.toFixed(2);
const mb = (bp: number) => `${bp / 1e6}`;

const T4 = R.tier4;
const [S6, S20, S40] = T4.scaling;
const H6 = S6.hotspots1m, W6 = S6.windows1m!, M6 = S6.mono;
const T3 = R.tier3;
const [r21, r17, r22] = T3.regions;
const C = R.where_to_cut;
const T1 = C.tier1;
const pick = (region: string, pggb: string) => T4.hard.find(h => h.region === region && h.pggb === pggb)!;
const P98_22 = pick("22q11.2", "-p 98"), P98_17 = pick("17q21.31", "-p 98"), P90_22 = pick("22q11.2", "-p 90");
// The fair comparison (-p 98) first, -p 90 as context.
const HARD_ORDER = [P98_22, pick("17q21.31", "-p 98"), P90_22, pick("17q21.31", "-p 90")];
const repDiff = Math.max(...T3.regions.flatMap(g => [g.mono, g.chunked, g.sdgroup].filter(Boolean).map(x => Math.abs(x!.f1[0] - x!.f1[1]))));

const WRITEUPS = [
  ["RECOMBINATION_CHUNKS.md", "Pilot: hotspot cuts on 1 Mb of chr21"],
  ["TIER1_CONTROLS.md", "Tier 1: shifted and random cuts, three replicates"],
  ["TIER2.md", "Tier 2: 40 haplotypes, duplicated regions, the seqwish bug"],
  ["TIER3.md", "Tier 3: scoring against a truth set, 10 Mb"],
  ["TIER4.md", "Tier 4: the chr21 arm, HG002 against GIAB"],
];

const AUTHORS = [
  ["Nhat Quang Tran", "pipeline architecture and graph stitching"],
  ["Dang Khoi Vu", "pipeline, variant comparison and validation"],
  ["Michael Olufemi", "pangenome graph construction"],
  ["Alexander W. Lehr", "documentation and presentation"],
  ["Ali Rizvi", "DNAnexus, integration and the website"],
];

export default function Home() {
  const span = T4.arm.end - T4.arm.start;
  const at = (bp: number) => `${((bp - T4.arm.start) / span) * 100}%`;
  const ruler = [15, 20, 25, 30, 35, 40];
  const maxChunkMin = Math.max(...T4.chunks_n6.map(c => c.wall_min));

  return (<>
    <a className="skip-link" href="#main">Skip to content</a>
    <nav className="topnav" aria-label="Sections">
      <div className="container topnav-inner">
        <a className="brand" href="#top">Chunked pangenome graphs</a>
        <ul>
          <li><a href="#results">Results</a></li>
          <li><a href="#hard">Hard regions</a></li>
          <li><a href="#cut">Where to cut</a></li>
          <li><a href="#explore">Explore</a></li>
          <li><a href="#methods">Methods</a></li>
          <li><a href="#team">Team</a></li>
          <li><a href={GH}>GitHub</a></li>
        </ul>
      </div>
    </nav>

    <header className="hero" id="top">
      <div className="container">
        <p className="locus">GRCh38 <span className="mono">chr21:{n(T4.arm.start)}-{n(T4.arm.end)}</span> <span>{T4.arm.mb} Mb, {S6.n} haplotypes</span></p>
        <h1>Partitioned pangenome graphs, measured against monolithic builds</h1>
        <p className="hero-lede">
          We cut a chromosome into windows on GRCh38 with no overlap, project the cuts onto every haplotype, build each window with PGGB
          on its own machine and join the paths end to end. Then we measure accuracy, memory and cost against one monolithic build.
          On {T4.arm.mb} Mb of chr21 both score F1 {f3(H6.f1_all!)} against the truth set; the partitioned build used {f2(H6.cpu_vs_mono!)} times
          the CPU, and its largest job needed {f2(H6.peak_gb!)} GB instead of {M6.peak_gb} GB.
        </p>
        <div className="hero-actions">
          <a className="btn btn-primary" href="#results">See the results</a>
          <a className="btn" href="#explore">Explore the graph</a>
        </div>

        <div className="tracks" role="img" aria-label={`Genome browser view of chr21 ${mb(T4.arm.start)} to ${mb(T4.arm.end)} Mb. The monolithic build is one job that peaked at ${M6.peak_gb} GB. The ${H6.builds} chunks each peaked at ${H6.peak_gb} GB or less.`}>
          <div className="track">
            <span className="track-name">chr21</span>
            <div className="track-body">
              <div className="ruler">
                {ruler.map(m => <span key={m} style={{ left: at(m * 1e6) }}>{m}</span>)}
                <span style={{ left: "100%" }}>{mb(T4.arm.end)} Mb</span>
              </div>
            </div>
          </div>
          <div className="track">
            <span className="track-name">{H6.builds} chunks</span>
            <div className="track-body">
              <div className="chunk-lane-arm">
                {T4.chunks_n6.map((c, i) => (
                  <span key={c.start} className={i % 2 ? "alt" : ""} style={{ left: at(c.start), width: `calc(${at(c.end + T4.arm.start - c.start)} - 1px)` }}
                    title={`chr21:${n(c.start)}-${n(c.end)}, ${c.wall_min} min, ${c.peak_gb} GB`} />
                ))}
              </div>
            </div>
          </div>
          <div className="track track-signal">
            <span className="track-name">Monolithic<br /><small>peak memory</small></span>
            <div className="track-body">
              <div className="signal">
                <div className="signal-limit"><span>{T4.runner.ram_gb} GB runner</span></div>
                <div className="signal-mono" style={{ height: `${(M6.peak_gb! / T4.runner.ram_gb) * 100}%` }}>
                  <span>one build: {M6.peak_gb} GB, {M6.longest_min} min</span>
                </div>
              </div>
            </div>
          </div>
          <div className="track track-signal">
            <span className="track-name">Chunked<br /><small>peak memory</small></span>
            <div className="track-body">
              <div className="signal">
                <div className="signal-limit" />
                {T4.chunks_n6.map(c => (
                  <span key={c.start} className="signal-chunk" style={{ left: at(c.start), width: `calc(${at(c.end + T4.arm.start - c.start)} - 2px)`, height: `${(c.peak_gb / T4.runner.ram_gb) * 100}%` }} />
                ))}
                <span className="signal-note">{H6.builds} builds: at most {f2(H6.peak_gb!)} GB and {maxChunkMin.toFixed(1)} min each</span>
              </div>
            </div>
          </div>
        </div>
        <dl className="readout">
          <div><dt>F1 vs truth</dt><dd>{f3(M6.f1_all!)} <span>monolithic</span> {f3(H6.f1_all!)} <span>chunked</span></dd></div>
          <div><dt>Seams</dt><dd>{H6.seams_pass}/{H6.seams} <span>passed</span></dd></div>
          <div><dt>Paths identical</dt><dd>{H6.paths_identical}/{H6.paths}</dd></div>
          <div><dt>CPU</dt><dd>{f2(H6.cpu_vs_mono!)}x <span>monolithic</span></dd></div>
        </dl>
      </div>
    </header>

    <main id="main">
      <section className="band band-white" aria-labelledby="how">
        <div className="container">
          <div className="band-head">
            <h2 id="how">How it works</h2>
            <p>PGGB aligns every haplotype against every other, so one build of a long region with many haplotypes needs a large machine.
              Reference-projected partitioning turns it into many small builds that fit on ordinary runners and run side by side.
              The idea is not new: impg&apos;s <code>partition</code> and <code>lace</code> do the same two jobs, and in Tier 2 <code>impg lace</code> gave
              the same graphs as our stitch. What this project adds is the measurement of how partitioned builds differ from monolithic ones.</p>
          </div>
          <ol className="steps-grid">
            <li className="step"><h3>Partition</h3><p>Split GRCh38 into positional windows, fixed or at haplotype-block boundaries, about 1 Mb apart, with zero overlap. Each cut is projected onto every haplotype through its alignment to GRCh38, so neighbouring chunks meet at the same base.</p></li>
            <li className="step"><h3>Build in parallel</h3><p>Run PGGB on each chunk on its own runner ({T4.runner.vcpu} vCPU, {T4.runner.ram_gb} GB). Chunks do not depend on each other.</p></li>
            <li className="step"><h3>Stitch</h3><p>Join the chunk graphs where each haplotype path crosses a cut. With no overlap there is nothing to deduplicate: on the chr21 arm {H6.seams_pass} of {H6.seams} seams passed and {H6.paths_identical} of {H6.paths} paths spell the same sequence as in the monolithic graph.</p></li>
            <li className="step"><h3>Score</h3><p>Compare with a monolithic build of the same input, and score both against variant calls made from the assemblies and, for HG002, against the GIAB Q100 benchmark.</p></li>
          </ol>
        </div>
      </section>

      <section className="band" id="results" aria-labelledby="results-h">
        <div className="container">
          <p className="locus">Tier 4 <span className="mono">chr21:{n(T4.arm.start)}-{n(T4.arm.end)}</span></p>
          <div className="band-head">
            <h2 id="results-h">A chromosome arm: same accuracy, far smaller jobs</h2>
            <p>At {S6.n} haplotypes the {H6.builds} chunks and the monolithic build score the same against the truth (F1 {f3(H6.f1_all!)} for all variants).
              The chunked build used {f2(H6.cpu_vs_mono!)} times the CPU, but no job ran longer than {H6.longest_min} minutes or used more than {f2(H6.peak_gb!)} GB,
              against {M6.longest_min} minutes and {M6.peak_gb} GB for the monolithic build. At {S20.n} and {S40.n} haplotypes the monolithic build did not
              finish on a {T4.runner.ram_gb} GB runner, while the chunks needed at most {f2(S20.hotspots1m.peak_gb!)} and {f2(S40.hotspots1m.peak_gb!)} GB.</p>
          </div>
          <div className="fig-grid">
            <Figure title="Peak memory of the largest job"
              legend={<Legend items={[{ cls: "k-chunked", label: "chunked, largest chunk" }, { cls: "k-mono", label: "monolithic" }, { cls: "k-dnf", label: "monolithic, did not finish" }]} />}
              caption={<>At {S20.n} haplotypes the kernel killed the monolithic build for running out of memory after {S20.mono.after_min} minutes, so it needed more than {T4.runner.ram_gb} GB. At {S40.n} the runner stopped responding after {S40.mono.after_min! / 60} hours inside PGGB and GitHub kept no log, so there is no measurement.</>}
              table={<ScaleTable metric="peak_gb" unit="GB" />}>
              <ScalingPlot xs={[6, 20, 40]} yMax={18} yTicks={[0, 4, 8, 12, 16]} unit="GB" fmt={v => v.toFixed(v >= 10 ? 1 : 2)}
                chunked={T4.scaling.map(s => ({ x: s.n, y: s.hotspots1m.peak_gb! }))}
                mono={[{ x: S6.n, y: M6.peak_gb! }]}
                unfinished={[{ x: S20.n, y: T4.runner.ram_gb, label: "killed: out of memory" }, { x: S40.n, y: 18, label: "runner lost, no data" }]}
                limit={{ y: T4.runner.ram_gb, label: `runner memory, ${T4.runner.ram_gb} GB` }} />
            </Figure>
            <Figure title="Longest single build"
              legend={<Legend items={[{ cls: "k-chunked", label: "chunked, slowest chunk" }, { cls: "k-mono", label: "monolithic" }, { cls: "k-dnf", label: "monolithic, did not finish" }]} />}
              caption={<>Wall time of the slowest job, which is how long a fully parallel build takes. Chunks ran up to {T4.max_concurrent_jobs} at a time. Stitching took {S6.hotspots1m.stitch_s}, {S20.hotspots1m.stitch_s} and {S40.hotspots1m.stitch_s} s.</>}
              table={<ScaleTable metric="longest_min" unit="min" />}>
              <ScalingPlot xs={[6, 20, 40]} yMax={66} yTicks={[0, 15, 30, 45, 60]} unit="min" fmt={v => v.toFixed(1)}
                chunked={T4.scaling.map(s => ({ x: s.n, y: s.hotspots1m.longest_min! }))}
                mono={[{ x: S6.n, y: M6.longest_min! }]}
                unfinished={[{ x: S20.n, y: S20.mono.after_min!, label: `killed after ${S20.mono.after_min} min` }, { x: S40.n, y: 66, label: `lost after ${S40.mono.after_min! / 60} h` }]} />
            </Figure>
          </div>

          <div className="panel table-panel">
            <h3 className="fig-title">Cost and accuracy on the arm</h3>
            <div className="table-scroll">
              <table className="data-table">
                <thead>
                  <tr><th scope="col">Haplotypes</th><th scope="col">Graph</th><th scope="col" className="num">Builds</th><th scope="col" className="num">CPU-h</th>
                    <th scope="col" className="num">CPU-s per Mb per haplotype</th><th scope="col" className="num">F1, all</th><th scope="col" className="num">F1, SVs</th>
                    <th scope="col" className="num">HG002 vs GIAB</th><th scope="col" className="num">Seams passed</th></tr>
                </thead>
                <tbody>
                  {T4.scaling.flatMap(s => {
                    const rows: [string, Build][] = [["monolithic", s.mono], [`${s.hotspots1m.builds} chunks at haplotype-block boundaries`, s.hotspots1m]];
                    if (s.windows1m) rows.push([`${s.windows1m.builds} fixed 1 Mb windows`, s.windows1m]);
                    return rows.map(([label, b]) => (
                      <tr key={`${s.n}${label}`}>
                        <th scope="row">{s.n}</th><td>{label}</td>
                        {b.status && b.status !== "done" ? <td colSpan={7} className="dnf">did not finish: {b.why}, after {b.status === "lost" ? `${b.after_min! / 60} h` : `${b.after_min} min`}</td> : <>
                          <td className="num">{b.builds}</td><td className="num">{f2(b.cpu_h!)}</td><td className="num">{b.cpu_s_per_mb_hap}</td>
                          <td className="num">{f3(b.f1_all!)}</td><td className="num">{f3(b.f1_sv!)}</td><td className="num">{b.hg002_f1!.map(f3).join(" / ")}</td>
                          <td className="num">{b.seams ? `${b.seams_pass}/${b.seams}` : ""}</td></>}
                      </tr>
                    ));
                  })}
                </tbody>
              </table>
            </div>
          </div>
          <div className="prose-grid">
            <p>Cost per Mb per haplotype stays between {S20.hotspots1m.cpu_s_per_mb_hap} and {W6.cpu_s_per_mb_hap} CPU-seconds up to {S20.n} haplotypes ({M6.cpu_s_per_mb_hap} for the monolithic build at {S6.n}) and rises to {S40.hotspots1m.cpu_s_per_mb_hap} at {S40.n}, as PGGB&apos;s all-to-all alignment, which grows with the square of the number of sequences, takes over.</p>
            <p>Against the curated GIAB HG002 Q100 v5.0q benchmark, HG002 scores {f3(M6.hg002_f1![0])} on both haplotypes in the monolithic graph and in both chunked graphs. The SV differences at {S6.n} haplotypes are {T4.arm_sv_diff_calls[0]} to {T4.arm_sv_diff_calls[1]} calls out of {n(T4.arm_truth_svs)} truth SVs and point in no consistent direction.</p>
            <p>F1 falls a little as haplotypes are added ({f3(H6.f1_all!)}, {f3(S20.hotspots1m.f1_all!)}, {f3(S40.hotspots1m.f1_all!)}). With no monolithic graph at {S20.n} or {S40.n} haplotypes, that drop cannot yet be split between chunking and graph size. <a href={doc("TIER4.md")}>Read TIER4.md</a></p>
          </div>
        </div>
      </section>

      <section className="band band-white" id="hard" aria-labelledby="hard-h">
        <div className="container">
          <p className="locus">Tiers 3 and 4 <span className="mono">{r22.locus}</span> <span className="mono">{r17.locus}</span></p>
          <div className="band-head">
            <h2 id="hard-h">Segmental duplications: better in 22q11.2, mixed in 17q21.31</h2>
            <p>With PGGB&apos;s human settings (-p 98) and a fair calling step (vcfbub -a 10000), the partitioned 22q11.2 graph scores F1 {f3(P98_22.chunked.f1_a10k[0])} against {f3(P98_22.mono.f1_a10k[0])} for the
              monolithic graph, and {f3(P98_22.chunked.f1_sd_a10k[0])} against {f3(P98_22.mono.f1_sd_a10k[0])} inside the duplications. Against GIAB, HG002 scores {f3(P98_22.chunked.hg002_f1![0])} in both graphs.
              17q21.31 is mixed: overall F1 is higher partitioned ({f3(P98_17.chunked.f1_a10k[0])} against {f3(P98_17.mono.f1_a10k[0])}), but in Tier 3 at -p 90 partitioning lowered SV F1 from {f2(r17.mono.sv_f1[0])} to {f2(r17.chunked.sv_f1[0])}, on only {T3.r17_truth_svs} truth SVs.</p>
            <p>The likely cause is paralog merging. wfmash aligns copies that are 90% or more identical and seqwish folds them into one node, so a haplotype path runs through the same nodes twice.
              Cuts that put the copies in different chunks stop this from happening.</p>
          </div>
          <div className="fig-grid">
            <Figure title="F1 against the truth, monolithic and chunked"
              legend={<Legend items={[{ cls: "k-chunked", label: "chunked" }, { cls: "k-mono", label: "monolithic" }]} />}
              caption={<>Tier 4: {S6.n} haplotypes (HG00673, HG00733, HG002), all variant sizes, vcfbub -a 10000, both replicates drawn. HG002 is scored against GIAB v5.0q, the others against calls from HPRC&apos;s assembly alignments. The axis starts at 0.6.</>}
              table={<table className="data-table"><thead><tr><th scope="col">Region</th><th scope="col">PGGB</th><th scope="col" className="num">Monolithic, rep 1 / 2</th><th scope="col" className="num">Chunked, rep 1 / 2</th></tr></thead>
                <tbody>{HARD_ORDER.map(h => <tr key={h.region + h.pggb}><th scope="row">{h.region}</th><td>{h.pggb}</td><td className="num">{h.mono.f1_a10k.map(f3).join(" / ")}</td><td className="num">{h.chunked.f1_a10k.map(f3).join(" / ")}</td></tr>)}</tbody></table>}>
              <DotRows lo={0.6} hi={1} ticks={[0.6, 0.7, 0.8, 0.9, 1]} fmt={v => v.toFixed(1)}
                rows={HARD_ORDER.map(h => {
                  const m = Math.min(...h.mono.f1_a10k), c = Math.max(...h.chunked.f1_a10k);
                  return {
                    label: <>{h.region} <span className="mono">{h.pggb}</span></>,
                    span: [m, c] as [number, number],
                    dots: [
                      ...h.mono.f1_a10k.map((v, i) => ({ v, cls: "d-mono", title: `monolithic, rep ${i + 1}: ${f3(v)}`, label: i === 0 ? f3(v) : undefined, side: "left" as const })),
                      ...h.chunked.f1_a10k.map((v, i) => ({ v, cls: "d-chunked", title: `chunked, rep ${i + 1}: ${f3(v)}`, label: i === 0 ? f3(v) : undefined })),
                    ],
                  };
                })} />
            </Figure>
            <Figure title="Share of each haplotype path that revisits a node"
              legend={<Legend items={[{ cls: "k-chunked", label: "chunked" }, { cls: "k-mono", label: "monolithic" }]} />}
              caption={<>Tier 3: 4 haplotypes, PGGB -p 90, length-weighted over all paths, replicate 1 (replicate 2 gives the same shares). A path can revisit a node only if the graph merged two copies of a sequence from the same haplotype.</>}
              table={<table className="data-table"><thead><tr><th scope="col">Region</th><th scope="col" className="num">Monolithic</th><th scope="col" className="num">Chunked</th><th scope="col" className="num">Duplication-aware groups</th></tr></thead>
                <tbody>{T3.regions.map(g => <tr key={g.key}><th scope="row">{g.name}</th><td className="num">{g.mono.revisit_pct.join(" / ")}%</td><td className="num">{g.chunked.revisit_pct.join(" / ")}%</td><td className="num">{g.sdgroup ? `${g.sdgroup.revisit_pct.join(" / ")}%` : ""}</td></tr>)}</tbody></table>}>
              <Bars max={14} unit="%" groups={T3.regions.map(g => ({
                label: g.name,
                bars: [{ v: g.mono.revisit_pct[0], cls: "b-mono", name: "monolithic" }, { v: g.chunked.revisit_pct[0], cls: "b-chunked", name: "chunked" }],
              }))} />
            </Figure>
          </div>
          <div className="prose-grid">
            <p>Tier 3, with 4 haplotypes at -p 90, gave a wider gap in 22q11.2 (F1 {f3(r22.chunked.f1[0])} against {f3(r22.mono.f1[0])}), but part of it is the vcfbub artefact described below:
              rescored with -a 10000 the monolithic graph reaches {f3(r22.mono.f1_a10k!)}, against {f3(r22.chunked.f1_a10k!)} partitioned. Replicates agree to within {f3(repDiff)}.</p>
            <p>Building the chunks that share a duplication in one job does the opposite of what was hoped. That job collapses the paralogs even harder: paths revisit {r22.sdgroup!.revisit_pct[0].toFixed(1)}% in 22q11.2 and {r17.sdgroup!.revisit_pct[0].toFixed(1)}% in 17q21.31, and F1 in 22q11.2 drops to {f3(r22.sdgroup!.f1[0])}.</p>
            <p>At -p 90, HG002&apos;s GIAB F1 in 22q11.2 is {f3(P90_22.chunked.hg002_f1![0])} partitioned against {f3(P90_22.mono.hg002_f1![0])} and {f3(P90_22.mono.hg002_f1![1])} monolithic for the two haplotypes.
              GIAB&apos;s benchmark regions cover about {T4.giab_22q_coverage_pct}% of 22q11.2 and none of 17q21.31, and they leave out much of the LCR22s, where the graphs differ most.</p>
          </div>
          <div className="panel trap">
            <div>
              <h3 className="fig-title">The vcfbub calling trap</h3>
              <p>With PGGB&apos;s recommended human settings (-p 98 -s 10000 -k 47) and the PGGB pipeline&apos;s <code>vcfbub -l 0 -a 100000</code>, the monolithic 22q11.2 VCF lost almost every call: F1 {f3(T3.vcfbub_trap[1].mono_a100k)} in Tier 3 and {f3(T4.hard[1].mono.f1_a100k[0])} in Tier 4.
                The reference path loops through a short top-level bubble that hides hundreds of kb of nested sites, and vcfbub drops them. The graph itself looks fine. Rescored with <code>-a 10000</code> the monolithic graphs recover, but still trail the chunked graphs.
                Before trusting a monolithic VCF of a duplicated region, check whether the reference path revisits nodes.</p>
            </div>
            <div className="table-scroll">
              <table className="data-table">
                <thead><tr><th scope="col">Tier 3, -p 98</th><th scope="col" className="num">Monolithic, -a 100000</th><th scope="col" className="num">Monolithic, -a 10000</th><th scope="col" className="num">Chunked, -a 10000</th></tr></thead>
                <tbody>{T3.vcfbub_trap.map(t => <tr key={t.region}><th scope="row">{t.region}</th><td className="num">{f3(t.mono_a100k)}</td><td className="num">{f3(t.mono_a10k)}</td><td className="num">{f3(t.chunked_a10k)}</td></tr>)}</tbody>
              </table>
            </div>
          </div>
          <p className="links"><a href={doc("TIER3.md")}>Read TIER3.md</a> <a href={doc("TIER4.md")}>Read TIER4.md</a></p>
        </div>
      </section>

      <section className="band" id="cut" aria-labelledby="cut-h">
        <div className="container">
          <p className="locus">Pilot, Tiers 1 and 2 <span className="mono">{r21.locus}</span></p>
          <div className="band-head">
            <h2 id="cut-h">Where to cut matters less than not overlapping</h2>
            <p>The pilot cut 1 Mb of chr21 at recombination hotspots and got a graph closer to the monolithic one: {C.pilot.hotspots_extra_bp} bp of extra node sequence,
              against {n(C.pilot.windows_extra_bp)} bp for fixed windows with {C.pilot.overlap_kb} kb of overlap. Tier 1 asked why. All {n(T1.near_cut_bp.windows[0])} bp sit within 10 kb of a seam,
              in all three replicates. Remove the overlap and the seams carry {T1.near_cut_bp.windows0[0]} bp wherever the cut falls.</p>
          </div>
          <div className="fig-grid fig-grid-wide">
            <Figure title="Agreement with the monolithic graph, against run-to-run noise"
              caption={<>truvari F1 of each chunked graph against a monolithic build of the same input, 1 Mb, 5 paths. Hotspot and fixed cuts: three replicates. Random cuts: one run each, median {T1.random_f1_median} (line). The grey band is the range between two monolithic builds.</>}
              table={<table className="data-table"><thead><tr><th scope="col">Chunking</th><th scope="col" className="num">truvari F1</th><th scope="col" className="num">Extra bp within 10 kb of a cut</th></tr></thead>
                <tbody>
                  <tr><th scope="row">Hotspot cuts</th><td className="num">{T1.f1.hotspots300.join(" / ")}</td><td className="num">{T1.near_cut_bp.hotspots300.join(" / ")}</td></tr>
                  <tr><th scope="row">Fixed cuts, no overlap</th><td className="num">{T1.f1.windows0.join(" / ")}</td><td className="num">{T1.near_cut_bp.windows0.join(" / ")}</td></tr>
                  <tr><th scope="row">Fixed cuts, {C.pilot.overlap_kb} kb overlap</th><td className="num">{T1.f1.windows.join(" / ")}</td><td className="num">{T1.near_cut_bp.windows.map(n).join(" / ")}</td></tr>
                  <tr><th scope="row">20 random cut sets</th><td className="num">{T1.random_f1.join(", ")}</td><td /></tr>
                  <tr><th scope="row">Two monolithic builds</th><td className="num">{T1.mono_vs_mono_f1.join(" / ")}</td><td /></tr>
                </tbody></table>}>
              <DotRows lo={0.9955} hi={1} ticks={[0.996, 0.997, 0.998, 0.999, 1]} fmt={v => v.toFixed(3)}
                rows={[
                  { label: "20 random cut sets", band: [T1.mono_vs_mono_f1[0], 1], mark: T1.random_f1_median,
                    dots: T1.random_f1.map(v => ({ v, cls: "d-noise", title: `random cuts: ${v}` })) },
                  { label: "Hotspot cuts", band: [T1.mono_vs_mono_f1[0], 1], dots: T1.f1.hotspots300.map((v, i) => ({ v, cls: "d-chunked", title: `hotspot cuts, rep ${i + 1}: ${v}`, })) },
                  { label: "Fixed cuts, no overlap", band: [T1.mono_vs_mono_f1[0], 1], dots: T1.f1.windows0.map((v, i) => ({ v, cls: "d-chunked", title: `fixed, no overlap, rep ${i + 1}: ${v}` })) },
                  { label: `Fixed cuts, ${C.pilot.overlap_kb} kb overlap`, band: [T1.mono_vs_mono_f1[0], 1], dots: T1.f1.windows.map((v, i) => ({ v, cls: "d-chunked d-hollow", title: `fixed, ${C.pilot.overlap_kb} kb overlap, rep ${i + 1}: ${v}` })) },
                  { label: "Two monolithic builds", band: [T1.mono_vs_mono_f1[0], 1], dots: T1.mono_vs_mono_f1.map(v => ({ v, cls: "d-mono", title: `monolithic against monolithic: ${v}` })) },
                ]} />
            </Figure>
            <div className="prose-stack">
              <p>Hotspot cuts land at the median of 20 random cut sets (truvari F1 {T1.f1.hotspots300[0]} against {T1.random_f1_median}).
                Fixed zero-overlap windows did better than hotspot cuts in all three replicates, and again at 10, 20 and 40 haplotypes in Tier 2.</p>
              <p>PGGB is not deterministic. Of three monolithic builds of the same input, two were identical and the third differed by {T1.mono_drift_nodes} nodes and {T1.mono_drift_bp} bp.
                So every comparison here is repeated and read against that noise, not against zero.</p>
              <p>At 40 haplotypes on the same 1 Mb, three chunks used about the same CPU as the monolithic build ({n(C.tier2.n40.chunk_cpu_s)} against {n(C.tier2.n40.mono_cpu_s)} s),
                the slowest chunk finished in {C.tier2.n40.longest_vs_mono_pct}% of its wall time, and no chunk used more than {C.tier2.n40.chunk_peak_mb} MB, against {C.tier2.n40.mono_peak_gb} GB.</p>
            </div>
          </div>
          <div className="panel bug">
            <h3 className="fig-title">A seqwish bug, found and fixed</h3>
            <p>The pinned PGGB v0.6.0 image passes its default seqwish <code>-B 10M</code> without converting the suffix. seqwish reads it as 1,000,000,
              and its batched transitive closure then loses merges across batch edges. Rerun on the same 22q11.2 alignments, seqwish made {n(C.tier2.seqwish_bug_bp)} bp of node sequence
              with <code>-B 10M</code> and {n(C.tier2.seqwish_fixed_bp)} bp with <code>-B 30000000</code>. Every runner in this repository now passes a numeric value.</p>
          </div>
          <p className="links"><a href={doc("RECOMBINATION_CHUNKS.md")}>Read RECOMBINATION_CHUNKS.md</a> <a href={doc("TIER1_CONTROLS.md")}>Read TIER1_CONTROLS.md</a> <a href={doc("TIER2.md")}>Read TIER2.md</a></p>
        </div>
      </section>

      <section className="band band-white explore-band" id="explore" aria-labelledby="explore-h">
        <div className="container">
          <p className="locus">Pilot graph <span className="mono">{r21.locus}</span></p>
          <div className="band-head">
            <h2 id="explore-h">Explore the pilot graph</h2>
            <p>The 1 Mb pilot with GRCh38, HG00673 and HG00733, built in one go and from three zero-overlap chunks. Pick a graph and a haplotype:
              green nodes are shared with other haplotypes, blue ones are unique. The viewer draws at most 4,000 nodes around the selected path.</p>
          </div>
          <Explorer />
        </div>
      </section>

      <section className="band" id="methods" aria-labelledby="methods-h">
        <div className="container">
          <div className="band-head"><h2 id="methods-h">Methods and how to reproduce them</h2></div>
          <div className="methods">
            <dl className="method-list">
              <div><dt>Data</dt><dd>HPRC release 2 assemblies and GRCh38 from Ensembl release 112. The contig holding each region is found from the release&apos;s minigraph-cactus chains and fetched with remote <code>samtools faidx</code>, so no whole assembly is downloaded.</dd></div>
              <div><dt>Cuts</dt><dd>Manifests from fixed windows or haploblocks.org block boundaries, projected onto each haplotype through its minimap2 alignment to GRCh38 (longest colinear chain), with zero overlap.</dd></div>
              <div><dt>Graphs</dt><dd>PGGB v0.6.0, pinned image, <code>-p 90 -s 5000 -K 19 -k 29</code> with a numeric seqwish <code>-B</code>, or <code>-p 98 -s 10000 -k 47</code> where marked. CPU time and peak memory from GNU <code>time -v</code>.</dd></div>
              <div><dt>Stitch</dt><dd>Chunk graphs joined where each path crosses a cut. <code>impg lace</code> gives the same graphs.</dd></div>
              <div><dt>Scoring</dt><dd><code>vg deconstruct</code>, <code>vcfbub</code>, <code>vcfwave</code>, <code>bcftools norm</code> and truvari 5.5.0. Truth calls from HPRC&apos;s winnowmap assembly alignments, dipcall style; GIAB v5.0q for HG002. <code>path_revisits.py</code> measures paralog merging.</dd></div>
            </dl>
            <div className="repro">
              <h3>Reproduce</h3>
              <p>Every Tier 4 number comes from a dispatch of <code>.github/workflows/experiment.yml</code>: stage <code>tier4-arm</code> with <code>n</code> = 6, 20 or 40, and <code>tier4-hard</code> with <code>pggb</code> = p90 or p98 and <code>rep</code> = 1 or 2.
                One job runs per chunk, and each run commits a summary JSON to <code>docs/tier4_runs/</code>. The <code>tier3-*</code> stages run the Tier 3 scripts on one runner.</p>
              <p>This page reads its numbers from <code>web/public/data/results.json</code>, which <code>scripts/build_site_results.py</code> builds from the files in <code>docs/</code>. A test fails when the two drift apart.</p>
              <pre><code>python3 scripts/build_site_results.py{"\n"}python3 -m pytest tests/test_site_results.py</code></pre>
              <p><a href={`${GH}/blob/main/.github/workflows/experiment.yml`}>experiment.yml</a> <a href={`${GH}/tree/main/pipeline/benchmark`}>pipeline/benchmark</a></p>
            </div>
          </div>
        </div>
      </section>

      <section className="band band-white" id="team" aria-labelledby="team-h">
        <div className="container team">
          <div>
            <h2 id="team-h">Team</h2>
            <p className="team-note">Built at the BCM SV Hackathon 2026. Paper authors:</p>
            <ul className="authors">
              {AUTHORS.map(([name, role]) => <li key={name}><strong>{name}</strong><span>{role}</span></li>)}
            </ul>
            <p className="team-note team-after">The hackathon team also included Lex, who ran the PGGB cost benchmarks.</p>
          </div>
          <div>
            <h2>Write-ups</h2>
            <ul className="writeups">
              {WRITEUPS.map(([file, what]) => <li key={file}><a href={doc(file)}><span className="mono">{file}</span></a><span>{what}</span></li>)}
              <li><a href={doc("paper/manuscript.md")}><span className="mono">paper/manuscript.md</span></a><span>Manuscript draft</span></li>
            </ul>
          </div>
        </div>
      </section>
    </main>
    <footer className="footer">
      <div className="container"><p>Code, data and write-ups: <a href={GH}>github.com/khoivd-dot/linear-scaling-pangenome-sv</a>. MIT licence.</p></div>
    </footer>
  </>);
}

function ScaleTable({ metric, unit }: { metric: "peak_gb" | "longest_min"; unit: string }) {
  return (
    <table className="data-table">
      <thead><tr><th scope="col">Haplotypes</th><th scope="col" className="num">Monolithic</th><th scope="col" className="num">Chunked, largest job</th></tr></thead>
      <tbody>{T4.scaling.map(s => (
        <tr key={s.n}><th scope="row">{s.n}</th>
          <td className="num">{s.mono.status === "done" ? `${s.mono[metric]} ${unit}` : `did not finish (${s.mono.status})`}</td>
          <td className="num">{s.hotspots1m[metric]} {unit}</td></tr>
      ))}</tbody>
    </table>
  );
}
