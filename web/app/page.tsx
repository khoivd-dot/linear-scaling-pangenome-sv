"use client";
import { useState, useEffect, useCallback, type KeyboardEvent } from "react";
import type { SampleGraph, NodeInfo, EdgeInfo } from "../types";
import GraphExplorer from "../components/GraphExplorer";
import PerfChart from "../components/PerfChart";

interface ManifestData {
  data_mode?: string;
  run?: { run_id?: string; data_mode?: string; git_sha?: string; scientific_status?: string };
  pipeline_status?: Record<string, string>;
  samples?: { sample: string; haplotypes: string[]; hap_labels?: Record<string, string> }[];
  graphs?: Record<string, { nodes: number; edges: number; paths: number; walks: number }>;
}
interface LatestData {
  data_mode: string;
  metrics?: Record<string, Record<string, number>>;
  stitch?: { status?: string };
  equivalence?: { verdict?: string };
  boundaries?: { boundary?: string; status?: string; haplotypes_joined?: number }[];
}
interface Finding {
  step: string; question: string; finding: string;
  status: "done" | "running"; writeup: { label: string; url: string };
}
interface OverviewData {
  chunks?: { chunk_id: string; reference_start: number; reference_end: number;
    core_start: number; core_end: number; overlap_left: number; overlap_right: number }[];
}

const FALLBACK_SAMPLES: { sample: string; haplotypes: string[]; hap_labels?: Record<string, string> }[] = [
  { sample: "GRCh38", haplotypes: ["0"], hap_labels: { "0": "reference" } },
  { sample: "HG00673", haplotypes: ["1", "2"], hap_labels: { "1": "paternal", "2": "maternal" } },
  { sample: "HG00733", haplotypes: ["1", "2"], hap_labels: { "1": "paternal", "2": "maternal" } },
];

export default function Home() {
  const [loading, setLoading] = useState(true);
  const [sg, setSG] = useState<SampleGraph | null>(null);
  const [nodeInfo, setNodeInfo] = useState<NodeInfo | null>(null);
  const [edgeInfo, setEdgeInfo] = useState<EdgeInfo | null>(null);
  const [tab, setTab] = useState("explore");
  const [manifest, setManifest] = useState<ManifestData | null>(null);
  const [latest, setLatest] = useState<LatestData | null>(null);
  const [overview, setOverview] = useState<OverviewData | null>(null);
  const [findings, setFindings] = useState<Finding[] | null>(null);
  const [graphMode, setGraphMode] = useState("baseline");
  const [selSample, setSelSample] = useState("GRCh38");
  const [selHap, setSelHap] = useState("0");
  const [selSampleB, setSelSampleB] = useState("HG00673");
  const [selHapB, setSelHapB] = useState("1");
  const [sgB, setSGB] = useState<SampleGraph | null>(null);

  const basePath = process.env.NEXT_PUBLIC_BASE_PATH ?? "";

  useEffect(() => {
    Promise.allSettled([
      fetch(`${basePath}/data/manifest.json`).then(r => (r.ok ? r.json() : Promise.reject())),
      fetch(`${basePath}/data/latest.json`).then(r => (r.ok ? r.json() : Promise.reject())),
      fetch(`${basePath}/data/overview.json`).then(r => (r.ok ? r.json() : Promise.reject())),
      fetch(`${basePath}/data/findings.json`).then(r => (r.ok ? r.json() : Promise.reject())),
    ]).then(([m, l, o, f]) => {
      if (m.status === "fulfilled") {
        setManifest(m.value);
        const s = m.value.samples;
        if (s && s.length) { setSelSample(s[0].sample); setSelHap(s[0].haplotypes?.[0] || "0"); }
      }
      if (l.status === "fulfilled") setLatest(l.value);
      if (o.status === "fulfilled") setOverview(o.value);
      if (f.status === "fulfilled") setFindings(f.value.findings);
      setLoading(false);
    });
  }, []);

  const loadGraph = useCallback(() => {
    fetch(`${basePath}/data/graphs/${graphMode}/${selSample}_${selHap}.json`)
      .then(r => (r.ok ? r.json() : Promise.reject()))
      .then((g: SampleGraph) => { setSG(g); setLoading(false); })
      .catch(() => { setSG(null); setLoading(false); });
  }, [basePath, graphMode, selSample, selHap]);

  useEffect(() => {
    const h = window.location.hash.slice(1);
    if (["explore", "dashboard", "findings", "chunks", "compare", "performance"].includes(h)) setTab(h);
  }, []);

  useEffect(() => { loadGraph(); }, [loadGraph]);
  const loadGraphB = useCallback(() => {
    fetch(`${basePath}/data/graphs/${graphMode}/${selSampleB}_${selHapB}.json`)
      .then(r => (r.ok ? r.json() : Promise.reject()))
      .then((g: SampleGraph) => setSGB(g))
      .catch(() => setSGB(null));
  }, [basePath, graphMode, selSampleB, selHapB]);
  useEffect(() => { loadGraphB(); }, [loadGraphB]);

  const samples: { sample: string; haplotypes: string[]; hap_labels?: Record<string, string> }[] =
    manifest?.samples?.length ? manifest.samples : FALLBACK_SAMPLES;
  const hapLabels: Record<string, string> = samples.find(s => s.sample === selSample)?.hap_labels || {};
  const graphNames = manifest?.graphs && Object.keys(manifest.graphs).length
    ? Object.keys(manifest.graphs) : ["baseline", "stitched"];
  const bm = latest?.metrics?.baseline || {};
  const mm = latest?.metrics?.merged || {};

  if (loading) return <div className="container loading">Loading…</div>;

  const renderExplore = () => (
    <div className="explorer-wrapper">
      <div className="sidebar">
        <div className="sidebar-section">
          <h3>What this graph shows</h3>
          <div className="what-graph">
            <p>
              A <strong>pangenome graph</strong> of a ~1 Mb slice of chromosome 21.
              Each colored circle is a piece of DNA. Lines show how pieces connect.
            </p>
            <div className="key-row"><span className="key-dot shared" aria-hidden="true" />Shared: same across people</div>
            <div className="key-row"><span className="key-dot unique" aria-hidden="true" />Unique: specific to this person</div>
            <div className="key-row"><span className="key-dot other" aria-hidden="true" />Nearby pieces</div>
            <div className="key-why">
              <strong>Why it matters:</strong> a pangenome represents many people at once.
              Where the graph <strong>forks</strong>, people differ. Where it <strong>merges</strong>, they match again.
            </div>
          </div>
        </div>

        <div className="sidebar-section">
          <h3>1. Choose a graph</h3>
          <p className="sidebar-note">
            Baseline = built all at once (comparison standard). Stitched = parallel chunks → merged.
          </p>
          {graphNames.map(name => (
            <button key={name} onClick={() => setGraphMode(name)} aria-pressed={graphMode === name}
              className="graph-mode-btn">
              {name === "stitched" ? "Stitched (parallel)" : name === "baseline" ? "Baseline" : name}
            </button>
          ))}
          <div className="graph-mode-explainer">
            {graphMode === "stitched" ? (
              <>
                <div className="gm-row"><span className="gm-chip" /> Cut the region into three chunks with no overlap.</div>
                <div className="gm-row"><span className="gm-chip" /> Build a graph for each chunk in parallel.</div>
                <div className="gm-row"><span className="gm-chip" /> Join each person&apos;s path across the cuts.</div>
              </>
            ) : (
              <>
                <div className="gm-row"><span className="gm-chip" /> Build the whole region at once with PGGB.</div>
                <div className="gm-row"><span className="gm-chip" /> Works for small regions, slow for whole genomes.</div>
                <div className="gm-row gm-fast"><span className="gm-chip" /> Comparison standard (not biological truth).</div>
              </>
            )}
          </div>
        </div>
        <div className="sidebar-section">
          <h3>2. Choose a person</h3>
          <p className="sidebar-note">
            Green = shared with others. Blue = unique to them.
          </p>
          {samples.map(s => (
            <div key={s.sample}>
              <button onClick={() => { setSelSample(s.sample); setSelHap(s.haplotypes?.[0] || "0"); }}
                aria-pressed={selSample === s.sample}
                className="sample-btn">{s.sample}</button>
              {selSample === s.sample && (
                <div className="hap-list">
                  {s.haplotypes.map(h => (
                    <button key={h} onClick={() => setSelHap(h)} aria-pressed={selHap === h}
                      className="hap-btn">{hapLabels[h] || h}</button>
                  ))}
                </div>
              )}
            </div>
          ))}
        </div>
      </div>
      <div>
        <GraphExplorer graph={sg}
          onNodeSelect={(i) => { setNodeInfo(i); setEdgeInfo(null); }}
          onEdgeSelect={(e) => { setEdgeInfo(e); setNodeInfo(null); }} />
        {sg && sg.nodes.length <= 3 && graphMode === "baseline" && (
          <div className="card" style={{ marginTop: 12, fontSize: 15, color: "var(--ink-2)" }}>
            <strong>Heads up:</strong> this view shows only {sg.nodes.length} node{sg.nodes.length===1?"":"s"} because
            {selSample === "GRCh38"
              ? " GRCh38 is stored as one large reference piece in this graph. Try HG00673 or HG00733 to see a richly branched path with shared and unique regions."
              : " this person's DNA is very similar to the others here — few differences to show."}
          </div>
        )}
      </div>
      <div className="inspector" aria-live="polite">
        <h3>Inspector</h3>
        {nodeInfo ? <div>
          <div className="inspector-row"><span>DNA piece ID</span><span>{nodeInfo.id}</span></div>
          <div className="inspector-row"><span>Length</span><span>{(nodeInfo.length/1000).toFixed(1)} Kb</span></div>
          <div className="inspector-row"><span>Connections</span><span>{nodeInfo.degree}</span></div>
          <div className="inspector-row"><span>On this person</span><span>{nodeInfo.on_selected_path ? "Yes" : "No"}</span></div>
          <div className="inspector-row"><span>Shared with others</span><span>{nodeInfo.is_shared ? "Yes" : "No"}</span></div>
        </div> : edgeInfo ? <div>
          <div className="inspector-row"><span>Connection</span><span>{edgeInfo.source} → {edgeInfo.target}</span></div>
          <div className="inspector-row"><span>Direction</span><span>{edgeInfo.source_orientation}→{edgeInfo.target_orientation}</span></div>
        </div> : <p className="inspector-empty">Click any circle or line to see details.</p>}
      </div>
    </div>
  );

  const metricRows: [string, string][] = [["Nodes", "nodes"], ["Edges", "edges"], ["Paths", "paths"], ["Node sequence (bp)", "sequence_bp"]];

  const renderDashboard = () => (
    <section className="section">
      <h2>Results</h2>
      <p className="section-lede">
        The whole region built in one go (baseline) against the three chunks stitched back together, same five haplotypes.
        PGGB is not deterministic, so two baseline builds also differ slightly; see Findings.
      </p>
      <div className="card results-table">
        <table className="data-table">
          <thead>
            <tr><th scope="col">Graph metric</th><th scope="col" className="num">Baseline</th><th scope="col" className="num">Stitched</th></tr>
          </thead>
          <tbody>
            {metricRows.filter(([, k]) => bm[k] != null || mm[k] != null).map(([label, k]) => (
              <tr key={k}><th scope="row">{label}</th><td className="num">{bm[k]?.toLocaleString("en-US") ?? "–"}</td><td className="num">{mm[k]?.toLocaleString("en-US") ?? "–"}</td></tr>
            ))}
          </tbody>
        </table>
      </div>
    </section>
  );

  const renderChunks = () => {
    const chunks = overview?.chunks || [];
    const lo = Math.min(...chunks.map(c => c.reference_start));
    const hi = Math.max(...chunks.map(c => c.reference_end));
    const pct = (bp: number) => `${((bp - lo) / Math.max(hi - lo, 1)) * 100}%`;
    const span = (a: number, b: number) => ({ left: pct(a), width: `${((b - a) / Math.max(hi - lo, 1)) * 100}%` });
    const hasOverlap = chunks.some(c => c.overlap_left > 0 || c.overlap_right > 0);
    const bp = (n: number) => n.toLocaleString("en-US");
    return (
      <section className="section">
        <h2>Chunks</h2>
        {chunks.length > 1 && (
          <p className="section-lede">
            {chunks.length} chunks{hasOverlap ? " with overlaps" : ", no overlap"}. Cuts at chr21:{chunks.slice(1).map(c => bp(c.reference_start)).join(" and ")}.
            Each chunk was built on its own, then each path was joined across the cuts.
          </p>
        )}
        {chunks.length > 0 && (
          <div className="card">
            {hasOverlap && (
              <div className="chunk-legend" aria-hidden="true">
                <span><span className="swatch swatch-span" /> chunk with overlaps</span>
                <span><span className="swatch swatch-core" /> core</span>
              </div>
            )}
            {chunks.map(c => (
              <div key={c.chunk_id} className="chunk-row">
                <strong className="chunk-id">{c.chunk_id}</strong>
                <div>
                  <div className="chunk-track" aria-hidden="true">
                    <div className="chunk-span" style={span(c.reference_start, c.reference_end)} />
                    <div className="chunk-core" style={span(c.core_start, c.core_end)} />
                  </div>
                  <span className="chunk-coords">
                    chr21:{bp(c.reference_start)}–{bp(c.reference_end)}{hasOverlap && ` (core ${bp(c.core_start)}–${bp(c.core_end)})`}
                  </span>
                </div>
              </div>
            ))}
          </div>
        )}
        {(latest?.boundaries || []).length > 0 && (
          <div className="card">
            {(latest?.boundaries || []).map(b => (
              <p key={b.boundary} className="boundary-row"><span>{b.boundary}: {b.status}</span><span>joined={b.haplotypes_joined}</span></p>
            ))}
          </div>
        )}
      </section>
    );
  };

  const renderCompare = () => {
    const hapLabelsB: Record<string, string> = samples.find(s => s.sample === selSampleB)?.hap_labels || {};
    const idA = new Set((sg?.nodes || []).map(n => n.id));
    const idB = new Set((sgB?.nodes || []).map(n => n.id));
    const shared = sg && sgB ? Array.from(idA).filter(id => idB.has(id)).length : null;
    const onlyA = sg && sgB ? Array.from(idA).filter(id => !idB.has(id)).length : null;
    const onlyB = sg && sgB ? Array.from(idB).filter(id => !idA.has(id)).length : null;
    const sameGraphSpace = sg && sgB && sg.graph === sgB.graph;
    const sharePct = shared != null && sg && sgB ? Math.round(shared / Math.max(idA.size, idB.size, 1) * 100) : null;

    return (
      <section className="section">
        <h2>Compare two people</h2>
        <p className="section-lede">
          Pick two people and see how much of their DNA is shared vs. unique.
        </p>
        <div className="compare-layout">
          <div className="compare-panel">
            <h3>Sample A</h3>
            <select aria-label="Sample A" className="select" value={selSample} onChange={e => { setSelSample(e.target.value); setSelHap(samples.find(s => s.sample === e.target.value)?.haplotypes[0] || "0"); }}>
              {samples.map(s => <option key={s.sample} value={s.sample}>{s.sample}</option>)}
            </select>
            <div className="hap-row">
              {(samples.find(s => s.sample === selSample)?.haplotypes || []).map(h => (
                <button key={h} onClick={() => setSelHap(h)} aria-pressed={selHap === h}
                  className="hap-btn">{hapLabels[h] || h}</button>
              ))}
            </div>
            {sg && <div>
              <div className="inspector-row"><span>Path length</span><span>{(sg.path?.length_bp ?? 0).toLocaleString()} bp</span></div>
              <div className="inspector-row"><span>Nodes</span><span>{sg.nodes?.length ?? 0}</span></div>
              <div className="inspector-row"><span>Edges</span><span>{sg.edges?.length ?? 0}</span></div>
            </div>}
          </div>
          <div className="compare-panel">
            <h3>Sample B</h3>
            <select aria-label="Sample B" className="select" value={selSampleB} onChange={e => { setSelSampleB(e.target.value); setSelHapB(samples.find(s => s.sample === e.target.value)?.haplotypes[0] || "0"); }}>
              {samples.map(s => <option key={s.sample} value={s.sample}>{s.sample}</option>)}
            </select>
            <div className="hap-row">
              {(samples.find(s => s.sample === selSampleB)?.haplotypes || []).map(h => (
                <button key={h} onClick={() => setSelHapB(h)} aria-pressed={selHapB === h}
                  className="hap-btn">{hapLabelsB[h] || h}</button>
              ))}
            </div>
            {sgB && <div>
              <div className="inspector-row"><span>Path length</span><span>{(sgB.path?.length_bp ?? 0).toLocaleString()} bp</span></div>
              <div className="inspector-row"><span>Nodes</span><span>{sgB.nodes?.length ?? 0}</span></div>
              <div className="inspector-row"><span>Edges</span><span>{sgB.edges?.length ?? 0}</span></div>
            </div>}
          </div>
        </div>
        {sg && sgB && (
          <div className="card compare-diff" style={{ marginTop: 20 }}>
            <h3>How they differ</h3>
            <p className="sidebar-note" style={{ marginBottom: 14 }}>
              Each bar = the DNA pieces that make up that person. Green = pieces they share.
            </p>
            {!sameGraphSpace && (
              <p className="callout-warn" style={{ marginBottom: 14 }}>
                These two were built in separate graphs, so the numbers are approximate.
              </p>
            )}
            <div style={{ marginBottom: 16 }}>
              <div className="share-head">
                <span>Shared</span><span className="mono">{shared ?? "—"} pieces ({sharePct ?? "—"}%)</span>
              </div>
              <div className="share-bar"><div className="share-bar-fill" style={{ transform: `scaleX(${(sharePct ?? 0) / 100})` }} /></div>
            </div>
            <div className="inspector-row"><span>Shared DNA pieces</span><span>{shared ?? "—"}</span></div>
            <div className="inspector-row"><span>Only in {selSample}</span><span>{onlyA ?? "—"}</span></div>
            <div className="inspector-row"><span>Only in {selSampleB}</span><span>{onlyB ?? "—"}</span></div>
          </div>
        )}
      </section>
    );
  };

  const renderFindings = () => (
    <section className="section">
      <h2>Findings</h2>
      <p className="section-lede">What each experiment asked and what it found. Each row links to its write-up on GitHub.</p>
      {findings?.length ? (
        <div className="card findings-card">
          <table className="data-table findings-table">
            <thead>
              <tr><th scope="col">Step</th><th scope="col">Question</th><th scope="col">Finding</th><th scope="col">Status</th></tr>
            </thead>
            <tbody>
              {findings.map(f => (
                <tr key={f.question}>
                  <th scope="row">{f.step}</th>
                  <td className="f-question">{f.question}</td>
                  <td>{f.finding} <a className="f-link" href={f.writeup.url}>{f.writeup.label}</a></td>
                  <td><span className={"status status-" + f.status}>{f.status}</span></td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      ) : <p className="inspector-empty">Findings could not be loaded.</p>}
    </section>
  );

  const renderPerformance = () => (
    <section className="section">
      <h2>Performance scaling</h2>
      <p className="section-lede">
        How long a monolithic build takes against the longest chunk of the same region, measured at three sizes.
      </p>
      <div className="chart-frame">
        <PerfChart />
      </div>
      <div className="perf-wrap">
        <div className="card">
          <h3>Measured cost</h3>
          <div className="table-scroll">
            <table className="data-table">
              <thead>
                <tr>
                  <th scope="col">Region (chr21)</th>
                  <th scope="col" className="c-baseline">Monolithic</th>
                  <th scope="col" className="c-stitched">Chunked (longest chunk)</th>
                </tr>
              </thead>
              <tbody>
                <tr><th scope="row">1 Mb, 5 haplotypes</th><td>122 s, about 970 MB</td><td>78 s, about 970 MB; 1.3x the CPU</td></tr>
                <tr><th scope="row">10 Mb, 5 haplotypes</th><td>582 s, 2.8 GB</td><td>109 s, 971 MB; 1.22x the CPU</td></tr>
                <tr><th scope="row">31 Mb, 6 haplotypes</th><td>48.8 min, 12.1 GB</td><td>3.4 min, 0.97 GB; 1.05x the CPU</td></tr>
                <tr><th scope="row">31 Mb, 20 haplotypes</th><td className="c-baseline">ran out of 16 GB</td><td className="c-stitched">8.2 min, 1.5 GB</td></tr>
              </tbody>
            </table>
          </div>
        </div>
        <p className="section-lede" style={{ marginTop: 12 }}>
          Chunks run on separate workers, so wall time is set by the longest chunk. Sources:
          TIER1_CONTROLS.md (1 Mb, median of 3 runs), TIER3.md E4 (10 Mb, 1 Mb chunks) and TIER4.md E6
          (31 Mb, 29 hotspot chunks on GitHub runners) in the repository&apos;s docs folder.
        </p>
      </div>
    </section>
  );

  const tabs = [
    { key: "explore", label: "Explore", hint: "look at the graph" },
    { key: "dashboard", label: "Results", hint: "graph numbers" },
    { key: "findings", label: "Findings", hint: "what the experiments found" },
    { key: "chunks", label: "Chunks", hint: "the parallel pieces" },
    { key: "compare", label: "Compare", hint: "two people side by side" },
    { key: "performance", label: "Speed", hint: "measured build cost" },
  ];

  const selectTab = (key: string) => {
    setTab(key);
    history.replaceState(null, "", `#${key}`);
  };
  // WAI-ARIA tabs: arrow keys, Home and End move between tabs.
  const onTabKey = (e: KeyboardEvent) => {
    const i = tabs.findIndex(t => t.key === tab);
    const next = { ArrowRight: i + 1, ArrowLeft: i - 1, Home: 0, End: tabs.length - 1 }[e.key];
    if (next === undefined) return;
    e.preventDefault();
    const key = tabs[(next + tabs.length) % tabs.length].key;
    selectTab(key);
    document.getElementById(`tab-${key}`)?.focus();
  };

  // The benchmark region, drawn as a genome-browser ruler with the real chunk layout under it.
  const REGION = [20_000_000, 21_000_000];
  const at = (bp: number) => `${((bp - REGION[0]) / (REGION[1] - REGION[0])) * 100}%`;
  const bar = (a: number, b: number) => ({ left: at(a), width: `${((b - a) / (REGION[1] - REGION[0])) * 100}%` });
  const heroChunks = overview?.chunks || [];

  return (<>
      <a className="skip-link" href="#explore-heading">Skip to the graph explorer</a>
      <header className="hero">
        <div className="container">
          <p className="locus">GRCh38 <span className="mono">chr21:20,000,000–21,000,000</span></p>
          <h1>Can We Build Pangenome Graphs in Parallel?</h1>
          <p className="hero-lede">Normally, building a pangenome graph requires comparing every genome against every other, which is very slow. We are testing a new approach: <strong>build small pieces in parallel</strong>, then stitch them back together.</p>
          <a className="hero-action" href="#explore-heading">Explore the graph</a>
          <div className="tracks" aria-hidden="true">
            <div className="track">
              <span className="track-name">chr21</span>
              <div className="track-body">
                <div className="ruler">
                  {[20, 20.2, 20.4, 20.6, 20.8, 21].map(mb => (
                    <span key={mb} style={{ left: at(mb * 1e6) }}>{mb.toFixed(1)}{mb === 21 ? " Mb" : ""}</span>
                  ))}
                </div>
              </div>
            </div>
            {heroChunks.length > 0 && (
              <div className="track">
                <span className="track-name">chunks</span>
                <div className="track-body chunk-lane">
                  {heroChunks.map(c => (
                    <div key={c.chunk_id}>
                      <div className="chunk-span" style={bar(c.reference_start, c.reference_end)} />
                      <div className="chunk-core" style={bar(c.core_start, c.core_end)} />
                    </div>
                  ))}
                </div>
              </div>
            )}
            <div className="track">
              <span className="track-name">graph</span>
              <div className="track-body"><HeroMotif /></div>
            </div>
          </div>
          <dl className="readout">
            <div><dt>Region</dt><dd>1 Mb <span>real HPRC chr21</span></dd></div>
            <div><dt>Monolithic</dt><dd>2.0 min</dd></div>
            <div><dt>3 chunks in parallel</dt><dd>1.6× sooner</dd></div>
          </dl>
          <p className="hero-note">
            <strong>Benchmark region:</strong> GRCh38 <span className="mono">chr21:20,000,000–21,000,000</span> (1 Mb, 5 haplotypes).
            The graph viewer draws at most 4,000 nodes around the selected path so it stays fast;
            the Results table counts the full graphs.
          </p>
        </div>
      </header>

      <main>
      <section className="band band-white">
        <div className="container">
          <div className="band-head"><h2>What are we actually doing?</h2></div>
          <div className="concepts">
            <div className="concept">
              <h3>Genomes</h3>
              <p>Each person has unique DNA. We compare 5 versions: one reference + 4 from real people.</p>
            </div>
            <div className="concept">
              <h3>Pangenome graph</h3>
              <p>Instead of comparing to one reference, a graph shows where genomes share the same route and where they take different paths.</p>
            </div>
            <div className="concept">
              <h3>Our idea</h3>
              <p>Split the genome into pieces, build each piece in parallel, then stitch them back together.</p>
            </div>
          </div>
        </div>
      </section>

      <section className="band">
        <div className="container">
          <div className="band-head">
            <h2>How it works, in five steps</h2>
          </div>
          <ol className="steps-grid">
            <PlainStep title="Pick a region" desc="We focus on a tiny slice of chromosome 21: just 1 million DNA letters. Small enough to test fast, big enough to matter." />
            <PlainStep title="Line everyone up" desc="We find where each person's DNA matches our reference, like aligning everyone to the same starting line." />
            <PlainStep title="Split and build" desc="We cut the region into three chunks of about 333 Kb, with no overlap, and build a little graph for each. All chunks run at the same time." />
            <PlainStep title="Stitch together" desc="Every person's path crosses each cut, so we join the chunk graphs where the paths cross, like connecting puzzle pieces." />
            <PlainStep title="Check our work" desc="We compare against building the whole thing in one go (no splitting). The Findings tab lists what each round of tests showed." />
          </ol>
        </div>
      </section>

      <section className="container explore-band" aria-labelledby="explore-heading">
        <div className="band-head">
          <h2 id="explore-heading" tabIndex={-1}>Explore the graph yourself</h2>
          <p>Use the tabs below to dig in.</p>
        </div>
        <div className="tabs" role="tablist" aria-label="Views" onKeyDown={onTabKey}>
          {tabs.map(t => (
            <button key={t.key} id={`tab-${t.key}`} role="tab" aria-selected={tab === t.key} aria-controls="tab-panel"
              tabIndex={tab === t.key ? 0 : -1} onClick={() => selectTab(t.key)} title={t.hint} className="tab">
              <span className="tab-label">{t.label}</span>
              <span className="tab-hint">{t.hint}</span>
            </button>
          ))}
        </div>
        <div id="tab-panel" role="tabpanel" aria-labelledby={`tab-${tab}`}>
          {tab === "explore" && renderExplore()}
          {tab === "dashboard" && renderDashboard()}
          {tab === "findings" && renderFindings()}
          {tab === "chunks" && renderChunks()}
          {tab === "compare" && renderCompare()}
          {tab === "performance" && renderPerformance()}
        </div>
      </section>
      </main>
  </>);
}

function PlainStep({ title, desc }: { title: string; desc: string }) {
  return (
    <li className="step">
      <h3>{title}</h3>
      <p>{desc}</p>
    </li>
  );
}

// Decorative sketch: a shared backbone (green) with a branch where one haplotype differs (blue).
// Drawn in percentages so it spans the track at any width without stretching the nodes.
function HeroMotif() {
  return (
    <svg className="hero-motif" aria-hidden="true" focusable="false">
      <svg x="36%" width="28%" height="64" viewBox="0 0 100 64" preserveAspectRatio="none">
        <path className="motif-unique" d="M0 46 C 25 46, 25 14, 50 14 S 75 46, 100 46" strokeWidth="3" vectorEffect="non-scaling-stroke" />
      </svg>
      <line className="motif-shared" x1="1%" y1="46" x2="99%" y2="46" strokeWidth="4" strokeLinecap="round" />
      {[2, 19, 36, 64, 81, 98].map(x => <circle key={x} className="motif-node-shared" cx={`${x}%`} cy="46" r="6" strokeWidth="2" />)}
      {[43, 50, 57].map(x => <circle key={x} className="motif-node-unique" cx={`${x}%`} cy={x === 50 ? 14 : 30} r="5.5" strokeWidth="2" />)}
      <text className="motif-label-shared" x="20%" y="32">shared</text>
      <text className="motif-label-unique" x="59%" y="14">unique</text>
    </svg>
  );
}
