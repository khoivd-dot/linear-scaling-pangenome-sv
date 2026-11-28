"use client";
import { useState, useEffect, useCallback, type KeyboardEvent } from "react";
import type { SampleGraph, NodeInfo, EdgeInfo } from "../types";
import GraphExplorer from "./GraphExplorer";

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
interface OverviewData {
  chunks?: { chunk_id: string; reference_start: number; reference_end: number;
    core_start: number; core_end: number; overlap_left: number; overlap_right: number }[];
}

const FALLBACK_SAMPLES: { sample: string; haplotypes: string[]; hap_labels?: Record<string, string> }[] = [
  { sample: "GRCh38", haplotypes: ["0"], hap_labels: { "0": "reference" } },
  { sample: "HG00673", haplotypes: ["1", "2"], hap_labels: { "1": "paternal", "2": "maternal" } },
  { sample: "HG00733", haplotypes: ["1", "2"], hap_labels: { "1": "paternal", "2": "maternal" } },
];

// The pilot graph (chr21:20-21 Mb, 5 haplotypes): explorer, two-person compare, chunk layout and graph counts.
export default function Explorer() {
  const [loading, setLoading] = useState(true);
  const [sg, setSG] = useState<SampleGraph | null>(null);
  const [nodeInfo, setNodeInfo] = useState<NodeInfo | null>(null);
  const [edgeInfo, setEdgeInfo] = useState<EdgeInfo | null>(null);
  const [tab, setTab] = useState("explore");
  const [manifest, setManifest] = useState<ManifestData | null>(null);
  const [latest, setLatest] = useState<LatestData | null>(null);
  const [overview, setOverview] = useState<OverviewData | null>(null);
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
    ]).then(([m, l, o]) => {
      if (m.status === "fulfilled") {
        setManifest(m.value);
        const s = m.value.samples;
        if (s && s.length) { setSelSample(s[0].sample); setSelHap(s[0].haplotypes?.[0] || "0"); }
      }
      if (l.status === "fulfilled") setLatest(l.value);
      if (o.status === "fulfilled") setOverview(o.value);
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
    if (["explore", "dashboard", "chunks", "compare"].includes(h)) setTab(h);
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

  if (loading) return <p className="loading">Loading the pilot graph…</p>;

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
      <h3>Graph numbers</h3>
      <p className="section-lede">
        The whole region built in one go (baseline) against the three chunks stitched back together, same five haplotypes.
        PGGB is not deterministic, so two baseline builds also differ by a few nodes. The viewer draws at most 4,000 nodes
        around the selected path; these counts are for the full graphs.
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
        <h3>Chunks</h3>
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
        <h3>Compare two people</h3>
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

  const tabs = [
    { key: "explore", label: "Explore", hint: "look at the graph" },
    { key: "compare", label: "Compare", hint: "two people side by side" },
    { key: "chunks", label: "Chunks", hint: "the parallel pieces" },
    { key: "dashboard", label: "Graph numbers", hint: "full-graph counts" },
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

  return (<>
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
      {tab === "chunks" && renderChunks()}
      {tab === "compare" && renderCompare()}
    </div>
  </>);
}
