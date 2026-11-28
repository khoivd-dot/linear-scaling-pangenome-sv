// Charts drawn from results.json as HTML and CSS, so labels keep their size at any width
// and the page reads without JavaScript. Colour: --chunked for chunked graphs, --mono for
// monolithic, --noise for controls. Each mark carries a title for hover; each chart has a table.
import type { CSSProperties, ReactNode } from "react";

const pct = (v: number, lo: number, hi: number) => `${((v - lo) / (hi - lo)) * 100}%`;

export function Legend({ items }: { items: { cls: string; label: string }[] }) {
  return (
    <ul className="legend">
      {items.map(i => <li key={i.label}><span className={`key ${i.cls}`} aria-hidden="true" />{i.label}</li>)}
    </ul>
  );
}

export function Figure({ title, caption, legend, children, table }: {
  title: string; caption?: ReactNode; legend?: ReactNode; children: ReactNode; table: ReactNode;
}) {
  return (
    <figure className="panel fig">
      <h3 className="fig-title">{title}</h3>
      {legend}
      <div className="fig-plot" aria-hidden="true">{children}</div>
      {caption && <figcaption>{caption}</figcaption>}
      <details className="fig-data">
        <summary>Show the numbers</summary>
        <div className="table-scroll">{table}</div>
      </details>
    </figure>
  );
}

// ---------- Line chart over haplotype count, with builds that did not finish marked ----------
export interface LinePoint { x: number; y: number }
export interface Unfinished { x: number; y: number; label: string }
export function ScalingPlot({ xs, yMax, yTicks, unit, chunked, mono, unfinished, limit, fmt = String }: {
  xs: number[]; yMax: number; yTicks: number[]; unit: string; fmt?: (v: number) => string;
  chunked: LinePoint[]; mono: LinePoint[]; unfinished: Unfinished[];
  limit?: { y: number; label: string };
}) {
  const X = (x: number) => pct(x, xs[0] - 3, xs[xs.length - 1] + 3);
  const Y = (y: number) => pct(y, 0, yMax);
  const num = (x: number) => (x - (xs[0] - 3)) / (xs[xs.length - 1] - xs[0] + 6) * 100;
  const end = (x: number) => (x === xs[xs.length - 1] ? " is-end" : "");
  return (
    <div className="lc">
      <div className="lc-plot">
        {yTicks.map(t => (
          <div key={t} className="lc-grid" style={{ bottom: Y(t) }}><span>{t}{t === yTicks[yTicks.length - 1] ? ` ${unit}` : ""}</span></div>
        ))}
        {limit && <div className="lc-limit" style={{ bottom: Y(limit.y) }}><span>{limit.label}</span></div>}
        <svg className="lc-lines" viewBox="0 0 100 100" preserveAspectRatio="none">
          <polyline className="line-chunked" vectorEffect="non-scaling-stroke"
            points={chunked.map(p => `${num(p.x)},${100 - (p.y / yMax) * 100}`).join(" ")} />
        </svg>
        {chunked.map(p => (
          <span key={`c${p.x}`} className={`pt pt-chunked${end(p.x)}`} style={{ left: X(p.x), bottom: Y(p.y) }} title={`chunked, ${p.x} haplotypes: ${fmt(p.y)} ${unit}`}>
            <span className="pt-label">{fmt(p.y)} {unit}</span>
          </span>
        ))}
        {mono.map(p => (
          <span key={`m${p.x}`} className={`pt pt-mono${end(p.x)}`} style={{ left: X(p.x), bottom: Y(p.y) }} title={`monolithic, ${p.x} haplotypes: ${fmt(p.y)} ${unit}`}>
            <span className="pt-label">{fmt(p.y)} {unit}</span>
          </span>
        ))}
        {unfinished.map(u => (
          <span key={`u${u.x}`} className={`pt pt-dnf${end(u.x)}`} style={{ left: X(u.x), bottom: Y(u.y) }} title={`monolithic, ${u.x} haplotypes: ${u.label}`}>
            <span className="pt-label">{u.label}</span>
          </span>
        ))}
      </div>
      <div className="lc-x">
        {xs.map(x => <span key={x} style={{ left: X(x) }}>{x}</span>)}
      </div>
      <p className="lc-xtitle">haplotypes in the graph</p>
    </div>
  );
}

// ---------- Dot rows on a shared axis (dumbbells, strips) ----------
export interface Dot { v: number; cls: string; title: string; label?: string; side?: "left" | "right" }
export function DotRows({ lo, hi, ticks, rows, fmt }: {
  lo: number; hi: number; ticks: number[]; fmt: (v: number) => string;
  rows: { label: ReactNode; dots: Dot[]; span?: [number, number]; band?: [number, number]; mark?: number }[];
}) {
  return (
    <div className="dots">
      {rows.map((row, i) => (
        <div className="dots-row" key={i}>
          <div className="dots-label">{row.label}</div>
          <div className="dots-track">
            {row.band && <span className="dots-band" style={{ left: pct(row.band[0], lo, hi), width: pct(row.band[1] - row.band[0] + lo, lo, hi) }} />}
            {row.span && <span className="dots-span" style={{ left: pct(row.span[0], lo, hi), width: pct(row.span[1] - row.span[0] + lo, lo, hi) }} />}
            {row.mark !== undefined && <span className="dots-mark" style={{ left: pct(row.mark, lo, hi) }} />}
            {row.dots.map((d, j) => (
              <span key={j} className={`dot ${d.cls}`} style={{ left: pct(d.v, lo, hi) }} title={d.title}>
                {d.label && <span className={`dot-label ${d.side === "left" ? "is-left" : ""}`}>{d.label}</span>}
              </span>
            ))}
          </div>
        </div>
      ))}
      <div className="dots-row dots-axis">
        <div className="dots-label" />
        <div className="dots-track">
          {ticks.map(t => <span key={t} style={{ left: pct(t, lo, hi) }}>{fmt(t)}</span>)}
        </div>
      </div>
    </div>
  );
}

// ---------- Horizontal bars ----------
export function Bars({ max, groups, unit }: {
  max: number; unit: string;
  groups: { label: string; bars: { v: number; cls: string; name: string }[] }[];
}) {
  return (
    <div className="bars">
      {groups.map(g => (
        <div className="bars-group" key={g.label}>
          <div className="bars-label">{g.label}</div>
          <div>
            {g.bars.map(b => (
              <div className="bars-row" key={b.name} title={`${g.label}, ${b.name}: ${b.v.toFixed(1)}${unit}`}>
                <span className="bars-name">{b.name}</span>
                <span className="bars-track">
                  <span className={`bar ${b.cls}`} style={{ width: pct(b.v, 0, max) } as CSSProperties} />
                  <span className="bars-value">{b.v.toFixed(1)}{unit}</span>
                </span>
              </div>
            ))}
          </div>
        </div>
      ))}
    </div>
  );
}
