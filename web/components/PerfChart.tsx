"use client";
import { useEffect, useRef } from "react";
import { token } from "./token";

export default function PerfChart() {
  const canvasRef = useRef<HTMLCanvasElement>(null);

  useEffect(() => {
    let live = true;
    // Wait for the self-hosted face so canvas labels match the page.
    document.fonts.ready.then(() => { if (live) draw(); });
    function draw() {
      const c = canvasRef.current;
      if (!c) return;
      const ctx = c.getContext("2d");
      if (!ctx) return;

      // Draw at device resolution; keep a 600px minimum so labels stay legible (the frame scrolls on phones).
      const W = Math.max(600, Math.min(740, (c.parentElement?.clientWidth || 720) - 24));
      const H = 420;
      const dpr = Math.min(window.devicePixelRatio || 1, 2);
      c.width = W * dpr; c.height = H * dpr;
      c.style.width = `${W}px`; c.style.height = `${H}px`;
      ctx.scale(dpr, dpr);
      const FONT = getComputedStyle(c).fontFamily;
      const col = { paper: token("--surface"), grid: token("--rule-soft"), tick: token("--rule-strong"), ink: token("--ink"),
        muted: token("--muted"), baseline: token("--baseline"), stitched: token("--shared"), stitchedInk: token("--shared-ink"), tint: token("--shared-tint") };
      const pad = { top: 30, right: 30, bottom: 76, left: 70 };
      const pw = W - pad.left - pad.right;
      const ph = H - pad.top - pad.bottom;

      // Measured, not modelled: monolithic wall time against the longest chunk (chunks run in parallel).
      // 1 Mb: chr21:20-21 Mb, 5 haplotypes, median of 3 runs, 3 zero-overlap chunks (docs/TIER1_CONTROLS.md).
      // 10 Mb: chr21:20-30 Mb, 5 haplotypes, 10 x 1 Mb chunks (docs/TIER3.md, E4).
      // 31 Mb: chr21:15-46 Mb, 6 haplotypes, 29 hotspot chunks of ~1 Mb, GitHub runner (docs/TIER4.md, E6).
      const maxMb = 35;
      const xLabels = [1, 10, 20, 31];
      const xToPx = (mb: number) => pad.left + (mb / maxMb) * pw;

      const maxTimeS = 10000;
      const yLogMax = Math.log10(maxTimeS);
      const yToPx = (t: number) => pad.top + ph - ((Math.log10(Math.max(t, 1)) / yLogMax) * ph);

      const baselineTimes: [number, number][] = [[1, 122], [10, 582], [31, 2928]];
      const stitchedTimes: [number, number][] = [[1, 78], [10, 109], [31, 204]];

      // Canvas background
      ctx.fillStyle = col.paper;
      ctx.fillRect(0, 0, W, H);

      // Grid
      ctx.strokeStyle = col.grid;
      ctx.lineWidth = 1;
      // Gridlines at round clock times so the labels are exact.
      for (const [t, label] of [[1, "1 s"], [10, "10 s"], [60, "1 min"], [600, "10 min"], [3600, "1 h"]] as [number, string][]) {
        const y = yToPx(t);
        ctx.beginPath();
        ctx.moveTo(pad.left, y);
        ctx.lineTo(W - pad.right, y);
        ctx.stroke();

        ctx.fillStyle = col.muted;
        ctx.font = "12px " + FONT;
        ctx.textAlign = "right";
        ctx.fillText(label, pad.left - 8, y + 4);
      }

      // Y-axis label
      ctx.save();
      ctx.translate(16, pad.top + ph / 2);
      ctx.rotate(-Math.PI / 2);
      ctx.fillStyle = col.muted;
      ctx.font = "13px " + FONT;
      ctx.textAlign = "center";
      ctx.fillText("Wall time (log scale)", 0, 0);
      ctx.restore();

      // X-axis labels
      ctx.fillStyle = col.muted;
      ctx.font = "12px " + FONT;
      ctx.textAlign = "center";
      for (const mb of xLabels) {
        const x = xToPx(mb);
        ctx.fillText(mb === 0 ? "0" : `${mb}Mb`, x, pad.top + ph + 18);
        ctx.beginPath();
        ctx.moveTo(x, pad.top + ph);
        ctx.lineTo(x, pad.top + ph + 5);
        ctx.strokeStyle = col.tick;
        ctx.stroke();
      }
      ctx.fillText("Region size", pad.left + pw / 2, pad.top + ph + 38);

      // Draw lines helper
      const drawLine = (points: [number, number][], color: string, dash: number[], width: number) => {
        if (points.length < 2) return;
        ctx.save();
        ctx.strokeStyle = color;
        ctx.lineWidth = width;
        ctx.setLineDash(dash);
        ctx.lineJoin = "round";
        ctx.beginPath();
        const first = points[0];
        ctx.moveTo(xToPx(first[0]), yToPx(first[1]));
        let cliff = false;
        for (let i = 1; i < points.length; i++) {
          const pt = points[i];
          if (pt[1] >= 500000) {
            cliff = true;
            ctx.stroke();
            ctx.setLineDash([6, 8]);
            ctx.globalAlpha = 0.5;
            ctx.beginPath();
            const prev = points[i - 1];
            ctx.moveTo(xToPx(prev[0]), yToPx(prev[1]));
            ctx.lineTo(xToPx(pt[0]), pad.top + 2);
            ctx.stroke();
            break;
          }
          ctx.lineTo(xToPx(pt[0]), yToPx(pt[1]));
        }
        if (!cliff) ctx.stroke();
        ctx.restore();
      };

      // Monolithic wall time
      drawLine(baselineTimes, col.baseline, [], 2.5);

      // Stitched line (green, near-flat)
      drawLine(stitchedTimes, col.stitched, [], 2.5);

      // Point markers, so three measurements never read as a modelled curve
      for (const [pts, color] of [[baselineTimes, col.baseline], [stitchedTimes, col.stitched]] as [[number, number][], string][]) {
        ctx.fillStyle = color;
        for (const [mb, t] of pts) {
          ctx.beginPath();
          ctx.arc(xToPx(mb), yToPx(t), 4.5, 0, Math.PI * 2);
          ctx.fill();
        }
      }
      ctx.textAlign = "left";
      ctx.fillStyle = col.baseline;

      // Name the baseline line directly so the two series never rely on colour alone.
      const last = baselineTimes[baselineTimes.length - 1];
      ctx.font = "bold 13px " + FONT;
      ctx.fillText("Monolithic", xToPx(last[0]) - 92, yToPx(last[1]) - 10);
      const lastS = stitchedTimes[stitchedTimes.length - 1];
      ctx.fillStyle = col.stitchedInk;
      ctx.fillText("Longest chunk", xToPx(lastS[0]) - 110, yToPx(lastS[1]) + 22);

      // Speedup callout
      const callX = xToPx(12);
      const callY = yToPx(12); // empty band below both series
      ctx.fillStyle = col.tint;
      ctx.beginPath();
      ctx.roundRect(callX - 10, callY - 10, 300, 34, 6);
      ctx.fill();
      ctx.fillStyle = col.stitchedInk;
      ctx.font = "bold 13px " + FONT;
      ctx.textAlign = "left";
      ctx.fillText("Chunks in parallel finish 1.6x, 5.3x and 14x", callX, callY + 6);
      ctx.fillText("sooner at 1, 10 and 31 Mb", callX, callY + 24);

      // Title
      ctx.fillStyle = col.ink;
      ctx.font = "bold 15px " + FONT;
      ctx.textAlign = "center";
      ctx.fillText("Wall time: monolithic build vs longest chunk (measured)", W / 2, pad.top - 8);

      ctx.fillStyle = col.muted;
      ctx.font = "11px " + FONT;
      ctx.fillText(
        "PGGB commit 4225c6c (pinned image), chr21; 5 haplotypes at 1 and 10 Mb, 6 at 31 Mb. At 20 haplotypes the 31 Mb monolithic build ran out of 16 GB",
        W / 2, H - 8
      );
    }
    return () => { live = false; };
  }, []);

  return (
    <canvas ref={canvasRef} role="img"
      aria-label="Measured wall time on a log scale at 1, 10 and 31 Mb of chr21: the monolithic build takes 122 s, 582 s and 2,928 s; the longest chunk takes 78 s, 109 s and 204 s."
      style={{ display: "block" }} />
  );
}