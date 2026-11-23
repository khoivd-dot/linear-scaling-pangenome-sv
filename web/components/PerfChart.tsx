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

      const maxMb = 100;
      const xLabels = [0, 10, 25, 50, 100]; // 1 and 5 Mb collided with 0; the 5 Mb limit has its own marker
      const xToPx = (mb: number) => pad.left + (mb / maxMb) * pw;

      const maxTimeS = 360000;
      const yLogMax = Math.log10(maxTimeS);
      const yToPx = (t: number) => pad.top + ph - ((Math.log10(Math.max(t, 1)) / yLogMax) * ph);

      // Baseline: O(L) linear, hits ~5Mb memory wall
      const baselineTimes: [number, number][] = [];
      for (let mb = 0.1; mb <= 100; mb *= 1.15) {
        const t = mb <= 5 ? mb * 100 : Infinity;
        if (t <= 500000) baselineTimes.push([mb, t]);
      }
      baselineTimes.push([5.5, 550]);
      baselineTimes.push([6, 600]);

      // Stitched: O(L/k) wall time, k chunks of ~400Kb
      const stitchedTimes: [number, number][] = [];
      const chunkSizeKb = 400;
      const timePerKb = 0.25;
      const stitchOverhead = 5;
      for (let mb = 0.1; mb <= 100; mb *= 1.15) {
        const kb = mb * 1000;
        const numChunks = Math.max(1, Math.ceil(kb / chunkSizeKb));
        const perChunkTime = timePerKb * chunkSizeKb;
        const wallTime = perChunkTime + stitchOverhead + (numChunks > 1 ? 10 : 0);
        stitchedTimes.push([mb, wallTime]);
      }

      // Canvas background
      ctx.fillStyle = col.paper;
      ctx.fillRect(0, 0, W, H);

      // Grid
      ctx.strokeStyle = col.grid;
      ctx.lineWidth = 1;
      for (let pow = 0; pow <= 5; pow++) {
        const t = Math.pow(10, pow);
        const y = yToPx(t);
        ctx.beginPath();
        ctx.moveTo(pad.left, y);
        ctx.lineTo(W - pad.right, y);
        ctx.stroke();

        ctx.fillStyle = col.muted;
        ctx.font = "12px " + FONT;
        ctx.textAlign = "right";
        const label = t >= 3600 ? `${(t / 3600).toFixed(0)}h`
          : t >= 60 ? `${(t / 60).toFixed(0)}min`
          : `${t.toFixed(0)}s`;
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

      // Baseline line (red, slopes up then hits wall)
      drawLine(baselineTimes, col.baseline, [], 2.5);

      // Stitched line (green, near-flat)
      drawLine(stitchedTimes, col.stitched, [], 2.5);

      // OOM cliff vertical marker at ~5Mb
      ctx.strokeStyle = col.baseline;
      ctx.globalAlpha = 0.45;
      ctx.lineWidth = 1;
      ctx.setLineDash([3, 5]);
      const wallX = xToPx(5);
      ctx.beginPath();
      ctx.moveTo(wallX, pad.top);
      ctx.lineTo(wallX, pad.top + ph);
      ctx.stroke();
      ctx.setLineDash([]);
      ctx.globalAlpha = 1;
      ctx.fillStyle = col.baseline;
      ctx.font = "12px " + FONT;
      ctx.textAlign = "center";
      ctx.textAlign = "left";
      ctx.fillText("single-machine", wallX + 6, pad.top + ph - 20);
      ctx.fillText("memory limit ~5 Mb", wallX + 6, pad.top + ph - 6);

      // Name the baseline line directly so the two series never rely on colour alone.
      const last = baselineTimes[baselineTimes.length - 1];
      ctx.font = "bold 13px " + FONT;
      ctx.fillText("Baseline", xToPx(last[0]) + 8, yToPx(last[1]) + 4);

      // Speedup callout
      const callX = xToPx(50);
      const callY = yToPx(stitchedTimes[stitchedTimes.length - 1][1]) - 60;
      ctx.fillStyle = col.tint;
      ctx.beginPath();
      ctx.roundRect(callX - 60, callY - 10, 230, 34, 6);
      ctx.fill();
      ctx.fillStyle = col.stitchedInk;
      ctx.font = "bold 13px " + FONT;
      ctx.textAlign = "center";
      ctx.fillText("Stitched: near-constant wall time", callX + 55, callY + 6);
      ctx.fillText("regardless of region size", callX + 55, callY + 24);

      // Title
      ctx.fillStyle = col.ink;
      ctx.font = "bold 15px " + FONT;
      ctx.textAlign = "center";
      ctx.fillText("Time Complexity: Baseline vs Parallel Stitched", W / 2, pad.top - 8);

      ctx.fillStyle = col.muted;
      ctx.font = "11px " + FONT;
      ctx.fillText(
        "5 haplotypes, PGGB v0.6.0, all-pairs alignment O(n\u00b2\u00d7L) \u2192 O(L) for n=5, 400 Kb chunks",
        W / 2, H - 8
      );
    }
    return () => { live = false; };
  }, []);

  return (
    <canvas ref={canvasRef} role="img"
      aria-label="Wall time versus region size on a log scale: baseline grows linearly up to a single-machine memory limit near 5 Mb; stitched stays near-constant."
      style={{ display: "block" }} />
  );
}