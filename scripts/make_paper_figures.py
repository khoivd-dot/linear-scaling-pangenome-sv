#!/usr/bin/env python3
"""Render paper figures 1-5 (png 300 dpi + svg) from docs/*.json and the path-revisit TSVs.

Every plotted number is read from those files. Run:
    /home/user/work/venv/bin/python scripts/make_paper_figures.py
"""
import csv
import glob
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
DOCS = ROOT / "docs"
OUT = DOCS / "paper" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

# Okabe-Ito colour-blind-safe palette, fixed entity -> colour
C = {
    "mono": "#222222", "windows": "#0072B2", "windows0": "#0072B2",
    "hotspots300": "#D55E00", "hotspots50": "#CC79A7", "windows1m": "#009E73",
    "windows2500k": "#E69F00", "sdgroup": "#56B4E9", "hotspots1m": "#D55E00",
    "random": "#999999", "rep": "#444444",
}
LAB = {"mono": "monolithic", "windows0": "windows0", "windows": "windows (overlap)",
       "hotspots300": "hotspots300", "hotspots50": "hotspots50", "windows1m": "windows1m",
       "windows2500k": "windows2500k", "sdgroup": "sdgroup", "hotspots1m": "hotspots1m"}

plt.rcParams.update({
    "font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
    "legend.fontsize": 7, "axes.spines.top": False, "axes.spines.right": False,
    "axes.edgecolor": "#555555", "axes.linewidth": 0.7, "svg.fonttype": "none",
    "axes.grid": True, "grid.color": "#e6e6e6", "grid.linewidth": 0.5, "axes.axisbelow": True,
})


def load(name):
    return json.load(open(DOCS / name))


T1, T2, T3 = load("tier1_controls.json"), load("tier2.json"), load("tier3.json")
R3 = T3["runs"]


def letter(ax, s):
    ax.text(-0.02, 1.04, s, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom", ha="right")


def save(fig, n):
    fig.savefig(OUT / f"fig{n}.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / f"fig{n}.svg", bbox_inches="tight")
    plt.close(fig)


def revisit_shares():
    """pooled share per GFA = sum(revisit_bp)/sum(path_bp), all paths incl. GRCh38."""
    out = {}
    for f in ["rep1", "rep2", "p98"]:
        num, den = {}, {}
        with open(DOCS / f"tier3_path_revisits_{f}.tsv") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                num[r["gfa"]] = num.get(r["gfa"], 0) + int(r["revisit_bp"])
                den[r["gfa"]] = den.get(r["gfa"], 0) + int(r["path_bp"])
        for g in num:
            out[g] = num[g] / den[g]
    return out


# ---------------------------------------------------------------- Figure 1
def fig1():
    reps = T1["replicates"]
    fig = plt.figure(figsize=(7.2, 6.4))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1], hspace=0.5, wspace=0.3)
    ax = fig.add_subplot(gs[0, :])
    kinds = ["windows", "windows0", "hotspots300", "hotspots50"]
    w = 0.25
    for i, k in enumerate(kinds):
        for j, r in enumerate(["rep1", "rep2", "rep3"]):
            d = reps[r][k]
            x = i + (j - 1) * (w + 0.03)
            near, els = d["excess_bp_within_10kb_of_cut"], d["excess_bp_elsewhere"]
            ax.bar(x, near, w, color=C["windows0"], label="within 10 kb of a cut" if (i, j) == (0, 0) else None)
            ax.bar(x, els, w, bottom=max(near, 0), color="#E69F00", edgecolor="white", linewidth=0.6,
                   label="elsewhere" if (i, j) == (0, 0) else None)
            ax.text(x, -0.012, str(j + 1), transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=6, color="#555")
            if near > 1000:
                ax.text(x, near + 80, f"{near:,}", ha="center", fontsize=6.5)
    ax.set_xticks(range(4))
    ax.set_xticklabels([LAB[k] for k in kinds])
    ax.tick_params(axis="x", pad=10)
    ax.set_ylabel("Extra node sequence vs monolithic (bp)")
    ax.legend(frameon=False, loc="upper right")
    ax.text(0.0, -0.2, "numbers under bars: replicate", transform=ax.transAxes, fontsize=6.5, color="#555")
    letter(ax, "(a)")

    rnd = [f"random_{i:02d}" for i in range(1, 21)]
    r1 = reps["rep1"]
    mono_f1 = float(np.mean(list(T1["monolithic_vs_monolithic_truvari_f1"].values())))
    marks = [("hotspots300 (rep1)", "hotspots300"), ("windows0 (rep1)", "windows0"),
             ("shift50", "shift50"), ("shift100", "shift100")]
    mcol = ["#D55E00", "#0072B2", "#009E73", "#CC79A7"]
    mshape = ["D", "s", "^", "v"]
    for pidx, (key, xl, ref) in enumerate([("truvari_f1", "truvari F1 vs monolithic", mono_f1),
                                            ("excess_bp", "Extra node sequence vs monolithic (bp)", None)]):
        a = fig.add_subplot(gs[1, pidx])
        vals = np.array([r1[k][key] for k in rnd])
        rng = np.random.default_rng(1)
        a.scatter(vals, rng.uniform(-0.12, 0.12, len(vals)), s=14, color=C["random"], alpha=0.85,
                  edgecolor="white", linewidth=0.4, label="random cut sets (n = 20)", zorder=3)
        for (lab, k), c, m, yy in zip(marks, mcol, mshape, [0.45, 0.3, 0.15, 0.0]):
            a.scatter(r1[k][key], 0.7 + 0.0 * yy - yy * 0.9, s=34, color=c, marker=m, zorder=4, label=lab,
                      edgecolor="white", linewidth=0.5)
        if ref is not None:
            a.axvline(ref, color="#222", ls="--", lw=0.9, label="monolithic vs monolithic (mean)")
        else:
            a.axvline(0, color="#222", ls="--", lw=0.9, label="no extra sequence")
        a.set_yticks([])
        a.set_ylim(-0.4, 1.0)
        a.set_xlabel(xl)
        a.grid(axis="y", visible=False)
        a.spines["left"].set_visible(False)
        letter(a, "(b)" if pidx == 0 else "(c)")
        if pidx == 0:
            a.legend(frameon=False, loc="upper center", bbox_to_anchor=(1.12, -0.28), ncol=3, fontsize=6.8)
    save(fig, 1)


# ---------------------------------------------------------------- Figure 2
def fig2():
    res = T2["map_robustness"]["results"]
    maps = ["HapMapII", "PyrhoYRI", "PyrhoCEU", "PyrhoCHB", "DeCodeSexAveraged"]
    nice = {"HapMapII": "HapMap II", "PyrhoYRI": "Pyrho YRI", "PyrhoCEU": "Pyrho CEU",
            "PyrhoCHB": "Pyrho CHB", "DeCodeSexAveraged": "deCODE"}

    def mapplot(ax, sfx, ylab=True):
        x = np.arange(len(maps))
        obs = [res[f"{m}_{sfx}"]["haploblock_cuts_within_5kb"] * 100 for m in maps]
        rnd = [res[f"{m}_{sfx}"]["random_cuts_within_5kb_mean"] * 100 for m in maps]
        ax.bar(x - 0.19, obs, 0.36, color="#0072B2", label="haploblocks.org cuts")
        ax.bar(x + 0.19, rnd, 0.36, color="#999999", label="random cuts (mean)")
        for i, m in enumerate(maps):
            p = res[f"{m}_{sfx}"]["p_value"]
            ax.text(i - 0.19, obs[i] + 1.5, "p < 0.001" if p < 0.001 else f"p = {p:.3g}", ha="center", fontsize=6.3)
        ax.set_xticks(x)
        ax.set_xticklabels([nice[m] for m in maps], rotation=20, ha="right")
        ax.set_ylabel("Cuts within 5 kb of a hotspot (%)")
        ax.set_ylim(0, 100)
        ax.legend(frameon=False, loc="upper right", fontsize=6.8)

    fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.2), gridspec_kw={"wspace": 0.3})
    mapplot(axs[0], "x5")
    letter(axs[0], "(a)")
    ax = axs[1]
    S = T2["block_compression"]["summary"]
    Ns = [10, 20, 40]
    ax.fill_between(Ns, [S[str(n)]["random_min"] * 100 for n in Ns], [S[str(n)]["random_max"] * 100 for n in Ns],
                    color="#999999", alpha=0.3, lw=0, label="random layouts, min to max")
    ax.plot(Ns, [S[str(n)]["random_median"] * 100 for n in Ns], "-o", color="#666666", ms=4, label="random layouts, median")
    ax.plot(Ns, [S[str(n)]["hotspot_blocks"] * 100 for n in Ns], "-o", color="#D55E00", ms=5, label="hotspot blocks")
    ax.set_xticks(Ns)
    ax.set_xlabel("Haplotypes N")
    ax.set_ylabel("Duplicate fraction (%)")
    ax.set_ylim(0, None)
    ax.legend(frameon=False, loc="upper left")
    letter(ax, "(b)")
    save(fig, 2)
    f2, a2 = plt.subplots(figsize=(3.6, 3.2))
    mapplot(a2, "x10")
    f2.savefig(OUT / "fig2_supp_x10.png", dpi=300, bbox_inches="tight")
    f2.savefig(OUT / "fig2_supp_x10.svg", bbox_inches="tight")
    plt.close(f2)


# ---------------------------------------------------------------- Figure 3
def fig3():
    kinds = ["windows0", "hotspots300", "hotspots50"]
    Ns = [10, 20, 40]
    sc = T2["scaling_rep1"]
    reps1 = T1["replicates"]
    series = {}  # metric -> kind -> list of (N, value, is_rep2)
    metr = {"cpu": ("cpu", "cpu_sum", "cpu_s", "cpu_sum_s", "cpu_s"),
            "wall": ("wall", "max_wall", "wall_s", "max_chunk_wall_s", "max_wall_s"),
            "rss": ("rss", "max_rss", "peak_rss_mb", "max_peak_rss_mb", "peak_rss_mb")}
    pts = {m: {k: [] for k in ["mono"] + kinds} for m in metr}
    for m, (mono2, ch2, mono1, ch1, ch3) in metr.items():
        for r in reps1.values():
            pts[m]["mono"].append((4, r["monolithic"][{"cpu": "cpu_s", "wall": "wall_s", "rss": "peak_rss_mb"}[m]], 0))
        for k in kinds:
            for r in reps1.values():
                if k in r:
                    pts[m][k].append((4, r[k][ch1], 0))
        for n in Ns:
            d = sc[f"N{n}/rep1"]
            pts[m]["mono"].append((n, d["mono"][mono2], 0))
            for k in kinds:
                pts[m][k].append((n, d[k][ch2], 0))
            e = R3[f"scaling/N{n}/rep2"]
            pts[m]["mono"].append((n, e["monolithic"][{"cpu": "cpu_s", "wall": "wall_s", "rss": "peak_rss_mb"}[m]], 1))
            for k in kinds:
                pts[m][k].append((n, e[k][ch3], 1))
    fig = plt.figure(figsize=(7.2, 6.4))
    gs = fig.add_gridspec(2, 3, hspace=0.55, wspace=0.42)
    labs = {"cpu": "CPU time, total (s)", "wall": "Wall time, longest job (s)", "rss": "Peak memory, largest job (MB)"}
    for i, m in enumerate(["cpu", "wall", "rss"]):
        ax = fig.add_subplot(gs[0, i])
        for k in ["mono"] + kinds:
            p = pts[m][k]
            xs = np.array([a for a, b, c in p])
            ys = np.array([b for a, b, c in p])
            rp = np.array([c for a, b, c in p])
            ax.scatter(xs[rp == 0], ys[rp == 0], s=18, color=C[k], zorder=3, label=LAB[k] if i == 0 else None)
            ax.scatter(xs[rp == 1], ys[rp == 1], s=22, facecolor="white", edgecolor=C[k], linewidth=1, zorder=3)
            med = {x: np.mean(ys[(xs == x) & (rp == 0)]) for x in sorted(set(xs))}
            ax.plot(list(med), list(med.values()), color=C[k], lw=1.1, ls="-" if k == "mono" else "--")
        ax.set_xscale("log")
        ax.set_xticks([4, 10, 20, 40])
        ax.set_xticklabels(["4", "10", "20", "40"])
        ax.minorticks_off()
        ax.set_xlabel("Haplotypes N")
        ax.set_ylabel(labs[m])
        if m == "rss":
            ax.set_ylim(0, None)
        letter(ax, "(a)" if i == 0 else "(b)" if i == 1 else "(c)")
    h, l = fig.axes[0].get_legend_handles_labels()
    h += [plt.Line2D([], [], marker="o", ls="", mfc="white", mec="#444", label="second replicate (N = 10, 20, 40)")]
    fig.legend(h, l + [h[-1].get_label()], loc="upper center", bbox_to_anchor=(0.5, 0.01), ncol=5, frameon=False, fontsize=6.8)
    # (d) region length, 10 Mb
    rx = R3["r21x/rep1"]
    ks = ["monolithic", "windows0", "windows1m", "windows2500k"]
    key = {"monolithic": "mono", "windows0": "windows0", "windows1m": "windows1m", "windows2500k": "windows2500k"}
    lab = {"monolithic": "monolithic", "windows0": "25 x 400 kb", "windows1m": "1 Mb windows", "windows2500k": "2.5 Mb windows"}
    cpu0 = rx["monolithic"]["cpu_s"]
    vals = [
        [rx[k]["cpu_s"] / cpu0 for k in ks],
        [(rx[k]["wall_s"] if k == "monolithic" else rx[k]["max_wall_s"]) / 60 for k in ks],
        [rx[k]["peak_rss_mb"] for k in ks],
    ]
    yl = ["CPU time relative to monolithic", "Wall time, longest job (min)", "Peak memory, largest job (MB)"]
    for i in range(3):
        ax = fig.add_subplot(gs[1, i])
        ax.bar(range(4), vals[i], 0.6, color=[C[key[k]] for k in ks])
        ax.set_xticks(range(4))
        ax.set_xticklabels([lab[k] for k in ks], rotation=35, ha="right")
        ax.set_ylabel(yl[i])
        for x, v in enumerate(vals[i]):
            ax.text(x, v, f"{v:.2f}" if i < 2 else f"{v:,.0f}", ha="center", va="bottom", fontsize=6.3)
        ax.set_ylim(0, max(vals[i]) * 1.15)
        letter(ax, "(d)" if i == 0 else "(e)" if i == 1 else "(f)")
    save(fig, 3)


# ---------------------------------------------------------------- Figure 4
def fig4():
    regs = [("r17", "17q21.31"), ("r22", "22q11.2"), ("r21", "chr21 (control)")]
    gk = ["monolithic", "windows0", "sdgroup"]
    gc = {"monolithic": C["mono"], "windows0": C["windows0"], "sdgroup": C["sdgroup"]}
    gl = {"monolithic": "monolithic", "windows0": "windows0", "sdgroup": "sdgroup"}
    fig = plt.figure(figsize=(7.2, 6.8))
    gs = fig.add_gridspec(2, 2, hspace=0.6, wspace=0.28)

    def getv(reg, rep, g, f):
        d = R3[f"{reg}/{rep}"].get(g)
        return None if d is None else f(d)

    # (a) F1 inside SDs and elsewhere
    ax = fig.add_subplot(gs[0, 0])
    groups = [(r, part) for r, _ in regs for part in ["sd", "other"]]
    pos, ticks = 0, []
    for ri, (reg, rl) in enumerate(regs):
        for part in ["sd", "other"]:
            if reg == "r21" and part == "sd":
                continue  # chr21 control window has no SD truth calls
            for gi, g in enumerate(gk):
                x = pos + (gi - 1) * 0.25
                v1 = getv(reg, "rep1", g, lambda d: d["truth"]["all"][part]["f1"])
                if v1 is None:
                    continue
                v2 = getv(reg, "rep2", g, lambda d: d["truth"]["all"][part]["f1"])
                ax.bar(x, v1, 0.22, color=gc[g], label=gl[g] if (ri, part, gi) == (0, "sd", 0) or (ri == 0 and part == "sd") and gi else None)
                if v2 is not None:
                    ax.scatter(x, v2, s=14, facecolor="white", edgecolor="#222", linewidth=0.9, zorder=4)
            ticks.append((pos, "in SD" if part == "sd" else "elsewhere"))
            pos += 1
        pos += 0.5
    ax.set_xticks([t[0] for t in ticks])
    ax.set_xticklabels([t[1] for t in ticks], rotation=35, ha="right")
    for i, (reg, rl) in enumerate(regs):
        ax.text(i * 2.5 + 0.5, 1.02, rl, transform=ax.get_xaxis_transform(), ha="center", fontsize=7)
    ax.set_ylabel("Truth F1, all variants")
    ax.set_ylim(0, 1.0)
    h, l = ax.get_legend_handles_labels()
    seen = dict(zip(l, h))
    seen["replicate 2"] = plt.Line2D([], [], marker="o", ls="", mfc="white", mec="#222")
    fig.legend(seen.values(), seen.keys(), frameon=False, loc="upper center", bbox_to_anchor=(0.5, 0.01), fontsize=7, ncol=4)
    letter(ax, "(a)")
    ax.text(-0.02, 1.16, "", transform=ax.transAxes)

    # (b) SV F1 with truth SV count
    ax = fig.add_subplot(gs[0, 1])
    for ri, (reg, rl) in enumerate(regs):
        for gi, g in enumerate(gk):
            x = ri + (gi - 1) * 0.27
            v1 = getv(reg, "rep1", g, lambda d: d["truth"]["sv"]["total"]["f1"])
            if v1 is None:
                continue
            ax.bar(x, v1, 0.24, color=gc[g])
            v2 = getv(reg, "rep2", g, lambda d: d["truth"]["sv"]["total"]["f1"])
            if v2 is not None:
                ax.scatter(x, v2, s=14, facecolor="white", edgecolor="#222", linewidth=0.9, zorder=4)
        tot = R3[f"{reg}/rep1"]["monolithic"]["truth"]["sv"]["total"]
        ax.text(ri, 1.02, f"{rl}\n{tot['tp_base'] + tot['fn']} truth SVs", transform=ax.get_xaxis_transform(), ha="center", fontsize=7, va="bottom")
    ax.set_xticks([])
    ax.set_ylabel("Truth F1, SVs (50 bp or more)")
    ax.set_ylim(0, 1.0)
    letter(ax, "(b)")

    # (c) path revisit share from TSV
    sh = revisit_shares()
    ax = fig.add_subplot(gs[1, 0])
    tsv = {"monolithic": "monolithic/final.gfa", "windows0": "windows0/merged.gfa", "sdgroup": "sdgroup/laced.gfa"}
    for ri, (reg, rl) in enumerate(regs):
        for gi, g in enumerate(gk):
            x = ri + (gi - 1) * 0.27
            k1 = f"runs/{reg}/rep1/{tsv[g]}"
            if k1 not in sh:
                continue
            ax.bar(x, sh[k1] * 100, 0.24, color=gc[g])
            k2 = f"runs/{reg}/rep2/{tsv[g]}"
            if k2 in sh:
                ax.scatter(x, sh[k2] * 100, s=14, facecolor="white", edgecolor="#222", linewidth=0.9, zorder=4)
    ax.set_xticks(range(3))
    ax.set_xticklabels([r[1] for r in regs])
    ax.set_ylabel("Path revisit share (% of path length)")
    letter(ax, "(c)")

    # (d) -p 90 vs -p 98, both vcfbub cuts
    ax = fig.add_subplot(gs[1, 1])
    combos = []
    for reg, rl in regs[:2]:
        for g in ["monolithic", "windows0"]:
            for cut, tk in [("-a 100k", "truth"), ("-a 10k", "truth_a10k")]:
                combos.append((reg, rl, g, cut, tk))
    xs, tick = [], []
    pos = 0
    for reg, rl in regs[:2]:
        for g in ["monolithic", "windows0"]:
            for ci, (cut, tk) in enumerate([("-a 100k", "truth"), ("-a 10k", "truth_a10k")]):
                v90 = R3[f"{reg}/rep1"][g][tk]["all"]["total"]["f1"]
                v98 = R3[f"p98/{reg}/rep1"][g][tk]["all"]["total"]["f1"]
                ax.plot([pos, pos], [v90, v98], color="#bbbbbb", lw=1, zorder=2)
                ax.scatter(pos, v90, s=26, color="#0072B2", zorder=3, label="-p 90" if pos == 0 else None)
                ax.scatter(pos, v98, s=26, color="#D55E00", marker="s", zorder=3, label="-p 98" if pos == 0 else None)
                tick.append((pos, f"{gl[g]}\nvcfbub {cut}"))
                pos += 1
            pos += 0.5
        pos += 0.6
    ax.set_xticks([t[0] for t in tick])
    ax.set_xticklabels([t[1] for t in tick], rotation=60, ha="right", fontsize=6.3)
    ax.set_ylabel("Truth F1, all variants")
    ax.set_ylim(0, 1.0)
    ax.legend(frameon=False, loc="lower right")
    ax.text(0.25, 1.02, regs[0][1], transform=ax.transAxes, fontsize=7, ha="center")
    ax.text(0.75, 1.02, regs[1][1], transform=ax.transAxes, fontsize=7, ha="center")
    letter(ax, "(d)")
    save(fig, 4)


# ---------------------------------------------------------------- Figure 5
def fig5():
    runs = {}
    for f in glob.glob(str(DOCS / "tier4_runs" / "tier4_arm_n*_p90_rep1.json")):
        d = json.load(open(f))
        runs[d["n_haplotypes"]] = d
    Ns = sorted(runs)
    MB = 31.0  # chr21:15-46 Mb (TIER4.md: CPU divided by 31 Mb and by N)
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.3), gridspec_kw={"wspace": 0.3})
    ax = axs[0]
    for g, st in [("arm.hotspots1m", "hotspots1m"), ("arm.windows1m", "windows1m")]:
        xs = [n for n in Ns if g in runs[n]["graphs"]]
        ys = [runs[n]["graphs"][g]["cost"]["cpu_s"] / MB / n for n in xs]
        ax.plot(xs, ys, "-o", color=C[st], ms=5, label=LAB[st] + " (chunked)")
    m = runs[6]["graphs"]["arm.mono"]["cost"]
    ax.scatter([6], [m["cpu_s"] / MB / 6], s=45, color=C["mono"], marker="D", zorder=4, label="monolithic")
    for n in Ns:
        mg = runs[n]["graphs"]["arm.mono"]
        if mg["missing"]:
            pass
    ax.set_xticks(Ns)
    ax.set_xlabel("Haplotypes N")
    ax.set_ylabel("CPU-seconds per Mb per haplotype")
    ax.set_ylim(0, None)
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.45, 1.0))
    ax.text(0.98, 0.04, "monolithic: no result at N = 20, 40", transform=ax.transAxes, ha="right", fontsize=6.5, color="#555")
    letter(ax, "(a)")
    ax = axs[1]
    for g, st in [("arm.hotspots1m", "hotspots1m"), ("arm.windows1m", "windows1m")]:
        xs = [n for n in Ns if g in runs[n]["graphs"]]
        ys = [runs[n]["graphs"][g]["cost"]["max_peak_rss_mb"] / 1000 for n in xs]
        ax.plot(xs, ys, "-o", color=C[st], ms=5, label=LAB[st] + " (chunked)")
    ax.scatter([6], [m["max_peak_rss_mb"] / 1000], s=45, color=C["mono"], marker="D", zorder=4, label="monolithic")
    ax.axhline(16, color="#B22222", lw=1.2, ls="--")
    ax.text(40, 16.3, "16 GB runner limit", ha="right", va="bottom", fontsize=7, color="#B22222")
    for n in Ns:
        if runs[n]["graphs"]["arm.mono"]["missing"]:
            ax.scatter([n], [16], marker="x", s=40, color=C["mono"], zorder=5,
                       label="monolithic, no result" if n == 20 else None)
    ax.set_xticks(Ns)
    ax.set_xlabel("Haplotypes N")
    ax.set_ylabel("Peak memory, largest job (GB)")
    ax.set_ylim(0, 19.5)
    ax.legend(frameon=False, loc="center right", fontsize=6.8, bbox_to_anchor=(1.0, 0.55))
    letter(ax, "(b)")
    save(fig, 5)


if __name__ == "__main__":
    fig1(); fig2(); fig3(); fig4(); fig5()
    print("wrote", sorted(p.name for p in OUT.iterdir()))
