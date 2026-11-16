#!/usr/bin/env python3
"""生成正文四图 + 两张补充图。全部只用存档数据，不重算模型。"""
import glob

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy import stats

D = "b_data"
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8, "axes.linewidth": 0.8,
    "axes.spines.top": False, "axes.spines.right": False,
    "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "pdf.fonttype": 42, "savefig.bbox": "tight", "savefig.dpi": 300,
})
C = {"matched": "#111111", "stim": "#C0392B", "unstim": "#2C6FA6",
     "fit": "#999999", "hl": "#E67E22"}


def panel(ax, letter):
    # above the axes frame, left-aligned: at (-0.16, 1.06) the letter collides with the
    # top y tick label of narrow panels
    ax.text(0.0, 1.03, letter, transform=ax.transAxes, fontsize=10, fontweight="bold",
            va="bottom", ha="left")


# ---------------------------------------------------------------- Figure 1
def fig1():
    R = pd.read_csv(f"{D}/calibration_final.csv")
    fig, ax = plt.subplots(figsize=(4.4, 3.8))
    for kind, mk, lab in [("unstim", "o", "unstimulated"), ("stim", "^", "stimulated"),
                          ("matched", "*", "matched (IFN-beta)")]:
        d = R[R.kind == kind]
        ax.errorbar(d.measured, d.model, xerr=[d.measured - d.lo, d.hi - d.measured],
                    fmt=mk, ms=7 if kind == "matched" else 5, mfc=C[kind], mec=C[kind],
                    ecolor=C[kind], elinewidth=0.7, capsize=1.6, alpha=0.85,
                    label=lab, zorder=3, linestyle="none")
    sl, ic = np.polyfit(R.measured, R.model, 1)
    xs = np.linspace(-0.24, 1.05, 50)
    ax.plot(xs, sl * xs + ic, "--", color=C["fit"], lw=1, zorder=1,
            label=f"OLS: slope {sl:.2f}")
    ax.axhline(0, color="0.86", lw=0.6, zorder=0)
    ax.axvline(0, color="0.86", lw=0.6, zorder=0)

    # --- 标注 1：三个 T/NK 点各自加圈，虚线框住，文字放下方空白
    tnk = R[(R.kind == "stim") & (R.cond.isin(["CD4 T cells", "CD8 T cells", "NK cells"]))]
    x0, x1 = tnk.measured.min() - 0.026, tnk.measured.max() + 0.026
    y0, y1 = tnk.model.min() - 0.034, tnk.model.max() + 0.034
    ax.add_patch(FancyBboxPatch((x0, y0), x1 - x0, y1 - y0,
                                boxstyle="round,pad=0.006,rounding_size=0.012",
                                fill=False, edgecolor=C["stim"], lw=0.7,
                                ls=(0, (2.2, 1.6)), zorder=2))
    ax.plot([(x0 + x1) / 2, 0.178], [y0, -0.20], color=C["stim"], lw=0.6, zorder=2)
    ax.text(0.295, -0.235,
            "T/NK stimulated, 3 conditions:\n"
            "systematic underestimation (p = 0.029)",
            fontsize=6.0, color=C["stim"], ha="center", va="top", linespacing=1.35)

    # --- 标注 2：CD8 / B 同一实测重叠、模型差 0.23
    c8 = R[(R.cond == "CD8 T cells") & (R.kind == "stim")].iloc[0]
    bb = R[(R.cond == "B cells") & (R.kind == "stim")].iloc[0]
    ax.annotate("", xy=(c8.measured, c8.model), xytext=(bb.measured, bb.model),
                arrowprops=dict(arrowstyle="<->", color=C["hl"], lw=1.1,
                                shrinkA=4, shrinkB=4), zorder=5)
    ax.plot([bb.measured + 0.012, 0.335], [0.115, 0.155], color=C["hl"], lw=0.6, zorder=2)
    ax.text(0.34, 0.155, "delta = 0.23 at the same\nmeasured overlap",
            fontsize=6.0, color=C["hl"], ha="left", va="center", linespacing=1.35)

    # --- 标注 3：巨核，文字置左上空白，引线从点向上接
    mg = R[(R.cond == "Megakaryocytes") & (R.kind == "stim")].iloc[0]
    ax.scatter([mg.measured], [mg.model], s=120, facecolors="none",
               edgecolors="0.30", lw=0.9, zorder=4)
    ax.plot([mg.measured - 0.008, -0.045], [mg.model + 0.016, 0.545],
            color="0.45", lw=0.6, zorder=2)
    ax.text(-0.245, 0.665, "megakaryocyte:\nlargest residual (0.15),\n134 cells (widest CI)",
            fontsize=6.0, color="0.30", ha="left", va="top", linespacing=1.35)

    ax.set_xlabel("measured overlap with the IFN-$\\beta$ response (Spearman $\\rho$)",
                  fontsize=7.5)
    ax.set_ylabel("model overlap (Spearman $\\rho$)", fontsize=7.5)
    ax.tick_params(labelsize=7)
    ax.set_xlim(-0.30, 1.10)
    ax.set_ylim(-0.52, 0.72)
    ax.legend(loc="lower right", fontsize=6.2, frameon=False, handletextpad=0.5,
              borderaxespad=0.2, labelspacing=0.35)
    fig.tight_layout(pad=0.3)
    fig.savefig("fig1_calibration.png", dpi=600)
    fig.savefig("fig1_calibration.pdf")
    plt.close(fig)


# ---------------------------------------------------------------- Figure 2
def fig2():
    p = pd.read_csv(f"{D}/b_pred_disjoint_1115_115447.csv")
    base = p[(p[["ctrl", "stim"]].max(axis=1) > 1) & p.obs.notna()]
    pred = np.log2((base.stim + 1) / (base.ctrl + 1))
    r = stats.spearmanr(pred, base.obs).correlation
    ho = pd.read_csv(f"{D}/holdout_pred_1115_121546.csv").set_index("ensg").loc[base.ensg]
    hl = np.log2((ho.stim + 1) / (ho.ctrl + 1))
    r_ho = stats.spearmanr(hl, base.obs).correlation
    r_two = stats.spearmanr(hl, pred).correlation

    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.6),
                             gridspec_kw={"width_ratios": [1.5, 0.7, 0.7]})
    ax = axes[0]
    ax.scatter(base.obs, pred, s=10, c=C["unstim"], alpha=0.75, lw=0)
    lim = [min(base.obs.min(), pred.min()) - 0.5, max(base.obs.max(), pred.max()) + 0.5]
    ax.plot(lim, lim, "--", color="0.7", lw=0.8)
    ax.set_xlabel("measured log$_2$FC (IFN-β)")
    ax.set_ylabel("predicted log$_2$FC")
    ax.text(0.03, 0.95, f"ρ = {r:.3f}\nn = {len(base)} genes", transform=ax.transAxes,
            va="top", fontsize=7)
    panel(ax, "a")

    ax = axes[1]
    ax.bar([0, 1], [r, r_ho], color=[C["unstim"], C["hl"]], width=0.6)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["full", "gene\nhold-out"], fontsize=7)
    ax.set_ylim(0, 0.8)
    ax.set_ylabel("ρ vs measured")
    for i, v in enumerate([r, r_ho]):
        ax.text(i, v + 0.02, f"{v:.3f}", ha="center", fontsize=7)
    panel(ax, "b")

    ax = axes[2]
    ax.scatter(pred, hl, s=10, c="0.3", alpha=0.7, lw=0)
    ax.plot(lim, lim, "--", color="0.7", lw=0.8)
    ax.set_xlabel("prediction, full input")
    ax.set_ylabel("prediction, held-out")
    ax.text(0.03, 0.95, f"ρ = {r_two:.3f}", transform=ax.transAxes, va="top", fontsize=7)
    panel(ax, "c")
    fig.tight_layout()
    fig.savefig("fig2_expression.png", dpi=600)
    fig.savefig("fig2_expression.pdf")
    plt.close(fig)


# ---------------------------------------------------------------- Figure 3
def fig3():
    gates = [("max(A,B) > 0.10", 0.0019, 12, 18), ("min(A,B) > 0.05", 0.0168, 9, 18),
             ("min(A,B) > 0.10", 0.4430, 5, 13), ("min(A,B) > 0.15", 0.1429, 2, 6)]
    fig = plt.figure(figsize=(7.2, 2.8))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 0.85, 1.0])

    ax = fig.add_subplot(gs[0])
    x = np.arange(len(gates))
    ax.bar(x, [g[1] for g in gates], color=["#C0392B", "#E67E22", "0.6", "0.6"], width=0.62)
    ax.axhline(0.05, color="0.4", ls=":", lw=0.8)
    ax.set_yscale("log")
    ax.set_ylim(1e-3, 1)
    ax.set_xticks(x)
    ax.set_xticklabels([g[0] for g in gates], rotation=30, ha="right", fontsize=6.5)
    ax.set_ylabel("Mann-Whitney p")
    for i, g in enumerate(gates):
        ax.text(i, g[1] * 1.35, f"{g[1]:.3g}\n{g[2]}/{g[3]}", ha="center", fontsize=6)
    panel(ax, "a")

    ax = fig.add_subplot(gs[1])
    ax.bar([0, 1], [0.0037, 0.1442], color=[C["stim"], "0.55"], width=0.55)
    ax.axhline(0.05, color="0.4", ls=":", lw=0.8)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["denominator\n(resting)", "numerator\n(stimulated)"], fontsize=7)
    ax.set_ylabel("p (act vs rest)")
    for i, v in enumerate([0.0037, 0.1442]):
        ax.text(i, v + 0.006, f"{v:.4f}", ha="center", fontsize=7)
    ax.set_ylim(0, 0.18)
    panel(ax, "b")

    ax = fig.add_subplot(gs[2])
    d = pd.read_csv(f"{D}/flu_visible_sites.csv")
    d = d[d.group.isin(["act_specific", "rest_specific"])]
    A, B = d.A.abs(), d.B.abs()
    pos = {"act_specific": [0.82, 1.18], "rest_specific": [1.82, 2.18]}
    for grp, col in [("act_specific", C["stim"]), ("rest_specific", C["unstim"])]:
        m = d.group == grp
        for k, (v, xpos) in enumerate(zip([A[m], B[m]], pos[grp])):
            ax.boxplot(v, positions=[xpos], widths=0.26, showfliers=False,
                       medianprops=dict(color=col, lw=1.2),
                       boxprops=dict(color=col, lw=0.8),
                       whiskerprops=dict(color=col, lw=0.8),
                       capprops=dict(color=col, lw=0.8))
            ax.scatter(np.random.default_rng(0).normal(xpos, 0.035, len(v)), v,
                       s=5, color=col, alpha=0.55, lw=0, zorder=3)
    ax.axhline(0.1, color="0.5", ls=":", lw=0.8)
    ax.set_xticks([1, 2])
    ax.set_xticklabels(["act-specific", "rest-specific"], fontsize=7)
    ax.set_ylabel("|effect| (predicted log$_2$ ratio)")
    ax.set_yscale("log")
    h = [Line2D([], [], color="0.35", lw=0.8), Line2D([], [], color="0.35", lw=0.8, ls="")]
    ax.legend([plt.Line2D([], [], color="0.35", lw=1.4),
               plt.Line2D([], [], color="0.35", lw=1.4, alpha=0.4)],
              ["numerator A", "denominator B"], fontsize=6, frameon=False, loc="lower left")
    panel(ax, "c")
    fig.tight_layout()
    fig.savefig("fig3_site_layer.png", dpi=600)
    fig.savefig("fig3_site_layer.pdf")
    plt.close(fig)


# ---------------------------------------------------------------- Figure 4
def fig4():
    rows = [
        ("Fairfax", {"|z$_{other}$| ≥ 2": 1077, "untested in the\nother condition": 18 + 2, "direction": 4}, 825),
        ("Quach", {"|z$_{other}$| ≥ 2": 2754, "untested in the\nother condition": 421 + 34, "direction": 0}, 2170),
    ]
    fig, ax = plt.subplots(figsize=(3.4, 2.6))
    for i, (name, cats, kept) in enumerate(rows):
        total = sum(cats.values()) + kept
        bottom = 0
        for lab, v in cats.items():
            ax.bar(i, v / total * 100, bottom=bottom, width=0.55,
                   color={"|z$_{other}$| ≥ 2": "#C0392B", "untested in the\nother condition": "#E67E22",
                          "direction": "0.65"}[lab], label=lab if i == 0 else None)
            bottom += v / total * 100
        ax.bar(i, kept / total * 100, bottom=bottom, width=0.55, color="#2C6FA6",
               label="retained" if i == 0 else None)
        ax.text(i, 101, f"n = {total:,}", ha="center", fontsize=6.5)
        ax.text(i, bottom + kept / total * 50, f"{kept/total*100:.1f}%", ha="center",
                fontsize=6.5, color="white")
    ax.set_xticks([0, 1])
    ax.set_xticklabels([r[0] for r in rows])
    ax.set_ylabel("% of variant–gene–stimulus rows")
    ax.set_ylim(0, 108)
    ax.legend(fontsize=6, frameon=False, loc="lower center", bbox_to_anchor=(0.5, -0.42), ncol=2)
    fig.savefig("fig4_attrition.png", dpi=600)
    fig.savefig("fig4_attrition.pdf")
    plt.close(fig)


# ---------------------------------------------------------------- Figure S1
def figS1():
    th = pd.read_csv(f"{D}/depth_thinning.csv").rho_shallow.values
    fig, ax = plt.subplots(figsize=(3.2, 2.6))
    vals = [0.486, th.mean(), 0.6155]
    errs = [0, th.std(ddof=1), 0]
    ax.bar(range(3), vals, yerr=errs, capsize=3, width=0.55,
           color=["0.65", "#E67E22", "#2C6FA6"],
           error_kw=dict(lw=0.8, ecolor="0.3"))
    ax.set_xticks(range(3))
    ax.set_xticklabels(["circular,\nshallow", "disjoint,\nshallow", "disjoint,\nfull depth"],
                       fontsize=7)
    ax.set_ylabel("ρ (prediction vs measured)")
    ax.set_ylim(0, 0.75)
    for i, v in enumerate(vals):
        lab = f"{v:.3f}" if i != 1 else f"{v:.3f}\n(SD {errs[i]:.3f})"
        ax.text(i, max(v + errs[i] + 0.02, 0.52 if i == 1 else 0),
                lab, ha="center", fontsize=7, va="bottom")
    fig.savefig("figS1_depth.png", dpi=600)
    fig.savefig("figS1_depth.pdf")
    plt.close(fig)


# ---------------------------------------------------------------- Figure S2
def figS2():
    ns = pd.read_csv(f"{D}/noise_splits_rho_newobs.csv").rho.values
    nul = pd.read_csv(f"{D}/null_all_disjoint.csv")
    OBS = 0.6155
    k = int((ns >= OBS).sum())
    pval = (1 + k) / (1 + len(ns))
    fig, axes = plt.subplots(1, 2, figsize=(6.4, 2.5))
    ax = axes[0]
    ax.hist(ns, bins=14, color="0.7", edgecolor="white", lw=0.5, range=(-0.45, 0.65))
    ax.axvline(OBS, color=C["stim"], lw=1.2, zorder=3)
    ax.annotate(f"observed 0.616\n({k}/{len(ns)} draws reach it;\n"
                f"one-sided p = {pval:.3f})",
                xy=(OBS, ax.get_ylim()[1] * 0.55), xytext=(0.30, 0.62),
                textcoords="axes fraction", color=C["stim"], fontsize=6.0,
                ha="left", va="center", linespacing=1.35,
                arrowprops=dict(arrowstyle="-", color=C["stim"], lw=0.6))
    ax.set_xlabel("ρ from split-half conditions")
    ax.set_ylabel("count (50 draws)")
    ax.text(0.02, 0.96, f"mean {ns.mean():+.3f}\nSD {ns.std(ddof=1):.3f}",
            transform=ax.transAxes, va="top", fontsize=6.5, linespacing=1.35)
    panel(ax, "a")

    ax = axes[1]
    EN = {"CD4 T": "CD4 T", "CD8 T": "CD8 T", "NK": "NK", "B 细胞": "B",
          "树突细胞": "dendritic", "巨核细胞": "megakaryocyte", "FCGR3A+ 单核": "FCGR3A+ mono"}
    nul = nul.assign(cond=nul.cond.map(lambda v: EN.get(v, v)))
    d = nul[nul.tag != "real"].sort_values("rho")
    ax.scatter(d.rho, d.mag, s=22, c="#E67E22", lw=0)
    ax.scatter(nul[nul.tag == "real"].rho, nul[nul.tag == "real"].mag, s=55,
               marker="*", c=C["matched"])
    OFF = {"CD4 T": (-16, 2), "CD8 T": (5, -8), "NK": (2, 6), "B": (4, 4),
           "dendritic": (5, 0), "megakaryocyte": (5, 0), "FCGR3A+ mono": (-6, 7)}
    for _, r in d.iterrows():
        dx, dy = OFF.get(r.cond, (4, 3))
        ha = "right" if dx < 0 else "left"
        ax.annotate(r.cond if len(r.cond) < 14 else r.cond[:12], (r.rho, r.mag),
                    textcoords="offset points", xytext=(dx, dy), fontsize=5.5,
                    color="0.35", ha=ha)
    ax.annotate("matched", (nul[nul.tag == "real"].rho.iloc[0],
                            nul[nul.tag == "real"].mag.iloc[0]),
                textcoords="offset points", xytext=(-30, 6), fontsize=6.5)
    ax.set_xlabel("ρ with measured IFN-β response")
    ax.set_ylabel("median |predicted log$_2$FC|")
    ax.axvline(0, color="0.85", lw=0.6)
    ax.set_xlim(-0.26, 0.66)
    ax.set_ylim(0.16, 1.09)
    panel(ax, "b")
    fig.tight_layout()
    fig.savefig("figS2_nulls.png", dpi=600)
    fig.savefig("figS2_nulls.pdf")
    plt.close(fig)


if __name__ == "__main__":
    for f in (fig1, fig2, fig3, fig4, figS1, figS2):
        f()
        print("ok", f.__name__, flush=True)
