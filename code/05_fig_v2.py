#!/usr/bin/env python3
"""出图：v2 结果（散点 + 逐轨道一致率带 95% CI + 前列腺/其他组织对比）。"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                     "axes.spines.right": False, "font.sans-serif": ["DejaVu Sans"]})

V, S, SM = "variants_v2.csv", "scores_v2.csv", "summary_v2.csv"
v = pd.read_csv(V, dtype={"variant_id": str})
s = pd.read_csv(S, dtype={"input_id": str}, low_memory=False)
df = s.merge(v[["variant_id", "effect", "significant"]],
             left_on="input_id", right_on="variant_id", suffixes=("_ag", ""))
sm = pd.read_csv(SM)
tr = sm[sm["level"].isna()].copy()

fig = plt.figure(figsize=(13.2, 4.6))
gs = fig.add_gridspec(1, 3, width_ratios=[1, 1.18, 0.8], wspace=0.72)

# A 散点
ax = fig.add_subplot(gs[0, 0])
sub = df[(df.output_type == "DNASE") & (df.biosample_name == "LNCaP clone FGC")
         & (df.significant == 1)]
x, y = sub.effect.values, sub.raw_score.values
ax.axhline(0, color="#ccc", lw=.8); ax.axvline(0, color="#ccc", lw=.8)
ax.scatter(x, y, s=9, alpha=.45, color="#2f6f9f", edgecolor="none")
b = np.polyfit(x, y, 1); xs = np.linspace(x.min(), x.max(), 50)
ax.plot(xs, np.polyval(b, xs), color="#c0392b", lw=1.4)
row = tr[(tr.output_type == "DNASE") & (tr.biosample_name == "LNCaP clone FGC")].iloc[0]
ax.set_xlabel("experimental ALT fraction − 0.5")
ax.set_ylabel("predicted effect (DNase, LNCaP)")
ax.set_title(f"DNase / LNCaP clone FGC\n"
             f"ρ={row.spearman:.2f} [{row.rho_lo:.2f},{row.rho_hi:.2f}], "
             f"agree={row['方向一致率']:.0%}, n={int(row.n)}", fontsize=8.5)

# B 逐轨道一致率 + CI
ax = fig.add_subplot(gs[0, 1])
top = tr.sort_values("方向一致率", ascending=False).head(16).iloc[::-1]
lab = []
ABBR = {"CHIP_TF": "TF", "CHIP_HISTONE": "HIST", "DNASE": "DNase", "ATAC": "ATAC"}
for r in top.itertuples():
    mark = r.histone_mark if isinstance(r.histone_mark, str) else (
        r.transcription_factor if isinstance(r.transcription_factor, str) else "")
    name = str(r.biosample_name).replace("epithelial cell of prostate", "prostate epi.")
    mark = {"POLR2AphosphoS5": "POLR2A-pS5"}.get(mark, mark)
    lab.append(f"{ABBR.get(r.output_type, r.output_type)} | {name[:16]}"
               + (f" | {mark[:13]}" if isinstance(mark, str) else ""))
yy = np.arange(len(top))
ax.errorbar(top["方向一致率"], yy,
            xerr=[top["方向一致率"] - top["一致率_lo"], top["一致率_hi"] - top["方向一致率"]],
            fmt="o", color="#2f6f9f", ecolor="#2f6f9f", ms=4, lw=1, capsize=2)
ax.axvline(.5, color="#c0392b", ls="--", lw=1)
ax.set_yticks(yy); ax.set_yticklabels(lab, fontsize=6.6)
ax.set_xlim(.45, .9)
ax.set_xlabel("direction agreement (95% CI), chance = 0.5")
ax.set_title(f"per-track, {int(tr['n'].max())} PRAD as-aQTLs", fontsize=8.5)

# C 阴性对照：前列腺 vs 其他组织
ax = fig.add_subplot(gs[0, 2])
try:
    g = pd.read_csv("negctrl_biosample.csv")
    data = [g[g.prostate == 1]["agree"], g[g.prostate == 0]["agree"]]
    try:
        bp = ax.boxplot(data, tick_labels=["prostate\n(%d tracks)" % len(data[0]),
                                           "other tissues\n(%d tracks)" % len(data[1])],
                        widths=.55, patch_artist=True, showfliers=False)
    except TypeError:
        bp = ax.boxplot(data, widths=.55, patch_artist=True, showfliers=False)
        ax.set_xticks([1, 2])
        ax.set_xticklabels(["prostate\n(%d tracks)" % len(data[0]),
                            "other tissues\n(%d tracks)" % len(data[1])])
    for patch, c in zip(bp["boxes"], ["#2f6f9f", "#b0b7bd"]):
        patch.set_facecolor(c); patch.set_alpha(.55); patch.set_edgecolor("#444")
    for med in bp["medians"]:
        med.set_color("#c0392b"); med.set_linewidth(1.4)
    for i, d in enumerate(data, 1):
        ax.scatter(np.random.normal(i, .045, len(d)), d, s=6, alpha=.35,
                   color="#2f6f9f" if i == 1 else "#7f8c8d", edgecolor="none")
    ax.axhline(.5, color="#c0392b", ls="--", lw=1)
    ax.set_ylabel("direction agreement")
    ax.set_title("same variants, all 493 biosamples\n(model, not experimental, control)", fontsize=8.5)
except Exception as err:
    ax.text(.5, .5, f"control not available\n{err}", ha="center", va="center", fontsize=7)
    ax.axis("off")

fig.savefig("fig_v2.png", dpi=195, bbox_inches="tight")
print("saved fig_v2.png")
