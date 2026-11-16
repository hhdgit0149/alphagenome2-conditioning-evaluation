#!/usr/bin/env python3
"""头敏感性分析（审计第 3、4 条，零 API 调用）。

1) 主对比在四个头上并排：dn_lfc / rna_lfc（对数比）与 dn_act / rna_act（活性比）
   每个头都：可见阈值 0.1，秩和检验，激活特异 vs 静息特异。
2) 可见性门控：同一批 200 个抽样位点，RNA_SEQ 门控 vs DENOISED 门控的可见率。
3) 用 DENOISED 门控判定的"已用位点"重跑一次去除敏感性。

用法：python 42_head_sensitivity.py b_data/fairfax_long_fixed.csv
"""
import glob
import sys

import numpy as np
import pandas as pd
from scipy import stats

THR = 0.1
path = sys.argv[1] if len(sys.argv) > 1 else "b_data/fairfax_long_fixed.csv"
d = pd.read_csv(path)
ROLE = {"刺激": "stim", "对照": "ctrl", "对照对半A": "hA", "对照对半B": "hB"}
d["role"] = d.cond_role.map(ROLE)
assert d.role.notna().all()

key = ["variant", "gene", "stim", "group"]
piv = {}
for col in ["dn_lfc", "rna_lfc", "dn_act", "rna_act"]:
    w = d.pivot_table(index=key, columns="role", values=col, aggfunc="first")
    w.columns = [f"{col}__{c}" for c in w.columns]
    piv[col] = w.reset_index()
wide = piv["dn_lfc"]
for col in ["rna_lfc", "dn_act", "rna_act"]:
    wide = wide.merge(piv[col], on=key, how="outer")

print(f"长表 {path}：{len(wide)} 个（变异×基因×刺激）单元\n")
print("== 四个头并排（按刺激分别算，可见阈值 0.1，秩和双侧）==")
rows = []
for head, unsigned in [("dn_lfc", False), ("rna_lfc", False),
                       ("dn_act", True), ("rna_act", True)]:
    a, b = f"{head}__stim", f"{head}__ctrl"
    if a not in wide.columns or b not in wide.columns:
        print(f"{head}: 缺列，跳过")
        continue
    for stim in sorted(wide.stim.unique()):
        x = wide[wide.stim == stim].dropna(subset=[a, b]).copy()
        x["A"] = x[a] if unsigned else x[a].abs()
        x["B"] = x[b] if unsigned else x[b].abs()
        x = x[(x[["A", "B"]].max(axis=1) > THR) & (x.A > 0) & (x.B > 0)].copy()
        if len(x) < 6:
            continue
        x["lr"] = np.log2(x.A / x.B)
        A = x.loc[x.group == "act_specific", "lr"]
        R = x.loc[x.group == "rest_specific", "lr"]
        if len(A) < 3 or len(R) < 3:
            print(f"{head:8s} {stim:10s} 位点不足（{len(A)}/{len(R)}）")
            continue
        p = stats.mannwhitneyu(A, R, alternative="two-sided").pvalue
        print(f"{head:8s} {stim:10s} n={len(A):3d}/{len(R):3d}  "
              f"中位数 {A.median():+.3f} / {R.median():+.3f}  秩和 p={p:.4f}")
        rows.append(dict(head=head, stim=stim, n_act=len(A), n_rest=len(R),
                         med_act=A.median(), med_rest=R.median(), p_rank=p))

print("\n== 可见性门控：同一批 200 个抽样位点，两种头 ==")
vis = pd.concat([pd.read_csv(f) for f in glob.glob("b_data/visibility_fairfax_*.csv")])
sc = pd.concat([pd.read_csv("b_data/fairfax_v2_dedup.csv")] +
               [pd.read_csv(f) for f in glob.glob("b_data/fairfax_v2b_*.csv")],
               ignore_index=True).drop_duplicates(["variant", "gene", "cond_role"])
z = (sc[sc.cond_role.isin(["stim", "ctrl"])]
     .pivot_table(index=["variant", "gene"], columns="cond_role",
                  values="dn_lfc", aggfunc="first")
     .rename(columns={"stim": "dn_stim", "ctrl": "dn_ctrl"}).reset_index())
j = vis.merge(z, on=["variant", "gene"], how="left")
print(f"抽样位点 {len(j)}，DENOISED 打分覆盖 {int(j.dn_stim.notna().sum())}")
j["vis_rna"] = j[["ctrl", "stim"]].abs().max(axis=1) > THR
j["vis_dn"] = j[["dn_stim", "dn_ctrl"]].abs().max(axis=1) > THR
print(j.groupby("group")[["vis_rna", "vis_dn"]].mean().round(3).to_string())
print(f"总体：RNA_SEQ 门控 {j.vis_rna.mean():.1%}，DENOISED 门控 {j.vis_dn.mean():.1%}（阈值 0.1）")

print("\n== 用 DENOISED 门控重跑'去除已用位点'敏感性 ==")
x = wide.dropna(subset=["dn_lfc__stim", "dn_lfc__ctrl"]).copy()
x["A"], x["B"] = x["dn_lfc__stim"].abs(), x["dn_lfc__ctrl"].abs()
x = x[(x[["A", "B"]].max(axis=1) > THR) & (x.A > 0) & (x.B > 0)].copy()
x["lr"] = np.log2(x.A / x.B)
used_dn = set(j.loc[j.vis_dn, "variant"].astype(str))
A = x[(x.group == "act_specific") & ~x.variant.astype(str).isin(used_dn)]["lr"]
R = x[(x.group == "rest_specific") & ~x.variant.astype(str).isin(used_dn)]["lr"]
if len(A) > 2 and len(R) > 2:
    print(f"去掉 DENOISED 门控下判可见的抽样位点：n={len(A)}/{len(R)}  "
          f"中位数 {A.median():+.3f} / {R.median():+.3f}  "
          f"秩和 p={stats.mannwhitneyu(A, R, alternative='two-sided').pvalue:.4f}")
pd.DataFrame(rows).to_csv("b_data/head_sensitivity.csv", index=False)
print("\n写入 b_data/head_sensitivity.csv")
