#!/usr/bin/env python3
"""Flu 臂的敏感性合集：阈值 n、去重叠、门控、表达混淆（相关 + 调整）。"""
import glob

import numpy as np
import pandas as pd
from scipy import stats

d = pd.read_csv("b_data/quach_long_strict.csv")
f = d[d.stim == "Influenza_6h"].copy()
w = f.pivot_table(index=["variant", "gene", "group"], columns="cond_role",
                  values=["dn_lfc", "dn_act"], aggfunc="first").reset_index()
w.columns = ["variant", "gene", "group"] + [f"{a}_{b}" for a, b in w.columns[3:]]
print("列:", [c for c in w.columns])

obs = pd.concat([pd.read_csv(x) for x in glob.glob("b_data/obs_cpm_pooled_*.csv")]).drop_duplicates("ensg")
obs["d_expr"] = np.log2((obs.cpm_stim_full + 1) / (obs.cpm_ctrl_full + 1))
w = w.merge(obs[["ensg", "d_expr", "cpm_ctrl_full"]], left_on="gene", right_on="ensg", how="left")


def prep(x, thr=0.1):
    x = x.dropna(subset=["dn_lfc_刺激", "dn_lfc_对照"]).copy()
    x["A"], x["B"] = x.dn_lfc_刺激.abs(), x.dn_lfc_对照.abs()
    x = x[(x[["A", "B"]].max(axis=1) > thr) & (x.A > 0) & (x.B > 0)].copy()
    x["lr"] = np.log2(x.A / x.B)
    return x


print("\n== 各阈值下的 n 与 p ==")
for thr in (0.05, 0.10, 0.15, 0.20):
    x = prep(w, thr)
    a, r = x[x.group == "act_specific"].lr, x[x.group == "rest_specific"].lr
    print(f"阈值 {thr:.2f}: n={len(a)}/{len(r)}  p={stats.mannwhitneyu(a, r, alternative='two-sided').pvalue:.4f}  "
          f"中位 {a.median():+.3f}/{r.median():+.3f}")

x = prep(w, 0.1)

print("\n== 去重叠子集（排除与 Fairfax 重叠的变异-基因）==")
lab = pd.read_csv("b_data/fairfax_labels_fixed.csv")
ov = set(zip(lab.variant, lab.gene))
x["overlap"] = [(r.variant, r.gene) in ov for r in x.itertuples(index=False)]
sub = x[~x.overlap]
a, r = sub[sub.group == "act_specific"].lr, sub[sub.group == "rest_specific"].lr
print(f"去掉 {int(x.overlap.sum())}/{len(x)} → n={len(a)}/{len(r)}  "
      f"p={stats.mannwhitneyu(a, r, alternative='two-sided').pvalue:.4f}  "
      f"中位 {a.median():+.3f}/{r.median():+.3f}")

print("\n== 门控敏感性（DENOISED 门控下判可见的抽样位点，去掉后）==")
vis = pd.concat([pd.read_csv(p) for p in glob.glob("b_data/visibility_fairfax_*.csv")])
used = set(vis.variant.astype(str))
res = x[~x.variant.astype(str).isin(used)]
a, r = res[res.group == "act_specific"].lr, res[res.group == "rest_specific"].lr
print(f"去掉已用位点 → n={len(a)}/{len(r)}  "
      f"p={stats.mannwhitneyu(a, r, alternative='two-sided').pvalue:.4f}")

print("\n== 表达混淆：相关 ==")
q = x.dropna(subset=["d_expr"])
rho, p = stats.spearmanr(q.lr, q.d_expr)
print(f"lr vs 观测表达变化 rho={rho:+.3f} p={p:.4f}（n={len(q)}）")

print("\n== 表达混淆：调整 ==")
q = x.dropna(subset=["d_expr", "dn_act_刺激", "dn_act_对照"]).copy()
q["d_act"] = np.log2((q.dn_act_刺激.abs() + 1e-9) / (q.dn_act_对照.abs() + 1e-9))
q["is_act"] = (q.group == "act_specific").astype(int)
def ols(y, X):
    X = np.column_stack([np.ones(len(X)), np.asarray(X, float)])
    beta, *_ = np.linalg.lstsq(X, np.asarray(y, float), rcond=None)
    resid = np.asarray(y, float) - X @ beta
    dof = len(y) - X.shape[1]
    s2 = resid @ resid / dof
    cov = s2 * np.linalg.inv(X.T @ X)
    se = np.sqrt(np.diag(cov))
    t = beta / se
    p = 2 * stats.t.sf(np.abs(t), dof)
    r2 = 1 - (resid @ resid) / ((y - y.mean()) ** 2).sum()
    return beta, p, r2

for cols, label in [(["is_act"], "只放分组"),
                    (["is_act", "d_act"], "分组 + 预测活性变化"),
                    (["is_act", "d_act", "d_expr"], "分组 + 活性变化 + 观测表达变化")]:
    beta, pv, r2 = ols(q.lr.values, q[cols].values)
    print(f"{label}: 分组系数 {beta[1]:+.3f}  p={pv[1]:.4f}  R²={r2:.3f}  n={len(q)}")

print("\n== 单侧检查（act 中位 > rest）==")
a, r = x[x.group == "act_specific"].lr, x[x.group == "rest_specific"].lr
print(f"单侧 p={stats.mannwhitneyu(a, r, alternative='greater').pvalue:.4f}")
