#!/usr/bin/env python3
"""Fairfax 主分析的修正版（按外部复算意见）。

修正点：
  1. p 值统一为双侧（置换 + Mann-Whitney 各报一次，n_perm=20000）
  2. 报告极端值：找出贡献最大的位点，给去掉后的结果；再给一个事后加的有界指标
  3. 描述改为「三组在刺激条件下的效应中位数」+ 趋势检验
  4. 标签用按刺激独立重建的版本（29_rebuild_labels.py）
"""
import glob

import numpy as np
import pandas as pd
from scipy import stats

RNG = np.random.default_rng(20261114)
THR = 0.1
EXPR_THR = 0.850


def load_scores():
    base = pd.read_csv("b_data/fairfax_v2_dedup.csv")
    add = pd.concat([pd.read_csv(f) for f in sorted(glob.glob("b_data/fairfax_v2b_*.csv"))])
    sc = pd.concat([base, add], ignore_index=True)
    sc = sc.drop_duplicates(["variant", "gene", "cond_role"])
    assert not sc.duplicated(["variant", "gene", "cond_role"]).any()
    print(f"打分 {len(sc)} 行｜变异-基因 {sc[['variant','gene']].drop_duplicates().shape[0]}")
    return sc


def wide(sc):
    w = sc.pivot_table(index=["variant", "gene"], columns="cond_role",
                       values=["dn_lfc", "rna_lfc", "dn_act"], aggfunc="first")
    w.columns = [f"{a}_{b}" for a, b in w.columns]
    return w.reset_index()


def contrast(x, a_col, r_col, n=20000):
    x = x.copy()
    x["a"], x["r"] = x[a_col].abs(), x[r_col].abs()
    x["lr"] = np.log2((x.a + 1e-6) / (x.r + 1e-6))
    A = x[x.group == "act_specific"].lr
    R = x[x.group == "rest_specific"].lr
    if len(A) < 5 or len(R) < 5:
        return None
    obs = A.mean() - R.mean()
    lab = x.group.values
    hits = 0
    for _ in range(n):
        p = RNG.permutation(lab)
        hits += abs(x.lr[p == "act_specific"].mean() -
                    x.lr[p == "rest_specific"].mean()) >= abs(obs)
    mwu = stats.mannwhitneyu(A, R, alternative="two-sided").pvalue
    return dict(n_act=len(A), n_rest=len(R), diff=round(float(obs), 3),
                p_perm_2s=round((hits + 1) / (n + 1), 4), p_mwu_2s=round(float(mwu), 4))


def main():
    sc = load_scores()
    w = wide(sc)
    obs = pd.concat([pd.read_csv(f) for f in glob.glob("b_data/obs_cpm_pooled_*.csv")])
    obs = obs.drop_duplicates("ensg")
    w = w.merge(obs[["ensg", "cpm_ctrl_full", "cpm_stim_full"]],
                left_on="gene", right_on="ensg", how="left").drop(columns=["ensg"])
    lab = pd.read_csv("b_data/fairfax_labels_fixed.csv")

    print("\n" + "=" * 70)
    print("IFNg_24h 主分析（标签按刺激独立重建）")
    l = lab[lab.stim == "IFNg_24h"].drop_duplicates(["variant", "gene", "group"])
    d = w.merge(l, on=["variant", "gene"], how="inner")
    d["vis"] = d[["dn_lfc_stim", "dn_lfc_ctrl"]].abs().max(axis=1) > THR
    v = d[d.vis].copy()
    print(f"合并 {len(d)}｜可见 {len(v)}｜分组 {v.group.value_counts().to_dict()}")

    # 3. 三组在刺激条件下的效应中位数（正确的描述方式）
    print("\n各组的预测效应（|dn_lfc|）中位数：")
    for g in ["act_specific", "shared", "rest_specific"]:
        x = v[v.group == g]
        if len(x) == 0:
            continue
        print(f"  {g:14s} n={len(x):3d}｜刺激 {x.dn_lfc_stim.abs().median():.3f}"
              f"｜对照 {x.dn_lfc_ctrl.abs().median():.3f}"
              f"｜刺激更大比例 {(x.dn_lfc_stim.abs() > x.dn_lfc_ctrl.abs()).mean():.2f}")
    order = {"act_specific": 0, "shared": 1, "rest_specific": 2}
    v["ord"] = v.group.map(order)
    v["val"] = v.dn_lfc_stim.abs()
    tr = stats.spearmanr(v.ord, v.val)
    print(f"  趋势检验（组序 vs 刺激条件效应）: rho={tr.correlation:+.3f}, p={tr.pvalue:.4f}")

    # 1. 主对比（双侧）
    r = contrast(v, "dn_lfc_stim", "dn_lfc_ctrl")
    print(f"\n主对比（双侧）：{r}")

    # 2. 极端值
    v = v.assign(lr=np.log2((v.dn_lfc_stim.abs() + 1e-6) / (v.dn_lfc_ctrl.abs() + 1e-6)))
    top = v.reindex(v.lr.abs().sort_values(ascending=False).index).head(3)
    print("\n贡献最大的位点：")
    print(top[["variant", "gene", "group", "dn_lfc_stim", "dn_lfc_ctrl", "lr"]]
          .round(3).to_string(index=False))
    worst = top.iloc[0]["variant"]
    r2 = contrast(v[v.variant != worst], "dn_lfc_stim", "dn_lfc_ctrl")
    print(f"去掉 {worst} 后：{r2}")

    # 有界指标（事后）：用 |a|-|r| 的符号化差值，避免分母趋零
    b = v.copy()
    b["bounded"] = np.sign(b.dn_lfc_stim) * b.dn_lfc_stim.abs() - np.sign(b.dn_lfc_ctrl) * b.dn_lfc_ctrl.abs()
    A = b[b.group == "act_specific"]["bounded"]
    R = b[b.group == "rest_specific"]["bounded"]
    print(f"有界指标（事后）：差值 {A.mean() - R.mean():+.4f}, "
          f"MWU 双侧 p={stats.mannwhitneyu(A, R).pvalue:.4f}")

    # 安慰剂
    p = d.copy()
    p["vis_p"] = p[["dn_lfc_halfA", "dn_lfc_halfB"]].abs().max(axis=1) > THR
    print(f"\n安慰剂（单独判定 n={int(p.vis_p.sum())}）："
          f"{contrast(p[p.vis_p], 'dn_lfc_halfA', 'dn_lfc_halfB')}")

    # 替代解释
    v2 = v.dropna(subset=["dn_act_stim", "dn_act_ctrl"])
    chg = np.abs(v2.dn_act_stim - v2.dn_act_ctrl)
    rho = stats.spearmanr(v2.lr, chg)
    print(f"\n效应比 vs 预测表达变化：rho={rho.correlation:+.3f}, p={rho.pvalue:.2f}")
    for g in ["act_specific", "shared", "rest_specific"]:
        x = v[v.group == g]
        if len(x) == 0:
            continue
        print(f"  {g:14s} 观测表达变化中位（log2 CPM）"
              f"={np.median(np.log2((x.cpm_stim_full + 1) / (x.cpm_ctrl_full + 1))):+.2f}")

    # 4. 三刺激比较
    print("\n" + "=" * 70)
    print("三种刺激（修正标签后）")
    rows = []
    for s in ["IFNg_24h", "LPS_24h", "LPS_2h"]:
        ls = lab[lab.stim == s].drop_duplicates(["variant", "gene", "group"])
        dd = w.merge(ls, on=["variant", "gene"], how="inner")
        dd["vis"] = dd[["dn_lfc_stim", "dn_lfc_ctrl"]].abs().max(axis=1) > THR
        rr = contrast(dd[dd.vis], "dn_lfc_stim", "dn_lfc_ctrl")
        if rr:
            rows.append(dict(stim=s, **rr))
            print(f"  {s:10s} {rr}")
    pd.DataFrame(rows).to_csv("b_data/fairfax_three_stim_fixed.csv", index=False)

    # 探索性：同一 IFN 条件下，IFN 特异 vs LPS_2h 特异
    print("\n" + "=" * 70)
    lp = lab[lab.stim == "LPS_2h"].drop_duplicates(["variant", "gene", "group"])
    s2 = lp[lp.group == "act_specific"][["variant", "gene"]]
    a2 = d[d.group == "act_specific"][["variant", "gene", "dn_lfc_stim"]].merge(
        s2, on=["variant", "gene"], how="inner")
    r3 = d[(d.group == "act_specific") & (~d.variant.isin(set(s2.variant)))][
        ["variant", "gene", "dn_lfc_stim"]]
    if len(a2) > 5 and len(r3) > 5:
        print(f"探索性（IFN 条件内）：LPS_2h 特异 n={len(a2)} 中位 "
              f"{a2.dn_lfc_stim.abs().median():.3f}｜其余 IFN 特异 n={len(r3)} 中位 "
              f"{r3.dn_lfc_stim.abs().median():.3f}｜"
              f"MWU p={stats.mannwhitneyu(a2.dn_lfc_stim.abs(), r3.dn_lfc_stim.abs()).pvalue:.4f}")


if __name__ == "__main__":
    main()
