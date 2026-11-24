#!/usr/bin/env python3
"""表达量门槛的校准（在组间分析之前跑）。

把单细胞去噪头的 `_ACTIVE` 值与观测 CPM（全部单核细胞合并）做单调对应，
找出相当于 1 CPM 的 `_ACTIVE` 值。只做技术校准，不看分组。
"""
import glob

import numpy as np
import pandas as pd
from scipy import stats

SC = "b_data/fairfax_v2_dedup.csv"
OBS = sorted(glob.glob("b_data/obs_cpm_pooled_*.csv"))[-1]


def main():
    d = pd.read_csv(SC)
    obs = pd.read_csv(OBS)
    print(f"打分 {len(d)} 行，覆盖 {d.variant.nunique()} 个变异、{d.gene.nunique()} 个基因")
    print(f"dn_act 空值 {d.dn_act.isna().sum()}｜rna_act 空值 {d.rna_act.isna().sum()}")

    m = d.merge(obs[["ensg", "cpm_ctrl_full", "cpm_stim_full"]],
                left_on="gene", right_on="ensg", how="left")
    m["cpm_obs"] = np.where(m.cond_role.isin(["stim", "halfA"]),
                            m.cpm_stim_full, m.cpm_ctrl_full)
    m.attrs = {}
    m.to_csv("b_data/calib_input.csv", index=False)

    x = m.dropna(subset=["dn_act", "cpm_obs"])
    print(f"\n可用 (dn_act, 观测 CPM) 对：{len(x)}｜涉及基因 {x.gene.nunique()}")
    rho = stats.spearmanr(x.dn_act, x.cpm_obs)
    print(f"秩相关（单调性检查）：rho={rho.correlation:.3f} p={rho.pvalue:.2g}")

    # 单调对应：按 dn_act 排序后取分位映射
    x = x.sort_values("dn_act")
    q_obs = x.cpm_obs.rank(pct=True)
    x = x.assign(q=q_obs)
    target = x[x.cpm_obs >= 1]
    if len(target) == 0:
        print("没有基因达到 1 CPM，无法校准")
        return
    # 1 CPM 对应观测分位 → 找同分位的 dn_act
    p1 = (x.cpm_obs < 1).mean()
    thr = float(x.dn_act.quantile(p1))
    print(f"\n观测 CPM ≥1 的比例：{1-p1:.1%}")
    print(f"对应分位 {p1:.3f} 的 dn_act 门槛 = {thr:.3f}")
    print(f"（门槛下保留 {int((x.dn_act>=thr).sum())}/{len(x)} 行）")

    pd.DataFrame([dict(metric="dn_act", threshold_1cpm=thr,
                       spearman_vs_obs_cpm=round(float(rho.correlation), 3),
                       n_pairs=len(x), n_genes=int(x.gene.nunique()),
                       obs_file=OBS, score_file=SC)]).to_csv(
        "b_data/expr_threshold.csv", index=False)
    print("\n写出 b_data/expr_threshold.csv")

    # 备用：转录组头同样校准
    y = m.dropna(subset=["rna_act", "cpm_obs"])
    if len(y):
        p1y = (y.cpm_obs < 1).mean()
        print(f"（备用）rna_act 门槛 = {float(y.rna_act.quantile(p1y)):.4f}")


if __name__ == "__main__":
    main()
