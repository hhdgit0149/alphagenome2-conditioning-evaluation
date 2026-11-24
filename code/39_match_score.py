#!/usr/bin/env python3
"""匹配度计算（PREREG_v2 补充C §5）：Quach 每种刺激 vs 静息的表达变化，
与 Kang IFN-β 条件表达变化的 Spearman 相关。

数据源（无额外网络请求）：
  naive median_tpm    b_data/strict_z_quach.csv（35 号：other=naive 行，act 基因）
                      + b_data/strict_z_own_quach.csv（38 号：other=naive 行，rest/shared 基因）
  刺激 median_tpm     b_data/strict_z_own_quach.csv（other=刺激 行，act 基因）
                      + b_data/strict_z_quach.csv（other=刺激 行，rest 基因）
  Δ_IFNβ             b_data/obs_cpm_pooled_*.csv（Kang GSE96583 全部 CD14 细胞）
基因集：各刺激的 act+rest 标签基因中，两条件 median_tpm 均 >0 者。
输出 b_data/match_scores.csv。只算 Quach（Fairfax microarray 的 median_tpm 整列 NA）。
"""
import glob

import numpy as np
import pandas as pd
from scipy import stats

SEED = 20261114
STIM = ["LPS_6h", "Pam3CSK4_6h", "R848_6h", "Influenza_6h"]

lab = pd.read_csv("b_data/Quach_2016_labels_fixed.csv")
zo = pd.read_csv("b_data/strict_z_quach.csv")            # 另一条件（35 号）
ow = pd.read_csv("b_data/strict_z_own_quach.csv")        # 本条件（38 号）
for df in (zo, ow):
    assert not (df.n_pair_rows < 0).any(), "查询输出有失败行"


def tpm_map(df, cond):
    """(cond, gene) -> median_tpm。同一基因多行时值必须一致。"""
    d = df[(df.other_cond == cond) & (df.n_pair_rows > 0)].copy()
    d["tpm"] = pd.to_numeric(d.tpm, errors="coerce")
    g = d.groupby("gene").tpm.agg(["min", "max", "count"])
    assert (g["min"] == g["max"]).all(), f"{cond} 同基因 tpm 不一致"
    return g["min"].to_dict()


obs = (pd.concat([pd.read_csv(f) for f in glob.glob("b_data/obs_cpm_pooled_*.csv")])
       .drop_duplicates("ensg"))
d_ifn = np.log2((obs.cpm_stim_full + 1) / (obs.cpm_ctrl_full + 1))
d_ifn.index = obs.ensg.values

out = []
for s in STIM:
    genes = set(lab[(lab.stim == s) & (lab.group.isin(["act_specific", "rest_specific"]))].gene)
    naive = tpm_map(zo, "naive")
    naive.update({k: v for k, v in tpm_map(ow, "naive").items() if k not in naive})
    stim = tpm_map(ow, s)
    stim.update({k: v for k, v in tpm_map(zo, s).items() if k not in stim})
    use = [g for g in genes if g in naive and g in stim and naive[g] > 0 and stim[g] > 0
           and g in d_ifn.index and pd.notna(d_ifn[g])]
    ds = np.log2([stim[g] / naive[g] for g in use])
    di = d_ifn.loc[list(use)].values
    rho, p = stats.spearmanr(ds, di)
    rng = np.random.default_rng(SEED)
    boots = []
    ds_a, di_a = np.array(ds), np.array(di)
    for _ in range(1000):
        ix = rng.integers(0, len(use), len(use))
        if np.std(ds_a[ix]) == 0 or np.std(di_a[ix]) == 0:
            continue
        boots.append(stats.spearmanr(ds_a[ix], di_a[ix])[0])
    lo, hi = np.percentile(boots, [2.5, 97.5])
    clearly_neg = (rho < 0) and (hi < 0)
    out.append(dict(stim=s, n_genes=len(use), spearman=rho, p=p,
                    ci95_lo=lo, ci95_hi=hi, clearly_negative=bool(clearly_neg)))
    print(f"{s:14s} n={len(use):4d}  r={rho:+.3f} (p={p:.4f})  "
          f"bootstrap95%CI [{lo:+.3f}, {hi:+.3f}]  明显为负={clearly_neg}")

res = pd.DataFrame(out)
res.to_csv("b_data/match_scores.csv", index=False)
print("\n写入 b_data/match_scores.csv")
print("预注册判读：r_S<0 且 CI 上界<0 的刺激，反向不计入跨刺激方向不一致；"
      "近零或为正的刺激反向才计入。")
