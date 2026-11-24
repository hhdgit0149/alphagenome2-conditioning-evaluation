#!/usr/bin/env python3
"""Fairfax 主分析（按预注册第 10–13 条）。

流程：
  1. 按刺激拆表，标签用相对该刺激的 group（同一刺激内同变异+基因两个标签才算冲突）
  2. 主指标：dn_lfc，可见 = 任一条件 |dn_lfc| > 0.1
  3. 主对比：激活特异 vs 静息特异，效应比 log2(|stim|/|ctrl|)，标签置换
  4. 安慰剂：halfA vs halfB（可见位点单独判定，另报沿用主集合的版本）
  5. 三刺激排序：各自的差值 + 重采样 CI
  6. 敏感性：表达量门槛、|z|<2、表达量变化 <2 倍、rna_lfc 头、含/不含 237 个已用位点
"""
import glob

import numpy as np
import pandas as pd
from scipy import stats

RNG = np.random.default_rng(20261113)
STIMS = ["IFNg_24h", "LPS_24h", "LPS_2h"]
THR = 0.1


def load():
    sc = pd.read_csv("b_data/fairfax_v2_dedup.csv")
    lab = pd.concat([pd.read_csv(f) for f in sorted(glob.glob("b_data/fairfax_labels_*.csv"))])
    assert not sc.duplicated(["variant", "gene", "cond_role"]).any(), "打分表主键重复"
    w = sc.pivot_table(index=["variant", "gene"], columns="cond_role",
                       values=["dn_lfc", "rna_lfc", "dn_act", "rna_act"],
                       aggfunc="first")
    w.columns = [f"{a}_{b}" for a, b in w.columns]
    w = w.reset_index()
    print(f"宽表 {len(w)} 个变异-基因")
    return w, lab


def prep(w, lab, stim, metric="dn_lfc"):
    """取出某刺激的标签，与打分合并。"""
    l = lab[lab.condition == stim].drop_duplicates(["variant", "gene", "group"])
    n0 = len(l)
    conf = l.groupby(["variant", "gene"]).group.nunique()
    n_conf = int((conf > 1).sum())
    if n_conf:
        keep = conf[conf == 1].index
        l = l.set_index(["variant", "gene"]).loc[keep].reset_index()
        print(f"  {stim}：标签冲突 {n_conf} 个，已剔除")
    d = w.merge(l, on=["variant", "gene"], how="inner")
    d["a"] = d[f"{metric}_stim"].abs()
    d["r"] = d[f"{metric}_ctrl"].abs()
    d["lr"] = np.log2((d.a + 1e-6) / (d.r + 1e-6))
    d["vis"] = d[["a", "r"]].max(axis=1) > THR
    print(f"  {stim}：标签 {n0} → 合并后 {len(d)}｜可见 {int(d.vis.sum())}")
    return d


def perm_diff(x, n=5000):
    obs = x[x.group == "act_specific"].lr.mean() - x[x.group == "rest_specific"].lr.mean()
    lab = x.group.values
    hits = 0
    for _ in range(n):
        p = RNG.permutation(lab)
        d = x.lr[p == "act_specific"].mean() - x.lr[p == "rest_specific"].mean()
        hits += abs(d) >= abs(obs)
    return obs, (hits + 1) / (n + 1)


def boot_diff(x, n=2000):
    vals = []
    for _ in range(n):
        idx = RNG.integers(0, len(x), len(x))
        s = x.iloc[idx]
        a = s[s.group == "act_specific"].lr
        r = s[s.group == "rest_specific"].lr
        if len(a) > 2 and len(r) > 2:
            vals.append(a.mean() - r.mean())
    return np.percentile(vals, [2.5, 97.5])


def main():
    w, lab = load()
    summary = []
    for stim in STIMS:
        d = prep(w, lab, stim)
        v = d[d.vis & ~d.group.eq("shared")].copy()
        v = v[v.group.isin(["act_specific", "rest_specific"])]
        if len(v) < 10:
            print(f"  {stim}：可见特异位点 {len(v)}，跳过检验")
            continue
        rows = []
        for g, x in v.groupby("group"):
            rows.append(dict(组=g, n=len(x), act大于rest=(x.a > x.r).mean(),
                             几何均值=2 ** x.lr.mean()))
        obs, p = perm_diff(v)
        ci = boot_diff(v)
        print(f"\n=== {stim}（可见特异 n={len(v)}）===")
        print(pd.DataFrame(rows).round(3).to_string(index=False))
        print(f"  差值（log2 均值）= {obs:+.3f}，重采样 95% CI [{ci[0]:+.3f}, {ci[1]:+.3f}]，"
              f"置换 p={p:.4f}")
        summary.append(dict(stim=stim, n_act=int((v.group == "act_specific").sum()),
                            n_rest=int((v.group == "rest_specific").sum()),
                            diff=round(obs, 3), ci_lo=round(ci[0], 3),
                            ci_hi=round(ci[1], 3), p=round(p, 4)))

    s = pd.DataFrame(summary)
    print("\n=== 三刺激差值排序（预期 IFNg > LPS_24h > LPS_2h）===")
    print(s.to_string(index=False))
    if len(s) == 3:
        ifng = s[s.stim == "IFNg_24h"].iloc[0]
        lps2 = s[s.stim == "LPS_2h"].iloc[0]
        if lps2["diff"] >= ifng["diff"]:
            print("→ LPS_2h 差值不低于 IFNg：按预注册第 8 条，**不支持**")
        elif lps2["ci_hi"] < ifng["ci_lo"]:
            print("→ IFNg 明显更大且区间不重叠：支持排序")
        else:
            print("→ IFNg 更大但区间重叠：只能写「方向符合预期」")
    s.to_csv("b_data/fairfax_summary.csv", index=False)


if __name__ == "__main__":
    main()
