#!/usr/bin/env python3
"""A 组分析：三类位点在「激活 / 静息」条件下的预测效应强度对比。

判据（事先定死）：
  激活特异  → |pred_act| > |pred_rest|
  静息特异  → |pred_act| < |pred_rest|
  共享      → 无差异
三组同时成立才算模型区分了状态；只看激活特异不够。
"""
import sys

import numpy as np
import pandas as pd
from scipy import stats

SC = "scores_tcell.csv"
NEG = "scores_tcell_negcond.csv"
RNG = np.random.default_rng(20261113)


def main():
    d = pd.read_csv(SC).dropna(subset=["pred_act", "pred_rest"])
    d["a"], d["r"] = d.pred_act.abs(), d.pred_rest.abs()
    d["hi"] = d.a > d.r
    print(f"位点 {len(d)}｜分组 {d.group.value_counts().to_dict()}")
    print(f"研究 {d.study.value_counts().to_dict()}\n")

    print("=== 三组主表（主指标 = |预测效应|）===")
    print(f"{'组':14s} {'n':>4s} {'|act|中位':>10s} {'|rest|中位':>10s} "
          f"{'act>rest':>9s} {'Wilcoxon p':>11s}")
    for g, x in d.groupby("group"):
        p = stats.wilcoxon(x.a, x.r).pvalue
        print(f"{g:14s} {len(x):>4d} {x.a.median():>10.5f} {x.r.median():>10.5f} "
              f"{x.hi.mean():>9.3f} {p:>11.4f}")

    print("\n=== 按研究分层（act>rest 比例，n>=15 才列）===")
    for (g, s), x in d.groupby(["group", "study"]):
        if len(x) < 15:
            continue
        print(f"  {g:14s} {s:16s} n={len(x):>3d} act>rest={x.hi.mean():.3f} "
              f"p={stats.wilcoxon(x.a, x.r).pvalue:.3f}")

    # 置换：在每项研究内打乱组标签，看主统计量是否超出随机
    print("\n=== 置换检验（研究内打乱组标签，1000 次）===")
    for g in ["act_specific", "rest_specific"]:
        obs = d[d.group == g].hi.mean()
        null = []
        for _ in range(1000):
            perm = []
            for _, x in d.groupby("study"):
                lab = RNG.permutation(x.group.values)
                perm.append(x.assign(grp=lab))
            dd = pd.concat(perm)
            sel = dd[dd.grp == g]
            if len(sel):
                null.append(sel.hi.mean())
        null = np.array(null)
        p = (null >= obs).mean() if g == "act_specific" else (null <= obs).mean()
        print(f"  {g:14s} 观测 {obs:.3f}｜零分布中位 {np.median(null):.3f}｜单侧 p={p:.3f}")

    # 阴性条件
    try:
        n = pd.read_csv(NEG)
        m = d.merge(n, on="variant", how="inner")
        if len(m) >= 30:
            m["g"] = m.pred_neg.abs()
            print(f"\n=== 阴性条件（肝细胞）配对 n={len(m)} ===")
            for grp, x in m.groupby("group"):
                if len(x) < 10:
                    continue
                print(f"  {grp:14s} n={len(x):>3d} act>rest={(x.a > x.r).mean():.3f} "
                      f"| act>neg={(x.a > x.g).mean():.3f} "
                      f"| rest>neg={(x.r > x.g).mean():.3f}")
            print("  若 act>neg 与 act>rest 接近，说明差异只是「换了条件」，不是激活本身。")
    except FileNotFoundError:
        print("\n（无阴性条件文件）")


if __name__ == "__main__":
    main()
