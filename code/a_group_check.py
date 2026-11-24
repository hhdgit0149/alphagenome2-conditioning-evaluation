#!/usr/bin/env python3
"""A 组复核：激活特异 vs 静息特异位点，在两种条件下预测效应强度的对比。

用法：
    python a_group_check.py tcell_per_variant_long.csv
    python a_group_check.py tcell_per_variant_long.csv --min-expr 0.1   # pred_ref_expr 补上之后再用

输入是长表，每行一个「位点 x 条件」，需要列：
    variant, gene, study, cell, group, cond_role(act/rest), pred_raw_score
    可选：pred_ref_expr（预测的参考等位基因表达量，CPM）

主检验（事先固定，不要看了结果再改）：
    只取模型看得见的位点（任一条件下 |效应| > 阈值），
    按「变异-基因」去重并剔除分组标签互相矛盾的重复项，
    比较 act_specific 与 rest_specific 两组的 log2(|act| / |rest|)，
    用研究内置换检验。shared 组只用来估计基线，基线不一定是 0.5。
"""
import argparse
import sys

import numpy as np
import pandas as pd
from scipy import stats

KEY = ["variant", "gene", "study", "cell", "group"]


def load(path, thr, min_expr):
    df = pd.read_csv(path)
    w = df.pivot_table(index=KEY, columns="cond_role", values="pred_raw_score",
                       aggfunc="first").reset_index().dropna(subset=["act", "rest"])
    w["absA"], w["absR"] = w["act"].abs(), w["rest"].abs()
    n_all = len(w)
    if min_expr is not None:
        if df["pred_ref_expr"].notna().sum() == 0:
            sys.exit("pred_ref_expr 这一列是空的，不能加表达量门槛。")
        e = df.pivot_table(index=KEY, columns="cond_role", values="pred_ref_expr",
                           aggfunc="first").reset_index()
        e = e.rename(columns={"act": "exprA", "rest": "exprR"})
        w = w.merge(e, on=KEY)
        w = w[(w["exprA"] >= min_expr) & (w["exprR"] >= min_expr)]
        print(f"表达量门槛：两种条件下都 >= {min_expr}，保留 {len(w)}/{n_all}")
    w["vis"] = w[["absA", "absR"]].max(axis=1) > thr
    return w


def dedup(v):
    """同一个变异-基因在多个细胞类型或研究里出现时预测值相同，只算一次。"""
    n_groups = v.groupby(["variant", "gene"])["group"].nunique()
    ok = n_groups[n_groups == 1].index
    u = v.set_index(["variant", "gene"]).loc[ok].reset_index()
    u = u.drop_duplicates(["variant", "gene"]).copy()
    u["lr"] = np.log2(u["absA"] / u["absR"])
    u["a_gt_r"] = u["absA"] > u["absR"]
    return u, int((n_groups > 1).sum())


def perm_within_study(x, n_perm, seed=0):
    def diff(labels):
        return x.loc[labels == "act_specific", "lr"].mean() - \
               x.loc[labels == "rest_specific", "lr"].mean()
    obs = diff(x["group"])
    rng = np.random.default_rng(seed)
    hits = 0
    for _ in range(n_perm):
        lab = x.groupby("study")["group"].transform(lambda s: rng.permutation(s.values))
        hits += abs(diff(lab)) >= abs(obs)
    return obs, (hits + 1) / (n_perm + 1)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("table")
    ap.add_argument("--thr", type=float, default=0.1, help="可见位点的效应阈值")
    ap.add_argument("--min-expr", type=float, default=None,
                    help="两种条件下预测表达量的下限（需要 pred_ref_expr）")
    ap.add_argument("--n-perm", type=int, default=20000)
    args = ap.parse_args()

    w = load(args.table, args.thr, args.min_expr)
    print(f"有预测值的位点：{len(w)}；看得见的：{int(w['vis'].sum())}"
          f"（{w['vis'].mean():.1%}）")
    print("\n各组可见比例：")
    print(w.groupby("group")["vis"].agg(比例="mean", 可见="sum", 总数="count")
          .round(3).to_string())

    v = w[w["vis"]]
    print(f"\n两种条件下预测值的相关系数（可见位点）："
          f"{np.corrcoef(v['act'], v['rest'])[0, 1]:.3f}；"
          f"方向相反的位点数：{int((np.sign(v['act']) != np.sign(v['rest'])).sum())}")

    u, n_conflict = dedup(v)
    print(f"\n去重：{len(v)} 行 -> {len(u)} 个唯一的变异-基因"
          f"（剔除分组标签矛盾的 {n_conflict} 个）")

    print("\n各组结果（去重后）。注意基线看 shared 组，不是 0.5：")
    rows = []
    for g, d in u.groupby("group"):
        rows.append({"组": g, "n": len(d), "act>rest 比例": round(d["a_gt_r"].mean(), 3),
                     "|act|/|rest| 几何均值": round(2 ** d["lr"].mean(), 2),
                     "中位数": round(2 ** d["lr"].median(), 2)})
    print(pd.DataFrame(rows).to_string(index=False))

    def contrast(a, b):
        da, db = u[u["group"] == a], u[u["group"] == b]
        if len(da) < 3 or len(db) < 3:
            return
        tab = [[int(da["a_gt_r"].sum()), int((~da["a_gt_r"]).sum())],
               [int(db["a_gt_r"].sum()), int((~db["a_gt_r"]).sum())]]
        print(f"  {a} vs {b}: Fisher p={stats.fisher_exact(tab)[1]:.4f}，"
              f"Mann-Whitney p={stats.mannwhitneyu(da['lr'], db['lr']).pvalue:.4f}")

    print("\n组间对比：")
    contrast("act_specific", "rest_specific")
    contrast("act_specific", "shared")
    contrast("rest_specific", "shared")

    x = u[u["group"].isin(["act_specific", "rest_specific"])].reset_index(drop=True)
    obs, p = perm_within_study(x, args.n_perm)
    print(f"\n主检验（研究内置换）：两组 log2 比值的均值差 = {obs:.3f}，p = {p:.4f}")

    print("\n分研究看：")
    for s, d in x.groupby("study"):
        a = d[d["group"] == "act_specific"]["lr"]
        r = d[d["group"] == "rest_specific"]["lr"]
        p_s = stats.mannwhitneyu(a, r).pvalue if len(a) > 2 and len(r) > 2 else float("nan")
        print(f"  {s}: 激活特异 n={len(a)} 中位 {a.median():+.3f}；"
              f"静息特异 n={len(r)} 中位 {r.median():+.3f}；p={p_s:.4f}")

    print("\n阈值敏感性（去重后 act>rest 比例）：")
    for thr in (0.02, 0.05, 0.1, 0.2):
        uu, _ = dedup(w[w[["absA", "absR"]].max(axis=1) > thr])
        parts = [f"{g} {int(d['a_gt_r'].sum())}/{len(d)}={d['a_gt_r'].mean():.2f}"
                 for g, d in uu.groupby("group")]
        print(f"  阈值 {thr}: " + " | ".join(parts))


if __name__ == "__main__":
    main()
