#!/usr/bin/env python3
"""刺激特异位点的条件化检验：按刺激分别复算主检验、安慰剂和稳健性。

用法：
    python stim_site_check.py fairfax_long_table.csv
    python stim_site_check.py 长表.csv --primary-stim IFNg_24h --thr 0.1

长表每行一个「变异 x 基因 x 刺激 x 条件」，需要列：
    variant, gene, stim, group, cond_role, dn_lfc
    group 取 act_specific / rest_specific / shared，必须是相对该行 stim 的标签
    cond_role 取 刺激 / 对照 / 对照对半A / 对照对半B（或 stim / ctrl / hA / hB）
可选列：rna_lfc, dn_act, cpm_ctrl_full, cpm_stim_full, used_in_visibility_check

所有 p 值都是双侧。
主检验用秩和检验（Mann-Whitney）：它只看排序，对数比和有界指标给出同一个 p，
也不会被单个极端值左右。均值差的置换检验只作参考。
"""
import argparse
import sys

import numpy as np
import pandas as pd
from scipy import stats

ROLE = {"刺激": "stim", "对照": "ctrl", "对照对半A": "hA", "对照对半B": "hB",
        "stim": "stim", "ctrl": "ctrl", "hA": "hA", "hB": "hB"}
ACT, REST, SHARED = "act_specific", "rest_specific", "shared"
KEY = ["variant", "gene", "stim", "group"]


def load(path):
    df = pd.read_csv(path)
    df["role"] = df["cond_role"].map(ROLE)
    if df["role"].isna().any():
        sys.exit(f"cond_role 有不认识的取值：{sorted(df.loc[df['role'].isna(), 'cond_role'].unique())}")
    dup = df.duplicated(["variant", "gene", "stim", "role"]).sum()
    if dup:
        sys.exit(f"主键（变异+基因+刺激+条件）有 {dup} 行重复，先查清再算。")
    parts = []
    for col in ("dn_lfc", "rna_lfc", "dn_act"):
        if col in df.columns:
            w = df.pivot_table(index=KEY, columns="role", values=col, aggfunc="first")
            w.columns = [f"{col}_{c}" for c in w.columns]
            parts.append(w)
    wide = pd.concat(parts, axis=1).reset_index()
    extra = [c for c in ("cpm_ctrl_full", "cpm_stim_full", "used_in_visibility_check")
             if c in df.columns]
    if extra:
        meta = df.drop_duplicates(["variant", "gene", "stim"])[["variant", "gene", "stim"] + extra]
        wide = wide.merge(meta, on=["variant", "gene", "stim"], how="left")
    return wide


def label_checks(wide):
    """专门用来抓「跨刺激去重」这类建表错误。"""
    print("== 分组标签检查 ==")
    n = wide.groupby(["variant", "gene", "stim"])["group"].nunique()
    print(f"同一刺激内同一变异-基因有多个标签的：{int((n > 1).sum())}（应为 0）")
    print("各刺激各组位点数：")
    print(pd.crosstab(wide["stim"], wide["group"]).to_string())
    stims = sorted(wide["stim"].unique())
    for g in (REST, SHARED):
        sets = {s: set(map(tuple, wide[(wide["stim"] == s) & (wide["group"] == g)]
                           [["variant", "gene"]].values)) for s in stims}
        pairs = [f"{a}∩{b}={len(sets[a] & sets[b])}"
                 for i, a in enumerate(stims) for b in stims[i + 1:]]
        print(f"{g} 在不同刺激间的重叠：" + "，".join(pairs))
    print("提示：静息特异和共享位点在不同刺激间重叠为 0 是不正常的，"
          "多半是建表时跨刺激去重了。\n")


def prep(d, a, b, thr, head="dn_lfc"):
    ca, cb = f"{head}_{a}", f"{head}_{b}"
    if ca not in d.columns or cb not in d.columns:
        return d.iloc[0:0].assign(lr=[], idx=[])
    x = d.dropna(subset=[ca, cb]).copy()
    x["A"], x["B"] = x[ca].abs(), x[cb].abs()
    x = x[x[["A", "B"]].max(axis=1) > thr].copy()
    x["lr"] = np.log2(x["A"] / x["B"])                    # 对数比，分母接近 0 时会很大
    x["idx"] = (x["A"] - x["B"]) / (x["A"] + x["B"])      # 有界指标，范围 -1 到 1
    return x


def perm_mean(a, r, rng, n_perm):
    obs = a.mean() - r.mean()
    pool, na, hits = np.concatenate([a, r]), len(a), 0
    for _ in range(n_perm):
        p = rng.permutation(pool)
        hits += abs(p[:na].mean() - p[na:].mean()) >= abs(obs)
    return obs, (hits + 1) / (n_perm + 1)


def contrast(x, label, rng, n_perm=10000, col="lr"):
    a = x.loc[x["group"] == ACT, col].values
    r = x.loc[x["group"] == REST, col].values
    if len(a) < 3 or len(r) < 3:
        print(f"{label:34s} 位点太少（{len(a)}/{len(r)}），不检验")
        return None
    mw = stats.mannwhitneyu(a, r, alternative="two-sided").pvalue
    diff, p_perm = perm_mean(a, r, rng, n_perm)
    print(f"{label:34s} n={len(a):3d}/{len(r):3d}  中位数 {np.median(a):+.3f} / {np.median(r):+.3f}"
          f"  秩和 p={mw:.4f} | 均值差 {diff:+.3f} 置换 p={p_perm:.4f}")
    return mw


def influence(x):
    """逐个去掉位点，看秩和检验 p 最差变成多少，并列出影响最大的位点。"""
    sub = x[x["group"].isin([ACT, REST])].reset_index(drop=True)
    base = stats.mannwhitneyu(sub.loc[sub.group == ACT, "lr"], sub.loc[sub.group == REST, "lr"]).pvalue
    rows = []
    for i in range(len(sub)):
        d = sub.drop(index=i)
        a, r = d.loc[d.group == ACT, "lr"], d.loc[d.group == REST, "lr"]
        rows.append((sub.loc[i, "variant"], sub.loc[i, "gene"], sub.loc[i, "group"],
                     sub.loc[i, "lr"], stats.mannwhitneyu(a, r).pvalue, a.mean() - r.mean()))
    out = pd.DataFrame(rows, columns=["variant", "gene", "group", "lr", "去掉后秩和p", "去掉后均值差"])
    print(f"逐一剔除：秩和 p 原值 {base:.4f}，最差 {out['去掉后秩和p'].max():.4f}，"
          f"最好 {out['去掉后秩和p'].min():.4f}")
    print("对数比绝对值最大的 5 个位点：")
    print(out.reindex(out["lr"].abs().sort_values(ascending=False).index).head(5)
          .round(4).to_string(index=False))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("table")
    ap.add_argument("--primary-stim", default="IFNg_24h")
    ap.add_argument("--thr", type=float, default=0.1)
    ap.add_argument("--n-perm", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)

    wide = load(args.table)
    label_checks(wide)
    stims = [args.primary_stim] + [s for s in sorted(wide["stim"].unique()) if s != args.primary_stim]

    print(f"== 主检验：刺激 vs 对照，可见阈值 {args.thr}，每种刺激单独算 ==")
    vis = {}
    for s in stims:
        vis[s] = prep(wide[wide["stim"] == s], "stim", "ctrl", args.thr)
        contrast(vis[s], s, rng, args.n_perm)

    x = vis[args.primary_stim]
    print(f"\n== {args.primary_stim}：三组概况（基线看 shared，不是 0）==")
    print(x.groupby("group").agg(n=("lr", "size"), 对数比中位数=("lr", "median"),
                                 有界指标中位数=("idx", "median"),
                                 刺激更大的比例=("lr", lambda v: (v > 0).mean())).round(3).to_string())
    order = x["group"].map({REST: 0, SHARED: 1, ACT: 2})
    rho, p = stats.spearmanr(order, x["lr"])
    print(f"趋势检验（静息特异 < 共享 < 激活特异）：rho={rho:.3f}，p={p:.4f}")
    for a, b in ((ACT, SHARED), (REST, SHARED)):
        u, v = x.loc[x.group == a, "lr"], x.loc[x.group == b, "lr"]
        if len(u) > 2 and len(v) > 2:
            print(f"  {a} vs {b}：秩和 p={stats.mannwhitneyu(u, v).pvalue:.4f}")

    print(f"\n== {args.primary_stim}：单个位点的影响 ==")
    influence(x)

    print("\n== 安慰剂：对照对半 A vs B ==")
    for s in stims:
        contrast(prep(wide[wide["stim"] == s], "hA", "hB", args.thr), s + "（各自判可见）", rng, args.n_perm)
    m = x[["variant", "gene"]].merge(wide[wide["stim"] == args.primary_stim], on=["variant", "gene"])
    contrast(prep(m, "hA", "hB", 0.0), args.primary_stim + "（沿用主可见集合）", rng, args.n_perm)

    print(f"\n== {args.primary_stim}：稳健性（方向和 p 都要看，样本量会变小）==")
    d0 = wide[wide["stim"] == args.primary_stim]
    for thr in (0.05, 0.15, 0.2):
        contrast(prep(d0, "stim", "ctrl", thr), f"可见阈值 {thr}", rng, args.n_perm)
    contrast(x[(x.A > 0.05) & (x.B > 0.05)], "两种条件下都 >0.05", rng, args.n_perm)
    if "rna_lfc_stim" in wide.columns:
        contrast(prep(d0, "stim", "ctrl", args.thr, head="rna_lfc"), "换普通转录组头", rng, args.n_perm)
    if "used_in_visibility_check" in x.columns:
        contrast(x[~x["used_in_visibility_check"].astype(bool)], "去掉已用过的位点", rng, args.n_perm)
    if {"cpm_ctrl_full", "cpm_stim_full"} <= set(x.columns):
        contrast(x[(x.cpm_ctrl_full >= 1) & (x.cpm_stim_full >= 1)], "观测表达量都 >=1 CPM", rng, args.n_perm)
        fc = np.log2((x.cpm_stim_full + 1) / (x.cpm_ctrl_full + 1))
        contrast(x[fc.abs() < 1], "观测表达变化 <2 倍", rng, args.n_perm)
        ok = fc.notna()
        rho, p = stats.spearmanr(x.loc[ok, "lr"], fc[ok])
        print(f"效应比与目标基因表达变化的相关：rho={rho:.3f}，p={p:.4f}")

    others = [s for s in stims if s != args.primary_stim]
    if others:
        print(f"\n== 探索性：同在该条件下，{args.primary_stim} 特异位点 vs 其他刺激的特异位点 ==")
        a = x[x.group == ACT]
        mine = set(map(tuple, a[["variant", "gene"]].values))
        for s in others:
            c = vis[s][vis[s].group == ACT]
            c = c[~c.set_index(["variant", "gene"]).index.isin(mine)]
            if len(a) > 2 and len(c) > 2:
                p = stats.mannwhitneyu(a["lr"], c["lr"]).pvalue
                print(f"  vs {s}: n={len(a)}/{len(c)}  中位数 {a.lr.median():+.3f} / {c.lr.median():+.3f}  秩和 p={p:.4f}")

    print("\n== 可见比例 ==")
    w = wide.dropna(subset=["dn_lfc_stim", "dn_lfc_ctrl"]).copy()
    w["vis"] = w[["dn_lfc_stim", "dn_lfc_ctrl"]].abs().max(axis=1) > args.thr
    print(w.groupby(["stim", "group"])["vis"].agg(比例="mean", 可见="sum", 总数="count").round(3).to_string())


if __name__ == "__main__":
    main()
