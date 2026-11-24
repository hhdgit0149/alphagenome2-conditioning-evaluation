#!/usr/bin/env python3
"""安慰剂 + 敏感性分析（对应预注册第 3、7、12 条与三项敏感性）。

安慰剂：halfA 当作「刺激」、halfB 当作「对照」——同一种细胞、无刺激。
       可见位点单独判定（规则 7），另报沿用主集合的版本。
"""
import glob

import numpy as np
import pandas as pd
from scipy import stats

RNG = np.random.default_rng(7)
THR = 0.1
EXPR_THR = 0.850          # dn_act，由 26_calibrate_expr.py 得出（≈1 CPM）


def load():
    sc = pd.read_csv("b_data/fairfax_v2_dedup.csv")
    lab = pd.concat([pd.read_csv(f) for f in sorted(glob.glob("b_data/fairfax_labels_*.csv"))])
    w = sc.pivot_table(index=["variant", "gene"], columns="cond_role",
                       values=["dn_lfc", "rna_lfc", "dn_act"], aggfunc="first")
    w.columns = [f"{a}_{b}" for a, b in w.columns]
    return w.reset_index(), lab


def contrast(d, a_col, r_col, n=5000):
    x = d.copy()
    x["a"] = x[a_col].abs()
    x["r"] = x[r_col].abs()
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
    return dict(n_act=len(A), n_rest=len(R), diff=round(float(obs), 3),
                p=round((hits + 1) / (n + 1), 4),
                geo_a=round(float(2 ** A.mean()), 2), geo_r=round(float(2 ** R.mean()), 2))


def main():
    w, lab = load()
    l = lab[lab.condition == "IFNg_24h"].drop_duplicates(["variant", "gene", "group"])
    d = w.merge(l, on=["variant", "gene"], how="inner")
    print(f"IFNg_24h：合并后 {len(d)} 个变异-基因")

    out = []

    # 主分析（复现）
    m = d.copy()
    m["vis"] = m[["dn_lfc_stim", "dn_lfc_ctrl"]].abs().max(axis=1) > THR
    r = contrast(m[m.vis], "dn_lfc_stim", "dn_lfc_ctrl")
    print(f"\n[主分析] {r}")
    out.append(dict(分析="主分析", **r))

    # 安慰剂（可见位点单独判定）
    p = d.copy()
    p["vis"] = p[["dn_lfc_halfA", "dn_lfc_halfB"]].abs().max(axis=1) > THR
    r = contrast(p[p.vis], "dn_lfc_halfA", "dn_lfc_halfB")
    print(f"[安慰剂·单独判定] 可见 {int(p.vis.sum())}｜{r}")
    out.append(dict(分析="安慰剂(单独判定)", **r))

    # 安慰剂（沿用主集合）
    r = contrast(m[m.vis], "dn_lfc_halfA", "dn_lfc_halfB")
    print(f"[安慰剂·沿用主集合] {r}")
    out.append(dict(分析="安慰剂(沿用主集合)", **r))

    # 敏感性 1：表达量门槛（两条件都 ≥ 门槛）
    s = m[m.vis & (m.dn_act_stim >= EXPR_THR) & (m.dn_act_ctrl >= EXPR_THR)]
    r = contrast(s, "dn_lfc_stim", "dn_lfc_ctrl")
    print(f"\n[敏感性1 表达量门槛] 保留 {len(s)}｜{r}")
    out.append(dict(分析="敏感性1_表达量门槛", **r))

    # 敏感性 3：排除表达量本身变化 >2 倍
    chg = (m.dn_act_stim - m.dn_act_ctrl).abs()
    s = m[m.vis & (chg < 1.0)]
    r = contrast(s, "dn_lfc_stim", "dn_lfc_ctrl")
    print(f"[敏感性3 表达量变化<2倍] 保留 {len(s)}｜{r}")
    out.append(dict(分析="敏感性3_表达量稳定", **r))

    # 敏感性 4：换普通转录组头
    m2 = m.copy()
    m2["vis"] = m2[["rna_lfc_stim", "rna_lfc_ctrl"]].abs().max(axis=1) > THR
    r = contrast(m2[m2.vis], "rna_lfc_stim", "rna_lfc_ctrl")
    print(f"[敏感性4 转录组头] 可见 {int(m2.vis.sum())}｜{r}")
    out.append(dict(分析="敏感性4_转录组头", **r))

    # 敏感性 5：去掉上一轮可见率测试用过的位点
    vis_used = set(pd.concat([pd.read_csv(x) for x in
                              glob.glob("b_data/visibility_fairfax_*.csv")]).variant.astype(str))
    s = m[m.vis & ~m.variant.isin(vis_used)]
    r = contrast(s, "dn_lfc_stim", "dn_lfc_ctrl")
    print(f"[敏感性5 去掉已用位点] 保留 {len(s)}｜{r}")
    out.append(dict(分析="敏感性5_去已用位点", **r))

    df = pd.DataFrame(out)
    df.to_csv("b_data/fairfax_placebo_sens.csv", index=False)
    print("\n" + df.to_string(index=False))


if __name__ == "__main__":
    main()
