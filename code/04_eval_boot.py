#!/usr/bin/env python3
"""bootstrap 版评估：方向一致率、Spearman、AUROC 都给 95% 置信区间。

用法: python 04_eval_boot.py variants_v2.csv scores_v2.csv summary_v2.csv
"""
import sys

import numpy as np
import pandas as pd
from scipy import stats

V = sys.argv[1] if len(sys.argv) > 1 else "variants_v2.csv"
S = sys.argv[2] if len(sys.argv) > 2 else "scores_v2.csv"
OUT = sys.argv[3] if len(sys.argv) > 3 else "summary_v2.csv"
B = 2000
RNG = np.random.default_rng(20261113)


def boot_ci(fn, n, b=B):
    vals = []
    for _ in range(b):
        idx = RNG.integers(0, n, n)
        try:
            vals.append(fn(idx))
        except Exception:
            pass
    vals = np.array([v for v in vals if np.isfinite(v)])
    if len(vals) < 50:
        return (np.nan, np.nan)
    return (float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5)))


def track_stats(sub, min_n=15):
    pos = sub[sub["significant"] == 1]
    neg = sub[sub["significant"] == 0]
    use = pos[(pos["effect"] != 0) & (pos["raw_score"] != 0)]
    if len(use) < min_n:
        return None
    e = use["effect"].to_numpy()
    r = use["raw_score"].to_numpy()
    agree = float((np.sign(r) == np.sign(e)).mean())
    rho = float(stats.spearmanr(r, e).correlation)
    n = len(use)

    ci_agree = boot_ci(lambda idx: float((np.sign(r[idx]) == np.sign(e[idx])).mean()), n)
    ci_rho = boot_ci(lambda idx: float(stats.spearmanr(r[idx], e[idx]).correlation), n)

    rec = dict(n=n, 方向一致率=round(agree, 3),
               一致率_lo=round(ci_agree[0], 3), 一致率_hi=round(ci_agree[1], 3),
               spearman=round(rho, 3),
               rho_lo=round(ci_rho[0], 3), rho_hi=round(ci_rho[1], 3),
               二项检验p=float(stats.binomtest(int(round(agree * n)), n, 0.5).pvalue))

    if len(neg) >= min_n:
        pr = pos["raw_score"].abs().to_numpy()
        nr = neg["raw_score"].abs().to_numpy()

        def auroc(pi, ni):
            u = stats.mannwhitneyu(pr[pi], nr[ni], alternative="two-sided").statistic
            return float(u / (len(pi) * len(ni)))

        rec["auroc"] = round(auroc(np.arange(len(pr)), np.arange(len(nr))), 3)
        vals = []
        for _ in range(500):
            pi = RNG.integers(0, len(pr), len(pr))
            ni = RNG.integers(0, len(nr), len(nr))
            vals.append(auroc(pi, ni))
        rec["auroc_lo"] = round(float(np.percentile(vals, 2.5)), 3)
        rec["auroc_hi"] = round(float(np.percentile(vals, 97.5)), 3)
    return rec


def main():
    v = pd.read_csv(V, dtype={"variant_id": str})
    s = pd.read_csv(S, dtype={"input_id": str}, low_memory=False)
    df = s.merge(v[["variant_id", "effect", "significant"]],
                 left_on="input_id", right_on="variant_id", suffixes=("_ag", ""))
    print(f"位点 {df['input_id'].nunique()}（病例 {v[v.significant == 1]['variant_id'].nunique()}，"
          f"对照 {v[v.significant == 0]['variant_id'].nunique()}），评分行 {len(df)}")

    groups = ["output_type", "biosample_name", "transcription_factor", "histone_mark"]
    rows = []
    for keys, g in df.groupby(groups, dropna=False):
        rec = track_stats(g)
        if rec:
            rows.append(dict(zip(groups, keys), **rec))
    per = pd.DataFrame(rows).sort_values("方向一致率", ascending=False)

    pooled = []
    for ot, g in df.groupby("output_type"):
        rec = track_stats(g)
        if rec:
            pooled.append(dict(output_type=ot, level="按输出类型汇总", **rec))
    pooled = pd.DataFrame(pooled)

    out = pd.concat([per, pooled], ignore_index=True)
    out.to_csv(OUT, index=False)
    cols = ["output_type", "biosample_name", "transcription_factor", "histone_mark",
            "n", "方向一致率", "一致率_lo", "一致率_hi", "spearman", "rho_lo", "rho_hi",
            "auroc", "auroc_lo", "auroc_hi", "level"]
    print("\n" + out[[c for c in cols if c in out.columns]].head(25).to_string(index=False))
    print(f"\n写出 {OUT}")


if __name__ == "__main__":
    main()
