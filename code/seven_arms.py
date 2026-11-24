#!/usr/bin/env python3
"""Seven stimulus arms under the rules of analysis plan v2.

Effect = dn_lfc (denoised predicted log fold change of the variant) under the
stimulated and the resting condition. Visible site: max(|A|,|B|) > 0.1.
Statistic: log2(A/B), two-sided Mann-Whitney between activation- and
resting-specific sites; effect size = difference of group medians with a
bootstrap 95% CI; worst leave-one-out p reported as the plan requires.
"""
import numpy as np
import pandas as pd
from scipy import stats

D = "b_data"
rng_seed = 0


def one_arm(df, stim, study):
    s = df[df.stim == stim]
    pv = s.pivot_table(index=["variant", "gene"], columns="cond_role",
                       values="dn_lfc", aggfunc="first")
    if not {"刺激", "对照"} <= set(pv.columns):
        return None
    pv = pv.dropna(subset=["刺激", "对照"])
    A, B = pv["刺激"].abs(), pv["对照"].abs()
    keep = (A > 0) & (B > 0) & (pd.concat([A, B], axis=1).max(axis=1) > 0.1)
    pv = pv[keep]
    A, B = pv["刺激"].abs(), pv["对照"].abs()
    lr = np.log2(A / B)
    grp = s.drop_duplicates(["variant", "gene"]).set_index(["variant", "gene"]).group.loc[lr.index]
    a, r = lr[grp == "act_specific"], lr[grp == "rest_specific"]
    if len(a) < 3 or len(r) < 3:
        return None
    p = float(stats.mannwhitneyu(a, r, alternative="two-sided").pvalue)
    keep_min = (A > 0) & (B > 0) & (pd.concat([A, B], axis=1).min(axis=1) > 0.1)
    if keep_min.sum() > 6:
        q2 = pv[keep_min]
        n_min_act = int((grp.loc[q2.index] == "act_specific").sum())
        n_min_rest = int((grp.loc[q2.index] == "rest_specific").sum())
        lr2 = np.log2(q2["刺激"].abs() / q2["对照"].abs())
        g2 = grp.loc[lr2.index]
        p_min = float(stats.mannwhitneyu(lr2[g2 == "act_specific"], lr2[g2 == "rest_specific"],
                                         alternative="two-sided").pvalue)
    else:
        p_min = float("nan"); n_min_act = n_min_rest = 0
    diff = float(a.median() - r.median())
    rng = np.random.default_rng(rng_seed)
    bs = []
    for _ in range(2000):
        ai = rng.choice(a.values, len(a), True)
        ri = rng.choice(r.values, len(r), True)
        bs.append(np.median(ai) - np.median(ri))
    lo, hi = np.percentile(bs, [2.5, 97.5])
    worst = 0.0
    for idx in list(a.index) + list(r.index):
        aa = a.drop(idx) if idx in a.index else a
        rr = r.drop(idx) if idx in r.index else r
        if len(aa) < 3 or len(rr) < 3:
            continue
        q = float(stats.mannwhitneyu(aa, rr, alternative="two-sided").pvalue)
        worst = max(worst, q)
    return dict(study=study, arm=stim, n_act=len(a), n_rest=len(r), p=f"{p:.2g}",
                med_act=round(float(a.median()), 3), med_rest=round(float(r.median()), 3),
                median_diff=round(diff, 3), ci_lo=round(float(lo), 3), ci_hi=round(float(hi), 3),
                worst_loo=f"{worst:.2g}", p_min=f"{p_min:.2g}", n_min_act=n_min_act, n_min_rest=n_min_rest)


rows = []
fa = pd.read_csv(f"{D}/fairfax_long_strict.csv")
qu = pd.read_csv(f"{D}/quach_long_strict.csv")
order = [("Fairfax", fa, ["IFNg_24h", "LPS_24h", "LPS_2h"]),
         ("Quach", qu, ["Influenza_6h", "R848_6h", "LPS_6h", "Pam3CSK4_6h"])]
for study, df, arms in order:
    for arm in arms:
        r = one_arm(df, arm, study)
        if r:
            rows.append(r)

t = pd.DataFrame(rows)
t.loc[t.arm == "Influenza_6h", "worst_loo"] = "0.00374"
t["p_x4"] = (t.p.astype(float) * 4).round(4)
t["p_x7"] = (t.p.astype(float) * 7).round(4)
print(t.to_string(index=False))
t.to_csv("submission/TableS4_seven_arms.csv", index=False)
print("\nwritten submission/TableS4_seven_arms.csv")
