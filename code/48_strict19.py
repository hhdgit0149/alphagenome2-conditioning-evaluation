#!/usr/bin/env python3
"""|z|<1.9 阈值敏感性（补充E 承诺）：重建 Quach 严格标签并重跑 Flu 臂。"""
import numpy as np
import pandas as pd
from scipy import stats

THR = 1.9
j = pd.read_csv("b_data/cs_join_quach.csv")
z = pd.read_csv("b_data/strict_z_quach.csv")
rows = []
for r in j.itertuples(index=False):
    if r.group == "shared":
        g = "shared" if (r.cs_stim_beta > 0) == (r.cs_naive_beta > 0) else "excluded"
    else:
        other = "naive" if r.group == "act_specific" else r.stim
        trait = r.cs_stim_trait if r.group == "act_specific" else r.cs_naive_trait
        pr = z[(z.other_cond == other) & (z.variant == r.variant) & (z.gene == r.gene)]
        if pr.empty or int(pr.iloc[0].n_pair_rows) == 0:
            g = "excluded"
        else:
            hit = pr[pr.trait == trait]
            g = ("excluded" if hit.empty
                 else (r.group if abs(float(hit.iloc[0].beta) / float(hit.iloc[0].se)) < THR
                       else "excluded"))
    rows.append(dict(stim=r.stim, group=g, variant=r.variant, gene=r.gene))
lab = pd.DataFrame(rows)
print(lab[lab.stim == "Influenza_6h"].group.value_counts().to_dict())

sc = pd.read_csv("b_data/quach_v2_1114_005902.csv")
sc = sc[sc.cond_role.isin(["stim", "ctrl"])]
w = (sc.pivot_table(index=["variant", "gene"], columns="cond_role",
                    values="dn_lfc", aggfunc="first")
     .rename(columns={"stim": "lfc_stim", "ctrl": "lfc_ctrl"}).reset_index())
m = lab[(lab.stim == "Influenza_6h") & lab.group.isin(["act_specific", "rest_specific"])]
m = m.merge(w, on=["variant", "gene"], how="inner")
m["A"], m["B"] = m["lfc_stim"].abs(), m["lfc_ctrl"].abs()
m = m[(m[["A", "B"]].max(axis=1) > 0.1) & (m.A > 0) & (m.B > 0)]
m["lr"] = np.log2(m.A / m.B)
a, r_ = m[m.group == "act_specific"].lr, m[m.group == "rest_specific"].lr
print(f"阈值 |z|<{THR}：n={len(a)}/{len(r_)}  中位 {a.median():+.3f}/{r_.median():+.3f}  "
      f"秩和 p={stats.mannwhitneyu(a, r_, alternative='two-sided').pvalue:.4f}")
print("参照（|z|<2.0）：n=12/18  中位 +1.504/−0.388  p=0.0019")
