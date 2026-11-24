#!/usr/bin/env python3
"""构建严格/宽松两套长表（Fairfax 与 Quach），供 stim_site_check.py 复算。

输出：
  b_data/fairfax_long_strict.csv   Fairfax 严格标签（主口径）
  b_data/quach_long_fixed.csv      Quach 宽松标签（敏感性）
  b_data/quach_long_strict.csv     Quach 严格标签（主口径）
列结构与 fairfax_long_fixed.csv 完全一致。
"""
import glob

import pandas as pd

ROLE_CN = {"stim": "刺激", "ctrl": "对照", "halfA": "对照对半A", "halfB": "对照对半B"}
COLS = ["variant", "gene", "study", "cell", "stim", "group", "cond_role",
        "dn_lfc", "rna_lfc", "dn_act", "rna_act", "n_tracks",
        "exp_beta", "exp_z", "exp_pip", "used_in_visibility_check",
        "cpm_ctrl_full", "cpm_stim_full"]

# 打分来源
fa_sc = pd.concat([pd.read_csv("b_data/fairfax_v2_dedup.csv")] +
                  [pd.read_csv(f) for f in sorted(glob.glob("b_data/fairfax_v2b_*.csv"))],
                  ignore_index=True).drop_duplicates(["variant", "gene", "cond_role"])
qu_files = sorted(glob.glob("b_data/quach_v2_*.csv"))
qu_sc = pd.concat([pd.read_csv(f) for f in qu_files],
                  ignore_index=True).drop_duplicates(["variant", "gene", "cond_role"])
assert not qu_sc.duplicated(["variant", "gene", "cond_role"]).any(), "Quach 打分主键重复"

obs = (pd.concat([pd.read_csv(f) for f in glob.glob("b_data/obs_cpm_pooled_*.csv")])
       .drop_duplicates("ensg"))
used = set()
for f in glob.glob("b_data/visibility_fairfax_*.csv"):
    used |= set(pd.read_csv(f).variant.astype(str))

for labp, sc, outp, tag in [
    ("b_data/fairfax_labels_strict.csv", fa_sc, "b_data/fairfax_long_strict.csv", "fairfax-strict"),
    ("b_data/Quach_2016_labels_fixed.csv", qu_sc, "b_data/quach_long_fixed.csv", "quach-loose"),
    ("b_data/quach_labels_strict.csv", qu_sc, "b_data/quach_long_strict.csv", "quach-strict"),
]:
    lab = pd.read_csv(labp)
    if "group_strict" in lab.columns:
        lab = lab.drop(columns=["group"]).rename(columns={"group_strict": "group"})
    lab = lab[lab.group.isin(["act_specific", "rest_specific", "shared"])]
    long = lab.merge(sc, on=["variant", "gene"], how="left")
    long["used_in_visibility_check"] = long.variant.astype(str).isin(used)
    long = long.merge(obs[["ensg", "cpm_ctrl_full", "cpm_stim_full"]],
                      left_on="gene", right_on="ensg", how="left").drop(columns=["ensg"])
    long = long.rename(columns={"beta": "exp_beta", "z": "exp_z", "pip": "exp_pip"})
    long["cond_role"] = long.cond_role.map(ROLE_CN)
    assert long.cond_role.notna().all(), f"[{tag}] cond_role 映射失败"
    long = long[[c for c in COLS if c in long.columns]]
    dup = int(long.duplicated(["variant", "gene", "stim", "cond_role"]).sum())
    assert dup == 0, f"[{tag}] 长表主键重复 {dup}"
    long.to_csv(outp, index=False)
    nmiss = int(long.dn_lfc.isna().sum())
    print(f"[{tag}] {outp}：{len(long)} 行，打分缺失 {nmiss}，"
          f"重叠已用位点 {int(long.used_in_visibility_check.sum())}")
    print(pd.crosstab(long.stim, long.group).to_string(), "\n")

# Quach 与 Fairfax 位点是否完全不重叠（PREREG_v2 第 4 节要求）
fa_lab = pd.read_csv("b_data/fairfax_labels_fixed.csv")
qu_lab = pd.read_csv("b_data/Quach_2016_labels_fixed.csv")
ov = len(set(zip(fa_lab.variant, fa_lab.gene)) & set(zip(qu_lab.variant, qu_lab.gene)))
print(f"Fairfax 与 Quach 重叠的变异-基因：{ov}（应为 0）")
