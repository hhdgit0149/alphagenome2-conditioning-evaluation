#!/usr/bin/env python3
"""生成用修正标签的长表，供 stim_site_check.py 复算。"""
import glob

import numpy as np
import pandas as pd

base = pd.read_csv("b_data/fairfax_v2_dedup.csv")
add = pd.concat([pd.read_csv(f) for f in sorted(glob.glob("b_data/fairfax_v2b_*.csv"))])
sc = pd.concat([base, add], ignore_index=True).drop_duplicates(["variant", "gene", "cond_role"])
print(f"打分 {len(sc)} 行，主键重复 {int(sc.duplicated(['variant','gene','cond_role']).sum())}")

lab = pd.read_csv("b_data/fairfax_labels_fixed.csv")
obs = pd.concat([pd.read_csv(f) for f in glob.glob("b_data/obs_cpm_pooled_*.csv")]).drop_duplicates("ensg")
used = set()
for f in glob.glob("b_data/visibility_fairfax_*.csv"):
    used |= set(pd.read_csv(f).variant.astype(str))

name = {"stim": "刺激", "ctrl": "对照", "halfA": "对照对半A", "halfB": "对照对半B"}
sc["cond_role"] = sc.cond_role.map(name)

long = lab.merge(sc, on=["variant", "gene"], how="left")
long["used_in_visibility_check"] = long.variant.isin(used)
long = long.merge(obs[["ensg", "cpm_ctrl_full", "cpm_stim_full"]],
                  left_on="gene", right_on="ensg", how="left").drop(columns=["ensg"])
long = long.rename(columns={"beta": "exp_beta", "z": "exp_z", "pip": "exp_pip"})
cols = ["variant", "gene", "study", "cell", "stim", "group", "cond_role",
        "dn_lfc", "rna_lfc", "dn_act", "rna_act", "n_tracks",
        "exp_beta", "exp_z", "exp_pip", "used_in_visibility_check",
        "cpm_ctrl_full", "cpm_stim_full"]
long = long[[c for c in cols if c in long.columns]]
long.to_csv("b_data/fairfax_long_fixed.csv", index=False)
print(f"\n写出 b_data/fairfax_long_fixed.csv：{len(long)} 行")
print("主键重复（变异+基因+刺激+条件）:",
      int(long.duplicated(["variant", "gene", "stim", "cond_role"]).sum()))
print(pd.crosstab(long.stim, long.group).to_string())
print("\n空值：\n" + long[["dn_lfc", "dn_act", "rna_lfc", "cpm_stim_full"]].isna().sum().to_string())
