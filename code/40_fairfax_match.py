#!/usr/bin/env python3
"""Fairfax 匹配度（PREREG_v2 补充D §5.1 对称口径）。

E-MTAB-2232 原始表达矩阵，按供体配对算每基因 Δ（刺激 − naive），
再与 Kang IFN-β 的 Δ（基因级）求 Spearman。
探针→基因：本地 Fairfax 四个可信集文件的 (molecular_trait_id, gene_id)；
一基因多探针时基因 Δ 取探针 Δ 的中位数。
缺哪个条件的文件，该臂记「不可评估」并跳过（不默认匹配）。
"""
import glob
import os

import numpy as np
import pandas as pd
from scipy import stats

SEED = 20261114
D = "mpra"
COND_FILES = {
    "naive": "CD14.47231.414.b.txt",
    "IFNg_24h": "IFN.47231.367.b.txt",
    "LPS_2h": "LPS2.47231.261.b.txt",
    "LPS_24h": "LPS24.47231.322.b.txt",
}


def load_probes(path):
    d = pd.read_csv(path, sep="\t", low_memory=False)
    d = d.rename(columns={"PROBE_ID": "probe"}).set_index("probe")
    d.columns = [str(c).strip() for c in d.columns]
    return d.apply(pd.to_numeric, errors="coerce")


def probe_gene_map():
    maps = []
    for qtd in ["QTD000081", "QTD000082", "QTD000083", "QTD000084"]:
        p = f"eqtl_cache/{qtd}.credible_sets.tsv.gz"
        if os.path.exists(p):
            cs = pd.read_csv(p, sep="\t", usecols=["molecular_trait_id", "gene_id"])
            maps.append(cs.drop_duplicates())
    m = pd.concat(maps).drop_duplicates()
    m.columns = ["probe", "gene"]
    return m


def main():
    for f in COND_FILES.values():
        if not os.path.exists(f"{D}/{f}"):
            print(f"缺 {f}（该臂将记为不可评估）")
    a = load_probes(f"{D}/{COND_FILES['naive']}")
    print(f"naive {a.shape}", flush=True)
    pm = probe_gene_map()
    print(f"探针→基因映射 {len(pm)} 对 / {pm.gene.nunique()} 基因")

    obs = (pd.concat([pd.read_csv(f) for f in glob.glob("b_data/obs_cpm_pooled_*.csv")])
           .drop_duplicates("ensg"))
    d_ifn = np.log2((obs.cpm_stim_full + 1) / (obs.cpm_ctrl_full + 1))
    d_ifn.index = obs.ensg.values
    lab = pd.read_csv("b_data/fairfax_labels_fixed.csv")

    out = []
    for stim, fname in COND_FILES.items():
        if stim == "naive" or not os.path.exists(f"{D}/{fname}"):
            continue
        b = load_probes(f"{D}/{fname}")
        donors = sorted(set(a.columns) & set(b.columns),
                        key=lambda x: int(x) if str(x).isdigit() else 0)
        probes = a.index.intersection(b.index)
        print(f"{stim}: {b.shape}，共同供体 {len(donors)}，共同探针 {len(probes)}",
              flush=True)
        delta = (b.loc[probes, donors] - a.loc[probes, donors]).median(axis=1)
        dd = delta.reset_index()
        dd.columns = ["probe", "d_mean"]
        dd = dd.merge(pm, on="probe", how="inner")
        g = dd.groupby("gene").d_mean.median()

        genes = set(lab[lab.stim == stim].gene)
        use = [x for x in genes if x in g.index and x in d_ifn.index
               and pd.notna(g[x]) and pd.notna(d_ifn[x])]
        if len(use) < 30:
            print(f"  可用基因 {len(use)}，跳过")
            continue
        ds, di = g.loc[use].values, d_ifn.loc[use].values
        rho, p = stats.spearmanr(ds, di)
        rng = np.random.default_rng(SEED)
        boots = []
        for _ in range(1000):
            ix = rng.integers(0, len(use), len(use))
            if np.std(ds[ix]) == 0 or np.std(di[ix]) == 0:
                continue
            boots.append(stats.spearmanr(ds[ix], di[ix])[0])
        lo, hi = np.percentile(boots, [2.5, 97.5])
        neg = bool(rho < 0 and hi < 0)
        out.append(dict(study="Fairfax_2014", stim=stim, n_genes=len(use),
                        spearman=rho, p=p, ci95_lo=lo, ci95_hi=hi,
                        clearly_negative=neg))
        print(f"  n={len(use):4d}  r={rho:+.3f} (p={p:.4f})  "
              f"CI [{lo:+.3f}, {hi:+.3f}]  明显为负={neg}")

    pd.DataFrame(out).to_csv("b_data/fairfax_match_scores.csv", index=False)
    print("\n写入 b_data/fairfax_match_scores.csv")


if __name__ == "__main__":
    main()
