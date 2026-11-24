#!/usr/bin/env python3
"""把 eQTL Catalogue 的配对可信集整理成三类位点：激活特异 / 静息特异 / 共享。

只取 PIP>0.9 的单碱基变异（每可信集取 PIP 最高者）。
  act_specific  : 刺激条件里有，且该基因在静息条件没有任何可信集
  rest_specific : 静息条件里有，且该基因在刺激条件没有任何可信集
  shared        : 同一 (基因, 变异) 两边都在

输出 stim_variants_groups.csv
"""
import os

import pandas as pd

META = "eqtl_cache/dataset_metadata_r8.tsv"
CS = "eqtl_cache/{d}.credible_sets.tsv.gz"
PIP = 0.9
BASE = {"naive"}
# 只保留有内置条件可用的研究；Fairfax 是 microarray，之前被 quant_method 筛掉了
STUDIES = {
    "Cytoimmgen": "T cell",
    "Schmiedel_2018": "T cell",
    "Fairfax_2014": "monocyte",
    "Quach_2016": "monocyte",
}


def load(d):
    p = CS.format(d=d)
    if not os.path.exists(p) or os.path.getsize(p) < 200:
        return None
    df = pd.read_csv(p, sep="\t")
    v = df["variant"].astype(str).str.split("_")
    df["snv"] = v.str.len().eq(4) & v.str[-2].str.len().eq(1) & v.str[-1].str.len().eq(1)
    return df


def hi(df):
    top = df.sort_values("pip", ascending=False).drop_duplicates("cs_id")
    return top[(top.pip > PIP) & top.snv]


def parse_var(s):
    p = str(s).split("_")
    return p[0], int(p[1]), p[2], p[3]


def main():
    meta = pd.read_csv(META, sep="\t")
    rows = []
    for (study, cell), x in meta.groupby(["study_label", "tissue_label"]):
        if study not in STUDIES:
            continue
        b = x[x.condition_label.astype(str).str.lower().isin(BASE)]
        st = x[~x.condition_label.astype(str).str.lower().isin(BASE)]
        if b.empty or st.empty:
            continue
        for sid_, scond in zip(st.dataset_id, st.condition_label):
            bcs = load(b.iloc[0].dataset_id)
            scs = load(sid_)
            if bcs is None or scs is None:
                continue
            bh, sh = hi(bcs), hi(scs)
            bkeys = set(zip(bh.gene_id, bh.variant))
            skeys = set(zip(sh.gene_id, sh.variant))
            bgene = set(bh.gene_id)
            sgene = set(sh.gene_id)
            for tag, sub in [("act_specific", sh[~sh.gene_id.isin(bgene)]),
                             ("rest_specific", bh[~bh.gene_id.isin(sgene)]),
                             ("shared", sh[[(g, v) in bkeys for g, v in
                                            zip(sh.gene_id, sh.variant)]])]:
                for r in sub.itertuples(index=False):
                    rows.append(dict(group=tag, study=study, cell=cell,
                                     condition=scond, gene=r.gene_id,
                                     variant=r.variant, pip=float(r.pip),
                                     beta=float(getattr(r, "beta", float("nan"))),
                                     z=float(getattr(r, "z", float("nan")))))
            print(f"  {study:16s} {cell[:16]:16s} {scond[:18]:18s} "
                  f"刺激 {len(sh):4d} 静息 {len(bh):4d}", flush=True)

    df = pd.DataFrame(rows)
    if df.empty:
        raise SystemExit("没有可用位点")
    # 同一研究同一位点只留一次（取 |z| 最大）
    df["absz"] = df["z"].abs()
    df = (df.sort_values("absz", ascending=False)
            .drop_duplicates(["group", "study", "variant"]))
    v = df["variant"].str.split("_", expand=True)
    df["chrom"], df["pos"], df["ref"], df["alt"] = v[0], v[1].astype(int), v[2], v[3]
    df.to_csv("stim_variants_groups.csv", index=False)
    print("\n== 三组位点数（去重后）==")
    print(df.groupby(["group", "study"]).size().to_string())
    print("\n合计:\n" + df.groupby("group").size().to_string())
    print("\n纯 T 细胞（Cytoimmgen + Schmiedel）:")
    tc = df[df.study.isin(["Cytoimmgen", "Schmiedel_2018"])]
    print(tc.groupby("group").size().to_string())


if __name__ == "__main__":
    main()
