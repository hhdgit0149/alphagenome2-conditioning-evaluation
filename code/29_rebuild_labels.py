#!/usr/bin/env python3
"""重建 Fairfax 标签：按每种刺激独立定义，不做跨刺激去重。

原 09_build_stim_variants.py 的 bug：
    df.sort_values("absz").drop_duplicates(["group","study","variant"])
会把同一 (组, 变异) 在多种刺激下的记录只留 |z| 最大的一条，
于是「共享」和「静息特异」被随机归给了某一种刺激，
造成 shared 数量按 IFNg 141 / LPS24h 53 / LPS2h 19 递减的假象。

修正：只在 (刺激, 组, 变异, 基因) 内去重。
"""
import os

import pandas as pd

META = "eqtl_cache/dataset_metadata_r8.tsv"
CS = "eqtl_cache/{d}.credible_sets.tsv.gz"
PIP = 0.9


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


def main():
    meta = pd.read_csv(META, sep="\t")
    rows = []
    for (study, cell), x in meta.groupby(["study_label", "tissue_label"]):
        if study != "Fairfax_2014":
            continue
        base = x[x.condition_label == "naive"]
        stim = x[x.condition_label != "naive"]
        if base.empty or stim.empty:
            continue
        bcs = load(base.iloc[0].dataset_id)
        if bcs is None:
            continue
        bh = hi(bcs)
        bkeys = set(zip(bh.gene_id, bh.variant))
        bgene = set(bh.gene_id)
        for did, cond in zip(stim.dataset_id, stim.condition_label):
            scs = load(did)
            if scs is None:
                continue
            sh = hi(scs)
            sgene = set(sh.gene_id)
            parts = [
                ("act_specific", sh[~sh.gene_id.isin(bgene)]),
                ("rest_specific", bh[~bh.gene_id.isin(sgene)]),
                ("shared", sh[[(g, v) in bkeys for g, v in
                               zip(sh.gene_id, sh.variant)]]),
            ]
            for tag, sub in parts:
                for r in sub.itertuples(index=False):
                    rows.append(dict(stim=cond, group=tag, study=study, cell=cell,
                                     variant=r.variant, gene=r.gene_id,
                                     beta=float(r.beta), z=float(r.z),
                                     pip=float(r.pip)))
            print(f"  {cond:10s} 刺激可信集 {len(sh):4d} 静息 {len(bh):4d}", flush=True)

    df = pd.DataFrame(rows)
    # 只在 (刺激, 组, 变异, 基因) 内去重
    before = len(df)
    df = df.drop_duplicates(["stim", "group", "variant", "gene"])
    v = df["variant"].str.split("_", expand=True)
    df["chrom"], df["pos"], df["ref"], df["alt"] = v[0], v[1].astype(int), v[2], v[3]
    df.to_csv("b_data/fairfax_labels_fixed.csv", index=False)
    print(f"\n重建：{before} 行 → 去重后 {len(df)} 行")
    print(pd.crosstab(df.group, df.stim).to_string())
    print("\n每种刺激的唯一变异-基因数：")
    print(df.groupby("stim").apply(lambda x: x[["variant", "gene"]].drop_duplicates().shape[0],
                                   include_groups=False).to_string())
    # 跨刺激一致性检查
    both = df.groupby(["variant", "gene"]).stim.nunique()
    print(f"\n出现在多种刺激里的变异-基因：{int((both>1).sum())}")
    x = df.groupby(["variant", "gene"]).group.nunique()
    print(f"在多种刺激下标签不同的：{int((x>1).sum())}")


if __name__ == "__main__":
    main()
