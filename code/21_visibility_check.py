#!/usr/bin/env python3
"""第三步：用 B 已建好的两个单核细胞条件，实测 Fairfax 位点的可见比例。

只看「模型是否认为这个位点有效应」（|效应|>0.1），不看组间差异。
这一步决定位点层面的主测试有没有可能做。
"""
import glob
import os
import time

import anndata as ad
import numpy as np
import pandas as pd
from alphagenome.data import genome
from alphagenome.models.v2 import dna_client as c2, variant_scorers as vs2

N = 100
STAMP = time.strftime("%m%d_%H%M%S")


def main():
    g = pd.read_csv("stim_variants_groups.csv")
    f = g[(g.study == "Fairfax_2014")].drop_duplicates(["group", "variant"])
    act = f[f.group == "act_specific"].sample(min(N, (f.group == "act_specific").sum()),
                                              random_state=7)
    rest = f[f.group == "rest_specific"].sample(min(N, (f.group == "rest_specific").sum()),
                                                random_state=7)
    tgt = pd.concat([act, rest])
    print(f"抽样：激活特异 {len(act)}，静息特异 {len(rest)}")

    ac = ad.read_h5ad(sorted(glob.glob("b_data/mono_ctrl_*.h5ad"))[-1])
    as_ = ad.read_h5ad(sorted(glob.glob("b_data/mono_stim_*.h5ad"))[-1])
    m = c2.create(os.environ["ALPHAGENOME_API_KEY"], timeout=600)
    sc = [vs2.RECOMMENDED_VARIANT_SCORERS["RNA_SEQ"]]

    rows = []
    for i, r in enumerate(tgt.itertuples(index=False), 1):
        try:
            v = genome.Variant(r.chrom, int(r.pos), r.ref, r.alt)
            iv = v.reference_interval.resize(c2.SEQUENCE_LENGTH_1MB)
        except Exception:
            continue
        rec = dict(variant=r.variant, gene=r.gene, group=r.group)
        for tag, cond in [("ctrl", ac), ("stim", as_)]:
            try:
                d = vs2.tidy_scores(m.score_variant(
                    interval=iv, variant=v, variant_scorers=sc,
                    cell_conditioning=cond, cell_condition_name=tag))
                h = d[d.gene_id == r.gene]
                rec[tag] = float(h.raw_score.mean()) if len(h) else np.nan
            except Exception as e:
                print(f"  {r.variant}/{tag} {type(e).__name__}", flush=True)
                rec[tag] = np.nan
        rows.append(rec)
        if i % 20 == 0:
            print(f"  {i}/{len(tgt)}", flush=True)

    p = pd.DataFrame(rows)
    p["vis"] = p[["ctrl", "stim"]].abs().max(axis=1) > 0.1
    p.to_csv(f"b_data/visibility_fairfax_{STAMP}.csv", index=False)
    print(f"\n打分 {len(p)} 个位点")
    print(p.groupby("group")["vis"].agg(可见率="mean", 可见="sum", 总数="count").round(3).to_string())
    v = p.vis.mean()
    print(f"\n总体可见率 {v:.1%}（T 细胞条件下的估计是 11%）")
    print(f"若按此比例，Fairfax 的 419+403 个特异位点约有 {int((419+403)*v)} 个可用")
    print("判据：<50 则位点层面主测试不做")


if __name__ == "__main__":
    main()
