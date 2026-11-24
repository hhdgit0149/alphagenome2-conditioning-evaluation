#!/usr/bin/env python3
"""补 5 个未刺激条件缺失的 9 个基因（落在 3 个窗口内）。

背景：52_null_tb.py 打分时有 3 个窗口对 5 个条件静默失败（异常被吞），
导致 null_tb 只有 375 个基因，base 119 里缺 9 个 → 那 5 个点只在 110 基因上算。

做法：定位这 9 个基因所在的 3 个窗口，只重打这 3 个窗口 × 5 个条件 = 15 次调用，
合并进 null_tb。
"""
import glob
import os

import anndata as ad
import numpy as np
import pandas as pd
from alphagenome.data import genome
from alphagenome.models.v2 import dna_client as c2, dna_output as do, interval_scorers as isc

D = "b_data"
CONDS = {"tcell": "CD4 T cells", "bcell": "B cells", "cd8": "CD8 T cells",
         "nk": "NK cells", "dc": "Dendritic cells"}


def main():
    pred = pd.read_csv(f"{D}/b_pred_disjoint_1115_115447.csv")
    base = pred[(pred[["ctrl", "stim"]].max(axis=1) > 1) & pred.obs.notna()]
    tb = pd.read_csv(f"{D}/null_tb.csv")
    have = set(tb.ensg)
    missing = [g for g in base.ensg if g not in have]
    print(f"base {len(base)} 基因；null_tb 已有 {len(have)}；缺 {len(missing)} 个")

    w = base.set_index("ensg").win
    wins = sorted({w[g] for g in missing})
    print(f"缺的基因落在 {len(wins)} 个窗口：{wins}")

    conds = {}
    for tag in CONDS:
        f = sorted(glob.glob(f"{D}/mm_{tag}_*.h5ad"))
        conds[tag] = ad.read_h5ad(f[-1])
        print(f"  {tag}: {conds[tag].shape}")

    m = c2.create(os.environ["ALPHAGENOME_API_KEY"], timeout=600)
    sc = [isc.GeneMaskScorer(requested_output=do.OutputType.SC_GENE_EXPRESSION_DENOISED,
                             width=200001, aggregation_type=isc.IntervalAggregationType.MEAN)]

    rows = []
    for wv in wins:
        chrom, b = wv.split(":")
        mid = int(b) * 1048576 + 524288
        iv = genome.Interval(chrom, mid - 524288, mid + 524288)
        d = {}
        for tag, cond in conds.items():
            try:
                a = m.score_interval(interval=iv, interval_scorers=sc,
                                     cell_conditioning=cond, cell_condition_name=tag)[0]
                ids = [str(x).split(".")[0] for x in a.obs.gene_id]
                for g, v in zip(ids, [float(np.nanmean(x)) for x in np.asarray(a.X)]):
                    d[(tag, g)] = v
                print(f"  {wv} {tag}: {len(ids)} 基因")
            except Exception as e:
                print(f"  ! {wv} {tag}: {type(e).__name__}: {e}")
        for g in missing:
            if g in {k[1] for k in d}:
                r = dict(ensg=g)
                for tag in conds:
                    r[tag] = d.get((tag, g), np.nan)
                rows.append(r)

    new = pd.DataFrame(rows)
    print(f"\n补到 {len(new)} 行 × {new.shape[1]-1} 条件")
    tb2 = tb.merge(new, on="ensg", how="outer", suffixes=("", "_new"))
    for tag in CONDS:
        if f"{tag}_new" in tb2:
            tb2[tag] = tb2[tag].fillna(tb2[f"{tag}_new"])
            tb2 = tb2.drop(columns=[f"{tag}_new"])
    tb2.to_csv(f"{D}/null_tb_filled.csv", index=False)
    print(f"写出 {D}/null_tb_filled.csv：{len(tb2)} 基因，"
          f"base 覆盖 {len(set(tb2.ensg) & set(base.ensg))}/{len(base)}")


if __name__ == "__main__":
    main()
