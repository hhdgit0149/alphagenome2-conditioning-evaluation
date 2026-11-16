#!/usr/bin/env python3
"""近邻对照：FCGR3A+ 单核（同为单核、未刺激）与巨核细胞，跑表达层同一流程。

pred_null = log2((近邻条件+1)/(CD14 ctrl 条件+1))，与不相交实测侧的 IFN-β 响应求相关。
FCGR3A+ 是关键一条：若它也给 ρ≈0，说明条件化分辨的是细胞类型而非刺激状态；
若给 0.3，则特异性不成立。
"""
import glob
import os
from concurrent.futures import ThreadPoolExecutor

import anndata as ad
import numpy as np
import pandas as pd
from alphagenome.data import genome
from alphagenome.models.v2 import dna_client as c2, dna_output as do, interval_scorers as isc
from scipy import stats

NBIN = 40
WORKERS = 5


def main():
    import importlib.util
    spec = importlib.util.spec_from_file_location("b18", "18_B_predict.py")
    b18 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(b18)

    de = pd.read_csv(sorted(glob.glob("b_data/de_pseudobulk_*.csv"))[-1])
    expr = de[de.expressed].copy()
    coords = b18.gene_coords_batch(expr.ensg.head(1500))
    dfc = pd.DataFrame([(e, c[0], c[1], c[2]) for e, c in coords.items()],
                       columns=["ensg", "chrom", "start", "end"])
    dfc["bin"] = (dfc.start + dfc.end) // 2 // 1048576
    top = dfc.groupby(["chrom", "bin"]).size().sort_values(ascending=False).head(NBIN)
    windows = list(top.index)
    print(f"{len(windows)} 个窗口")

    conds = {}
    for tag in ("fcgr3a", "mega"):
        f = sorted(glob.glob(f"b_data/mm_{tag}_*.h5ad"))
        conds[tag] = ad.read_h5ad(f[-1])
        print(f"{tag}: {conds[tag].shape}")

    ctrl = ad.read_h5ad(sorted(glob.glob("b_data/mono_ctrl_*.h5ad"))[-1])
    m = c2.create(os.environ["ALPHAGENOME_API_KEY"], timeout=600)
    sc = [isc.GeneMaskScorer(requested_output=do.OutputType.SC_GENE_EXPRESSION_DENOISED,
                             width=200001, aggregation_type=isc.IntervalAggregationType.MEAN)]

    def one(item):
        i, ((chrom, b), _) = item
        mid = int(b) * 1048576 + 524288
        iv = genome.Interval(chrom, mid - 524288, mid + 524288)
        d = dict(chrom=chrom, mid=mid)
        for tag, cond in list(conds.items()) + [("ctrl", ctrl)]:
            try:
                a = m.score_interval(interval=iv, interval_scorers=sc, cell_conditioning=cond,
                                     cell_condition_name=tag)[0]
                ids = [str(x).split(".")[0] for x in a.obs.gene_id]
                for g, v in zip(ids, [float(np.nanmean(x)) for x in np.asarray(a.X)]):
                    d[(tag, g)] = v
            except Exception as e:
                print(f"  {chrom}:{mid} {tag} 失败 {type(e).__name__}", flush=True)
        return d

    items = list(enumerate(top.items(), 1))
    with ThreadPoolExecutor(WORKERS) as ex:
        res = list(ex.map(one, items))
    print(f"完成 {len(res)} 个窗口")

    rows = []
    for d in res:
        genes = {k[1] for k in d if isinstance(k, tuple) and k[0] == "ctrl"}
        for g in genes:
            r = dict(ensg=g)
            for tag in ("fcgr3a", "mega", "ctrl"):
                r[tag] = d.get((tag, g), np.nan)
            rows.append(r)
    df = pd.DataFrame(rows).drop_duplicates("ensg")
    print(f"基因 {len(df)}")
    df.to_csv("b_data/null_near_1115.csv", index=False)

    # 与不相交实测侧合并
    pred = pd.read_csv(sorted(glob.glob("b_data/b_pred_disjoint_*.csv"))[-1])
    j = df.merge(pred[["ensg", "obs", "pred1", "win"]], on="ensg", how="inner")
    print(f"\n合并 {len(j)} 基因（base 119 的交集）")
    print(f"{'条件':10s} {'n':>5s} {'rho(新obs)':>11s} {'中位|pred|':>11s}")
    real = stats.spearmanr(j.pred1, j.obs).correlation
    print(f"{'真实IFN-β':10s} {len(j):5d} {real:+11.3f} {j.pred1.abs().median():11.3f}")
    for tag in ("fcgr3a", "mega"):
        j[f"pred_{tag}"] = np.log2((j[tag] + 1) / (j.ctrl + 1))
        m2 = j.dropna(subset=[f"pred_{tag}", "obs"])
        r = stats.spearmanr(m2[f"pred_{tag}"], m2.obs).correlation
        print(f"{tag:10s} {len(m2):5d} {r:+11.3f} {m2[f'pred_{tag}'].abs().median():11.3f}")


if __name__ == "__main__":
    main()
