#!/usr/bin/env python3
"""量级匹配零分布：用 T 细胞 / B 细胞条件跑表达层同一流程。

拆半零分布的预测 log2FC 中位只有 0.0011（真实 0.0775），尺度不匹配。
这里用真实细胞类型条件（预测有真实量级）作零分布：
pred_null = log2((T+1)/(mono_ctrl+1))，与 IFN-β 实测变化的相关应为零。
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
OUT = "b_data/null_tb.csv"


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

    conds = {}
    for tag, pat in [("tcell", "b_data/mm_tcell_*.h5ad"), ("bcell", "b_data/mm_bcell_*.h5ad"),
                     ("cd8", "b_data/mm_cd8_*.h5ad"), ("nk", "b_data/mm_nk_*.h5ad"),
                     ("dc", "b_data/mm_dc_*.h5ad")]:
        f = sorted(glob.glob(pat))
        if f:
            conds[tag] = ad.read_h5ad(f[-1])
    print("零分布条件:", list(conds))
    m = c2.create(os.environ["ALPHAGENOME_API_KEY"], timeout=600)
    sc = [isc.GeneMaskScorer(requested_output=do.OutputType.SC_GENE_EXPRESSION_DENOISED,
                             width=200001, aggregation_type=isc.IntervalAggregationType.MEAN)]

    def one(t):
        (chrom, b_), tag = t
        mid = int(b_) * 1048576 + 524288
        iv = genome.Interval(chrom, mid - 524288, mid + 524288)
        try:
            arr = m.score_interval(interval=iv, interval_scorers=sc,
                                   cell_conditioning=conds[tag],
                                   cell_condition_name=tag)[0]
            ids = [str(x).split(".")[0] for x in arr.obs.gene_id]
            return tag, dict(zip(ids, [float(np.nanmean(v)) for v in np.asarray(arr.X)]))
        except Exception as e:
            print(f"  {chrom}:{mid}/{tag} {type(e).__name__}", flush=True)
            return tag, {}

    res = {k: {} for k in conds}
    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        for tag, d in ex.map(one, [(w, t) for w in windows for t in conds]):
            res[tag].update(d)
    print(f"tcell {len(res['tcell'])} 基因，bcell {len(res['bcell'])} 基因")

    old = pd.read_csv(sorted(glob.glob("b_data/b_pred_*.csv"))[-1]).drop_duplicates("ensg")
    old["obs"] = np.log2((old.cpm_stim + 1) / (old.cpm_ctrl + 1))
    common = set.intersection(*[set(v) for v in res.values()]) if res else set()
    df = pd.DataFrame({"ensg": sorted(common)})
    for k, v in res.items():
        df[k] = df.ensg.map(v)
    df = df.merge(old[["ensg", "ctrl", "stim", "obs"]], on="ensg", how="left")
    for k in res:
        df[f"pred_null_{k}"] = np.log2((df[k] + 1) / (df.ctrl + 1))
    df["pred_real"] = np.log2((df.stim + 1) / (df.ctrl + 1))
    df.to_csv(OUT, index=False)

    x = df[df.obs.notna() & ((df.ctrl > 1) | (df.stim > 1))]
    print(f"\n可用基因 {len(x)}")
    cols = [c for c in df.columns if c.startswith("pred_null_")] + ["pred_real"]
    for c in cols:
        r = stats.spearmanr(x[c], x.obs)
        print(f"{c:12s} 中位|pred|={x[c].abs().median():.4f}  "
              f"rho(obs)={r.correlation:+.3f}  p={r.pvalue:.2g}")
    print(f"\n真实 pred 中位|log2FC|={x.pred_real.abs().median():.4f}；"
          f"拆半零分布中位 0.0011")
    print(f"输出 {OUT}")


if __name__ == "__main__":
    main()
