#!/usr/bin/env python3
"""伪散装条件 pilot：把 Kang 单细胞压成 4 个类 bulk 轮廓造条件，
看它能否复现单细胞条件给出的预测（真值 = 现有 0.486 那套）。

判读（事先写死）：
  - 伪散装 pred 与单细胞 pred 的 Spearman ≥ 0.7 → bulk 造条件可行（A 路）
  - 0.4–0.7 → 部分可行，需更大窗口复核
  - < 0.4 → bulk 造条件不可行，走 B 路
同时报伪散装 pred 与实测 log2FC 的 Spearman（应与 0.486 同量级）。
"""
import glob
import os
import time
from concurrent.futures import ThreadPoolExecutor

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from alphagenome.data import genome
from alphagenome.models.v2 import dna_client as c2, dna_output as do, interval_scorers as isc
from scipy import stats

STAMP = time.strftime("%m%d_%H%M%S")
OUT = f"b_data/pseudobulk_pilot_{STAMP}.csv"
NBIN = 40
NPB = 4          # 每个条件压成 4 个伪散装轮廓
WORKERS = 6


def pseudobulk(a, k=NPB, seed=11):
    x = sp.csr_matrix(a.X)
    n = x.shape[0]
    rng = np.random.default_rng(seed)
    idx = rng.permutation(n)
    groups = np.array_split(idx, k)
    rows = []
    for g in groups:
        m = np.asarray(x[g].mean(axis=0)).ravel()
        m = m / m.sum() * 1e6 if m.sum() > 0 else m
        rows.append(m)
    out = ad.AnnData(X=sp.csr_matrix(np.vstack(rows)),
                     obs=pd.DataFrame({"grp": [f"pb{i}" for i in range(k)]}),
                     var=pd.DataFrame(index=a.var.index.astype(str)))
    return out


def main():
    a = ad.read_h5ad(sorted(glob.glob("b_data/mono_ctrl_*.h5ad"))[-1])
    b = ad.read_h5ad(sorted(glob.glob("b_data/mono_stim_*.h5ad"))[-1])
    pbc, pbs = pseudobulk(a), pseudobulk(b)
    print(f"伪散装：ctrl {pbc.shape}，stim {pbs.shape}", flush=True)

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
    print(f"窗口 {len(top)}", flush=True)

    m = c2.create(os.environ["ALPHAGENOME_API_KEY"], timeout=600)
    sc = [isc.GeneMaskScorer(requested_output=do.OutputType.SC_GENE_EXPRESSION_DENOISED,
                             width=200001, aggregation_type=isc.IntervalAggregationType.MEAN)]
    conds = {"ctrl_pb": pbc, "stim_pb": pbs}

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
            print(f"  {chrom}:{mid}/{tag} {type(e).__name__}: {e}", flush=True)
            return tag, {}

    tasks = [(k, tag) for k in top.index for tag in conds]
    res = {"ctrl_pb": {}, "stim_pb": {}}
    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        for tag, d in ex.map(one, tasks):
            res[tag].update(d)
    print(f"拿到 ctrl {len(res['ctrl_pb'])} 基因 / stim {len(res['stim_pb'])} 基因")

    genes = sorted(set(res["ctrl_pb"]) & set(res["stim_pb"]))
    pb = pd.DataFrame({"ensg": genes,
                       "ctrl_pb": [res["ctrl_pb"][g] for g in genes],
                       "stim_pb": [res["stim_pb"][g] for g in genes]})
    old = pd.read_csv(sorted(glob.glob("b_data/b_pred_*.csv"))[-1]).drop_duplicates("ensg")
    old = old[["ensg", "ctrl", "stim", "log2fc", "cpm_ctrl", "cpm_stim"]]
    old = old.rename(columns={"ctrl": "ctrl_sc", "stim": "stim_sc"})
    mrg = pb.merge(old, on="ensg", how="inner")
    mrg["pred_pb"] = np.log2((mrg.stim_pb + 1) / (mrg.ctrl_pb + 1))
    mrg["pred_sc"] = np.log2((mrg.stim_sc + 1) / (mrg.ctrl_sc + 1))
    mrg["obs"] = np.log2((mrg.cpm_stim + 1) / (mrg.cpm_ctrl + 1))
    mrg.to_csv(OUT, index=False)
    print(f"合并 {len(mrg)} 基因 → {OUT}")

    x = mrg.dropna(subset=["obs"])
    if len(x) > 10:
        r1 = stats.spearmanr(x.pred_pb, x.pred_sc)
        r2 = stats.spearmanr(x.pred_pb, x.obs)
        r3 = stats.spearmanr(x.pred_sc, x.obs)
        print(f"\n伪散装 vs 单细胞预测：rho={r1.correlation:+.3f} (n={len(x)})")
        print(f"伪散装 pred vs 实测：rho={r2.correlation:+.3f} (p={r2.pvalue:.2g})")
        print(f"单细胞 pred vs 实测：rho={r3.correlation:+.3f}")
        verdict = ("A 路可行" if r1.correlation >= 0.7 else
                   "部分可行（需复核）" if r1.correlation >= 0.4 else "B 路")
        print(f"判读：{verdict}")


if __name__ == "__main__":
    main()
