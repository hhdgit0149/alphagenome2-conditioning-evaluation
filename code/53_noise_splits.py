#!/usr/bin/env python3
"""只做拆半噪声基线（跳过梯度，续跑，并发 10）。目标 50 次拆分。"""
import glob
import os
import time
from concurrent.futures import ThreadPoolExecutor

import anndata as ad
import numpy as np
import pandas as pd
from alphagenome.data import genome
from alphagenome.models.v2 import dna_client as c2, dna_output as do, interval_scorers as isc
from scipy import stats

NBIN = 40
WORKERS = 10
N_NOISE = 50
OUT = "b_data/noise_splits_50.csv"


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

    ctrl = ad.read_h5ad(sorted(glob.glob("b_data/mono_ctrl_*.h5ad"))[-1])
    m = c2.create(os.environ["ALPHAGENOME_API_KEY"], timeout=600)
    sc = [isc.GeneMaskScorer(requested_output=do.OutputType.SC_GENE_EXPRESSION_DENOISED,
                             width=200001, aggregation_type=isc.IntervalAggregationType.MEAN)]

    def score_cond(cond, tag):
        def one(w):
            chrom, b_ = w
            mid = int(b_) * 1048576 + 524288
            iv = genome.Interval(chrom, mid - 524288, mid + 524288)
            try:
                arr = m.score_interval(interval=iv, interval_scorers=sc,
                                       cell_conditioning=cond, cell_condition_name=tag)[0]
                ids = [str(x).split(".")[0] for x in arr.obs.gene_id]
                return dict(zip(ids, [float(np.nanmean(v)) for v in np.asarray(arr.X)]))
            except Exception:
                return {}
        out = {}
        with ThreadPoolExecutor(max_workers=WORKERS) as ex:
            for d in ex.map(one, windows):
                out.update(d)
        return out

    done = set()
    for fp in glob.glob("b_data/noise_splits_*.csv"):
        if os.path.getsize(fp) > 200:
            done |= set(pd.read_csv(fp, usecols=["split"])["split"].unique().tolist())
    print(f"已完成拆分 {sorted(done)}", flush=True)

    rows = []
    if os.path.exists(OUT):
        rows = pd.read_csv(OUT).to_dict("records")
    for s in range(N_NOISE):
        if s in done:
            continue
        rng = np.random.default_rng(1000 + s)
        idx = rng.permutation(ctrl.shape[0])
        h1, h2 = ctrl[idx[:len(idx) // 2]].copy(), ctrl[idx[len(idx) // 2:]].copy()
        v1 = score_cond(h1, f"s{s}_h1")
        v2 = score_cond(h2, f"s{s}_h2")
        genes = sorted(set(v1) & set(v2))
        for g in genes:
            rows.append(dict(split=s, ensg=g,
                             noise=float(np.log2((v1[g] + 1) / (v2[g] + 1)))))
        pd.DataFrame(rows).to_csv(OUT, index=False)
        print(f"拆分 {s+1}/{N_NOISE}：{len(genes)} 基因（累计 {len(rows)} 行）", flush=True)

    nz = pd.read_csv(OUT)
    old = pd.read_csv(sorted(glob.glob("b_data/b_pred_*.csv"))[-1]).drop_duplicates("ensg")
    old["obs"] = np.log2((old.cpm_stim + 1) / (old.cpm_ctrl + 1))
    base = old[((old.ctrl > 1) | (old.stim > 1)) & old.obs.notna()]
    rs = []
    for s, g in nz.groupby("split"):
        x = g.merge(base[["ensg", "obs"]], on="ensg", how="inner")
        if len(x) >= 20:
            rs.append(stats.spearmanr(x.noise, x.obs).correlation)
    print(f"\n{len(rs)} 次拆分：均值 {np.mean(rs):+.4f}，SD {np.std(rs, ddof=1):.4f}，"
          f"范围 [{min(rs):+.3f}, {max(rs):+.3f}]")
    print(f"经验单侧 p = 1/(1+{len(rs)}) = {1/(1+len(rs)):.4f}")


if __name__ == "__main__":
    main()
