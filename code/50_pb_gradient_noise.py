#!/usr/bin/env python3
"""伪散装梯度 + 噪声基线多次拆分（同一套窗口、同一打分器）。

梯度：把单细胞压成 k=4/16/64 个伪散装轮廓，看 pred 与单细胞 pred 的一致性是否随 k 变化。
噪声：反复把对照细胞随机二分，各自造条件后跑同一流程，统计 noise_pred 与实测的
      Spearman 分布（报均值±SD），替代原先的"单次拆分"。
输出逐步落盘，可中断续跑。
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
NBIN = 40
WORKERS = 10
LEVELS = [4, 16, 64]
N_NOISE = 50         # 噪声基线的随机拆分次数
OUT_GRAD = f"b_data/pb_gradient_{STAMP}.csv"
OUT_NOISE = f"b_data/noise_splits_{STAMP}.csv"


def pseudobulk(a, k, seed):
    x = sp.csr_matrix(a.X)
    rng = np.random.default_rng(seed)
    idx = rng.permutation(x.shape[0])
    rows = []
    for g in np.array_split(idx, k):
        m = np.asarray(x[g].mean(axis=0)).ravel()
        m = m / m.sum() * 1e6 if m.sum() > 0 else m
        rows.append(m)
    return ad.AnnData(X=sp.csr_matrix(np.vstack(rows)),
                      obs=pd.DataFrame({"grp": [f"pb{i}" for i in range(k)]}),
                      var=pd.DataFrame(index=a.var.index.astype(str)))


def subset(a, idx):
    sub = a[idx].copy()
    return sub


def main():
    ctrl = ad.read_h5ad(sorted(glob.glob("b_data/mono_ctrl_*.h5ad"))[-1])
    stim = ad.read_h5ad(sorted(glob.glob("b_data/mono_stim_*.h5ad"))[-1])
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

    def score_cond(cond, tag, windows):
        def one(w):
            chrom, b_ = w
            mid = int(b_) * 1048576 + 524288
            iv = genome.Interval(chrom, mid - 524288, mid + 524288)
            try:
                arr = m.score_interval(interval=iv, interval_scorers=sc,
                                       cell_conditioning=cond, cell_condition_name=tag)[0]
                ids = [str(x).split(".")[0] for x in arr.obs.gene_id]
                return dict(zip(ids, [float(np.nanmean(v)) for v in np.asarray(arr.X)]))
            except Exception as e:
                print(f"  {chrom}:{mid}/{tag} {type(e).__name__}", flush=True)
                return {}
        out = {}
        with ThreadPoolExecutor(max_workers=WORKERS) as ex:
            for d in ex.map(one, windows):
                out.update(d)
        return out

    windows = list(top.index)

    # ---------- 梯度 ----------
    rows = []
    for k in LEVELS:
        pc = pseudobulk(ctrl, k, 11)
        ps = pseudobulk(stim, k, 23)
        cvals = score_cond(pc, f"k{k}_ctrl", windows)
        svals = score_cond(ps, f"k{k}_stim", windows)
        genes = sorted(set(cvals) & set(svals))
        for g in genes:
            rows.append(dict(level=k, ensg=g, ctrl_pb=cvals[g], stim_pb=svals[g]))
        pd.DataFrame(rows).to_csv(OUT_GRAD, index=False)
        print(f"水平 k={k}：{len(genes)} 基因", flush=True)

    # ---------- 噪声：多次随机二分 ----------
    nrows = []
    done = set()
    for fp in glob.glob("b_data/noise_splits_*.csv"):
        if os.path.getsize(fp) > 200:
            done |= set(pd.read_csv(fp, usecols=["split"])["split"].unique().tolist())
    print(f"已有拆分 {sorted(done)}", flush=True)
    for s in range(N_NOISE):
        if s in done:
            continue
        rng = np.random.default_rng(1000 + s)
        idx = rng.permutation(ctrl.shape[0])
        h1, h2 = subset(ctrl, idx[:len(idx) // 2]), subset(ctrl, idx[len(idx) // 2:])
        v1 = score_cond(h1, f"s{s}_h1", windows)
        v2 = score_cond(h2, f"s{s}_h2", windows)
        genes = sorted(set(v1) & set(v2))
        d = pd.DataFrame({"ensg": genes,
                          "h1": [v1[g] for g in genes],
                          "h2": [v2[g] for g in genes]})
        d["noise"] = np.log2((d.h1 + 1) / (d.h2 + 1))
        nrows.extend([dict(split=s, ensg=r.ensg, noise=r.noise) for r in d.itertuples(index=False)])
        pd.DataFrame(nrows).to_csv(OUT_NOISE, index=False)
        print(f"噪声拆分 {s+1}/{N_NOISE}：{len(genes)} 基因", flush=True)

    # ---------- 汇总 ----------
    old = pd.read_csv(sorted(glob.glob("b_data/b_pred_*.csv"))[-1]).drop_duplicates("ensg")
    old["obs"] = np.log2((old.cpm_stim + 1) / (old.cpm_ctrl + 1))
    old["pred_sc"] = np.log2((old.stim + 1) / (old.ctrl + 1))
    base = old[((old.ctrl > 1) | (old.stim > 1)) & old.obs.notna()]

    grad = pd.read_csv(OUT_GRAD)
    grad["pred_pb"] = np.log2((grad.stim_pb + 1) / (grad.ctrl_pb + 1))
    print("\n== 梯度（与单细胞 pred 的一致性、与实测的相关）==")
    for k, g in grad.groupby("level"):
        x = g.merge(base[["ensg", "pred_sc", "obs"]], on="ensg", how="inner")
        if len(x) < 20:
            continue
        print(f"k={k:3d} n={len(x):3d}  rho(pb,sc)={stats.spearmanr(x.pred_pb, x.pred_sc).correlation:+.3f}  "
              f"rho(pb,obs)={stats.spearmanr(x.pred_pb, x.obs).correlation:+.3f}")

    nz = pd.read_csv(OUT_NOISE)
    print("\n== 噪声基线（多次拆分）==")
    for s, g in nz.groupby("split"):
        x = g.merge(base[["ensg", "obs"]], on="ensg", how="inner")
        if len(x) < 20:
            continue
        r = stats.spearmanr(x.noise, x.obs).correlation
        print(f"split {s}: n={len(x)} rho={r:+.4f}")
    rs = []
    for s, g in nz.groupby("split"):
        x = g.merge(base[["ensg", "obs"]], on="ensg", how="inner")
        if len(x) >= 20:
            rs.append(stats.spearmanr(x.noise, x.obs).correlation)
    print(f"汇总：均值 {np.mean(rs):+.4f}，SD {np.std(rs, ddof=1):.4f}，范围 [{min(rs):+.4f}, {max(rs):+.4f}]（{len(rs)} 次）")


if __name__ == "__main__":
    main()
