#!/usr/bin/env python3
"""全基因留出：把 119 个目标基因从条件输入里全部剔除，重建条件，重打分。

问题：条件向量由表达数据构造，标准答案也是表达。若模型能从条件里直接读回目标基因，
      0.616 就含回流成分。此前"基因留出"只移除了输入的 1.2%，证明不了。

设计：
  - 条件 = 256 细胞 × 35635 基因，但 119 个目标基因的列**整列删除**（不是置零，避免
    "未表达"语义污染）
  - 窗口与打分器与主分析完全一致（40 个 1Mb 窗口、200kb 均值掩码、DENOISED）
  - ρ 与不相交实测侧比较；保持 → 从序列/相关基因推断；塌掉 → 输入回流

用法：python 58_gene_holdout.py --test   # 单窗口冒烟
      python 58_gene_holdout.py          # 全量 80 次调用
"""
import glob
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from alphagenome.data import genome
from alphagenome.models.v2 import dna_client as c2, dna_output as do, interval_scorers as isc
from scipy import stats

NBIN = 40
WORKERS = 5
STAMP = time.strftime("%m%d_%H%M%S")


def windows(b18):
    de = pd.read_csv(sorted(glob.glob("b_data/de_pseudobulk_*.csv"))[-1])
    expr = de[de.expressed].copy()
    coords = b18.gene_coords_batch(expr.ensg.head(1500))
    dfc = pd.DataFrame([(e, c[0], c[1], c[2]) for e, c in coords.items()],
                       columns=["ensg", "chrom", "start", "end"])
    dfc["bin"] = (dfc.start + dfc.end) // 2 // 1048576
    top = dfc.groupby(["chrom", "bin"]).size().sort_values(ascending=False).head(NBIN)
    return list(top.items())


def drop_genes(a, targets):
    keep = ~a.var.index.astype(str).isin(targets)
    b = ad.AnnData(X=a.X[:, np.where(keep)[0]], obs=a.obs.copy(),
                   var=a.var.loc[keep].copy())
    return b


def main():
    test = "--test" in sys.argv
    import importlib.util
    spec = importlib.util.spec_from_file_location("b18", "18_B_predict.py")
    b18 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(b18)

    pred = pd.read_csv(sorted(glob.glob("b_data/b_pred_disjoint_*.csv"))[-1])
    base = pred[(pred[["ctrl", "stim"]].max(axis=1) > 1) & np.isfinite(pred.obs)]
    targets = set(base.ensg)
    print(f"目标基因 {len(targets)}")

    ctrl = ad.read_h5ad(sorted(glob.glob("b_data/mono_ctrl_*.h5ad"))[-1])
    stim = ad.read_h5ad(sorted(glob.glob("b_data/mono_stim_*.h5ad"))[-1])
    hit_c = int(pd.Series(ctrl.var.index.astype(str)).isin(targets).sum())
    hit_s = int(pd.Series(stim.var.index.astype(str)).isin(targets).sum())
    print(f"条件矩阵里命中：ctrl {hit_c}/{ctrl.shape[1]}，stim {hit_s}/{stim.shape[1]}")

    c2d, s2d = drop_genes(ctrl, targets), drop_genes(stim, targets)
    print(f"删除后：ctrl {c2d.shape}，stim {s2d.shape}")

    wins = windows(b18)
    if test:
        wins = wins[:1]
    print(f"窗口 {len(wins)}")

    m = c2.create(os.environ["ALPHAGENOME_API_KEY"], timeout=600)
    sc = [isc.GeneMaskScorer(requested_output=do.OutputType.SC_GENE_EXPRESSION_DENOISED,
                             width=200001, aggregation_type=isc.IntervalAggregationType.MEAN)]

    def one(item):
        i, ((chrom, b), _) = item
        mid = int(b) * 1048576 + 524288
        iv = genome.Interval(chrom, mid - 524288, mid + 524288)
        out = dict(chrom=chrom, mid=mid)
        for tag, cond in (("ctrl", c2d), ("stim", s2d)):
            try:
                a = m.score_interval(interval=iv, interval_scorers=sc,
                                     cell_conditioning=cond, cell_condition_name=tag)[0]
                ids = [str(x).split(".")[0] for x in a.obs.gene_id]
                for g, v in zip(ids, [float(np.nanmean(x)) for x in np.asarray(a.X)]):
                    out[(tag, g)] = v
            except Exception as e:
                out[f"err_{tag}"] = f"{type(e).__name__}: {e}"
        return out

    items = list(enumerate(wins, 1))
    with ThreadPoolExecutor(WORKERS) as ex:
        res = list(ex.map(one, items))

    errs = [r for r in res if any(isinstance(k, str) and k.startswith("err_") for k in r)]
    for e in errs[:3]:
        print(" 错误:", {k: v for k, v in e.items() if isinstance(k, str) and k.startswith("err_")})
    if errs:
        print(f"{len(errs)} 个窗口报错")

    rows = []
    for r in res:
        genes = {k[1] for k in r if isinstance(k, tuple) and k[0] == "ctrl"}
        for g in genes:
            rows.append(dict(ensg=g, ctrl=r.get(("ctrl", g), np.nan),
                             stim=r.get(("stim", g), np.nan)))
    df = pd.DataFrame(rows).drop_duplicates("ensg")
    print(f"基因 {len(df)}（其中目标基因 {int(df.ensg.isin(targets).sum())}）")

    if test:
        print("冒烟通过：接口接受删列后的 AnnData")
        return

    df.to_csv(f"b_data/holdout_pred_{STAMP}.csv", index=False)
    j = df.merge(base[["ensg", "obs", "pred1"]], on="ensg", how="inner")
    j["pred_ho"] = np.log2((j.stim + 1) / (j.ctrl + 1))
    mm = j.dropna(subset=["pred_ho", "obs"])
    r_ho = stats.spearmanr(mm.pred_ho, mm.obs).correlation
    r_full = stats.spearmanr(mm.pred1, mm.obs).correlation
    r_two = stats.spearmanr(mm.pred_ho, mm.pred1).correlation
    print(f"\n== 全基因留出 ==")
    print(f"n={len(mm)}  留出后 rho(pred, obs) = {r_ho:+.4f}")
    print(f"          原条件 rho(pred, obs) = {r_full:+.4f}")
    print(f"          留出后与原预测的相关 = {r_two:+.4f}")
    print(f"写出 b_data/holdout_pred_{STAMP}.csv")


if __name__ == "__main__":
    main()
