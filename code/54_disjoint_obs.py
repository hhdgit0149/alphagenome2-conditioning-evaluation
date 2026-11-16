#!/usr/bin/env python3
"""用与条件向量不相交的细胞重算实测 log2FC。

动机：17_build_B.py 的 bulk() 直接对「造条件用的那两个 AnnData」求 pseudobulk，
即标准答案与条件向量来自同一批 256×2 个细胞（完全循环）。

本脚本：条件侧保持那 256 个细胞（预测值已在盘上，不重跑），
实测侧改用**排除这 512 个细胞后**的全部 CD14+ 单核细胞（batch 2）。
预测值不变 → 零 API 调用。

输出：b_data/b_pred_disjoint_<stamp>.csv
"""
import glob
import time

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
from scipy import stats

D = "b_data"
STAMP = time.strftime("%m%d_%H%M%S")


def load_genes():
    g = pd.read_csv(f"{D}/GSE96583_batch2.genes.tsv.gz", sep="\t", header=None,
                    names=["ensg", "symbol"])
    assert len(g) == 35635
    return g


def load_ann():
    return pd.read_csv(f"{D}/GSE96583_batch2.total.tsne.df.tsv.gz", sep="\t")


def cells_used():
    """条件侧实际用到的条形码（256 ctrl + 256 stim）。"""
    used = {}
    for tag in ("ctrl", "stim"):
        f = sorted(glob.glob(f"{D}/mono_{tag}_*.h5ad"))[-1]
        a = ad.read_h5ad(f)
        used[tag] = set(a.obs.index.astype(str))
        print(f"条件侧 {tag}: {a.shape[0]} 个细胞（{f}）")
    return used


def pseudobulk(tag, fn, ann, exclude):
    bc = pd.read_csv(f"{D}/{fn}.barcodes.tsv.gz", header=None)[0].astype(str).tolist()
    m = scipy.io.mmread(f"{D}/GSM2560248_2.1.mtx.gz" if tag == "ctrl"
                        else f"{D}/GSM2560249_2.2.mtx.gz").tocsr()
    assert m.shape == (35635, len(bc))
    sel = ann.reindex(bc)
    mono = (sel["cell"] == "CD14+ Monocytes") & (sel["multiplets"] == "singlet")
    keep = np.array([mono.iloc[i] and bc[i] not in exclude for i in range(len(bc))])
    print(f"{tag}: CD14 singlet {int(mono.sum())}，排除条件侧 {sum(1 for i in range(len(bc)) if mono.iloc[i] and bc[i] in exclude)}，"
          f"剩余 {int(keep.sum())}")
    sub = m[:, np.where(keep)[0]]
    tot = np.asarray(sub.sum(0)).ravel()
    counts = np.asarray(sub.sum(1)).ravel()
    print(f"  文库总量 {tot.sum():,.0f} counts；每细胞 UMI 中位 {int(np.median(tot))}")
    return counts, float(tot.sum()), int(keep.sum())


def main():
    genes, ann = load_genes(), load_ann()
    used = cells_used()
    out = {}
    for tag, fn in (("ctrl", "b2.1"), ("stim", "b2.2")):
        out[tag] = pseudobulk(tag, fn, ann, used[tag])

    cc, tc, nc = out["ctrl"]
    cs, ts, ns = out["stim"]
    pc, ps = cc / tc * 1e6, cs / ts * 1e6
    df = pd.DataFrame({"ensg": genes.ensg.values, "cpm_ctrl": pc, "cpm_stim": ps,
                       "cnt_ctrl": cc, "cnt_stim": cs})
    df["log2fc"] = np.log2((ps + 1e-9) / (pc + 1e-9))
    df["obs"] = np.log2((ps + 1) / (pc + 1))
    df.to_csv(f"{D}/de_disjoint_{STAMP}.csv", index=False)
    print(f"\n写出 {D}/de_disjoint_{STAMP}.csv（{len(df)} 基因）")

    # 与预测合并（预测不变）
    p = pd.read_csv(sorted(glob.glob(f"{D}/b_pred_*.csv"))[-1]).drop_duplicates("ensg")
    p["pred1"] = np.log2((p.stim + 1) / (p.ctrl + 1))
    win = pd.read_csv(f"{D}/b_recomputed.csv")[["ensg", "win"]]
    p = p.merge(win, on="ensg", how="left")
    j = p.merge(df[["ensg", "cpm_ctrl", "cpm_stim", "cnt_ctrl", "cnt_stim", "obs", "log2fc"]],
                on="ensg", how="inner", suffixes=("_old", ""))
    print(f"与预测合并 {len(j)} 个基因")
    j.to_csv(f"{D}/b_pred_disjoint_{STAMP}.csv", index=False)

    # 主口径：base = 预测 max>1（沿用原筛选）
    base = j[(j[["ctrl", "stim"]].max(axis=1) > 1) & np.isfinite(j.obs)].copy()
    print(f"base（max(预测)>1 且 obs 有限）: {len(base)}")

    def rep(d, lab):
        if len(d) < 8:
            print(f"{lab:26s} n={len(d)} 太少")
            return None
        r = stats.spearmanr(d.pred1, d.obs).correlation
        ro = stats.spearmanr(d.pred1, d.log2fc).correlation if (d.log2fc.abs() < 1e8).all() else np.nan
        print(f"{lab:26s} n={len(d):4d}  rho(+1)={r:+.4f}  rho(无伪计数)={ro:+.4f}")
        return r

    print("\n== 不相交细胞上的结果 ==")
    rep(base, "全部 base")
    for thr in (1, 5, 6, 20):
        rep(base[(base.cnt_ctrl >= thr) & (base.cnt_stim >= thr)], f"两条件都 ≥{thr} counts")

    print("\n== 窗口聚类自助（base，2000 次）==")
    rng = np.random.default_rng(0)
    byw = {w: g for w, g in base.groupby("win")} if "win" in base.columns else None
    if byw:
        wl = list(byw)
        bs = []
        for _ in range(2000):
            s = pd.concat([byw[w] for w in rng.choice(wl, len(wl), True)])
            if s.pred1.nunique() > 3:
                bs.append(stats.spearmanr(s.pred1, s.obs).correlation)
        print(f"窗口数 {len(wl)}，SE {np.std(bs, ddof=1):.4f}，95% [{np.percentile(bs,2.5):.3f}, {np.percentile(bs,97.5):.3f}]")
    else:
        print("b_pred 无 win 列，跳过")

    # 与旧（循环）版本对照
    old = pd.read_csv(f"{D}/b_recomputed.csv")
    m = old[["ensg", "obs", "pred"]].merge(j[["ensg", "obs", "cnt_ctrl", "cnt_stim"]],
                                           on="ensg", how="inner", suffixes=("_old", "_new"))
    print(f"\n== 新旧实测侧对照（{len(m)} 基因）==")
    print(f"rho(旧obs, 新obs) = {stats.spearmanr(m.obs_old, m.obs_new).correlation:.4f}")
    print(f"旧侧 rho(pred,obs) = {stats.spearmanr(m.pred, m.obs_old).correlation:+.4f}")
    print(f"新侧 rho(pred,obs) = {stats.spearmanr(m.pred, m.obs_new).correlation:+.4f}")


if __name__ == "__main__":
    main()
