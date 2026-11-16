#!/usr/bin/env python3
"""校准曲线：各细胞类型与 CD14 的「实测表达重叠」vs 模型给的重叠。

争议点：FCGR3A+ 条件给 +0.296，是模型只分辨细胞类型，还是模型算对了
（CD16+ 单核本身就是弱干扰素预激态）？
判据：实测 log2(CD16ctrl/CD14ctrl) 与实测 IFN 响应的相关。若 ≈0.30 且模型 0.296 落在其 CI 内，
则模型正确复现了真实的生物学重叠。

全部本地计算，不调 API。所有细胞类型均排除造条件用的 256 个细胞。
"""
import glob

import numpy as np
import pandas as pd
import scipy.io
from scipy import stats

D = "b_data"


def main():
    genes = pd.read_csv(f"{D}/GSE96583_batch2.genes.tsv.gz", sep="\t", header=None,
                        names=["ensg", "symbol"])
    ann = pd.read_csv(f"{D}/GSE96583_batch2.total.tsne.df.tsv.gz", sep="\t")
    bc = pd.read_csv(f"{D}/b2.1.barcodes.tsv.gz", header=None)[0].astype(str).tolist()
    m = scipy.io.mmread(f"{D}/GSM2560248_2.1.mtx.gz").tocsr()      # ctrl 文库
    print(f"ctrl 文库 {m.shape}")

    # 造条件用掉的条形码（各类型各自的 256）——统一全部排除，保证比较干净
    used = set()
    for f in glob.glob(f"{D}/mono_*.h5ad") + glob.glob(f"{D}/mm_*.h5ad"):
        import anndata as ad
        a = ad.read_h5ad(f)
        used |= set(a.obs.index.astype(str))
    print(f"条件侧共用掉 {len(used)} 个条形码")

    sel = ann.reindex(bc)
    base_ok = (sel["multiplets"] == "singlet").values
    types = ["CD14+ Monocytes", "FCGR3A+ Monocytes", "CD4 T cells", "CD8 T cells",
             "NK cells", "B cells", "Dendritic cells", "Megakaryocytes"]
    libs = {}
    for t in types:
        mask = base_ok & (sel["cell"] == t).values & np.array([b not in used for b in bc])
        if mask.sum() < 20:   # 细胞太少且已全被条件侧占用，退回用全部
            mask = base_ok & (sel["cell"] == t).values
            print(f"  {t}: 不相交细胞不足，改用全部 {int(mask.sum())} 个（含条件侧，仅作参考）")
        if mask.sum() < 20:
            print(f"  {t}: 只剩 {mask.sum()} 个，跳过")
            continue
        sub = m[:, np.where(mask)[0]]
        cnt = np.asarray(sub.sum(1)).ravel()
        libs[t] = cnt / cnt.sum() * 1e6
        print(f"  {t}: {int(mask.sum())} 细胞，文库 {cnt.sum():,.0f}")

    de = pd.read_csv(sorted(glob.glob(f"{D}/de_disjoint_*.csv"))[-1]).set_index("ensg")
    pred = pd.read_csv(sorted(glob.glob(f"{D}/b_pred_disjoint_*.csv"))[-1])
    base = pred[(pred[["ctrl", "stim"]].max(axis=1) > 1) & np.isfinite(pred.obs)].copy()
    keep = base.ensg.tolist()
    obs = base.set_index("ensg").obs
    win = base.set_index("ensg").win

    gi = pd.Series(range(len(genes)), index=genes.ensg.values)
    idx = [int(gi[g]) for g in keep]
    c14 = libs["CD14+ Monocytes"][idx]
    obs_v = obs.values

    print(f"\n{'条件':18s} {'细胞数':>7s} {'实测重叠':>9s} {'95% CI':>18s} {'模型':>7s}")
    mod = pd.read_csv(f"{D}/null_all_disjoint.csv")
    mod_map = {"CD4 T cells": "CD4 T", "CD8 T cells": "CD8 T", "NK cells": "NK",
               "B cells": "B 细胞", "Dendritic cells": "树突细胞",
               "Megakaryocytes": "巨核细胞", "FCGR3A+ Monocytes": "FCGR3A+ 单核"}
    rng = np.random.default_rng(0)
    byw = {w: np.where(win.values == w)[0] for w in win.unique()}
    rows = []
    for t in types:
        if t not in libs or t == "CD14+ Monocytes":
            continue
        mt = libs[t][idx]
        ov = np.log2((mt + 1) / (c14 + 1))
        r = stats.spearmanr(ov, obs_v).correlation
        bs = []
        for _ in range(2000):
            p = rng.choice(list(byw), len(byw), True)
            ii = np.concatenate([byw[w] for w in p])
            bs.append(stats.spearmanr(ov[ii], obs_v[ii]).correlation)
        lo, hi = np.percentile(bs, [2.5, 97.5])
        mv = mod.loc[mod.cond == mod_map[t], "rho"]
        mv = float(mv.iloc[0]) if len(mv) else np.nan
        print(f"{t:18s} {'':>7s} {r:+9.4f}   [{lo:+.3f}, {hi:+.3f}]  {mv:+7.3f}")
        rows.append(dict(cond=t, rho_measured=r, lo=lo, hi=hi, rho_model=mv))
    rows.append(dict(cond="matched IFN-β", rho_measured=1.0, lo=np.nan, hi=np.nan,
                     rho_model=float(mod.loc[mod.tag == "real", "rho"].iloc[0])))
    out = pd.DataFrame(rows)
    out.to_csv(f"{D}/calibration_curve.csv", index=False)
    print(f"\n写出 {D}/calibration_curve.csv")

    # 最小二乘 + 去掉匹配点后的秩相关
    x, y = out.rho_measured.values, out.rho_model.values
    A = np.column_stack([np.ones(len(x)), x])
    beta, *_ = np.linalg.lstsq(A, y, rcond=None)
    r = stats.pearsonr(x, y)
    print(f"全 8 点：模型 = {beta[1]:.3f} × 实测 + {beta[0]:+.3f}，Pearson {r.statistic:.3f}（p={r.pvalue:.2g}）")
    sub = out[out.cond != "matched IFN-β"]
    s7 = stats.spearmanr(sub.rho_measured, sub.rho_model)
    print(f"去掉匹配点（7 点）：Spearman {s7.statistic:+.3f}（p={s7.pvalue:.2f}）")

    # 标记基因
    print("\n标记基因：")
    sym = {"IFIT3": "IFIT3", "OAS1": "OAS1", "IFIT1": "IFIT1", "MX1": "MX1",
           "ISG15": "ISG15", "IRF7": "IRF7", "FCGR3A": "FCGR3A", "CD14": "CD14", "LYZ": "LYZ"}
    s2g = genes.drop_duplicates("symbol").set_index("symbol").ensg
    gi2 = pd.Series(range(len(genes)), index=genes.ensg.values)
    mt = libs["FCGR3A+ Monocytes"]
    print(f"{'gene':8s} {'log2(CD16/CD14)':>16s} {'log2(IFN响应)':>14s}")
    for s_ in sym:
        if s_ not in s2g.index:
            continue
        e = s2g[s_]
        if e not in de.index:
            continue
        if e not in gi2.index:
            continue
        i = int(gi2[e])
        c14i = libs["CD14+ Monocytes"][i]
        ov = np.log2((mt[i] + 1) / (c14i + 1))
        print(f"{s_:8s} {ov:+16.2f} {de.loc[e, 'obs']:+14.2f}")


if __name__ == "__main__":
    main()
