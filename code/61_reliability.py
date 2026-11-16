#!/usr/bin/env python3
"""校准斜率的衰减校正：估计 x 轴（实测重叠）的可靠性系数。

问题：12 点回归里 x 轴本身是估计值（细胞抽样），带测量误差会把 OLS 斜率向零压，
      所以 0.71 是保真度的下界。
做法：把每个细胞类型的细胞随机二分，各自算一份实测重叠，得到两条平行的 x 向量；
      两者跨条件的相关即奇半可靠度 r_half，Spearman-Brown 校正得全样本可靠度
      r_full = 2 r_half / (1 + r_half)。
      衰减校正后的斜率 = 观测斜率 / r_full。
重复 20 次取均值。全部本地，零 API。
"""
import glob

import numpy as np
import pandas as pd
import scipy.io
from scipy import stats

D = "b_data"
UNSTIM = ["CD4 T cells", "CD8 T cells", "NK cells", "B cells", "Dendritic cells",
          "Megakaryocytes", "FCGR3A+ Monocytes"]
STIM = ["CD4 T cells", "CD8 T cells", "NK cells", "B cells", "Dendritic cells",
        "Megakaryocytes", "FCGR3A+ Monocytes"]


def main():
    genes = pd.read_csv(f"{D}/GSE96583_batch2.genes.tsv.gz", sep="\t", header=None,
                        names=["ensg", "symbol"])
    ann = pd.read_csv(f"{D}/GSE96583_batch2.total.tsne.df.tsv.gz", sep="\t")
    libs = {}
    for tag, bcfile, mtx in [("ctrl", "b2.1", "GSM2560248_2.1"), ("stim", "b2.2", "GSM2560249_2.2")]:
        bc = pd.read_csv(f"{D}/{bcfile}.barcodes.tsv.gz", header=None)[0].astype(str).tolist()
        s = ann.reindex(bc)
        m = scipy.io.mmread(f"{D}/{mtx}.mtx.gz").tocsr()
        libs[tag] = dict(bc=bc, sel=s, m=m,
                         ok=(s["multiplets"] == "singlet").values,
                         cell=s["cell"].values)

    pred = pd.read_csv(sorted(glob.glob(f"{D}/b_pred_disjoint_*.csv"))[-1])
    base = pred[(pred[["ctrl", "stim"]].max(axis=1) > 1) & np.isfinite(pred.obs)]
    keep = base.ensg.tolist()
    obs = base.set_index("ensg").obs.values
    gi = pd.Series(range(len(genes)), index=genes.ensg.values)
    idx = [int(gi[g]) for g in keep]

    def counts(tag, celltype, cells):
        L = libs[tag]
        sub = L["m"][:, cells]
        return np.asarray(sub.sum(1)).ravel()[idx]

    # 基线 CD14（未刺激）细胞池
    c14_pool = np.where(libs["ctrl"]["ok"] & (libs["ctrl"]["cell"] == "CD14+ Monocytes"))[0]

    rng = np.random.default_rng(11)
    rel_half, rel_full, slopes = [], [], []
    vec = []
    for rep in range(20):
        p = rng.permutation(c14_pool)
        c14_a, c14_b = p[: len(p) // 2], p[len(p) // 2:]
        cnt_a, cnt_b = counts("ctrl", "CD14", c14_a), counts("ctrl", "CD14", c14_b)
        # 基线总量用两半各自的 sum（等价于两个独立的 baseline 估计）
        tot_a = libs["ctrl"]["m"][:, c14_a].sum()
        tot_b = libs["ctrl"]["m"][:, c14_b].sum()
        cpm_a = cnt_a / tot_a * 1e6
        cpm_b = cnt_b / tot_b * 1e6

        xa, xb = [], []
        for t in UNSTIM:
            pool = np.where(libs["ctrl"]["ok"] & (libs["ctrl"]["cell"] == t))[0]
            q = rng.permutation(pool)
            ta, tb = q[: len(q) // 2], q[len(q) // 2:]
            ca, cb = counts("ctrl", t, ta), counts("ctrl", t, tb)
            va = np.log2((ca / libs["ctrl"]["m"][:, ta].sum() * 1e6 + 1) / (cpm_a + 1))
            vb = np.log2((cb / libs["ctrl"]["m"][:, tb].sum() * 1e6 + 1) / (cpm_b + 1))
            xa.append(stats.spearmanr(va, obs).correlation)
            xb.append(stats.spearmanr(vb, obs).correlation)
        for t in STIM:
            pool = np.where(libs["stim"]["ok"] & (libs["stim"]["cell"] == t))[0]
            q = rng.permutation(pool)
            ta, tb = q[: len(q) // 2], q[len(q) // 2:]
            ca, cb = counts("stim", t, ta), counts("stim", t, tb)
            va = np.log2((ca / libs["stim"]["m"][:, ta].sum() * 1e6 + 1) / (cpm_a + 1))
            vb = np.log2((cb / libs["stim"]["m"][:, tb].sum() * 1e6 + 1) / (cpm_b + 1))
            xa.append(stats.spearmanr(va, obs).correlation)
            xb.append(stats.spearmanr(vb, obs).correlation)

        r = stats.pearsonr(xa, xb).statistic
        vec.append(dict(rep=rep, **{f"A_{i}": v for i, v in enumerate(xa)},
                        **{f"B_{i}": v for i, v in enumerate(xb)}))
        rel_half.append(r)
        rel_full.append(2 * r / (1 + r))

    rel_half = np.array(rel_half)
    rel_full = np.array(rel_full)
    print(f"奇半可靠度 r_half = {rel_half.mean():.3f}（SD {rel_half.std(ddof=1):.3f}）")
    print(f"Spearman-Brown 全样本可靠度 r_full = {rel_full.mean():.3f}")
    print(f"\n观测斜率 0.710 → 校正后 = {0.710 / rel_full.mean():.3f}")
    print(f"观测斜率 0.801（去匹配点）→ 校正后 = {0.801 / rel_full.mean():.3f}")
    pd.DataFrame(dict(r_half=rel_half, r_full=rel_full)).to_csv(f"{D}/reliability.csv", index=False)
    pd.DataFrame(vec).to_csv(f"{D}/reliability_vectors.csv", index=False)
    print(f"写出 {D}/reliability_vectors.csv（20 次 × 14 条件 × 2 半）")
    print(f"\n写出 {D}/reliability.csv")


if __name__ == "__main__":
    main()
