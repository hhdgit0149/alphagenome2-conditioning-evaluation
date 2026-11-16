#!/usr/bin/env python3
"""校准曲线中段：受刺激细胞类型与 CD14 IFN 响应的「实测重叠」。

先只算实测（本地，零 API）：若这些点的实测重叠落在 0.4–0.8，才值得花钱让模型去算。
实测重叠 = log2((该类型刺激后 CPM + 1) / (CD14 未刺激 CPM + 1))，与 CD14 的实测 IFN 响应求相关。
"""
import glob

import numpy as np
import pandas as pd
import scipy.io
from scipy import stats

D = "b_data"


def lib(fn, mask):
    m = scipy.io.mmread(fn).tocsr()
    sub = m[:, np.where(mask)[0]]
    cnt = np.asarray(sub.sum(1)).ravel()
    return cnt / cnt.sum() * 1e6, int(mask.sum())


def main():
    genes = pd.read_csv(f"{D}/GSE96583_batch2.genes.tsv.gz", sep="\t", header=None,
                        names=["ensg", "symbol"])
    ann = pd.read_csv(f"{D}/GSE96583_batch2.total.tsne.df.tsv.gz", sep="\t")
    b1 = pd.read_csv(f"{D}/b2.1.barcodes.tsv.gz", header=None)[0].astype(str).tolist()
    b2 = pd.read_csv(f"{D}/b2.2.barcodes.tsv.gz", header=None)[0].astype(str).tolist()
    s1, s2 = ann.reindex(b1), ann.reindex(b2)
    ok1 = (s1["multiplets"] == "singlet").values
    ok2 = (s2["multiplets"] == "singlet").values
    print(f"ctrl 文库 {len(b1)} 细胞，stim 文库 {len(b2)} 细胞")

    # CD14 未刺激基线（与不相交实测侧一致：用 b2.1 里全部 CD14 singlet）
    c14, n14 = lib(f"{D}/GSM2560248_2.1.mtx.gz", ok1 & (s1["cell"] == "CD14+ Monocytes").values)
    print(f"CD14 ctrl 基线 {n14} 细胞")

    de = pd.read_csv(sorted(glob.glob(f"{D}/de_disjoint_*.csv"))[-1]).set_index("ensg")
    pred = pd.read_csv(sorted(glob.glob(f"{D}/b_pred_disjoint_*.csv"))[-1])
    base = pred[(pred[["ctrl", "stim"]].max(axis=1) > 1) & np.isfinite(pred.obs)]
    keep = base.ensg.tolist()
    obs = base.set_index("ensg").obs
    win = base.set_index("ensg").win
    gi = pd.Series(range(len(genes)), index=genes.ensg.values)
    idx = [int(gi[g]) for g in keep]

    byw = {w: np.where(win.values == w)[0] for w in win.unique()}
    rng = np.random.default_rng(0)

    print(f"\n{'刺激后类型':18s} {'细胞':>6s} {'实测重叠':>9s} {'95% CI':>18s}")
    rows = []
    for t in ["CD14+ Monocytes", "CD4 T cells", "CD8 T cells", "NK cells", "B cells",
              "FCGR3A+ Monocytes", "Dendritic cells", "Megakaryocytes"]:
        mask = ok2 & (s2["cell"] == t).values
        if mask.sum() < 20:
            print(f"{t:18s} {int(mask.sum()):6d} 太少")
            continue
        cp, n = lib(f"{D}/GSM2560249_2.2.mtx.gz", mask)
        ov = np.log2((cp[idx] + 1) / (c14[idx] + 1))
        r = stats.spearmanr(ov, obs.values).correlation
        bs = []
        for _ in range(1500):
            ii = np.concatenate([byw[w] for w in rng.choice(list(byw), len(byw), True)])
            bs.append(stats.spearmanr(ov[ii], obs.values[ii]).correlation)
        lo, hi = np.percentile(bs, [2.5, 97.5])
        print(f"{t:18s} {n:6d} {r:+9.4f}   [{lo:+.3f}, {hi:+.3f}]")
        rows.append(dict(cond=t, n=n, rho_measured=r, lo=lo, hi=hi))
    pd.DataFrame(rows).to_csv(f"{D}/calib_middle_measured.csv", index=False)
    print(f"\n写出 {D}/calib_middle_measured.csv")
    print("判读：实测重叠在 0.4–0.8 的点越多，中段越值得让模型去算；接近 1.0 的点说明该类型刺激后几乎等同单核响应。")


if __name__ == "__main__":
    main()
