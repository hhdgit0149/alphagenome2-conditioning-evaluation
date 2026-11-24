#!/usr/bin/env python3
"""用全部 CD14 单核细胞（不是抽样的 256 个）算观测 CPM，供表达量门槛使用。"""
import time

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

D = "b_data"
STAMP = time.strftime("%m%d_%H%M%S")


def load(fn, bc):
    bar = pd.read_csv(f"{D}/{bc}", header=None)[0].astype(str).tolist()
    m = scipy.io.mmread(f"{D}/{fn}").tocsr()
    assert m.shape[1] == len(bar), f"{m.shape[1]} != {len(bar)}"
    return m, bar


def main():
    genes = pd.read_csv(f"{D}/GSE96583_batch2.genes.tsv.gz", sep="\t", header=None,
                        names=["ensg", "symbol"])
    ann = pd.read_csv(f"{D}/GSE96583_batch2.total.tsne.df.tsv.gz", sep="\t")

    out = {}
    for tag, mtx, bc in [("ctrl", "GSM2560248_2.1.mtx.gz", "b2.1.barcodes.tsv.gz"),
                         ("stim", "GSM2560249_2.2.mtx.gz", "b2.2.barcodes.tsv.gz")]:
        m, bar = load(mtx, bc)
        sel = ann.reindex(bar)
        mono = ((sel["cell"] == "CD14+ Monocytes") & (sel["multiplets"] == "singlet")
                & (sel["stim"] == ("ctrl" if tag == "ctrl" else "stim")))
        idx = np.where(mono.values)[0]
        sub = m[:, idx]
        tot = np.asarray(sub.sum(0)).ravel()
        gsum = np.asarray(sub.sum(1)).ravel()
        cpm = gsum / gsum.sum() * 1e6
        out[tag] = cpm
        print(f"{tag}: {len(idx)} 个单核细胞，总 UMI {int(tot.sum()):,}，"
              f"每细胞中位 {int(np.median(tot))}")

    df = pd.DataFrame({"ensg": genes.ensg.values, "symbol": genes.symbol.values,
                       "cpm_ctrl_full": out["ctrl"], "cpm_stim_full": out["stim"]})
    assert not df.ensg.duplicated().any(), "ensg 重复"
    df["both_ge1"] = (df.cpm_ctrl_full >= 1) & (df.cpm_stim_full >= 1)
    df["either_ge1"] = (df.cpm_ctrl_full >= 1) | (df.cpm_stim_full >= 1)
    p = f"{D}/obs_cpm_pooled_{STAMP}.csv"
    df.to_csv(p, index=False)
    print(f"\n写出 {p}")
    print(f"两条件都 ≥1 CPM 的基因：{int(df.both_ge1.sum())}/{len(df)}")
    print(f"任一条件 ≥1 CPM 的基因：{int(df.either_ge1.sum())}/{len(df)}")


if __name__ == "__main__":
    main()
