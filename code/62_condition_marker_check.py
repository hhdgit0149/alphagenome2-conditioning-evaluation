#!/usr/bin/env python3
"""校验七个条件（及两个单核条件）的细胞标签是否与标记基因一致。

两种可能的错都要能抓到：
  (1) 源注释里 cell 标签本身错位 → 全库按类型平均即可看出
  (2) 造条件时抽错细胞 / 文件名与内容不符 → 必须查实际打进 AnnData 的那 256 个细胞

做法：
  A. 按注释的 cell 标签分层，算标记基因的平均表达（CPM）
  B. 重建 60_calib_middle_score.py 里的抽样（rng seed 7、两次的 TARGETS 顺序），
     对实际被抽中的 256 个细胞再算一遍
  C. 逐个打开 b_data/mm_*.h5ad，用其 obs.index 回查注释标签，确认与文件名一致
"""
import glob

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io

D = "b_data"
MARKERS = {"CD8A": "CD8 T", "CD8B": "CD8 T", "CD3D": "T", "CD3E": "T", "IL7R": "CD4 T",
           "MS4A1": "B", "CD79A": "B", "NKG7": "NK", "GNLY": "NK",
           "FCGR3A": "CD16 mono", "CD14": "CD14 mono", "LYZ": "myeloid",
           "LILRA4": "DC", "CLEC4C": "DC", "PF4": "Mega", "PPBP": "Mega"}


def main():
    genes = pd.read_csv(f"{D}/GSE96583_batch2.genes.tsv.gz", sep="\t", header=None,
                        names=["ensg", "symbol"])
    sym2i = {}
    for i, s in enumerate(genes.symbol.values):
        if s in MARKERS and s not in sym2i:
            sym2i[s] = i
    missing = [s for s in MARKERS if s not in sym2i]
    print(f"标记基因命中 {len(sym2i)}/{len(MARKERS)}，缺：{missing}\n")

    ann = pd.read_csv(f"{D}/GSE96583_batch2.total.tsne.df.tsv.gz", sep="\t")
    lib = {}
    for tag, bcf, mtx in [("ctrl", "b2.1", "GSM2560248_2.1"), ("stim", "b2.2", "GSM2560249_2.2")]:
        bc = pd.read_csv(f"{D}/{bcf}.barcodes.tsv.gz", header=None)[0].astype(str).tolist()
        s = ann.reindex(bc)
        m = scipy.io.mmread(f"{D}/{mtx}.mtx.gz").tocsr()
        lib[tag] = dict(bc=bc, s=s, m=m,
                        ok=(s["multiplets"] == "singlet").values,
                        cell=s["cell"].values)

    rows = [s for s in MARKERS if s in sym2i]
    cols = [int(sym2i[s]) for s in rows]

    def profile(m, cells):
        sub = m[:, cells][cols, :]
        cnt = np.asarray(sub.sum(1)).ravel()
        tot = float(m[:, cells].sum())
        return cnt / tot * 1e6

    print("=" * 100)
    print("A. 按注释标签分层（stim 文库，全部 singlet）")
    print("=" * 100)
    L = lib["stim"]
    types = ["CD4 T cells", "CD8 T cells", "NK cells", "B cells", "CD14+ Monocytes",
             "FCGR3A+ Monocytes", "Dendritic cells", "Megakaryocytes"]
    print(f"{'类型':20s}" + "".join(f"{s:>8s}" for s in rows))
    for t in types:
        cells = np.where(L["ok"] & (L["cell"] == t))[0]
        p = profile(L["m"], cells)
        print(f"{t:20s}" + "".join(f"{v:8.1f}" for v in p))

    print("\n" + "=" * 100)
    print("B. 实际被抽中打分的 256 个细胞（重建 seed 7 抽样）")
    print("=" * 100)
    for tag, targets in [("ctrl", ["FCGR3A+ Monocytes", "Dendritic cells", "NK cells", "CD4 T cells"]),
                         ("stim", ["CD8 T cells", "B cells", "Megakaryocytes"]),
                         ("stim", ["FCGR3A+ Monocytes", "Dendritic cells", "NK cells", "CD4 T cells"])]:
        L = lib[tag]
        rng = np.random.default_rng(7)
        for t in targets:
            pool = np.where(L["ok"] & (L["cell"] == t))[0]
            if len(pool) == 0:
                continue
            pick = rng.choice(pool, min(256, len(pool)), replace=False)
            p = profile(L["m"], pick)
            print(f"[{tag}] {t:20s}" + "".join(f"{v:8.1f}" for v in p))

    print("\n" + "=" * 100)
    print("C. mm_*.h5ad 文件内容 vs 文件名")
    print("=" * 100)
    for f in sorted(glob.glob(f"{D}/mm_*.h5ad")):
        a = ad.read_h5ad(f)
        ids = a.obs.index.astype(str).tolist()
        lab = ann.reindex(ids)["cell"]
        vc = lab.value_counts()
        top = vc.index[0] if len(vc) else "?"
        frac = vc.iloc[0] / len(lab) if len(vc) else 0
        name = f.split("mm_")[1].split("_")[0]
        flag = "OK " if frac > 0.95 else "!! "
        print(f"{flag}{f.split(chr(92))[-1]:32s} n={len(ids):4d} 主导标签 {top:22s} {frac:.1%}  其余 {vc.drop(top).to_dict()}")


if __name__ == "__main__":
    main()
