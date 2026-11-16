#!/usr/bin/env python3
"""为臂 4 构建错配细胞条件 AnnData（同批 GSE96583 数据，非单核细胞）。

mm_tcell = CD4 T 细胞（ctrl），mm_bcell = B 细胞（ctrl）。
同一实验、同一 AnnData 通路，只换细胞身份，供 DENOISED 打分器使用。
"""
import time

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

D = "b_data"
STAMP = time.strftime("%m%d_%H%M%S")
MAXCELLS = 256
TARGETS = {"mm_tcell": "CD4 T cells", "mm_bcell": "B cells",
           "mm_cd8": "CD8 T cells", "mm_nk": "NK cells",
           "mm_dc": "Dendritic cells",
           "mm_fcgr3a": "FCGR3A+ Monocytes", "mm_mega": "Megakaryocytes"}


def main():
    genes = pd.read_csv(f"{D}/GSE96583_batch2.genes.tsv.gz", sep="\t", header=None,
                        names=["ensg", "symbol"])
    ann = pd.read_csv(f"{D}/GSE96583_batch2.total.tsne.df.tsv.gz", sep="\t")
    print(f"细胞类型：{ann['cell'].value_counts().to_dict()}")

    barcodes = pd.read_csv(f"{D}/b2.1.barcodes.tsv.gz", header=None)[0].astype(str).tolist()
    m = scipy.io.mmread(f"{D}/GSM2560248_2.1.mtx.gz").tocsr()
    assert m.shape[1] == len(barcodes)
    sel = ann.reindex(barcodes)

    for tag, celltype in TARGETS.items():
        hit = (sel["cell"] == celltype) & (sel["multiplets"] == "singlet") & (sel["stim"] == "ctrl")
        idx = np.where(hit.values)[0]
        print(f"{tag}: {celltype} 单细胞 {len(idx)}")
        rng = np.random.default_rng(7)
        pick = rng.choice(idx, min(MAXCELLS, len(idx)), replace=False)
        sub = m[:, pick].T.tocsc()
        obs = sel.iloc[pick][["stim", "ind"]].copy()
        obs.index = [barcodes[i] for i in pick]
        a = ad.AnnData(X=sp.csr_matrix(sub), obs=obs,
                       var=pd.DataFrame(index=genes.ensg.values.astype(str)))
        a.obs["condition"] = tag
        p = f"{D}/{tag}_{STAMP}.h5ad"
        a.write_h5ad(p)
        print(f"  写出 {p}：{a.shape[0]} 细胞 × {a.shape[1]} 基因")


if __name__ == "__main__":
    main()
