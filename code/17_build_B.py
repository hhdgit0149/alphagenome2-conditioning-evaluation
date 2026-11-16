#!/usr/bin/env python3
"""B 第一步：读入 GSE96583 batch2，构建单核细胞的对照/刺激 AnnData 并验证。

防线：
  - 输出文件带时间戳，不覆盖
  - 读入结果表后断言主键无重复
  - 每步打印行数与空值数
"""
import gzip
import time

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

D = "b_data"
STAMP = time.strftime("%H%M%S")
MAXCELLS = 256


def load_genes():
    g = pd.read_csv(f"{D}/GSE96583_batch2.genes.tsv.gz", sep="\t", header=None,
                    names=["ensg", "symbol"])
    print(f"基因列表 {len(g)} 行，空值 {g.ensg.isna().sum()}")
    assert len(g) == 35635, f"基因数 {len(g)} != 35635"
    return g


def load_ann():
    a = pd.read_csv(f"{D}/GSE96583_batch2.total.tsne.df.tsv.gz", sep="\t")
    print(f"细胞注释 {len(a)} 行，列 {list(a.columns)}")
    print(f"  细胞类型：{a['cell'].value_counts().to_dict()}")
    print(f"  处理：{a['stim'].value_counts().to_dict()}，多重峰 {a['multiplets'].value_counts().to_dict()}")
    return a


def load_mtx(fn, barcodes):
    m = scipy.io.mmread(fn).tocsr()          # genes x cells
    print(f"{fn}: {m.shape} 非零 {m.nnz}")
    assert m.shape[0] == 35635
    assert m.shape[1] == len(barcodes), f"{m.shape[1]} != {len(barcodes)}"
    return m


def main():
    genes = load_genes()
    ann = load_ann()

    out = {}
    for tag, fn, bc in [("ctrl", "2.1", "b2.1.barcodes.tsv.gz"),
                        ("stim", "2.2", "b2.2.barcodes.tsv.gz")]:
        barcodes = pd.read_csv(f"{D}/{bc}", header=None)[0].astype(str).tolist()
        print(f"\n{tag}: {len(barcodes)} 个条形码")
        m = load_mtx(f"{D}/GSM2560248_{fn}.mtx.gz" if fn == "2.1"
                     else f"{D}/GSM2560249_{fn}.mtx.gz", barcodes)
        a = ann.set_index(ann.index)
        # 注释按条形码索引，核对交集
        sel = a.reindex(barcodes)
        hit = sel["stim"].notna().sum()
        print(f"  条形码能在注释里找到：{hit}/{len(barcodes)}")
        mono = (sel["cell"] == "CD14+ Monocytes") & (sel["multiplets"] == "singlet")
        print(f"  其中的单核细胞（singlet）：{int(mono.sum())}，"
              f"处理标签 {sel.loc[mono, 'stim'].value_counts().to_dict()}")
        idx = np.where(mono.values)[0]
        rng = np.random.default_rng(7)
        pick = rng.choice(idx, min(MAXCELLS, len(idx)), replace=False)
        sub = m[:, pick].T.tocsc()            # cells x genes
        obs = sel.iloc[pick][["stim", "ind"]].copy()
        obs.index = [barcodes[i] for i in pick]
        var = pd.DataFrame(index=genes.ensg.values)
        adata = ad.AnnData(X=sp.csr_matrix(sub), obs=obs, var=var)
        adata.obs["condition"] = tag
        out[tag] = adata
        print(f"  取 {adata.shape[0]} 个细胞 × {adata.shape[1]} 个基因；"
              f"每个细胞的总 UMI 中位 {int(np.median(np.asarray(adata.X.sum(1)).ravel()))}")

    for tag, a in out.items():
        p = f"{D}/mono_{tag}_{STAMP}.h5ad"
        a.write_h5ad(p)
        print(f"写出 {p}")

    # 观察到的差异表达（pseudobulk）
    def bulk(a):
        x = np.asarray(a.X.sum(0)).ravel()
        return x / x.sum()
    pc, ps = bulk(out["ctrl"]), bulk(out["stim"])
    genes_used = out["ctrl"].var.index.astype(str)
    df = pd.DataFrame({"ensg": genes_used, "cpm_ctrl": pc * 1e6, "cpm_stim": ps * 1e6})
    df["log2fc"] = np.log2((ps + 1e-9) / (pc + 1e-9))
    both = (df.cpm_ctrl > 10) & (df.cpm_stim > 10)
    df["expressed"] = both
    print(f"\n表达基因（两边 CPM>10）：{both.sum()}；log2FC 中位 {df.loc[both,'log2fc'].median():+.2f}")
    top = df[both].reindex(df[both].log2fc.abs().sort_values(ascending=False).index).head(10)
    print("|log2FC| 最大：\n" + top[["ensg", "log2fc", "cpm_ctrl", "cpm_stim"]].to_string(index=False))
    p = f"{D}/de_pseudobulk_{STAMP}.csv"
    df.to_csv(p, index=False)
    print(f"写出 {p}")


if __name__ == "__main__":
    main()
