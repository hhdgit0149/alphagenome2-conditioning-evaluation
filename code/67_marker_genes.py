#!/usr/bin/env python3
"""十个 ISG + 谱系标记的逐基因 CPM 小表，供核 CD16 预激那组数与细胞类型标签。

覆盖全部基因（不只 40 窗口内的 246 个），因此 ISG15/IFIT1 等可用。
列：各细胞类型 × 两个文库的 CPM（全部 singlet 细胞，不排除条件侧），
    外加 log2(CD16+/CD14+) 与 log2(IFN 响应) 两列，便于直接对表。
"""
import numpy as np
import pandas as pd
import scipy.io

D = "b_data"
ISG = ["IFIT3", "OAS1", "IFIT1", "MX1", "ISG15", "IRF7", "IFI6", "IFI44L", "STAT1", "OAS2"]
LINEAGE = ["FCGR3A", "CD14", "LYZ", "CD3D", "CD8A", "MS4A1", "NKG7", "PF4", "LILRA4"]
TYPES = ["CD14+ Monocytes", "CD4 T cells", "CD8 T cells", "NK cells", "B cells",
         "Dendritic cells", "Megakaryocytes", "FCGR3A+ Monocytes"]


def main():
    genes = pd.read_csv(f"{D}/GSE96583_batch2.genes.tsv.gz", sep="\t", header=None,
                        names=["ensg", "symbol"])
    ann = pd.read_csv(f"{D}/GSE96583_batch2.total.tsne.df.tsv.gz", sep="\t")
    want = ISG + LINEAGE
    sel = genes[genes.symbol.isin(want)].drop_duplicates("symbol")
    miss = [g for g in want if g not in set(sel.symbol)]
    print(f"命中 {len(sel)}/{len(want)}；缺 {miss}")
    rows = sel.set_index("symbol").reindex(want).dropna(subset=["ensg"])
    idx = rows.ensg.map({e: i for i, e in enumerate(genes.ensg.values)}).astype(int).values

    out = pd.DataFrame({"symbol": rows.index, "ensg": rows.ensg.values})
    for tag, bcf, mtx in [("ctrl", "b2.1", "GSM2560248_2.1"), ("stim", "b2.2", "GSM2560249_2.2")]:
        bc = pd.read_csv(f"{D}/{bcf}.barcodes.tsv.gz", header=None)[0].astype(str).tolist()
        s = ann.reindex(bc)
        m = scipy.io.mmread(f"{D}/{mtx}.mtx.gz").tocsr()
        ok = (s["multiplets"] == "singlet").values
        for t in TYPES:
            cells = np.where(ok & (s["cell"] == t).values)[0]
            cnt = np.asarray(m[:, cells][idx, :].sum(1)).ravel()
            out[f"{tag}|{t}"] = cnt / float(m[:, cells].sum()) * 1e6

    out["log2_CD16_vs_CD14_ctrl"] = np.log2(
        (out["ctrl|FCGR3A+ Monocytes"] + 1) / (out["ctrl|CD14+ Monocytes"] + 1))
    obs = pd.read_csv(f"{D}/de_disjoint_1115_115447.csv").set_index("ensg").obs
    out["log2_IFN_response"] = out.ensg.map(obs)
    cols = ["symbol", "ensg", "log2_CD16_vs_CD14_ctrl", "log2_IFN_response"] + \
           [c for c in out.columns if "|" in c]
    out = out[cols]
    out.to_csv(f"{D}/for_audit_marker_genes.csv", index=False)
    print(out[["symbol", "log2_CD16_vs_CD14_ctrl", "log2_IFN_response"]].round(2).to_string(index=False))
    print(f"\n写出 for_audit_marker_genes.csv：{out.shape}")


if __name__ == "__main__":
    main()
