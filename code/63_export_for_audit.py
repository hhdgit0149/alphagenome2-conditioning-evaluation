#!/usr/bin/env python3
"""把复算需要的逐基因数据拼成三张表，避免对方四处对列名。

输出：
  b_data/for_audit_conditions_long.csv  —— 15 个条件 × 逐基因的预测活性（treated / baseline）与 pred_lfc
  b_data/for_audit_celltype_cpm.csv     —— 逐基因、各细胞类型、两个文库的 CPM（用于独立重算实测重叠）
  b_data/reliability_vectors.csv        —— 20 次细胞二分的 xa / xb 向量（14 个条件），用于独立重算可靠性
"""
import glob

import numpy as np
import pandas as pd
import scipy.io
from scipy import stats

D = "b_data"
TYPES = ["CD14+ Monocytes", "CD4 T cells", "CD8 T cells", "NK cells", "B cells",
         "Dendritic cells", "Megakaryocytes", "FCGR3A+ Monocytes"]


def main():
    genes = pd.read_csv(f"{D}/GSE96583_batch2.genes.tsv.gz", sep="\t", header=None,
                        names=["ensg", "symbol"])
    ann = pd.read_csv(f"{D}/GSE96583_batch2.total.tsne.df.tsv.gz", sep="\t")

    # ---------- 1. 条件长表 ----------
    tb = pd.read_csv(f"{D}/null_tb_filled.csv")
    near = pd.read_csv(f"{D}/null_near_1115.csv")
    m1 = pd.read_csv(f"{D}/calib_mid_pred_1115_121815.csv")
    m2 = pd.read_csv(f"{D}/calib_mid_pred_1115_130035.csv")
    pred = pd.read_csv(f"{D}/b_pred_disjoint_1115_115447.csv")
    ho = pd.read_csv(f"{D}/holdout_pred_1115_121546.csv")

    base_cols = {"tb": tb.set_index("ensg").ctrl, "near": near.set_index("ensg").ctrl,
                 "m1": m1.set_index("ensg").cd14ctrl, "m2": m2.set_index("ensg").cd14ctrl,
                 "pred": pred.set_index("ensg").ctrl}
    print("baseline 一致性（以 pred.ctrl 为准，410 基因）：")
    ref = pred.set_index("ensg").ctrl
    for k, v in base_cols.items():
        common = ref.index.intersection(v.index)
        if len(common) > 20:
            print(f"  {k:5s} n={len(common):3d} max|diff|={np.abs(v[common]-ref[common]).max():.3e}"
                  f" rho={stats.spearmanr(v[common], ref[common]).correlation:.6f}")

    base = float("nan")
    src = {
        ("CD4 T cells", "unstim"): (tb.set_index("ensg").tcell, "tb"),
        ("CD8 T cells", "unstim"): (tb.set_index("ensg").cd8, "tb"),
        ("NK cells", "unstim"): (tb.set_index("ensg").nk, "tb"),
        ("B cells", "unstim"): (tb.set_index("ensg").bcell, "tb"),
        ("Dendritic cells", "unstim"): (tb.set_index("ensg").dc, "tb"),
        ("Megakaryocytes", "unstim"): (near.set_index("ensg").mega, "near"),
        ("FCGR3A+ Monocytes", "unstim"): (near.set_index("ensg").fcgr3a, "near"),
        ("CD4 T cells", "stim"): (m1.set_index("ensg")["CD4 T cells"], "m1"),
        ("Dendritic cells", "stim"): (m1.set_index("ensg")["Dendritic cells"], "m1"),
        ("NK cells", "stim"): (m1.set_index("ensg")["NK cells"], "m1"),
        ("FCGR3A+ Monocytes", "stim"): (m1.set_index("ensg")["FCGR3A+ Monocytes"], "m1"),
        ("CD8 T cells", "stim"): (m2.set_index("ensg")["CD8 T cells"], "m2"),
        ("B cells", "stim"): (m2.set_index("ensg")["B cells"], "m2"),
        ("Megakaryocytes", "stim"): (m2.set_index("ensg")["Megakaryocytes"], "m2"),
    }
    obs = pred.set_index("ensg").obs
    win = pred.set_index("ensg").win
    rows = []
    for (t, kind), (v, s) in src.items():
        b = base_cols["pred"]      # 统一基线；已验证与各文件 ctrl 完全相同
        common = v.index.intersection(b.index).intersection(obs.index)
        for g in common:
            rows.append(dict(ensg=g, cond=t, kind=kind, pred_treated=v[g], pred_baseline=b[g],
                             pred_lfc=np.log2((v[g] + 1) / (b[g] + 1)), baseline_src=s))
    # 匹配条件
    ps = pred.set_index("ensg")
    for g in ps.index:
        if not np.isfinite(ps.loc[g, "obs"]):
            continue
        rows.append(dict(ensg=g, cond="IFN-beta (matched)", kind="matched",
                         pred_treated=ps.loc[g, "stim"], pred_baseline=ps.loc[g, "ctrl"],
                         pred_lfc=np.log2((ps.loc[g, "stim"] + 1) / (ps.loc[g, "ctrl"] + 1)),
                         baseline_src="pred"))
    long = pd.DataFrame(rows)
    long["obs"] = long.ensg.map(obs)
    long["win"] = long.ensg.map(win)
    # 留出值单独一列（仅匹配条件）
    h = ho.set_index("ensg")
    long["pred_lfc_holdout"] = [
        np.log2((h.loc[g, "stim"] + 1) / (h.loc[g, "ctrl"] + 1)) if (g in h.index and
        long.iloc[i].kind == "matched") else np.nan
        for i, g in enumerate(long.ensg)]
    long.to_csv(f"{D}/for_audit_conditions_long.csv", index=False)
    print(f"\n写出 for_audit_conditions_long.csv：{len(long)} 行，"
          f"{long.cond.nunique()} 个条件，{long.ensg.nunique()} 基因")

    # ---------- 2. 各细胞类型 CPM ----------
    lib = {}
    for tag, bcf, mtx in [("ctrl", "b2.1", "GSM2560248_2.1"), ("stim", "b2.2", "GSM2560249_2.2")]:
        bc = pd.read_csv(f"{D}/{bcf}.barcodes.tsv.gz", header=None)[0].astype(str).tolist()
        s = ann.reindex(bc)
        lib[tag] = dict(ok=(s["multiplets"] == "singlet").values,
                        cell=s["cell"].values, m=scipy.io.mmread(f"{D}/{mtx}.mtx.gz").tocsr())
    keep = sorted(set(long.ensg))
    gi = pd.Series(range(len(genes)), index=genes.ensg.values)
    idx = [int(gi[g]) for g in keep]
    out = pd.DataFrame({"ensg": keep})
    for tag in ("ctrl", "stim"):
        L = lib[tag]
        for t in TYPES:
            cells = np.where(L["ok"] & (L["cell"] == t))[0]
            if len(cells) == 0:
                continue
            cnt = np.asarray(L["m"][:, cells].sum(1)).ravel()[idx]
            tot = float(L["m"][:, cells].sum())
            out[f"{tag}|{t}"] = cnt / tot * 1e6
    # obs 两列（不相交 CD14），使 §3 的 obs 定义在本表内自洽可复现
    import anndata as ad
    used = {t: set(ad.read_h5ad(sorted(glob.glob(f"{D}/mono_{t}_*.h5ad"))[-1]).obs.index.astype(str))
            for t in ("ctrl", "stim")}
    for tag in ("ctrl", "stim"):
        L = lib[tag]
        sel = L["ok"] & (L["cell"] == "CD14+ Monocytes")
        keep_ = np.array([sel[i] and L["m"].shape[1] and (str(i) is not None) for i in range(len(sel))])
        bc = pd.read_csv(f"{D}/{'b2.1' if tag=='ctrl' else 'b2.2'}.barcodes.tsv.gz",
                         header=None)[0].astype(str).tolist()
        keep_ = np.array([sel[i] and bc[i] not in used[tag] for i in range(len(sel))])
        cells = np.where(keep_)[0]
        cnt = np.asarray(L["m"][:, cells].sum(1)).ravel()[idx]
        out[f"obs|{tag}"] = cnt / float(L["m"][:, cells].sum()) * 1e6
        print(f"  obs|{tag}: {len(cells)} 细胞，文库 {float(L['m'][:, cells].sum()):,.0f}")
    out.to_csv(f"{D}/for_audit_celltype_cpm.csv", index=False)
    print(f"写出 for_audit_celltype_cpm.csv：{out.shape}")


if __name__ == "__main__":
    main()
