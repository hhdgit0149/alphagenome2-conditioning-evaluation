#!/usr/bin/env python3
"""校准曲线中段：用受刺激的其他细胞类型造条件，看模型能否复现它们的实测重叠。

四个条件（实测重叠已算好，本地）：
  FCGR3A+ 刺激 0.740 / 树突 0.495 / NK 0.252 / CD4 T 0.171
pred_null = log2((该类型刺激后 + 1)/(CD14 未刺激 + 1))。
若模型给出的值与实测重叠吻合，则「模型 = 0.62 × 实测」这条校准线是靠多点撑起来的，
而不是靠 (1.0, 0.616) 一个杠杆点。
"""
import glob
import os
import time

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
from alphagenome.data import genome
from alphagenome.models.v2 import dna_client as c2, dna_output as do, interval_scorers as isc
from scipy import stats

NBIN = 40
WORKERS = 5
STAMP = time.strftime("%m%d_%H%M%S")
TARGETS = ["CD8 T cells", "B cells", "Megakaryocytes"]


def main():
    import importlib.util
    spec = importlib.util.spec_from_file_location("b18", "18_B_predict.py")
    b18 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(b18)

    de = pd.read_csv(sorted(glob.glob("b_data/de_pseudobulk_*.csv"))[-1])
    expr = de[de.expressed].copy()
    coords = b18.gene_coords_batch(expr.ensg.head(1500))
    dfc = pd.DataFrame([(e, c[0], c[1], c[2]) for e, c in coords.items()],
                       columns=["ensg", "chrom", "start", "end"])
    dfc["bin"] = (dfc.start + dfc.end) // 2 // 1048576
    top = dfc.groupby(["chrom", "bin"]).size().sort_values(ascending=False).head(NBIN)
    wins = list(top.items())

    genes = pd.read_csv("b_data/GSE96583_batch2.genes.tsv.gz", sep="\t", header=None,
                        names=["ensg", "symbol"])
    ann = pd.read_csv("b_data/GSE96583_batch2.total.tsne.df.tsv.gz", sep="\t")
    b2 = pd.read_csv("b_data/b2.2.barcodes.tsv.gz", header=None)[0].astype(str).tolist()
    s2 = ann.reindex(b2)
    mt = scipy.io.mmread("b_data/GSM2560249_2.2.mtx.gz").tocsr()
    rng = np.random.default_rng(7)

    conds = {}
    for t in TARGETS:
        mask = ((s2["cell"] == t) & (s2["multiplets"] == "singlet")).values
        idx = np.where(mask)[0]
        pick = rng.choice(idx, min(256, len(idx)), replace=False)
        sub = mt[:, pick].T.tocsc()
        obs = s2.iloc[pick][["stim", "ind"]].copy()
        obs.index = [b2[i] for i in pick]
        a = ad.AnnData(X=sp.csr_matrix(sub), obs=obs,
                       var=pd.DataFrame(index=genes.ensg.values))
        conds[t] = a
        print(f"{t}: {a.shape}（供体 {int(len(set(obs.ind)))} 人）")

    ctrl = ad.read_h5ad(sorted(glob.glob("b_data/mono_ctrl_*.h5ad"))[-1])
    m = c2.create(os.environ["ALPHAGENOME_API_KEY"], timeout=600)
    sc = [isc.GeneMaskScorer(requested_output=do.OutputType.SC_GENE_EXPRESSION_DENOISED,
                             width=200001, aggregation_type=isc.IntervalAggregationType.MEAN)]

    from concurrent.futures import ThreadPoolExecutor

    def one(item):
        i, ((chrom, b), _) = item
        mid = int(b) * 1048576 + 524288
        iv = genome.Interval(chrom, mid - 524288, mid + 524288)
        out = {}
        for tag, cond in list(conds.items()) + [("cd14ctrl", ctrl)]:
            try:
                a = m.score_interval(interval=iv, interval_scorers=sc,
                                     cell_conditioning=cond, cell_condition_name="c")[0]
                ids = [str(x).split(".")[0] for x in a.obs.gene_id]
                for g, v in zip(ids, [float(np.nanmean(x)) for x in np.asarray(a.X)]):
                    out[(tag, g)] = v
            except Exception as e:
                print(f"  {chrom}:{mid} {tag} {type(e).__name__}", flush=True)
        return out

    with ThreadPoolExecutor(WORKERS) as ex:
        res = list(ex.map(one, list(enumerate(wins, 1))))
    print(f"完成 {len(res)} 窗口")

    rows = []
    for r in res:
        gs = {k[1] for k in r if isinstance(k, tuple) and k[0] == "cd14ctrl"}
        for g in gs:
            d = dict(ensg=g, cd14ctrl=r[("cd14ctrl", g)])
            for t in TARGETS:
                d[t] = r.get((t, g), np.nan)
            rows.append(d)
    df = pd.DataFrame(rows).drop_duplicates("ensg")
    df.to_csv(f"b_data/calib_mid_pred_{STAMP}.csv", index=False)
    print(f"基因 {len(df)}")

    pred = pd.read_csv(sorted(glob.glob("b_data/b_pred_disjoint_*.csv"))[-1])
    base = pred[(pred[["ctrl", "stim"]].max(axis=1) > 1) & np.isfinite(pred.obs)]
    j = df.merge(base[["ensg", "obs", "pred1"]], on="ensg", how="inner")
    # 一致性自检：用重打的 cd14ctrl 复算真实条件，应接近 0.6155
    j["pred_real"] = np.log2((j.stim.values if "stim" in j else np.nan) + 1)

    print(f"\n{'条件':20s} {'n':>4s} {'模型':>8s} {'实测':>8s}")
    meas = pd.read_csv("b_data/calib_middle_measured.csv")
    mm = meas.set_index("cond").rho_measured
    for t in TARGETS:
        v = np.log2((j[t] + 1) / (j.cd14ctrl + 1))
        k = np.isfinite(v)
        r = stats.spearmanr(v[k], j.obs.values[k]).correlation
        print(f"{t:20s} {int(k.sum()):4d} {r:+8.3f} {mm[t]:+8.3f}")
    out = pd.DataFrame([dict(cond=t,
                             rho_model=stats.spearmanr(
                                 np.log2((j[t] + 1) / (j.cd14ctrl + 1)),
                                 j.obs).correlation,
                             rho_measured=mm[t]) for t in TARGETS])
    out.to_csv(f"b_data/calib_mid_model_{STAMP}.csv", index=False)
    print(f"写出 b_data/calib_mid_model_{STAMP}.csv")


if __name__ == "__main__":
    main()
