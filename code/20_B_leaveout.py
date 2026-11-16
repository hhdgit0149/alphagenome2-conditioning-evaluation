#!/usr/bin/env python3
"""B 加固：把被评测的基因从细胞编码器输入里删掉，再看预测是否仍成立。

逻辑：细胞编码器的输入含这些基因的表达量，所以原版只能说明"喂进去能吐出来"。
删掉被评测窗口内的基因后，模型只能靠其余基因反映的细胞状态 + DNA 序列，
这时相关性还在，才算模型真的从序列/细胞状态推出了这些基因的变化。
"""
import glob
import time

import anndata as ad
import numpy as np
import pandas as pd
from alphagenome.data import genome
from alphagenome.models.v2 import dna_client as c2, dna_output as do, interval_scorers as isc
from scipy import stats

STAMP = time.strftime("%m%d_%H%M%S")


def main():
    ac = ad.read_h5ad(sorted(glob.glob("b_data/mono_ctrl_*.h5ad"))[-1])
    as_ = ad.read_h5ad(sorted(glob.glob("b_data/mono_stim_*.h5ad"))[-1])
    pred = pd.read_csv(sorted(glob.glob("b_data/b_pred_*.csv"))[-1])
    rec = pd.read_csv("b_data/b_recomputed.csv")
    print(f"原始 ctrl {ac.shape}｜stim {as_.shape}｜预测基因 {len(pred)}｜已分析 {len(rec)}")

    # 要删的基因：被预测过的窗口内基因
    drop = set(pred.ensg.astype(str))
    keep_c = ~ac.var_names.astype(str).isin(drop)
    keep_s = ~as_.var_names.astype(str).isin(drop)
    ac2, as2 = ac[:, keep_c].copy(), as_[:, keep_s].copy()
    print(f"删掉 {len(drop)} 个基因后：ctrl {ac2.shape}｜stim {as2.shape}")
    assert ac2.shape[1] == as2.shape[1], "两侧基因数不一致"

    idx = np.arange(ac2.n_obs)
    rng = np.random.default_rng(11)
    rng.shuffle(idx)
    h1, h2 = ac2[idx[:128]].copy(), ac2[idx[128:]].copy()

    m = c2.create(__import__("os").environ["ALPHAGENOME_API_KEY"], timeout=600)
    sc = [isc.GeneMaskScorer(requested_output=do.OutputType.SC_GENE_EXPRESSION_DENOISED,
                             width=200001,
                             aggregation_type=isc.IntervalAggregationType.MEAN)]
    wins = sorted({(w.split(":")[0], int(w.split(":")[1])) for w in rec.win.unique()
                   if w != "NA"})
    print(f"窗口数 {len(wins)}")

    conds = {"ctrl": ac2, "stim": as2, "half1": h1, "half2": h2}
    rows = []
    for i, (chrom, b) in enumerate(wins, 1):
        mid = b * 1048576 + 524288
        iv = genome.Interval(chrom, mid - 524288, mid + 524288)
        vals = {}
        for tag, cond in conds.items():
            try:
                a = m.score_interval(interval=iv, interval_scorers=sc,
                                     cell_conditioning=cond, cell_condition_name=tag)[0]
                ids = [str(x).split(".")[0] for x in a.obs.gene_id]
                vals[tag] = dict(zip(ids, [float(np.nanmean(x)) for x in np.asarray(a.X)]))
            except Exception as e:
                print(f"  {chrom}:{b} {tag} {type(e).__name__}", flush=True)
                vals[tag] = {}
        for e in set(vals["ctrl"]) & set(vals["stim"]):
            rows.append(dict(ensg=e, ctrl=vals["ctrl"][e], stim=vals["stim"][e],
                             half1=vals["half1"].get(e, np.nan),
                             half2=vals["half2"].get(e, np.nan)))
        if i % 10 == 0:
            print(f"  窗口 {i}/{len(wins)}｜基因 {len(rows)}", flush=True)

    p = pd.DataFrame(rows).drop_duplicates("ensg")
    p["pred"] = np.log2((p.stim + 1) / (p.ctrl + 1))
    p["noise"] = np.log2((p.half1 + 1) / (p.half2 + 1))
    p = p.merge(rec[["ensg", "obs"]], on="ensg", how="inner")
    p.to_csv(f"b_data/b_leaveout_{STAMP}.csv", index=False)
    print(f"\n可用于比较的基因 {len(p)}")

    sp = stats.spearmanr(p.pred, p.obs)
    rn = stats.spearmanr(p.noise, p.obs).correlation
    print(f"留出版本：Spearman = {sp.correlation:.3f} (p={sp.pvalue:.2g})")
    print(f"噪声基线：Spearman = {rn:.3f}")
    print("判据：>0.4 通过，<0.2 不通过 → " +
          ("通过" if sp.correlation > 0.4 else
           ("不通过" if sp.correlation < 0.2 else "灰区")))
    print(f"写出 b_data/b_leaveout_{STAMP}.csv")


if __name__ == "__main__":
    main()
