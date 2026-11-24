#!/usr/bin/env python3
"""补 pred_ref_expr：每个「位点 × 条件」的预测表达量（参考等位基因）。

用 interval_scorers.RECOMMENDED_INTERVAL_SCORERS["RNA_SEQ"]
（GeneMaskScorer, aggregation=MEAN），在参考序列上直接给出每个基因的预测表达量，
对 12 条 RNA-seq 轨道取均值。条件与 A 组打分完全一致。
"""
import os
import time

import numpy as np
import pandas as pd

from alphagenome.data import genome
from alphagenome.models.v2 import dna_client as c2, interval_scorers as isc

OUT = "expr_tcell.csv"
CONDS = {"act": "UBERON:0002116_CL:0000896", "rest": "UBERON:0002116_CL:0000624"}


def main():
    d = pd.read_csv("scores_tcell.csv").dropna(subset=["pred_act", "pred_rest"])
    print(f"位点 {len(d)}", flush=True)
    done = set()
    if os.path.exists(OUT):
        old = pd.read_csv(OUT)
        done = set(zip(old.variant, old.cond_role))
        print(f"已完成 {len(done)} 条（位点×条件），续跑", flush=True)

    m = c2.create(os.environ["ALPHAGENOME_API_KEY"], timeout=300)
    sc = [isc.RECOMMENDED_INTERVAL_SCORERS["RNA_SEQ"]]
    n = 0
    for r in d.itertuples(index=False):
        todo = [k for k in CONDS if (r.variant, k) not in done]
        if not todo:
            continue
        try:
            v = genome.Variant(r.variant.split("_")[0], int(r.variant.split("_")[1]),
                               r.variant.split("_")[2], r.variant.split("_")[3])
            iv = v.reference_interval.resize(c2.SEQUENCE_LENGTH_1MB)
        except Exception as e:
            print(f"  {r.variant} 构造失败 {type(e).__name__}", flush=True)
            continue
        for role in todo:
            try:
                a = m.score_interval(interval=iv, interval_scorers=sc,
                                     condition=CONDS[role])[0]
                ids = [str(g).split(".")[0] for g in a.obs.gene_id]
                gi = [i for i, g in enumerate(ids) if g == r.gene]
                val = (float(np.nanmean(np.asarray(a.X)[gi[0]])) if gi else np.nan)
                pd.DataFrame([dict(variant=r.variant, gene=r.gene, cond_role=role,
                                   pred_ref_expr=val)]).to_csv(
                    OUT, mode="a", index=False, header=not os.path.exists(OUT))
                n += 1
            except Exception as e:
                print(f"  {r.variant}/{role} 失败 {type(e).__name__}: {str(e)[:70]}", flush=True)
                time.sleep(2)
        if n % 40 == 0:
            print(f"  {n} 条", flush=True)
    print(f"完成 {n} 条", flush=True)


if __name__ == "__main__":
    main()
