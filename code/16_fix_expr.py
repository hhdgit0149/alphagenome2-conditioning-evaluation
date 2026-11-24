#!/usr/bin/env python3
"""修复表达量表：去重（保留非空）+ 补齐缺失键。

背景：13_expression_pass.py 的第一版（全 NaN）和修复版曾同时写同一个文件，
造成重复键。这里按 (variant, cond_role) 去重、保留非空值，再补真正缺失的。
"""
import os
import time

import numpy as np
import pandas as pd

from alphagenome.data import genome
from alphagenome.models.v2 import dna_client as c2, interval_scorers as isc

RAW = "expr_tcell.csv"
CLEAN = "expr_tcell_clean.csv"
CONDS = {"act": "UBERON:0002116_CL:0000896", "rest": "UBERON:0002116_CL:0000624"}


def fetch(m, sc, variant, gene, role):
    v = genome.Variant(variant.split("_")[0], int(variant.split("_")[1]),
                       variant.split("_")[2], variant.split("_")[3])
    iv = v.reference_interval.resize(c2.SEQUENCE_LENGTH_1MB)
    a = m.score_interval(interval=iv, interval_scorers=sc, condition=CONDS[role])[0]
    ids = [str(g).split(".")[0] for g in a.obs.gene_id]
    gi = [i for i, g in enumerate(ids) if g == gene]
    return float(np.nanmean(np.asarray(a.X)[gi[0]])) if gi else np.nan


def main():
    raw = pd.read_csv(RAW)
    print(f"原始 {len(raw)} 行，非空 {raw.pred_ref_expr.notna().sum()}")
    ok = raw.dropna(subset=["pred_ref_expr"]).drop_duplicates(["variant", "cond_role"])
    print(f"去重保留非空：{len(ok)} 条，覆盖 "
          f"{ok.variant.nunique()} 个位点（应覆盖 {795} 个）")

    need = []
    have = set(zip(ok.variant, ok.cond_role))
    d = pd.read_csv("scores_tcell.csv").dropna(subset=["pred_act", "pred_rest"])
    for r in d.itertuples(index=False):
        for role in CONDS:
            if (r.variant, role) not in have:
                need.append((r.variant, r.gene, role))
    print(f"需补 {len(need)} 条", flush=True)

    if need:
        m = c2.create(os.environ["ALPHAGENOME_API_KEY"], timeout=300)
        sc = [isc.RECOMMENDED_INTERVAL_SCORERS["RNA_SEQ"]]
        add = []
        for i, (var, gene, role) in enumerate(need, 1):
            try:
                add.append(dict(variant=var, gene=gene, cond_role=role,
                                pred_ref_expr=fetch(m, sc, var, gene, role)))
            except Exception as e:
                print(f"  {var}/{role} {type(e).__name__}", flush=True)
                time.sleep(2)
            if i % 50 == 0:
                print(f"  {i}/{len(need)}", flush=True)
        add = pd.DataFrame(add)
        ok = pd.concat([ok, add[add.pred_ref_expr.notna()]], ignore_index=True)

    ok.to_csv(CLEAN, index=False)
    print(f"\n写出 {CLEAN}：{len(ok)} 条，非空 {ok.pred_ref_expr.notna().sum()}，"
          f"位点 {ok.variant.nunique()}")
    print(f"空值率 {ok.pred_ref_expr.isna().mean():.1%}（基因不在 ±100 kb 掩膜内的属结构性缺失）")


if __name__ == "__main__":
    main()
