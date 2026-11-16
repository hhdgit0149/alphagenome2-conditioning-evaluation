#!/usr/bin/env python3
"""A 组主测试：模型能否区分「激活特异 / 静息特异 / 共享」三类 eQTL。

对每个位点，分别在「激活 T 细胞」和「静息 T 细胞」两个内置条件下预测
它对目标基因表达的影响，比较效应强度。

事先定死的判据（写在代码里，避免事后调参）：
  - 主指标 = |预测效应|，不是方向
  - 激活特异：|pred_act| > |pred_rest|
  - 静息特异：|pred_act| < |pred_rest|
  - 共享：两者相近
  - 用配对 Wilcoxon 检验，按研究分层
"""
import os
import sys
import time

import numpy as np
import pandas as pd

from alphagenome.data import genome
from alphagenome.models.v2 import dna_client as c2, variant_scorers as vs2

OUT = "scores_tcell.csv"
# 回肠 CD4 激活 vs 静息（唯一有同名对照的 T 细胞内置条件之一）
CONDS = {"act": "UBERON:0002116_CL:0000896", "rest": "UBERON:0002116_CL:0000624"}
SCORER = "RNA_SEQ"


def main():
    limit = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    df = pd.read_csv("stim_variants_groups.csv")
    df = df[df.study.isin(["Cytoimmgen", "Schmiedel_2018"])].copy()
    if limit:
        df = df.groupby("group", group_keys=False).head(limit)
    done = set()
    if os.path.exists(OUT):
        done = set(pd.read_csv(OUT, usecols=["variant"], dtype=str)["variant"])
        print(f"已完成 {len(done)}，续跑", flush=True)
    todo = df[~df.variant.isin(done)]
    print(f"待打分 {len(todo)} 个位点（T 细胞，两条条件）", flush=True)
    if todo.empty:
        return

    m = c2.create(os.environ["ALPHAGENOME_API_KEY"], timeout=300)
    sc = [vs2.RECOMMENDED_VARIANT_SCORERS[SCORER]]
    n_ok = 0
    for i, r in enumerate(todo.itertuples(index=False), 1):
        try:
            v = genome.Variant(r.chrom, int(r.pos), r.ref, r.alt)
            iv = v.reference_interval.resize(c2.SEQUENCE_LENGTH_1MB)
            rec = dict(variant=r.variant, gene=r.gene, group=r.group,
                       study=r.study, cell=r.cell, beta=r.beta)
            for lab, cond in CONDS.items():
                d = vs2.tidy_scores(m.score_variant(
                    interval=iv, variant=v, condition=cond, variant_scorers=sc))
                h = d[d.gene_id == r.gene]
                rec[f"pred_{lab}"] = float(h.raw_score.mean()) if len(h) else np.nan
                rec[f"n_{lab}"] = len(h)
            pd.DataFrame([rec]).to_csv(OUT, mode="a", index=False,
                                       header=not os.path.exists(OUT))
            n_ok += 1
        except Exception as e:
            print(f"  {r.variant} 失败 {type(e).__name__}: {str(e)[:80]}", flush=True)
            time.sleep(2)
        if i % 25 == 0:
            print(f"  {i}/{len(todo)}", flush=True)
    print(f"完成 {n_ok} 个", flush=True)


if __name__ == "__main__":
    main()
