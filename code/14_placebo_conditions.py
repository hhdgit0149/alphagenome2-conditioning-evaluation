#!/usr/bin/env python3
"""Claude 建议的第 2 步：安慰剂条件 + 复现条件。

对同一批「可见位点」，换两组条件重新打分：
  placebo   : 心脏成纤维细胞 激活 vs 非激活（与 T 细胞无关）→ 组间差异应消失
  replicate : 横结肠 CD4 激活 vs 非激活（另一对 T 细胞条件）→ 差异应重现
"""
import os
import time

import numpy as np
import pandas as pd

from alphagenome.data import genome
from alphagenome.models.v2 import dna_client as c2, variant_scorers as vs2

OUT = "scores_tcell_altcond.csv"
PAIRS = {
    "placebo":   ("UBERON:0002078_NTR:0008000", "UBERON:0002078_CL:0000057"),
    "replicate": ("UBERON:0001157_CL:0000896", "UBERON:0001157_CL:0000624"),
}
OLDCONDS = {"act": "UBERON:0002116_CL:0000896", "rest": "UBERON:0002116_CL:0000624"}


def main():
    d = pd.read_csv("scores_tcell.csv").dropna(subset=["pred_act", "pred_rest"])
    d["vis"] = d[["pred_act", "pred_rest"]].abs().max(axis=1) > 0.1
    d = d[d.vis]
    print(f"可见位点 {len(d)}", flush=True)
    done = set()
    if os.path.exists(OUT):
        old = pd.read_csv(OUT)
        done = set(zip(old.variant, old.pair))
    todo = [r for r in d.itertuples(index=False) if (r.variant, "placebo") not in done
            or (r.variant, "replicate") not in done]
    print(f"待打分 {len(todo)} 个位点 × 2 组条件", flush=True)

    m = c2.create(os.environ["ALPHAGENOME_API_KEY"], timeout=300)
    sc = [vs2.RECOMMENDED_VARIANT_SCORERS["RNA_SEQ"]]
    n = 0
    for r in todo:
        try:
            v = genome.Variant(r.variant.split("_")[0], int(r.variant.split("_")[1]),
                               r.variant.split("_")[2], r.variant.split("_")[3])
            iv = v.reference_interval.resize(c2.SEQUENCE_LENGTH_1MB)
        except Exception:
            continue
        for pair, (ca, cr) in PAIRS.items():
            if (r.variant, pair) in done:
                continue
            rec = dict(variant=r.variant, gene=r.gene, group=r.group,
                       study=r.study, cell=r.cell, pair=pair)
            ok = True
            for role, cond in [("act", ca), ("rest", cr)]:
                try:
                    t = vs2.tidy_scores(m.score_variant(
                        interval=iv, variant=v, condition=cond, variant_scorers=sc))
                    h = t[t.gene_id == r.gene]
                    rec[f"pred_{role}"] = float(h.raw_score.mean()) if len(h) else np.nan
                except Exception as e:
                    print(f"  {r.variant}/{pair}/{role} {type(e).__name__}", flush=True)
                    ok = False
                    time.sleep(2)
                    break
            if ok:
                pd.DataFrame([rec]).to_csv(OUT, mode="a", index=False,
                                           header=not os.path.exists(OUT))
                n += 1
        if n % 20 == 0:
            print(f"  {n}", flush=True)
    print(f"完成 {n}", flush=True)


if __name__ == "__main__":
    main()
