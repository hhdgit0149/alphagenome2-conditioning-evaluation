#!/usr/bin/env python3
"""臂 4：二代 + 错配条件（PREREG 下一阶段对照阶梯）。

Flu 臂的 12 act + 18 rest 位点固定、打分头固定（dn_lfc），只把分子条件换成
两个与 IFN-β 单核无关的内置条件：
  mm1 = 回肠 CD4 激活 T 细胞（UBERON:0002116_CL:0000896）
  mm2 = 回肠 CD4 静息 T 细胞（UBERON:0002116_CL:0000624）
分母沿用主检验的 mono_ctrl 条件（quach_v2 的 dn_lfc/ctrl），故比值只换了分子。
判读：若错配条件下 act/rest 仍分开 → 信号与条件化无关，结论要重写。
"""
import glob
import os
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import anndata as ad
import pandas as pd
from alphagenome.data import genome
from alphagenome.models.v2 import dna_client as c2, variant_scorers as vs2
from scipy import stats

STAMP = time.strftime("%m%d_%H%M%S")
OUT = f"b_data/arm4_mismatch_{STAMP}.csv"
# 错配条件用同批数据的非单核细胞（DENOISED 打分器要求 AnnData）
MM_FILES = {"mm1_Tcell": "b_data/mm_tcell_*.h5ad",
            "mm2_Bcell": "b_data/mm_bcell_*.h5ad"}
WORKERS = 6


def main():
    v = pd.read_csv("b_data/flu_visible_sites.csv")
    s = v[v.group.isin(["act_specific", "rest_specific"])].copy()
    MM = {k: ad.read_h5ad(sorted(glob.glob(v))[-1]) for k, v in MM_FILES.items()}
    print(f"位点 {len(s)}（act {(s.group=='act_specific').sum()} / "
          f"rest {(s.group=='rest_specific').sum()}），错配条件 {list(MM)}")
    lab = pd.read_csv("b_data/Quach_2016_labels_fixed.csv").drop_duplicates(
        ["variant", "gene"])[["variant", "gene", "chrom", "pos", "ref", "alt"]]
    s = s.merge(lab, on=["variant", "gene"], how="left")
    assert s.chrom.notna().all(), "缺坐标"

    done = set()
    for fp in glob.glob("b_data/arm4_mismatch_*.csv"):
        if os.path.getsize(fp) < 200:      # 空/未写入的文件跳过
            continue
        d = pd.read_csv(fp, usecols=["variant", "gene", "mm"])
        done |= set(d.itertuples(index=False, name=None))
    tasks = [(r, mm) for r in s.itertuples(index=False) for mm in MM
             if (r.variant, r.gene, mm) not in done]
    print(f"待打分 {len(tasks)}", flush=True)
    m = c2.create(os.environ["ALPHAGENOME_API_KEY"], timeout=600)
    sc = [vs2.RECOMMENDED_VARIANT_SCORERS["SC_GENE_EXPRESSION_DENOISED"]]

    def one(t):
        r, mm = t
        try:
            var = genome.Variant(r.chrom, int(r.pos), r.ref, r.alt)
            iv = var.reference_interval.resize(c2.SEQUENCE_LENGTH_1MB)
            d = vs2.tidy_scores(m.score_variant(
                interval=iv, variant=var, variant_scorers=sc,
                cell_conditioning=MM[mm], cell_condition_name=mm))
            h = d[d.gene_id == r.gene]
            val = float(h.raw_score.mean()) if len(h) else np.nan
            return dict(variant=r.variant, gene=r.gene, group=r.group, mm=mm,
                        lfc=val, n_tracks=len(h))
        except Exception as e:
            print(f"  {r.variant}/{mm} {type(e).__name__}: {e}", flush=True)
            return None

    rows = []
    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        for rec in ex.map(one, tasks):
            if rec:
                rows.append(rec)
                if len(rows) % 20 == 0:
                    print(f"  {len(rows)}/{len(tasks)}", flush=True)
    pd.DataFrame(rows).to_csv(OUT, index=False)
    print(f"完成 {len(rows)} 条 → {OUT}")


if __name__ == "__main__":
    main()
