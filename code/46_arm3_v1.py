#!/usr/bin/env python3
"""臂 3：一代模型 + CD14 单核轨道（对照阶梯的"当时能拿到的最好"）。

30 个 Flu 位点（12 act / 18 rest），v1 GeneMaskLFC 在 RNA_SEQ 与 DNASE 上打分，
只留 CD14+ monocyte 轨道。v1 无条件输入 → 两条件比值结构上恒为 1，
故本臂检验的是「分组是否与基线序列效应强度混淆」：
若 act/rest 在一代的 |LFC| 上就分得开，主臂的分离可能来自序列特征。
"""
import os
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd
from alphagenome.data import genome
from alphagenome.models.v1 import dna_client as c1, dna_output as dout, variant_scorers as v1s
from scipy import stats

OUT = "b_data/arm3_v1_cd14.csv"
WORKERS = 4


def main():
    v = pd.read_csv("b_data/flu_visible_sites.csv")
    s = v[v.group.isin(["act_specific", "rest_specific"])].copy()
    lab = pd.read_csv("b_data/Quach_2016_labels_fixed.csv").drop_duplicates(
        ["variant", "gene"])[["variant", "gene", "chrom", "pos", "ref", "alt"]]
    s = s.merge(lab, on=["variant", "gene"], how="left")
    print(f"位点 {len(s)}", flush=True)

    m = c1.create(os.environ["ALPHAGENOME_API_KEY"], timeout=600)
    sc = [v1s.GeneMaskLFCScorer(requested_output=dout.OutputType.RNA_SEQ),
          v1s.GeneMaskLFCScorer(requested_output=dout.OutputType.DNASE)]
    ivlen = c1.SUPPORTED_SEQUENCE_LENGTHS["SEQUENCE_LENGTH_1MB"]

    def one(r):
        var = genome.Variant(r.chrom, int(r.pos), r.ref, r.alt)
        iv = var.reference_interval.resize(ivlen)
        d = v1s.tidy_scores(m.score_variant(interval=iv, variant=var,
                                            variant_scorers=sc))
        d = d[d.gene_id == r.gene]
        out = dict(variant=r.variant, gene=r.gene, group=r.group)
        for head, pat in [("rna", "RNA_SEQ"), ("dnase", "DNASE")]:
            h = d[d.output_type.astype(str).str.upper().str.contains(pat)]
            cd14 = h[h.biosample_name.astype(str).str.contains("CD14", na=False)]
            use = cd14 if len(cd14) else h
            out[f"{head}_lfc"] = float(np.nanmean(use.raw_score)) if len(use) else np.nan
            out[f"{head}_n_tracks"] = len(use)
        return out

    rows = []
    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        for rec in ex.map(one, s.itertuples(index=False)):
            rows.append(rec)
    dd = pd.DataFrame(rows)
    dd.to_csv(OUT, index=False)
    print(f"写出 {OUT}：{len(dd)} 行；空值 {dd.isna().sum().to_dict()}")

    for head in ("rna", "dnase"):
        a = dd[dd.group == "act_specific"][f"{head}_lfc"].abs().dropna()
        r = dd[dd.group == "rest_specific"][f"{head}_lfc"].abs().dropna()
        p = stats.mannwhitneyu(a, r, alternative="two-sided").pvalue
        print(f"{head}: |LFC| 中位 {a.median():.3f}/{r.median():.3f}  "
              f"n={len(a)}/{len(r)}  秩和 p={p:.4f}")


if __name__ == "__main__":
    main()
