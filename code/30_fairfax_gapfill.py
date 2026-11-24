#!/usr/bin/env python3
"""Fairfax 主测试 v2：四个打分器一次调用拿全，8 并发，单写者。

打分器：
  lfc_dn  = SC_GENE_EXPRESSION_DENOISED 的对数比（**主指标**，自定义条件训练时接的头）
  act_dn  = SC_GENE_EXPRESSION_DENOISED 的活性（用于表达量门槛）
  lfc_rna = RNA_SEQ 的对数比（敏感性分析）
  act_rna = RNA_SEQ 的活性（备用门槛）

分组标签相对「该位点所属刺激」定义，一个位点在多种刺激下可有不同标签，
因此打分只按 (variant, gene) 去重，标签在分析阶段展开。
"""
import glob
import os
import time
from concurrent.futures import ThreadPoolExecutor

import anndata as ad
import numpy as np
import pandas as pd
from alphagenome.data import genome
from alphagenome.models.v2 import dna_client as c2, dna_output as do, variant_scorers as vs2

STAMP = time.strftime("%m%d_%H%M%S")
OUT = f"b_data/fairfax_v2b_{STAMP}.csv"
WORKERS = 8


def scorers():
    return [
        vs2.RECOMMENDED_VARIANT_SCORERS["SC_GENE_EXPRESSION_DENOISED"],
        vs2.GeneMaskActiveScorer(requested_output=do.OutputType.SC_GENE_EXPRESSION_DENOISED),
        vs2.RECOMMENDED_VARIANT_SCORERS["RNA_SEQ"],
        vs2.RECOMMENDED_VARIANT_SCORERS["RNA_SEQ_ACTIVE"],
    ]


def tag_of(s):
    s = str(s)
    if "DENOISED" in s:
        return "dn_lfc" if "LFC" in s else "dn_act"
    return "rna_lfc" if "LFC" in s else "rna_act"


def build_conditions():
    ctrl = ad.read_h5ad(sorted(glob.glob("b_data/mono_ctrl_*.h5ad"))[-1])
    stim = ad.read_h5ad(sorted(glob.glob("b_data/mono_stim_*.h5ad"))[-1])
    idx = np.arange(ctrl.n_obs)
    rng = np.random.default_rng(11)
    rng.shuffle(idx)
    return {"stim": stim, "ctrl": ctrl,
            "halfA": ctrl[idx[:128]].copy(), "halfB": ctrl[idx[128:]].copy()}


def main():
    g = pd.read_csv("b_data/fairfax_labels_fixed.csv")
    f = g.copy()
    lab = f[["variant", "gene", "group", "stim"]].drop_duplicates()
    print(f"打分单元 {f.drop_duplicates(['variant', 'gene']).shape[0]} 个（variant+gene）")
    print(f"标签行 {len(lab)}（同一变异在不同刺激下可有多行）")

    vv = f.drop_duplicates(["variant", "gene"])[
        ["variant", "gene", "chrom", "pos", "ref", "alt"]]

    done = set()
    base = pd.read_csv("b_data/fairfax_v2_dedup.csv",
                       usecols=["variant", "gene", "cond_role"])
    done = set(base.itertuples(index=False, name=None))
    print(f"已有 (variant,gene,cond_role) {len(done)} 条", flush=True)
    print(f"已完成 {len(done)} 条", flush=True)

    conds = build_conditions()
    m = c2.create(os.environ["ALPHAGENOME_API_KEY"], timeout=600)
    sc = scorers()

    tasks = [(r, role) for r in vv.itertuples(index=False) for role in conds
             if (r.variant, r.gene, role) not in done]
    print(f"待跑 {len(tasks)} 条，{WORKERS} 并发", flush=True)

    def one(t):
        r, role = t
        try:
            v = genome.Variant(r.chrom, int(r.pos), r.ref, r.alt)
            iv = v.reference_interval.resize(c2.SEQUENCE_LENGTH_1MB)
            d = vs2.tidy_scores(m.score_variant(
                interval=iv, variant=v, variant_scorers=sc,
                cell_conditioning=conds[role], cell_condition_name=role))
            d = d[d.gene_id == r.gene].copy()
            if d.empty:
                return dict(variant=r.variant, gene=r.gene, cond_role=role,
                            dn_lfc=np.nan, dn_act=np.nan, rna_lfc=np.nan,
                            rna_act=np.nan, n_tracks=0)
            d["tag"] = d.variant_scorer.map(tag_of)
            g_ = d.groupby("tag").raw_score.mean()
            return dict(variant=r.variant, gene=r.gene, cond_role=role,
                        dn_lfc=float(g_.get("dn_lfc", np.nan)),
                        dn_act=float(g_.get("dn_act", np.nan)),
                        rna_lfc=float(g_.get("rna_lfc", np.nan)),
                        rna_act=float(g_.get("rna_act", np.nan)),
                        n_tracks=len(d))
        except Exception as e:
            print(f"  {r.variant}/{role} {type(e).__name__}", flush=True)
            return None

    buf, n, t0 = [], 0, time.time()
    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        for rec in ex.map(one, tasks):
            if rec:
                buf.append(rec)
                n += 1
            if len(buf) >= 50:
                pd.DataFrame(buf).to_csv(OUT, mode="a", index=False,
                                         header=not os.path.exists(OUT))
                buf = []
                el = time.time() - t0
                print(f"  {n}/{len(tasks)}｜{n/max(el,1)*60:.0f} 条/分钟", flush=True)
    if buf:
        pd.DataFrame(buf).to_csv(OUT, mode="a", index=False,
                                 header=not os.path.exists(OUT))
    lab.to_csv(f"b_data/fairfax_labels_{STAMP}.csv", index=False)
    print(f"完成 {n} 条 → {OUT}")


if __name__ == "__main__":
    main()
