#!/usr/bin/env python3
"""Fairfax 主测试打分：1,035 个位点 × 4 个自定义条件。

四个条件都来自同一批单核细胞：
  stim  = IFN-β 刺激 6h（2.2）
  ctrl  = 未刺激（2.1）
  halfA / halfB = 未刺激细胞随机对半（安慰剂）

输出长表，断点续跑，文件名带时间戳。
"""
import glob
import os
import time

import anndata as ad
import numpy as np
import pandas as pd
from alphagenome.data import genome
from alphagenome.models.v2 import dna_client as c2, variant_scorers as vs2

STAMP = time.strftime("%m%d_%H%M%S")
OUT = f"b_data/fairfax_main_{STAMP}.csv"


def build_conditions():
    ctrl = ad.read_h5ad(sorted(glob.glob("b_data/mono_ctrl_*.h5ad"))[-1])
    stim = ad.read_h5ad(sorted(glob.glob("b_data/mono_stim_*.h5ad"))[-1])
    idx = np.arange(ctrl.n_obs)
    rng = np.random.default_rng(11)
    rng.shuffle(idx)
    return {"stim": stim, "ctrl": ctrl,
            "halfA": ctrl[idx[:128]].copy(), "halfB": ctrl[idx[128:]].copy()}


def main():
    g = pd.read_csv("stim_variants_groups.csv")
    f = g[g.study == "Fairfax_2014"].drop_duplicates(["group", "variant"]).copy()
    # 上一轮可见率测试用过的位点打标记
    vis_files = sorted(glob.glob("b_data/visibility_fairfax_*.csv"))
    used = set()
    if vis_files:
        v = pd.concat([pd.read_csv(x) for x in vis_files])
        used = set(v.variant.astype(str))
    f["used_in_visibility_check"] = f.variant.astype(str).isin(used)
    f = f.rename(columns={"condition": "stim"})
    f = f[["variant", "gene", "study", "cell", "group", "stim", "beta", "z", "pip",
           "chrom", "pos", "ref", "alt", "used_in_visibility_check"]]
    print(f"位点 {len(f)}（激活特异 {int((f.group=='act_specific').sum())}，"
          f"静息特异 {int((f.group=='rest_specific').sum())}，"
          f"共享 {int((f.group=='shared').sum())}）")
    print("按刺激：" + str(f.drop_duplicates('variant').stim.value_counts().to_dict()))
    print(f"含上一轮用过 {int(f.used_in_visibility_check.sum())} 个")

    # 断点：按 (variant, cond_role) 精确判断，避免半成品位点被整体跳过
    done = set()
    for f_ in sorted(glob.glob("b_data/fairfax_main_*.csv")):
        d_ = pd.read_csv(f_, usecols=["variant", "cond_role"])
        done |= set(d_.itertuples(index=False, name=None))
    print(f"已完成 {len(done)} 条（位点×条件）", flush=True)
    conds = build_conditions()
    m = c2.create(os.environ["ALPHAGENOME_API_KEY"], timeout=600)
    sc = [vs2.RECOMMENDED_VARIANT_SCORERS["RNA_SEQ"]]

    from concurrent.futures import ThreadPoolExecutor

    def one(r, role):
        cond = conds[role]
        v = genome.Variant(r.chrom, int(r.pos), r.ref, r.alt)
        iv = v.reference_interval.resize(c2.SEQUENCE_LENGTH_1MB)
        d = vs2.tidy_scores(m.score_variant(
            interval=iv, variant=v, variant_scorers=sc,
            cell_conditioning=cond, cell_condition_name=role))
        h = d[d.gene_id == r.gene]
        return dict(variant=r.variant, gene=r.gene, study=r.study, cell=r.cell,
                    group=r.group, stim=r.stim, cond_role=role,
                    pred_raw_score=float(h.raw_score.mean()) if len(h) else np.nan,
                    n_tracks=len(h), exp_beta=r.beta, exp_z=r.z, exp_pip=r.pip,
                    used_in_visibility_check=bool(r.used_in_visibility_check))

    n = 0
    buf = []
    for i, r in enumerate(f.itertuples(index=False), 1):
        todo = [x for x in conds if (r.variant, x) not in done]
        if not todo:
            continue
        with ThreadPoolExecutor(max_workers=4) as ex:
            for rec in ex.map(lambda role: one(r, role), todo):
                buf.append(rec); n += 1
        if len(buf) >= 40:
            pd.DataFrame(buf).to_csv(OUT, mode="a", index=False,
                                     header=not os.path.exists(OUT))
            buf = []
        if i % 25 == 0:
            print(f"  位点 {i}/{len(f)}｜写入 {n} 行", flush=True)
    if buf:
        pd.DataFrame(buf).to_csv(OUT, mode="a", index=False,
                                 header=not os.path.exists(OUT))
    print(f"完成，写入 {n} 行 → {OUT}")


if __name__ == "__main__":
    main()
