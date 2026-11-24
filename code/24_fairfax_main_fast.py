#!/usr/bin/env python3
"""Fairfax 主测试（提速版）：单一线程池 8 并发，任务按 (variant, role) 扁平化。

原版每个位点新建一个 4 线程池，位点之间串行，实测只有 40 行/分钟；
单独测 4 并发可达 102 次/分钟、8 并发 165 次/分钟、12 并发因限流掉到 43。
所以固定 8 并发、任务扁平化、主线程单写。
"""
import glob
import os
import time
from concurrent.futures import ThreadPoolExecutor

import anndata as ad
import numpy as np
import pandas as pd
from alphagenome.data import genome
from alphagenome.models.v2 import dna_client as c2, variant_scorers as vs2

STAMP = time.strftime("%m%d_%H%M%S")
OUT = f"b_data/fairfax_main_{STAMP}.csv"
WORKERS = 8


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
    vis = sorted(glob.glob("b_data/visibility_fairfax_*.csv"))
    used = set(pd.concat([pd.read_csv(x) for x in vis]).variant.astype(str)) if vis else set()
    f["used_vis"] = f.variant.astype(str).isin(used)
    f = f.rename(columns={"condition": "stim"})

    done = set()
    for fp in sorted(glob.glob("b_data/fairfax_main_*.csv")):
        d = pd.read_csv(fp, usecols=["variant", "cond_role"])
        done |= set(d.itertuples(index=False, name=None))
    print(f"已完成 {len(done)} 条，位点 {len(f)}", flush=True)

    conds = build_conditions()
    m = c2.create(os.environ["ALPHAGENOME_API_KEY"], timeout=600)
    sc = [vs2.RECOMMENDED_VARIANT_SCORERS["RNA_SEQ"]]

    tasks = []
    for r in f.itertuples(index=False):
        for role in conds:
            if (r.variant, role) not in done:
                tasks.append((r, role))
    print(f"待跑 {len(tasks)} 条，{WORKERS} 并发", flush=True)

    def one(task):
        r, role = task
        try:
            v = genome.Variant(r.chrom, int(r.pos), r.ref, r.alt)
            iv = v.reference_interval.resize(c2.SEQUENCE_LENGTH_1MB)
            d = vs2.tidy_scores(m.score_variant(
                interval=iv, variant=v, variant_scorers=sc,
                cell_conditioning=conds[role], cell_condition_name=role))
            h = d[d.gene_id == r.gene]
            return dict(variant=r.variant, gene=r.gene, study=r.study, cell=r.cell,
                        group=r.group, stim=r.stim, cond_role=role,
                        pred_raw_score=float(h.raw_score.mean()) if len(h) else np.nan,
                        n_tracks=len(h), exp_beta=r.beta, exp_z=r.z, exp_pip=r.pip,
                        used_in_visibility_check=bool(r.used_vis))
        except Exception as e:
            print(f"  {r.variant}/{role} {type(e).__name__}", flush=True)
            return None

    buf, n, t0 = [], 0, time.time()
    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        for rec in ex.map(one, tasks):
            if rec:
                buf.append(rec); n += 1
            if len(buf) >= 50:
                pd.DataFrame(buf).to_csv(OUT, mode="a", index=False,
                                         header=not os.path.exists(OUT))
                buf = []
                el = time.time() - t0
                print(f"  {n}/{len(tasks)}｜{n/max(el,1)*60:.0f} 条/分钟", flush=True)
    if buf:
        pd.DataFrame(buf).to_csv(OUT, mode="a", index=False,
                                 header=not os.path.exists(OUT))
    print(f"完成 {n} 条 → {OUT}")


if __name__ == "__main__":
    main()
