#!/usr/bin/env python3
"""导出逐位点长表（每行 = 一个「位点 × 条件」），给外部复算用。

列：
  variant, gene, gene_name, study, cell, group
  cond_role(act/rest), cond_curie, cond_label
  pred_raw_score           模型预测的对数比（ALT vs REF，GeneMaskLFC）
  pred_ref_expr            预测的参考等位基因表达量 —— 需额外 predict_variant 调用，暂缺
  exp_beta, exp_z, exp_pip 实验侧在该条件下的效应与统计量（可信集里取该变异）
  in_cs                    该 (基因,变异) 是否出现在该条件的可信集
"""
import os

import pandas as pd

META = "eqtl_cache/dataset_metadata_r8.tsv"
SC = "scores_tcell.csv"
OUT = "tcell_per_variant_long.csv"
CONDS = {"act": ("UBERON:0002116_CL:0000896", "ileum activated CD4+ T cell"),
         "rest": ("UBERON:0002116_CL:0000624", "ileum CD4+ T cell")}
BASE = {"naive"}


def main():
    meta = pd.read_csv(META, sep="\t")
    d = pd.read_csv(SC)
    g = pd.read_csv("stim_variants_groups.csv")[["variant", "study", "condition"]]
    d = d.merge(g.drop_duplicates(["variant", "study"]), on=["variant", "study"], how="left")
    print(f"评分位点 {len(d)}，其中带条件标签 {d.condition.notna().sum()}")

    # 建立 (study, cell, cond) -> dataset_id
    ds = {}
    for (s, c), x in meta.groupby(["study_label", "tissue_label"]):
        for r in x.itertuples(index=False):
            ds[(s, c, str(r.condition_label))] = r.dataset_id

    cache = {}

    def load(did):
        if did not in cache:
            p = f"eqtl_cache/{did}.credible_sets.tsv.gz"
            cache[did] = (pd.read_csv(p, sep="\t")
                          if os.path.exists(p) and os.path.getsize(p) > 200 else None)
        return cache[did]

    rows = []
    for r in d.itertuples(index=False):
        # 该位点所属刺激条件（分组时的 condition 列）
        stim_cond = r.condition if isinstance(r.condition, str) else None
        for role, (curie, label) in CONDS.items():
            did = ds.get((r.study, r.cell,
                          stim_cond if role == "act" else "naive"))
            beta = z = pip = None
            in_cs = False
            if did:
                cs = load(did)
                if cs is not None:
                    h = cs[(cs.gene_id == r.gene) & (cs.variant == r.variant)]
                    if len(h):
                        in_cs = True
                        j = h.sort_values("pip", ascending=False).iloc[0]
                        beta, z, pip = float(j.beta), float(j.z), float(j.pip)
            rows.append(dict(
                variant=r.variant, gene=r.gene, study=r.study, cell=r.cell,
                group=r.group, cond_role=role, cond_curie=curie, cond_label=label,
                pred_raw_score=(r.pred_act if role == "act" else r.pred_rest),
                pred_ref_expr=None,
                n_tracks=(r.n_act if role == "act" else r.n_rest),
                exp_beta=beta, exp_z=z, exp_pip=pip, in_cs=in_cs))

    out = pd.DataFrame(rows)
    out.to_csv(OUT, index=False)
    print(f"写出 {OUT}：{len(out)} 行（{out.variant.nunique()} 个位点 × 2 条件）")
    print("\n各组的可用性：")
    print(out.groupby(["group", "cond_role"]).agg(
        n=("variant", "size"),
        有过实验效应=("exp_beta", lambda s: s.notna().sum()),
        模型效应中位=("pred_raw_score", lambda s: round(s.abs().median(), 5))
    ).to_string())


if __name__ == "__main__":
    main()
