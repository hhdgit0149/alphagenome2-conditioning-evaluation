#!/usr/bin/env python3
"""导出两代模型的原始清单，供独立核验「2799 个条件零个受处理状态」与「一代轨道零命中」。

输出：
  v1_all_tracks.csv       一代全部轨道（每行一条）+ 是否命中五组关键词
  v1_keyword_summary.csv  一代：每种输出类型的总数与各关键词命中数
  ag2_conditions.csv      二代全部内置条件 + 是否命中处理/激素/前列腺关键词
  ag2_keyword_summary.csv 二代：各关键词命中数
"""
import os
import re

import pandas as pd

KEY = {
    "glucocorticoid": r"glucocort|dexameth|hydrocortis|methylprednis|cortisol|prednisol|steroid",
    "cytokine_or_stimulus": r"ifn\b|interferon|\blps\b|pam3|r848|influenza|cytokine|"
                            r"\btnf\b|\bil-?[0-9]|dex\b|treated|treatment|stimulat",
    "activated_label": r"activated|activation",
    "hormone_treated": r"dihydrotestosterone|\bdht\b|testosteron|estradiol|androgen|tamoxifen|"
                       r"estrogen|progest",
    "prostate": r"prostate|lncap|pc-?3|vcapp?|du145",
}
V1_TYPES = ["chip_tf", "chip_histone", "dnase", "atac", "cage", "rna_seq", "procap"]


def col(d, name):
    return d[name].astype(str) if name in d.columns else pd.Series("", index=d.index)


def main():
    key = os.environ["ALPHAGENOME_API_KEY"]
    summary = []

    # ---------- 一代 ----------
    try:
        from alphagenome.models.v1 import dna_client as v1
    except ModuleNotFoundError:
        from alphagenome.models import dna_client as v1
    m1 = v1.create(key, timeout=600)
    meta = m1.output_metadata(v1.Organism.HOMO_SAPIENS)
    frames = []
    for name in V1_TYPES:
        df = getattr(meta, name, None)
        if df is None or len(df) == 0:
            continue
        d = df.copy()
        d.insert(0, "output_type", name.upper())
        txt = (col(d, "biosample_name") + " | " + col(d, "ontology_curie") + " | "
               + col(d, "transcription_factor") + " | " + col(d, "histone_mark") + " | "
               + col(d, "Assay title") + " | " + col(d, "assay_title"))
        d["_txt"] = txt
        for k, pat in KEY.items():
            d[f"hit_{k}"] = txt.str.contains(re.compile(pat, re.I), na=False)
        d["hit_any_treated"] = d[[f"hit_{k}" for k in ("glucocorticoid", "hormone_treated",
                                                       "cytokine_or_stimulus")]].any(axis=1)
        frames.append(d)
        summary.append(dict(model="v1", output_type=name.upper(), total=len(d),
                            **{k: int(d[f"hit_{k}"].sum()) for k in KEY},
                            any_treated=int(d["hit_any_treated"].sum())))
        print(f"v1 {name:13s} 全部 {len(d):>5}  处理类命中 {int(d['hit_any_treated'].sum())}", flush=True)
    v1all = pd.concat(frames, ignore_index=True)
    v1all.drop(columns=["_txt"]).to_csv("v1_all_tracks.csv", index=False)
    print(f"\n写出 v1_all_tracks.csv：{len(v1all)} 行，处理类命中合计 {int(v1all['hit_any_treated'].sum())}")

    # ---------- 二代 ----------
    from alphagenome.models.v2 import dna_client as v2
    m2 = v2.create(key, timeout=600)
    c = m2.supported_conditions()
    c = c.copy()
    txt2 = (col(c, "biosample_name") + " | " + col(c, "biosample_type") + " | "
            + col(c, "condition") + " | " + col(c, "ontology_curie"))
    for k, pat in KEY.items():
        c[f"hit_{k}"] = txt2.str.contains(re.compile(pat, re.I), na=False)
    c["hit_any_treated"] = c[[f"hit_{k}" for k in ("glucocorticoid", "hormone_treated",
                                                  "cytokine_or_stimulus")]].any(axis=1)
    c.to_csv("ag2_conditions.csv", index=False)
    print(f"写出 ag2_conditions.csv：{len(c)} 个条件，处理类命中 {int(c['hit_any_treated'].sum())}")
    for k in KEY:
        print(f"  {k:20s} {int(c[f'hit_{k}'].sum())}")
    summary.append(dict(model="v2", output_type="conditions", total=len(c),
                        **{k: int(c[f"hit_{k}"].sum()) for k in KEY},
                        any_treated=int(c["hit_any_treated"].sum())))
    pd.DataFrame(summary).to_csv("keyword_summary.csv", index=False)
    print("\n写出 keyword_summary.csv")


if __name__ == "__main__":
    main()
