#!/usr/bin/env python3
"""一代（v1）模型轨道清单盘点（为对照臂 2/3 的可行性判定）。

查三类：
  1) 单核细胞/CD14 轨道
  2) 任何"受刺激状态"轨道（IFN/LPS/Pam3/流感/IL/刺激/激活等关键词）
  3) 前列腺/LNCaP 轨道
输出 v1_track_inventory.csv + 控制台结论。
"""
import os
import re
import sys

import pandas as pd

KEY_MONO = re.compile(r"monocyte|cd14|macrophage|pbmc|myeloid", re.I)
KEY_STIM = re.compile(r"ifn|interferon|lps|pam3|r848|influenza|stimulat|"
                      r"activated|activation|il-?1|il-?4|il-?6|tnf|dengue|"
                      r"vaccin|treated|treat", re.I)
KEY_PROST = re.compile(r"prostate|lncap|pc-?3|vcapp?|du145|bph", re.I)
TYPES = ["chip_tf", "chip_histone", "dnase", "atac", "cage", "rna_seq", "procap"]


def main():
    try:
        from alphagenome.models.v1 import dna_client
    except ModuleNotFoundError:
        from alphagenome.models import dna_client
    key = os.environ["ALPHAGENOME_API_KEY"]
    m = dna_client.create(key, timeout=600)
    meta = m.output_metadata(dna_client.Organism.HOMO_SAPIENS)
    print("organism metadata:", type(meta).__name__, flush=True)
    print("可用字段:", [a for a in dir(meta) if not a.startswith("_")][:30])

    rows = []
    for name in TYPES:
        df = getattr(meta, name, None)
        if df is None or len(df) == 0:
            print(f"  {name}: 无")
            continue
        txt = (df.get("biosample_name", pd.Series("", index=df.index)).astype(str)
               + " | " + df.get("ontology_curie", pd.Series("", index=df.index)).astype(str)
               + " | " + df.get("transcription_factor", pd.Series("", index=df.index)).astype(str)
               + " | " + df.get("histone_mark", pd.Series("", index=df.index)).astype(str))
        for tag, pat in [("mono", KEY_MONO), ("stim", KEY_STIM), ("prost", KEY_PROST)]:
            if tag == "stim":
                # 刺激关键词只用于找"受刺激状态"，不用于找处理过的细胞系
                hit = txt.str.contains(KEY_STIM, na=False)
            else:
                hit = txt.str.contains(pat, na=False)
            sub = df[hit]
            print(f"  {name:<13} 全部 {len(df):>5}  {tag} {len(sub):>4}", flush=True)
            if len(sub):
                s = sub.copy()
                s.insert(0, "output_type", name.upper())
                s.insert(1, "match_category", tag)
                rows.append(s)

    out = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    out.to_csv("v1_track_inventory.csv", index=False)
    print(f"\n写出 v1_track_inventory.csv：{len(out)} 行")

    # 结论性核查
    def samples(tag, ot=None):
        d = out[out.match_category == tag]
        if ot:
            d = d[d.output_type == ot]
        return sorted(set(d.biosample_name.astype(str))) if len(d) else []

    mono_all = samples("mono")
    stim_all = samples("stim")
    prost_all = samples("prost")
    print("\n=== 结论 ===")
    print(f"单核/巨噬/PBMC 样本：{len(mono_all)} → {mono_all[:15]}")
    print(f"受'刺激/激活'命中的样本：{len(stim_all)} → {stim_all[:15]}")
    print(f"前列腺/LNCaP 样本：{len(prost_all)} → {prost_all[:15]}")
    if not stim_all:
        print("判定：v1 没有任何受刺激状态的轨道 → 无法在原理上表达刺激依赖性"
              "（可作为文章的结构性陈述）。")
    else:
        print("判定：存在刺激相关样本，需逐条核对是否真为'同一细胞类型的刺激前后'。")


if __name__ == "__main__":
    main()
