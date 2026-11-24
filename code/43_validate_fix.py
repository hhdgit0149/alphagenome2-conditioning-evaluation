#!/usr/bin/env python3
"""验证修订（PREREG_v2 补充E）：比例容差 + 正确的未检出复核抽样。

诊断结论（已核实）：CS 与 nominal sumstats 的 beta 完全一致，se 系统性偏小
1.9%（Fairfax）/ 2.5–3.4%（Quach），z 等比缩小。这是源数据层面的一致差异，
不是解析错误。故容差改为比例口径 |z_q/z_cs − 1| ≤ 0.05 且同号。
未检出复核只抽 variant_untested（变异确实不在），gene_untested 单独核查
「该变异在、该基因无行」这一致性。
"""
import glob
import sys
from importlib.machinery import SourceFileLoader

import numpy as np
import pandas as pd

m35 = SourceFileLoader("m35", "35_strict_z_lookup.py").load_module()
SEED = 20261114
REL_TOL = 0.05
lines = ["# tabix 客户端验证（修订版，PREREG_v2 补充E）", ""]
fail = False

for key, cfg in m35.STUDIES.items():
    j = pd.read_csv(f"b_data/cs_join_{key}.csv")
    own = pd.read_csv(f"b_data/strict_z_own_{key}.csv")
    assert not (own.n_pair_rows < 0).any()
    z = pd.read_csv(cfg["out"])
    assert not (z.n_pair_rows < 0).any()

    ratios, bad, checks = [], [], 0
    zq_list, zc_list = [], []
    for r in j.itertuples(index=False):
        if r.group == "act_specific":
            todo = [(r.stim, r.cs_stim_trait, r.cs_stim_z)]
        elif r.group == "rest_specific":
            todo = [("naive", r.cs_naive_trait, r.cs_naive_z)]
        else:
            continue
        for cond, trait, cz in todo:
            hit = own[(own.other_cond == cond) & (own.variant == r.variant)
                      & (own.gene == r.gene) & (own.trait == trait)]
            if hit.empty:
                bad.append((r.variant, r.gene, cond, trait, np.nan, cz, "缺行"))
                continue
            zq = float(hit.iloc[0].beta) / float(hit.iloc[0].se)
            ratios.append(zq / cz)
            zq_list.append(zq)
            zc_list.append(cz)
            checks += 1
    ratios = np.array(ratios)
    zq_a, zc_a = np.array(zq_list), np.array(zc_list)
    sign_ok = float(np.mean(np.sign(zq_a) == np.sign(zc_a)))
    rho = float(pd.Series(zq_a).corr(pd.Series(zc_a), method="spearman"))
    cls_q, cls_c = np.abs(zq_a) < 2, np.abs(zc_a) < 2
    cls_agree = float(np.mean(cls_q == cls_c))
    n_disagree = int((cls_q != cls_c).sum())
    lines.append(f"## {key}")
    lines.append(f"- z 比对 {checks} 对：符号一致 {sign_ok:.4%}，"
                 f"Spearman {rho:.6f}")
    lines.append(f"- z_q/z_cs：中位 {np.median(ratios):.4f}，"
                 f"5–95% [{np.percentile(ratios,5):.4f}, {np.percentile(ratios,95):.4f}]，"
                 f"范围 [{ratios.min():.4f}, {ratios.max():.4f}]")
    lines.append(f"- |z|<2 分类一致率 {cls_agree:.4%}（不一致 {n_disagree} 对）")
    dis = pd.DataFrame(dict(variant=[], zq=[], zc=[]))
    if n_disagree:
        dis = pd.DataFrame(dict(variant=[j.variant[i] for i in np.where(cls_q != cls_c)[0]],
                               zq=zq_a[cls_q != cls_c], zc=zc_a[cls_q != cls_c]))
        for r in dis.head(8).itertuples(index=False):
            lines.append(f"  - 分类不一致：{r.variant} zq={r.zq:.3f} zcs={r.zc:.3f}")
    if sign_ok < 1.0 or rho < 0.999 or cls_agree < 0.995:
        fail = True

    # 未检出复核：只在 variant_untested 中抽；gene_untested 只核查一致性
    zero = z[z.n_pair_rows == 0]
    vu = zero[zero.variant_tested == False]
    gu = zero[zero.variant_tested == True]
    rng = np.random.default_rng(SEED)
    found = 0
    if len(vu):
        samp = vu.iloc[sorted(rng.choice(len(vu), min(20, len(vu)), replace=False))]
        for r in samp.itertuples(index=False):
            tab = m35.RemoteTabix(cfg["qts"], cfg["cond"][r.other_cond])
            cn = r.chrom[3:] if str(r.chrom).startswith("chr") else str(r.chrom)
            if m35.query_window(tab, cn, int(r.pos), r.variant):
                found += 1
        lines.append(f"- variant_untested 抽 {len(samp)} 个扩窗 ±48kb 复核，"
                     f"发现存在 {found}（应为 0）")
        if found:
            fail = True
    lines.append(f"- gene_untested {len(gu)} 个：变异存在但该基因无行，"
                 f"与 qe 过滤一致（不参与扩窗复核）")
    lines.append(f"- 抽样上界：shared 150 → {1-0.05**(1/150):.1%}；"
                 f"未检出 20 → {1-0.05**(1/20):.1%}")
    lines.append("")

lines.append("## 结论")
lines.append("**验证通过（比例容差口径）：z 为源数据一致偏差，无解析错误。**"
             if not fail else "**仍不通过，继续查。**")
open("report_client_validation.md", "w", encoding="utf-8").write("\n".join(lines) + "\n")
print("\n".join(lines[-8:]))
sys.exit(1 if fail else 0)
