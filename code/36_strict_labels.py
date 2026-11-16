#!/usr/bin/env python3
"""按 PREREG_v2 补充C 构建严格标签（前置：38 号验证已通过）。

规则：
  act_specific  其他条件=naive：z 取 (variant, 本条件刺激侧 CS 探针) 精确匹配行
  rest_specific 其他条件=该刺激：z 取 (variant, 静息侧 CS 探针) 精确匹配行
  未检出三分并计数：variant_untested / gene_untested / probe_untested
  |z|>=2 → excluded_notstrict
  shared：两条件 CS（各取 pip 最高）beta 符号需相同，不同 → excluded_discordant
输出 b_data/{study}_labels_strict.csv，列 group_strict + reason。
"""
import sys

import pandas as pd

assert "验证通过" in open("report_client_validation.md", encoding="utf-8").read(), \
    "38 号验证未通过（或未运行），禁止构建严格标签"

JOBS = {
    "fairfax": ("b_data/cs_join_fairfax.csv", "b_data/strict_z_fairfax.csv",
                "b_data/fairfax_labels_strict.csv"),
    "quach": ("b_data/cs_join_quach.csv", "b_data/strict_z_quach.csv",
              "b_data/quach_labels_strict.csv"),
}


def pair_rows(z, cond, variant, gene):
    d = z[(z.other_cond == cond) & (z.variant == variant) & (z.gene == gene)]
    assert len(d) > 0, f"35 号输出缺行：{cond}/{variant}/{gene}"
    return d


for key, (joinp, zp, outp) in JOBS.items():
    join = pd.read_csv(joinp)
    z = pd.read_csv(zp)
    assert not (z.n_pair_rows < 0).any(), f"[{key}] 35 号输出有失败行"
    rows = []
    for r in join.itertuples(index=False):
        d = dict(stim=r.stim, group=r.group, study=r.study, cell=r.cell,
                 variant=r.variant, gene=r.gene, beta=r.beta, z=r.z, pip=r.pip,
                 chrom=r.chrom, pos=r.pos, z_other="", reason="")
        if r.group == "shared":
            if (r.cs_stim_beta > 0) != (r.cs_naive_beta > 0):
                d.update(group_strict="excluded_discordant",
                         reason="两条件可信集效应方向相反")
            else:
                d.update(group_strict="shared")
            rows.append(d)
            continue
        other = "naive" if r.group == "act_specific" else r.stim
        cs_trait = r.cs_stim_trait if r.group == "act_specific" else r.cs_naive_trait
        d["other_cond"] = other
        pr = pair_rows(z, other, r.variant, r.gene)
        if int(pr.iloc[0].n_pair_rows) == 0:
            if not bool(pr.iloc[0].variant_tested):
                d.update(group_strict="excluded_noverify",
                         reason="variant_untested（另一条件未测此变异）")
            else:
                d.update(group_strict="excluded_noverify",
                         reason="gene_untested（变异在但基因无行）")
        else:
            hit = pr[pr.trait == cs_trait]
            if hit.empty:
                d.update(group_strict="excluded_noverify",
                         reason="probe_untested（基因有行但 CS 探针不在）")
            else:
                zv = float(hit.iloc[0].beta) / float(hit.iloc[0].se)
                d["z_other"] = zv
                if abs(zv) < 2:
                    d.update(group_strict=r.group)
                else:
                    d.update(group_strict="excluded_notstrict",
                             reason=f"|z|={abs(zv):.2f}>=2")
        rows.append(d)

    out = pd.DataFrame(rows)
    assert not out.duplicated(["stim", "group", "variant", "gene"]).any()
    out.to_csv(outp, index=False)
    print(f"[{key}] {outp}：{len(out)} 行")
    print(pd.crosstab(out.stim, out.group_strict).to_string())
    print(out.reason.value_counts().to_string(), "\n")
