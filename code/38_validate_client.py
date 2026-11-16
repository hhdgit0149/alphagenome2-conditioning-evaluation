#!/usr/bin/env python3
"""验证手写 tabix 客户端（PREREG_v2 补充C §4），并收集匹配度所需表达中位数。

必须在 35 号脚本完成后运行。产出：
  b_data/strict_z_own_{study}.csv   own-condition 按性状逐行（z 验证 + 匹配度数据）
  b_data/cs_join_{study}.csv        标签行 × 本条件可信集（trait/z/beta，36 号脚本复用）
  report_client_validation.md       验证报告（全量比对 + 未检出复核 + 各数据集缺失率）

退出码非零 = 验证不通过，禁止构建严格标签。
"""
import os
import sys
import time
from importlib.machinery import SourceFileLoader

import numpy as np
import pandas as pd

m35 = SourceFileLoader("m35", "35_strict_z_lookup.py").load_module()

SEED = 20261114
SHARED_N = 150
NOTFOUND_N = 20
TOL = 0.01


def read_cs(qtd):
    df = pd.read_csv(f"eqtl_cache/{qtd}.credible_sets.tsv.gz", sep="\t")
    for c in ("molecular_trait_id", "variant", "gene_id", "pip", "z", "beta"):
        assert c in df.columns, f"{qtd} CS 缺列 {c}"
    return df


def top_cs(df, variant, gene):
    d = df[(df.variant == variant) & (df.gene_id == gene)]
    assert not d.empty, f"标签行不在可信集: {variant}/{gene}"
    return d.sort_values("pip", ascending=False).iloc[0]


def run_pairs(cfg, need, out):
    """按 (own_cond, variant, gene) 查询，按性状逐行写 out（断点续跑，单写者）。"""
    done = set()
    if os.path.exists(out):
        prev = m35.clean_resume(out)
        done = set(map(tuple, prev[["other_cond", "variant", "gene"]]
                       .drop_duplicates().values))
    tasks = [v for k, v in need.items() if k not in done]
    print(f"  own-condition 需查 {len(need)}，已完成 {len(done)}，待查 {len(tasks)}",
          flush=True)
    if not tasks:
        return
    by_cond = {}
    for t in tasks:
        by_cond.setdefault(t[0], []).append(t)
    header = not os.path.exists(out)
    t0, n = time.time(), 0
    for cond, items in by_cond.items():
        items.sort(key=lambda x: (x[3], x[4]))
        tab = m35.RemoteTabix(cfg["qts"], cfg["cond"][cond])
        for own_cond, variant, gene, chrom, pos in items:
            cname = chrom[3:] if chrom.startswith("chr") else chrom
            pos_rows = None
            for att in range(2):
                try:
                    pos_rows = tab.query_point(cname, pos)
                    break
                except Exception as e:
                    print(f"  {variant}/{cond} 第{att+1}次失败: {e}", flush=True)
                    time.sleep(15)
            if pos_rows is None:
                rec = [dict(other_cond=own_cond, variant=variant, gene=gene,
                            chrom=chrom, pos=pos, trait="", beta="", se="",
                            tpm="", n_pair_rows=-1, n_pos_rows=-1,
                            variant_tested="", row_id=0)]
            else:
                tested = any(r["variant"] == variant for r in pos_rows)
                pair = [r for r in pos_rows if r["variant"] == variant
                        and r["gene_id"] == gene]
                seen, uniq = set(), []
                for r in pair:
                    k2 = (r["trait"], r["beta"], r["se"])
                    if k2 not in seen:
                        seen.add(k2)
                        uniq.append(r)
                if uniq:
                    rec = [dict(other_cond=own_cond, variant=variant, gene=gene,
                                chrom=chrom, pos=pos, trait=r["trait"],
                                beta=r["beta"], se=r["se"], tpm=r["tpm"],
                                n_pair_rows=len(uniq), n_pos_rows=len(pos_rows),
                                variant_tested=tested, row_id=i)
                           for i, r in enumerate(sorted(
                               uniq, key=lambda x: (x["trait"], x["beta"])))]
                else:
                    rec = [dict(other_cond=own_cond, variant=variant, gene=gene,
                                chrom=chrom, pos=pos, trait="", beta="", se="",
                                tpm="", n_pair_rows=0, n_pos_rows=len(pos_rows),
                                variant_tested=tested, row_id=0)]
            pd.DataFrame(rec).to_csv(out, mode="a", index=False, header=header)
            header = False
            n += 1
            if n % 50 == 0:
                el = time.time() - t0
                print(f"  {n}/{len(tasks)}｜{n/max(el,1)*60:.0f} 查/分", flush=True)


def main():
    lines = ["# tabix 客户端验证报告（PREREG_v2 补充C §4）",
             f"日期 {time.strftime('%Y-%m-%d %H:%M')}", ""]
    fail = False
    for key, cfg in m35.STUDIES.items():
        lines.append(f"## {key}")
        lab = pd.read_csv(cfg["labels"])
        cs = {c: read_cs(q) for c, q in cfg["cond"].items()}

        # 1. 可信集连接（本条件 CS：act/shared 用刺激侧，rest/shared 用静息侧）
        join_rows = []
        for r in lab.itertuples(index=False):
            row = dict(stim=r.stim, group=r.group, variant=r.variant, gene=r.gene,
                       study=r.study, cell=r.cell, chrom=r.chrom, pos=r.pos,
                       beta=r.beta, z=r.z, pip=r.pip)
            if r.group in ("act_specific", "shared"):
                t = top_cs(cs[r.stim], r.variant, r.gene)
                row.update(cs_stim_trait=t.molecular_trait_id, cs_stim_z=float(t.z),
                           cs_stim_beta=float(t.beta))
            if r.group in ("rest_specific", "shared"):
                t = top_cs(cs["naive"], r.variant, r.gene)
                row.update(cs_naive_trait=t.molecular_trait_id, cs_naive_z=float(t.z),
                           cs_naive_beta=float(t.beta))
            join_rows.append(row)
        join = pd.DataFrame(join_rows)
        assert not join.duplicated(["stim", "group", "variant", "gene"]).any()
        join.to_csv(f"b_data/cs_join_{key}.csv", index=False)
        lines.append(f"- 可信集连接：{len(join)} 行全部命中（act/shared→刺激侧，"
                     f"rest/shared→静息侧，同对多可信集取 pip 最高）")

        # 2. own-condition 查询任务：act+rest 全量；shared 抽样 150 行×两侧
        need = {}
        for r in join.itertuples(index=False):
            if r.group == "act_specific":
                need[(r.stim, r.variant, r.gene)] = (r.stim, r.variant, r.gene,
                                                     r.chrom, int(r.pos))
            elif r.group == "rest_specific":
                need[("naive", r.variant, r.gene)] = ("naive", r.variant, r.gene,
                                                      r.chrom, int(r.pos))
        sh = join[join.group == "shared"].reset_index(drop=True)
        rng = np.random.default_rng(SEED)
        pick = sh.iloc[sorted(rng.choice(len(sh), min(SHARED_N, len(sh)),
                                          replace=False))]
        for r in pick.itertuples(index=False):
            need[("naive", r.variant, r.gene)] = ("naive", r.variant, r.gene,
                                                  r.chrom, int(r.pos))
            need[(r.stim, r.variant, r.gene)] = (r.stim, r.variant, r.gene,
                                                  r.chrom, int(r.pos))
        own_out = f"b_data/strict_z_own_{key}.csv"
        run_pairs(cfg, need, own_out)
        own = pd.read_csv(own_out)
        assert not (own.n_pair_rows < 0).any(), f"{key} own 查询有失败行，先补齐"

        # 3. z 逐一比对（own-condition，CS 探针匹配）
        checks, bad = 0, []
        for r in join.itertuples(index=False):
            todo = []
            if r.group == "act_specific":
                todo = [(r.stim, r.cs_stim_trait, r.cs_stim_z)]
            elif r.group == "rest_specific":
                todo = [("naive", r.cs_naive_trait, r.cs_naive_z)]
            else:                                   # shared：仅抽样行比对两侧
                if ((pick.variant == r.variant) & (pick.gene == r.gene)).any():
                    todo = [(r.stim, r.cs_stim_trait, r.cs_stim_z),
                            ("naive", r.cs_naive_trait, r.cs_naive_z)]
            for cond, trait, cs_z in todo:
                hit = own[(own.other_cond == cond) & (own.variant == r.variant)
                          & (own.gene == r.gene) & (own.trait == trait)]
                if hit.empty:
                    bad.append((r.variant, r.gene, r.stim, r.group, cond, trait,
                                "own行缺失", np.nan, cs_z))
                    continue
                zq = float(hit.iloc[0].beta) / float(hit.iloc[0].se)
                if abs(zq - cs_z) > TOL or np.sign(zq) != np.sign(cs_z):
                    bad.append((r.variant, r.gene, r.stim, r.group, cond, trait,
                                "数值不符", zq, cs_z))
                checks += 1
        lines.append(f"- z 全量比对：{checks} 对，容差 ±{TOL}，"
                     f"不符 {len(bad)} 对"
                     + (f"，最大偏差 {max((abs(b[7]-b[8]) for b in bad), default=0):.4f}"
                        if bad else ""))
        for b in bad[:20]:
            lines.append(f"  - 不符：{b[0]}/{b[1]} {b[2]} {b[3]} {b[4]} "
                         f"trait={b[5]} {b[6]} z查询={b[7]:.4f} z可信集={b[8]:.4f}")
        if bad:
            fail = True

        # 4. 未检出复核：35 号输出 n_pair_rows==0 的行抽 20 个扩窗重查
        zother = pd.read_csv(cfg["out"])
        assert not (zother.n_pair_rows < 0).any(), f"{key} 35 号查询有失败行，先补齐"
        zerohit = zother[zother.n_pair_rows == 0]
        rng2 = np.random.default_rng(SEED)
        if len(zerohit):
            sample = zerohit.iloc[sorted(rng2.choice(len(zerohit),
                                                     min(NOTFOUND_N, len(zerohit)),
                                                     replace=False))]
            n_found_window = 0
            for r in sample.itertuples(index=False):
                tab = m35.RemoteTabix(cfg["qts"], cfg["cond"][r.other_cond])
                cname = r.chrom[3:] if str(r.chrom).startswith("chr") else str(r.chrom)
                if m35.query_window(tab, cname, int(r.pos), r.variant):
                    n_found_window += 1
                    bad.append((r.variant, r.gene, "-", "-", r.other_cond, "-",
                                "扩窗复核发现该变异存在（程序漏读）", np.nan, np.nan))
            lines.append(f"- 未检出复核：{len(sample)} 行扩窗 ±48kb 重查，"
                         f"发现存在 {n_found_window}（应为 0）")
            if n_found_window:
                fail = True
        else:
            lines.append("- 未检出复核：无 zero-hit 行")

        # 5. 各数据集缺失率（35 号输出，按 other_cond）
        per = zother.groupby("other_cond").agg(
            查询对数=("variant", "count"),
            zero_hit=("n_pair_rows", lambda v: (v == 0).sum()),
            variant_untested=("variant_tested", lambda v: (v == False).sum()))
        lines.append("- 另一条件查询缺失率：")
        lines.append("```\n" + per.to_string() + "\n```")
        lines.append("")

    # 抽样错误率上界（零失败，Clopper-Pearson 单侧 95%）
    for n, what in [(SHARED_N, "shared 抽样"), (NOTFOUND_N, "未检出抽样")]:
        ub = 1 - 0.05 ** (1.0 / n)
        lines.append(f"- {what} n={n}，零失败对应的 95% 错误率上界 {ub:.1%}")
    lines.append("  说明：shared 抽样足以支撑「代码路径正确」；未检出抽样仅点检，"
                 "「确实不存在」的论断由细则一的分类计数承担，不依赖抽样。")
    lines.append("")

    lines.append("## 结论")
    lines.append("**验证不通过，禁止构建严格标签，先查清不符项。**" if fail
                 else "**验证通过：z 全量比对一致，未检出复核未发现漏读。**")
    open("report_client_validation.md", "w", encoding="utf-8").write(
        "\n".join(lines) + "\n")
    print("\n".join(lines[-2:]))
    sys.exit(1 if fail else 0)


if __name__ == "__main__":
    main()
