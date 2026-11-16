#!/usr/bin/env python3
"""AlphaGenome 前列腺等位基因不平衡小试验。

接口已按 alphagenome 客户端实测校正（v1 客户端，2026-11 版本）。

四步用法：
  0) 自检（10 秒，确认密钥和网络通）
       python ar_pilot.py selftest
  1) 看模型里有哪些前列腺相关的轨道（几秒钟）
       python ar_pilot.py tracks --out prostate_tracks.csv
  2) 给位点打分（可中断，重跑会自动续上）
       python ar_pilot.py score --variants variants.csv --out scores.csv
  3) 和实验结果比较
       python ar_pilot.py eval --variants variants.csv --scores scores.csv --out summary.csv

variants.csv 需要你从论文补充表整理成下面 7 列：
  variant_id   任意唯一编号，例如 rs 号
  chrom        带 chr 前缀，GRCh38 坐标，例如 chr8
  pos          1-based 位置（GRCh38；原始数据若是 hg19 要先做坐标转换）
  ref          参考基因组上的碱基
  alt          另一个等位基因
  effect       实验效应，正数 = alt 等位基因信号更强（例如 alt 读段比例 - 0.5）
  significant  1 = 实验判定有不平衡，0 = 测了但没有不平衡（对照，可以没有）

依赖：pip install -U alphagenome pandas scipy
运行前：export ALPHAGENOME_API_KEY="你的密钥"
    Windows PowerShell:  $env:ALPHAGENOME_API_KEY="..."
    Windows cmd:         set ALPHAGENOME_API_KEY=...
"""
import argparse
import os
import re
import sys
import time

import numpy as np
import pandas as pd

# 前列腺相关 biosample 的名字关键词（大小写不敏感）
KEYWORDS = ["lncap", "22rv1", "c4-2b", "c42b", "vcap", "prostate", "rwpe", "pc-3", "du145"]
KEYWORD_RE = "|".join(re.escape(k) for k in KEYWORDS)

# 与上面关键词对应的 ontology CURIE（biosample_name 缺失时靠它兜底）
PROSTATE_CURIES = {
    "CL:0002231",   # epithelial cell of prostate
    "UBERON:0002367",  # prostate gland
    "EFO:0002071",  # LNCAP
    "EFO:0002074",  # PC-3
    "EFO:0002095",  # 22Rv1
    "EFO:0002323",  # RWPE1
    "EFO:0005726",  # LNCaP clone FGC
    "EFO:0006365",  # C4-2B
    "EFO:0007610",  # RWPE2
    "EFO:0007752",  # VCaP
}

# 只取和调控元件活性直接相关的四类输出
SCORER_KEYS = ["CHIP_TF", "CHIP_HISTONE", "DNASE", "ATAC"]
TRACK_COLS = ["output_type", "biosample_name", "transcription_factor", "histone_mark"]


def is_prostate(df):
    """按 biosample 名字或 ontology CURIE 判定前列腺相关行。"""
    name_hit = df.get("biosample_name", pd.Series("", index=df.index))
    name_hit = name_hit.astype(str).str.lower().str.contains(KEYWORD_RE, regex=True, na=False)
    curie_hit = df.get("ontology_curie", pd.Series("", index=df.index))
    curie_hit = curie_hit.astype(str).isin(PROSTATE_CURIES)
    return name_hit | curie_hit


def get_model(timeout=600):
    """建客户端。v1 是目前稳定版；地址默认 gdmscience.googleapis.com。"""
    try:
        from alphagenome.models.v1 import dna_client
    except ModuleNotFoundError:
        from alphagenome.models import dna_client  # 老版本包结构
    key = os.environ.get("ALPHAGENOME_API_KEY")
    if not key:
        sys.exit("请先设置环境变量 ALPHAGENOME_API_KEY")
    return dna_client.create(key, timeout=timeout), dna_client


def _scorers(keys):
    from alphagenome.models.v1 import variant_scorers
    missing = [k for k in keys if k not in variant_scorers.RECOMMENDED_VARIANT_SCORERS]
    if missing:
        sys.exit(f"这些 scorer 名不存在：{missing}\n可选：{list(variant_scorers.RECOMMENDED_VARIANT_SCORERS)}")
    return [variant_scorers.RECOMMENDED_VARIANT_SCORERS[k] for k in keys]


def cmd_selftest(args):
    """跑一个已知位点，确认密钥有效、网络通、返回结构正常。"""
    from alphagenome.data import genome
    from alphagenome.models.v1 import variant_scorers

    model, dna_client = get_model(timeout=180)
    v = genome.Variant(chromosome="chr8", position=127401060,
                       reference_bases="T", alternate_bases="G")  # rs6983267
    iv = v.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)
    t0 = time.time()
    res = model.score_variant(interval=iv, variant=v,
                              variant_scorers=_scorers(["DNASE", "ATAC"]))
    df = variant_scorers.tidy_scores(res)
    df["variant_id"] = df["variant_id"].astype(str)
    n = int(is_prostate(df).sum())
    print(f"密钥有效。rs6983267 打分用时 {time.time() - t0:.1f}s，"
          f"返回 {len(df)} 行，其中前列腺相关 {n} 行。")
    print("可以开始第 1 步：python ar_pilot.py tracks")


def cmd_tracks(args):
    """列出模型里所有前列腺相关的轨道，确认有没有 AR 的 ChIP 轨道。"""
    model, dna_client = get_model()
    meta = model.output_metadata(dna_client.Organism.HOMO_SAPIENS)
    frames = []
    for name in ["chip_tf", "chip_histone", "dnase", "atac", "cage", "rna_seq", "procap"]:
        df = getattr(meta, name, None)
        if df is None or len(df) == 0:
            continue
        df = df.copy()
        hit = is_prostate(df)
        print(f"  {name:<13} 全部 {len(df):>5} 条，前列腺相关 {int(hit.sum()):>3} 条")
        if hit.any():
            df = df[hit]
            df.insert(0, "output_type", name.upper())
            frames.append(df)
    if not frames:
        sys.exit("没有找到任何前列腺相关轨道，请检查关键词或接口。")

    out = pd.concat(frames, ignore_index=True)
    out.to_csv(args.out, index=False)
    show = [c for c in ["output_type", "biosample_name", "ontology_curie",
                        "transcription_factor", "histone_mark"] if c in out.columns]
    print("\n" + out[show].fillna("").to_string(index=False))

    print(f"\n共 {len(out)} 条前列腺相关轨道，已保存到 {args.out}")
    print("\n=== 关键判定 ===")
    if "transcription_factor" in out.columns:
        tfs = sorted(set(out["transcription_factor"].dropna().astype(str)))
        has_ar = any(t.upper() == "AR" for t in tfs)
        print("前列腺轨道里出现的转录因子：" + (", ".join(tfs) if tfs else "无"))
        print("有没有 AR 的 ChIP 轨道：" + ("有" if has_ar else "没有"))
    lines = []
    for ot, g in out.groupby("output_type"):
        bs = sorted(set(g["biosample_name"].dropna().astype(str)))
        lines.append(f"  {ot:<13} {len(g):>3} 条 / {len(bs)} 个样本：{', '.join(bs)}")
    print("各类输出可用的样本：\n" + "\n".join(lines))
    print("\n结论用法：没有 AR 直接结合轨道，就只能用染色质开放性（DNASE/ATAC）"
          "和组蛋白修饰（尤其 H3K27ac）间接评估；课题说法要改成「调控元件活性」而不是「AR 结合」。")


def _read_variants(path):
    from alphagenome.data import genome
    variants = pd.read_csv(path, dtype={"variant_id": str, "chrom": str})
    need = {"variant_id", "chrom", "pos", "ref", "alt", "effect"}
    lack = need - set(variants.columns)
    if lack:
        sys.exit(f"variants 文件缺这些列：{sorted(lack)}")
    if not variants["chrom"].astype(str).str.startswith("chr").all():
        sys.exit("chrom 列必须带 chr 前缀，例如 chr8")
    variants["pos"] = variants["pos"].astype(int)
    variants["ref"] = variants["ref"].astype(str).str.upper()
    variants["alt"] = variants["alt"].astype(str).str.upper()
    variants["_obj"] = [
        genome.Variant(chromosome=c, position=int(p), reference_bases=r, alternate_bases=a)
        for c, p, r, a in zip(variants["chrom"], variants["pos"],
                              variants["ref"], variants["alt"])]
    return variants


def cmd_score(args):
    """批量调用 AlphaGenome 打分，只保留前列腺相关轨道，按批落盘、可续跑。"""
    from alphagenome.models.v1 import variant_scorers

    model, dna_client = get_model()
    variants = _read_variants(args.variants)
    scorers = _scorers([s.strip().upper() for s in args.scorers.split(",") if s.strip()])
    seq_len = dna_client.SUPPORTED_SEQUENCE_LENGTHS["SEQUENCE_LENGTH_1MB"]

    done = set()
    if os.path.exists(args.out):
        done = set(pd.read_csv(args.out, usecols=["input_id"], dtype=str)["input_id"])
        print(f"已有 {len(done)} 个位点的结果，跳过它们继续。")

    todo = variants[~variants["variant_id"].isin(done)].reset_index(drop=True)
    if todo.empty:
        print("全部位点都已打过分。")
        return

    n_ok, n_fail = 0, 0
    for start in range(0, len(todo), args.batch):
        chunk = todo.iloc[start:start + args.batch]
        objs = list(chunk["_obj"])
        ivs = [v.reference_interval.resize(seq_len) for v in objs]
        res = None
        for attempt in range(4):
            try:
                res = model.score_variants(intervals=ivs, variants=objs,
                                           variant_scorers=scorers,
                                           progress_bar=False,
                                           max_workers=args.workers)
                break
            except Exception as err:  # 限流或网络抖动时等一会再重试
                wait = 15 * 2 ** attempt
                print(f"  批 {start + 1}-{start + len(chunk)} 失败（{err}），{wait}s 后重试")
                time.sleep(wait)
        if res is None:
            print(f"  批 {start + 1}-{start + len(chunk)} 连续失败，跳过")
            n_fail += len(chunk)
            continue

        df = variant_scorers.tidy_scores(res)
        if df is None or df.empty:
            continue
        df["variant_id"] = df["variant_id"].astype(str)
        # 用 Variant 的字符串形式把结果映射回自己的 variant_id
        key2id = {str(v): vid for v, vid in zip(objs, chunk["variant_id"])}
        df["input_id"] = df["variant_id"].map(key2id)
        df = df[df["input_id"].notna()]
        if not args.all_tracks:
            df = df[is_prostate(df)]
        if df.empty:
            print(f"  批 {start + 1}-{start + len(chunk)} 没有可用结果")
            continue
        df.insert(0, "input_id", df.pop("input_id"))
        df.drop(columns=["variant_id"], errors="ignore").to_csv(
            args.out, mode="a", index=False,
            header=not os.path.exists(args.out))
        n_ok += df["input_id"].nunique()
        print(f"  已完成 {min(start + args.batch, len(todo))}/{len(todo)}，"
              f"累计 {n_ok} 个位点有结果")

    print(f"打分结束。成功 {n_ok} 个位点，失败 {n_fail} 个。")


def _eval_block(df, group_cols, min_n):
    from scipy import stats
    rows = []
    for keys, g in df.groupby(group_cols, dropna=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        pos = g[g["significant"] == 1]
        neg = g[g["significant"] == 0]
        use = pos[(pos["effect"] != 0) & (pos["raw_score"] != 0)]
        if len(use) < min_n:
            continue
        agree = int((np.sign(use["raw_score"]) == np.sign(use["effect"])).sum())
        rec = dict(zip(group_cols, keys))
        rec["n"] = len(use)
        rec["方向一致率"] = round(agree / len(use), 3)
        rec["二项检验p"] = round(stats.binomtest(agree, len(use), 0.5).pvalue, 4)
        rho = stats.spearmanr(use["raw_score"], use["effect"]).correlation
        rec["spearman"] = round(float(rho), 3) if np.isfinite(rho) else np.nan
        if len(neg) >= min_n:
            u = stats.mannwhitneyu(pos["raw_score"].abs(),
                                   neg["raw_score"].abs(),
                                   alternative="two-sided").statistic
            rec["auroc"] = round(u / (len(pos) * len(neg)), 3)
        rows.append(rec)
    return rows


def cmd_eval(args):
    """对每条轨道计算：方向一致率、相关系数、区分有无不平衡的 AUROC。"""
    variants = pd.read_csv(args.variants, dtype={"variant_id": str})
    scores = pd.read_csv(args.scores, dtype={"input_id": str}, low_memory=False)
    if "significant" not in variants.columns:
        variants["significant"] = 1
    df = scores.merge(variants[["variant_id", "effect", "significant"]],
                      left_on="input_id", right_on="variant_id",
                      suffixes=("_ag", ""))
    if df.empty:
        sys.exit("scores 和 variants 没有对上的位点，检查 input_id / variant_id。")
    print(f"参与评估：{df['input_id'].nunique()} 个位点 × {len(df)} 行")

    group_cols = [c for c in TRACK_COLS if c in df.columns]
    min_n = args.min_n
    per_track = _eval_block(df, group_cols, min_n)
    if not per_track and min_n > 3:
        print(f"按轨道分组时位点数不足 {min_n}，降到 3 再看看。")
        min_n = 3
        per_track = _eval_block(df, group_cols, min_n)

    pooled = _eval_block(df, ["output_type"], min_n)
    all_rows = per_track + [dict(r, level="按输出类型汇总") for r in pooled]

    if not all_rows:
        sys.exit("没有足够的位点可以评估，请检查输入。")
    out = pd.DataFrame(all_rows).sort_values("方向一致率", ascending=False)
    out.to_csv(args.out, index=False)
    cols = [c for c in list(out.columns) if c != "level"] + ["level"]
    print(out[[c for c in cols if c in out.columns]].head(30).to_string(index=False))
    print(f"\n完整结果已保存到 {args.out}。随机水平的方向一致率是 0.5，AUROC 是 0.5。")
    print("AlphaGenome 官方在 TF binding QTL 上报告的方向相关系数约 0.55，可以拿来做参照。")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("selftest", help="确认密钥和网络可用")
    p.set_defaults(func=cmd_selftest)

    p = sub.add_parser("tracks")
    p.add_argument("--out", default="prostate_tracks.csv")
    p.set_defaults(func=cmd_tracks)

    p = sub.add_parser("score")
    p.add_argument("--variants", required=True)
    p.add_argument("--out", default="scores.csv")
    p.add_argument("--scorers", default=",".join(SCORER_KEYS),
                   help="逗号分隔，可选 CHIP_TF,CHIP_HISTONE,DNASE,ATAC")
    p.add_argument("--batch", type=int, default=8, help="每批位点数")
    p.add_argument("--workers", type=int, default=5, help="批内并发数")
    p.add_argument("--all-tracks", action="store_true", dest="all_tracks",
                   help="保留全部轨道（做阴性对照用，文件会很大）")
    p.set_defaults(func=cmd_score)

    p = sub.add_parser("eval")
    p.add_argument("--variants", required=True)
    p.add_argument("--scores", required=True)
    p.add_argument("--out", default="summary.csv")
    p.add_argument("--min-n", type=int, default=10, dest="min_n")
    p.set_defaults(func=cmd_eval)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
