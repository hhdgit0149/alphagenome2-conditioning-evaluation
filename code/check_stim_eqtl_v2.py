#!/usr/bin/env python3
"""可行性检查（修正版）：流感激刺激条件下有多少高置信 eQTL 可用。

原脚本 check_stim_eqtl.py 的三处问题已修：
  1. 元数据地址：现为 dataset_metadata_r8.tsv（原文件名已不存在，会 404）
  2. FTP 路径：susie/{study_id}/{dataset_id}/{dataset_id}.credible_sets.tsv.gz
     （原脚本把 study_label 当成目录名，但目录是 QTS000036 这种编号）
  3. 配对条件：condition_label 实际取值是 naive / Influenza_6h

用法：
    python check_stim_eqtl_v2.py
    python check_stim_eqtl_v2.py --studies Randolph_2021 Randolph_2021_reannotated
"""
import argparse
import os
import sys

import pandas as pd
import requests

META_URL = ("https://raw.githubusercontent.com/eQTL-Catalogue/"
            "eQTL-Catalogue-resources/master/data_tables/dataset_metadata_{ver}.tsv")
CS_URL = ("https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/susie/"
          "{study_id}/{dataset_id}/{dataset_id}.credible_sets.tsv.gz")


def need(df, cols, what):
    missing = [c for c in cols if c not in df.columns]
    if missing:
        sys.exit(f"{what} 缺少列 {missing}。实际列名：{list(df.columns)}")


def load_meta(cache, ver):
    path = os.path.join(cache, f"dataset_metadata_{ver}.tsv")
    if not os.path.exists(path) or os.path.getsize(path) < 1000:
        r = requests.get(META_URL.format(ver=ver), timeout=120)
        r.raise_for_status()
        open(path, "wb").write(r.content)
    meta = pd.read_csv(path, sep="\t")
    need(meta, ["study_id", "dataset_id", "study_label", "tissue_label",
                "condition_label", "quant_method", "sample_size"], "数据集列表")
    return meta


def load_cs(study_id, dataset_id, cache):
    path = os.path.join(cache, f"{dataset_id}.credible_sets.tsv.gz")
    if not os.path.exists(path) or os.path.getsize(path) < 100:
        url = CS_URL.format(study_id=study_id, dataset_id=dataset_id)
        r = requests.get(url, timeout=600)
        if r.status_code != 200 or len(r.content) < 100:
            print(f"  下载失败 {r.status_code}: {url}")
            return None
        open(path, "wb").write(r.content)
    cs = pd.read_csv(path, sep="\t")
    need(cs, ["gene_id", "cs_id", "variant", "pip"], f"{dataset_id}")
    return cs


def is_snv(v):
    p = str(v).split("_")
    return len(p) >= 4 and len(p[-2]) == 1 and len(p[-1]) == 1


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--studies", nargs="+",
                    default=["Randolph_2021", "Randolph_2021_reannotated"])
    ap.add_argument("--pip", type=float, default=0.9)
    ap.add_argument("--cache", default="eqtl_cache")
    ap.add_argument("--ver", default="r8")
    args = ap.parse_args()
    os.makedirs(args.cache, exist_ok=True)

    meta = load_meta(args.cache, args.ver)
    lab = meta["study_label"].astype(str)
    sel = meta[lab.isin(args.studies) & (meta["quant_method"] == "ge")].copy()
    if sel.empty:
        sys.exit(f"没找到 {args.studies}。可选：{sorted(lab.unique())[:40]}")
    print(f"找到 {len(sel)} 个数据集，{sel.study_label.nunique()} 项研究，"
          f"{sel.tissue_label.nunique()} 种细胞类型", flush=True)

    tops, rows = {}, []
    for r in sel.itertuples(index=False):
        cs = load_cs(r.study_id, r.dataset_id, args.cache)
        if cs is None:
            continue
        cs["is_snv"] = cs["variant"].map(is_snv)
        top = cs.sort_values("pip", ascending=False).drop_duplicates("cs_id")
        tops[r.dataset_id] = (top, cs)
        hi = top[(top["pip"] > args.pip) & top["is_snv"]]
        rows.append({"study": r.study_label, "dataset_id": r.dataset_id,
                     "cell_type": r.tissue_label, "condition": r.condition_label,
                     "n_donors": r.sample_size,
                     "可信集数": int(cs["cs_id"].nunique()),
                     f"PIP>{args.pip}单碱基eQTL": int(len(hi))})
    ds = pd.DataFrame(rows)
    if ds.empty:
        sys.exit("没有任何可信集文件下载成功")
    ds.to_csv("stim_eqtl_datasets.csv", index=False)
    print("\n== 每个数据集 ==\n" + ds.to_string(index=False))

    # 同一研究同一细胞类型内配对 naive vs 刺激
    pairs = []
    for (study, cell), g in ds.groupby(["study", "cell_type"]):
        base = g[g["condition"] == "naive"]
        stim = g[g["condition"] != "naive"]
        if base.empty or stim.empty:
            continue
        bid = base.iloc[0]["dataset_id"]
        b_top, b_cs = tops[bid]
        b_keys = set(zip(b_cs["gene_id"], b_cs["variant"]))       # 基因-变异对
        b_genes = set(b_cs["gene_id"])                             # 该基因有任何可信集
        for s in stim.itertuples(index=False):
            s_top, _ = tops[s.dataset_id]
            hi = s_top[(s_top["pip"] > args.pip) & s_top["is_snv"]]
            shared = [(g_, v) in b_keys for g_, v in zip(hi["gene_id"], hi["variant"])]
            gene_abs = ~hi["gene_id"].isin(b_genes)
            pairs.append({
                "study": study, "cell_type": cell,
                "静息": base.iloc[0]["condition"], "刺激": s.condition,
                "刺激下高置信eQTL": len(hi),
                "共享(同基因同变异)": int(sum(shared)),
                "刺激特异候选(静息下该基因无可信集)": int(gene_abs.sum()),
            })
    if not pairs:
        sys.exit("没能配对，检查 condition_label 取值")
    p = pd.DataFrame(pairs)
    p.to_csv("stim_eqtl_pairs.csv", index=False)
    print("\n== 刺激 vs 静息 ==\n" + p.to_string(index=False))

    tot = int(p["刺激特异候选(静息下该基因无可信集)"].sum())
    tot_hi = int(p["刺激下高置信eQTL"].sum())
    print(f"\n刺激下高置信 eQTL 合计 {tot_hi}，其中刺激特异候选 {tot}"
          f"（跨细胞类型有重复）")
    print("标准：≥300 够做配对比较；100–300 勉强；<100 需加入其他研究。")

    # 去重后的独立变异数（真正决定样本量）
    uniq = []
    for (study, cell), g in ds.groupby(["study", "cell_type"]):
        base = g[g["condition"] == "naive"]
        stim = g[g["condition"] != "naive"]
        if base.empty or stim.empty:
            continue
        b_cs = tops[base.iloc[0]["dataset_id"]][1]
        for s in stim.itertuples(index=False):
            s_top = tops[s.dataset_id][0]
            hi = s_top[(s_top["pip"] > args.pip) & s_top["is_snv"]]
            hi = hi[~hi["gene_id"].isin(set(b_cs["gene_id"]))]
            uniq += list(zip(hi["gene_id"], hi["variant"]))
    print(f"去重后（基因-变异对）刺激特异候选：{len(set(uniq))}")


if __name__ == "__main__":
    main()
