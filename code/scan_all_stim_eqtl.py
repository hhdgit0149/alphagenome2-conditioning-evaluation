#!/usr/bin/env python3
"""扫描 eQTL Catalogue 里所有「静息 vs 刺激」配对，统计刺激特异的高置信 eQTL。

判定（事先写死，避免事后调参）：
  - 只取 PIP>0.9 的单碱基变异，且是该可信集里 PIP 最高的变异
  - 「共享」= 同一 (基因, 变异) 在静息条件的高置信集里也出现
  - 「刺激特异候选」= 该基因在静息条件下没有任何可信集
  - 输出面板：配对级、变异级去重、基因级去重
"""
import os
import sys

import pandas as pd
import requests

META = "eqtl_cache/dataset_metadata_r8.tsv"
CS = ("https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/susie/"
      "{s}/{d}/{d}.credible_sets.tsv.gz")
BASE = {"naive", "control", "untreated", "vehicle", "baseline", "0h", "mock"}
PIP = 0.9


def get_cs(s, d, cache="eqtl_cache"):
    p = os.path.join(cache, f"{d}.credible_sets.tsv.gz")
    if os.path.exists(p) and os.path.getsize(p) < 200:
        return None                      # 已知不存在（196 字节的错误页已删，这里兜底）
    if not os.path.exists(p) or os.path.getsize(p) < 200:
        r = None
        for attempt in range(4):          # FTP 偶尔返回假 404，重试
            try:
                r = requests.get(CS.format(s=s, d=d), timeout=600)
            except Exception:
                r = None
            if r is not None and r.status_code == 200 and len(r.content) > 200:
                break
            import time
            time.sleep(2 * (attempt + 1))
        if r is None or r.status_code != 200 or len(r.content) < 200:
            print(f"    [{d}] HTTP {getattr(r, 'status_code', 'ERR')} 放弃")
            return None
        open(p, "wb").write(r.content)
    try:
        df = pd.read_csv(p, sep="\t")
    except Exception as e:
        print(f"    [{d}] 读取失败 {type(e).__name__}: {str(e)[:80]}")
        return None
    if not {"gene_id", "cs_id", "variant", "pip"} <= set(df.columns):
        print(f"    [{d}] 列名异常 {list(df.columns)[:6]}")
        return None
    v = df["variant"].astype(str).str.split("_")
    df["snv"] = v.str.len().ge(4) & v.str[-2].str.len().eq(1) & v.str[-1].str.len().eq(1)
    return df


def hi_cs(df):
    top = df.sort_values("pip", ascending=False).drop_duplicates("cs_id")
    return top[(top.pip > PIP) & top.snv]


def main():
    m = pd.read_csv(META, sep="\t")
    g = m[m.quant_method == "ge"]
    pairs = []
    for (s, t), x in g.groupby(["study_label", "tissue_label"]):
        b = x[x.condition_label.astype(str).str.lower().isin(BASE)]
        st = x[~x.condition_label.astype(str).str.lower().isin(BASE)]
        if b.empty or st.empty:
            continue
        pairs.append((s, t, b.iloc[0].dataset_id, b.iloc[0].condition_label,
                      list(zip(st.dataset_id, st.condition_label))))
    print(f"配对 {len(pairs)} 组，涉及 {len({p[0] for p in pairs})} 项研究", flush=True)

    rows, all_spec, all_hi = [], [], []
    for s, t, bid, bcond, stim in pairs:
        bcs = get_cs(s, bid) if True else None
        if bcs is None:
            print(f"  跳过 {s}/{t}（静息无文件）")
            continue
        b_keys = set(zip(bcs.gene_id, bcs.variant))
        b_genes = set(bcs.gene_id)
        for did, cond in stim:
            scs = get_cs(s, did)
            if scs is None:
                continue
            hi = hi_cs(scs)
            if hi.empty:
                rows.append(dict(study=s, cell=t, base=bcond, stim=cond,
                                 n_hi=0, shared=0, spec=0))
                continue
            shared = sum((a, b) in b_keys for a, b in zip(hi.gene_id, hi.variant))
            spec = hi[~hi.gene_id.isin(b_genes)]
            rows.append(dict(study=s, cell=t, base=bcond, stim=cond,
                             n_hi=len(hi), shared=int(shared), spec=len(spec)))
            for r in spec.itertuples(index=False):
                all_spec.append((r.gene_id, r.variant))
            for r in hi.itertuples(index=False):
                all_hi.append((r.gene_id, r.variant))
            print(f"  {s:26s} {t[:18]:18s} {cond[:16]:16s} 高置信 {len(hi):4d} "
                  f"特异 {len(spec):4d}", flush=True)

    res = pd.DataFrame(rows)
    res.to_csv("stim_scan_pairs.csv", index=False)
    print("\n== 汇总 ==")
    print(res.groupby("study")[["n_hi", "shared", "spec"]].sum().to_string())
    print(f"\n全部高置信 eQTL（基因-变异对，去重）: {len(set(all_hi))}")
    print(f"刺激特异候选（基因-变异对，去重）: {len(set(all_spec))}")
    print(f"刺激特异候选（仅变异去重）: {len({v for _, v in all_spec})}")
    print(f"刺激特异候选（仅基因去重）: {len({g_ for g_, _ in all_spec})}")
    print("\n标准（事先定）：≥300 够做配对比较；100–300 勉强；<100 需换方向。")
    pd.DataFrame(sorted({v for _, v in all_spec}), columns=["variant"]).to_csv(
        "stim_specific_variants.csv", index=False)


if __name__ == "__main__":
    main()
