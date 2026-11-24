#!/usr/bin/env python3
"""扩样版本：全量 PRAD as-aQTL（898 个 rsID）+ 600 个对照，并抢救多等位位点。

多等位位点（Ensembl allele_string 有 3~4 个碱基）用 1000G ALL 频率挑出真正的第二个
等位基因：ref 固定为参考碱基，alt 取非 ref 中频率最高的那个。这样能把上一版
因为 "G/A/C" 之类被丢掉的一百多个位点救回来。

输入：sup/Supplementary Table 2.xlsx (PRAD 页) + stratAS_results/PRAD_06.txt
输出：variants_v2.csv
"""
import json
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor

import pandas as pd

TABLE = "sup/Supplementary Table 2.xlsx"
PRAD_TXT = "stratAS_results/PRAD_06.txt"
COLS = ["CHR", "POS", "RSID", "P0", "P1", "NAME", "CENTER", "N.HET", "N.READS",
        "ALL.AF", "ALL.BBINOM.P", "C0.AF", "C0.BBINOM.P", "C1.AF",
        "C1.BBINOM.P", "DIFF.BBINOM.P", "IND.C0", "IND.C0.REF", "IND.C0.ALT",
        "IND.C1", "IND.C1.REF", "IND.C1.ALT"]
OUT = "variants_v2.csv"
CACHE = "ensembl_cache_v2.json"
N_CTRL = 600
HDR = {"Content-Type": "application/json", "User-Agent": "ar-pilot/1.0"}


def post_json(url, payload, tries=3):
    body = json.dumps(payload).encode()
    ids = payload["ids"]
    for i in range(tries):
        try:
            req = urllib.request.Request(url, data=body, headers=HDR, method="POST")
            with urllib.request.urlopen(req, timeout=240) as r:
                return json.load(r)
        except urllib.error.HTTPError as err:
            if err.code == 400:
                return None                      # 有查不到的 id，交给上层拆分
            print(f"  HTTP {err.code}，{4 * (i + 1)}s 后重试", flush=True)
            time.sleep(4 * (i + 1))
        except Exception as err:
            print(f"  {type(err).__name__}，{4 * (i + 1)}s 后重试", flush=True)
            time.sleep(4 * (i + 1))
    return None


def fetch_ids(url, ids, depth=0):
    """失败就把批次对半拆开重试，直到单个 id。"""
    d = post_json(url, {"ids": ids})
    if d is not None:
        return d
    if len(ids) <= 6 or depth >= 5:
        return {}
    mid = len(ids) // 2
    a = fetch_ids(url, ids[:mid], depth + 1)
    b = fetch_ids(url, ids[mid:], depth + 1)
    return {**(a or {}), **(b or {})}


def get_json(url, tries=4):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "ar-pilot/1.0"})
            with urllib.request.urlopen(req, timeout=60) as r:
                return json.load(r)
        except Exception:
            time.sleep(2 * (i + 1))
    return None


def ucsc_base(chrom, pos1):
    d = get_json(f"https://api.genome.ucsc.edu/getData/sequence?genome=hg38"
                 f";chrom={chrom};start={pos1 - 1};end={pos1}")
    return ((d or {}).get("dna") or "").upper() or None


def freq_map(rec):
    """每个等位基因的最高群体频率（优先 1000G phase3 ALL）。"""
    pops = rec.get("populations") or []
    allf = {}
    anyf = {}
    for p in pops:
        a = str(p.get("allele", "")).upper()
        f = p.get("frequency")
        if not a or f is None:
            continue
        anyf[a] = max(anyf.get(a, 0.0), float(f))
        if p.get("population") == "1000GENOMES:phase_3:ALL":
            allf[a] = max(allf.get(a, 0.0), float(f))
    return allf or anyf


def resolve(rec):
    """返回 (chrom, pos, ref, alt, 来源) 或 None。"""
    maps = [m for m in rec.get("mappings", [])
            if m.get("assembly_name") == "GRCh38" and m.get("coord_system") == "chromosome"]
    if not maps:
        return None
    m = maps[0]
    alleles = [a.upper() for a in str(m.get("allele_string", "")).split("/")]
    if not alleles or any(len(a) != 1 for a in alleles):
        return None
    ref = alleles[0]
    if len(alleles) == 2:
        return "chr" + str(m["seq_region_name"]), int(m["start"]), ref, alleles[1], "biallelic"
    fm = freq_map(rec)
    cands = [(fm.get(a, 0.0), a) for a in alleles[1:]]
    cands = [c for c in cands if c[1] in fm]
    if not cands:
        return None
    cands.sort(reverse=True)
    return ("chr" + str(m["seq_region_name"]), int(m["start"]), ref,
            cands[0][1], f"multiallelic->{cands[0][1]}")


def annotate(rsids, tag):
    """两段式：先不取群体频率（响应小、稳），再只对多等位位点补取频率。"""
    try:
        cache = json.load(open(CACHE))
    except Exception:
        cache = {}
    try:
        pops = json.load(open(CACHE.replace(".json", "_pops.json")))
    except Exception:
        pops = {}
    cache = {k: v for k, v in cache.items() if v}   # 丢掉上次失败留下的空值，重新查
    pops = {k: v for k, v in pops.items() if v}

    todo = [r for r in rsids if r not in cache]
    print(f"{tag}: {len(rsids)} 个 rsID，需查询 {len(todo)}", flush=True)
    chunks = [todo[i:i + 40] for i in range(0, len(todo), 40)]
    done = {"n": 0}

    def work(chunk):
        d = fetch_ids("https://rest.ensembl.org/variation/human?content-type=application/json",
                      chunk)
        done["n"] += len(chunk)
        print(f"  基础注释 {done['n']} / {len(todo)}", flush=True)
        return chunk, d or {}

    with ThreadPoolExecutor(max_workers=6) as ex:
        for chunk, d in ex.map(work, chunks):
            cache.update(d)
            for r in chunk:
                cache.setdefault(r, None)
            json.dump(cache, open(CACHE, "w"))

    def nalt(r):
        rec = cache.get(r) or {}
        maps = [m for m in rec.get("mappings", [])
                if m.get("assembly_name") == "GRCh38"]
        if not maps:
            return 0
        return len(str(maps[0].get("allele_string", "")).split("/"))

    multi = [r for r in rsids if nalt(r) > 2 and r not in pops]
    print(f"  多等位位点 {len(multi)} 个，补取群体频率", flush=True)
    mchunks = [multi[i:i + 20] for i in range(0, len(multi), 20)]
    mdone = {"n": 0}

    def work2(chunk):
        d = fetch_ids("https://rest.ensembl.org/variation/human?pops=1;content-type=application/json",
                      chunk)
        mdone["n"] += len(chunk)
        if mdone["n"] % 100 < 20:
            print(f"    频率 {mdone['n']} / {len(multi)}", flush=True)
        return chunk, d or {}

    with ThreadPoolExecutor(max_workers=4) as ex:
        for chunk, d in ex.map(work2, mchunks):
            pops.update(d)
            for r in chunk:
                pops.setdefault(r, None)
            json.dump(pops, open(CACHE.replace(".json", "_pops.json"), "w"))

    for r in rsids:
        if r in pops and pops[r]:
            cache[r] = pops[r]
    return cache


def main():
    # ---------- 病例：表 2 PRAD 页全部 rsID ----------
    x = pd.read_excel(TABLE, sheet_name="PRAD")
    x["RSID"] = x["RSID"].astype(str)
    x = x[x["RSID"].str.startswith("rs")].copy()
    print(f"表 2 PRAD 页 {len(x)} 行，唯一 rsID {x['RSID'].nunique()}", flush=True)

    full = pd.read_csv(PRAD_TXT, sep="\t", names=COLS,
                       usecols=["RSID", "NAME", "C1.AF", "N.HET", "N.READS",
                                "C1.BBINOM.P", "POS", "CHR"],
                       dtype={"RSID": str, "NAME": str}, on_bad_lines="skip")
    full = full[full["RSID"].isin(set(x["RSID"]))]
    for c in ["C1.AF", "N.HET", "N.READS", "C1.BBINOM.P"]:
        full[c] = pd.to_numeric(full[c], errors="coerce")
    full["dev"] = (full["C1.AF"] - 0.5).abs()
    cases = (full.sort_values("dev", ascending=False)
                 .drop_duplicates("RSID"))
    print(f"病例候选 {len(cases)} 个（每个 rsID 取不平衡最强的一条）", flush=True)

    # ---------- 对照：测过但不显著 ----------
    df = pd.read_csv(PRAD_TXT, sep="\t", names=COLS,
                     usecols=["RSID", "NAME", "C1.AF", "N.HET", "N.READS",
                              "C1.BBINOM.P", "POS", "CHR"],
                     dtype={"RSID": str, "NAME": str}, on_bad_lines="skip")
    for c in ["C1.AF", "N.HET", "N.READS", "C1.BBINOM.P"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    ok = df[(df["RSID"].str.startswith("rs"))
            & (df["N.READS"] >= 200) & (df["N.HET"] >= 20)
            & (df["C1.BBINOM.P"] > 0.2)
            & ((df["C1.AF"] - 0.5).abs() < 0.04)
            & (~df["RSID"].isin(set(cases["RSID"])))]
    ctrl = (ok.sort_values("N.READS", ascending=False)
              .drop_duplicates("RSID").head(N_CTRL))
    print(f"对照取 {len(ctrl)} 个", flush=True)

    # ---------- 注释 ----------
    cache = annotate(list(cases["RSID"]) + list(ctrl["RSID"]), "全部位点")

    def build(src, sig, suffix):
        rows, drop = [], 0
        for _, r in src.iterrows():
            got = resolve(cache.get(r["RSID"]) or {})
            if not got:
                drop += 1
                continue
            chrom, pos, ref, alt, how = got
            rows.append(dict(variant_id=r["RSID"] + suffix, rsid=r["RSID"],
                             chrom=chrom, pos=pos, ref=ref, alt=alt,
                             C1AF=float(r["C1.AF"]), NAME=r["NAME"],
                             N_HET=int(r["N.HET"]), N_READS=int(r["N.READS"]),
                             P=float(r["C1.BBINOM.P"]), POS19=int(r["POS"]),
                             CHR19=str(r["CHR"]), how=how, significant=sig))
        print(f"  {suffix or 'case'}: 保留 {len(rows)}，丢弃 {drop}", flush=True)
        return pd.DataFrame(rows)

    cdf = build(cases, 1, "")
    kdf = build(ctrl, 0, "_ctrl")
    allv = pd.concat([cdf, kdf], ignore_index=True)

    print("校验 GRCh38 参考碱基...", flush=True)
    with ThreadPoolExecutor(max_workers=10) as ex:
        allv["hg38"] = list(ex.map(lambda r: ucsc_base(r["chrom"], r["pos"]),
                                   [r for _, r in allv.iterrows()]))
    bad = int((allv["hg38"] != allv["ref"]).sum())
    allv = allv[allv["hg38"] == allv["ref"]]
    print(f"参考碱基不一致 {bad} 个，已剔除；剩 {len(allv)}", flush=True)

    allv["effect"] = (allv["C1AF"] - 0.5).round(4)
    out = allv[["variant_id", "chrom", "pos", "ref", "alt", "effect", "significant",
                "rsid", "NAME", "N_HET", "N_READS", "P", "POS19", "CHR19", "how"]]
    out = out.rename(columns={"P": "DIFF.BBINOM.P", "POS19": "POS", "CHR19": "CHR",
                              "N_HET": "N.HET", "N_READS": "N.READS"})
    out.to_csv(OUT, index=False)
    print(f"\n写出 {OUT}：{len(out)} 行 "
          f"（病例 {(out.significant == 1).sum()}，对照 {(out.significant == 0).sum()}）", flush=True)
    print("\n多等位抢救情况：\n" + out["how"].value_counts().to_string(), flush=True)


if __name__ == "__main__":
    main()
