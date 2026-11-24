#!/usr/bin/env python3
"""从补充表 2 的 PRAD 页整理 variants.csv。

来源：Grishin & Gusev, Nat Genet 2022, Supplementary Table 2 / sheet PRAD
  （as-aQTL variants at 10% FDR，hg19 坐标）
  C1.AF = ALT 等位基因的读段比例（stratAS 源码 stratas.R: binom.test(ALT, REF+ALT)）
  => effect = C1.AF - 0.5，正数表示 ALT 等位基因染色质更开放

等位基因与 GRCh38 坐标来自 Ensembl variation REST（批量 POST）；
用 UCSC hg38 参考碱基校验 ref，避免给 AlphaGenome 传错碱基。
"""
import json
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor

import pandas as pd

TABLE = "sup/Supplementary Table 2.xlsx"
SHEET = "PRAD"
OUT = "variants.csv"
CACHE = "ensembl_cache.json"
HDR = {"Content-Type": "application/json", "User-Agent": "ar-pilot/1.0"}


def post_json(url, payload, tries=4):
    body = json.dumps(payload).encode()
    for i in range(tries):
        try:
            req = urllib.request.Request(url, data=body, headers=HDR, method="POST")
            with urllib.request.urlopen(req, timeout=180) as r:
                return json.load(r)
        except Exception as err:
            print(f"  重试 {i + 1}/{tries}: {err}", flush=True)
            time.sleep(5 * (i + 1))
    return None


def get_json(url, tries=4):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "ar-pilot/1.0"})
            with urllib.request.urlopen(req, timeout=60) as r:
                return json.load(r)
        except Exception as err:
            print(f"  重试 {i + 1}/{tries}: {err}", flush=True)
            time.sleep(3 * (i + 1))
    return None


def fetch_ensembl(rsids):
    out = {}
    for i in range(0, len(rsids), 150):
        chunk = rsids[i:i + 150]
        print(f"Ensembl 批量 {i + 1}-{i + len(chunk)} / {len(rsids)}", flush=True)
        d = post_json("https://rest.ensembl.org/variation/human?content-type=application/json",
                      {"ids": chunk})
        if d:
            out.update(d)
        time.sleep(2)
    return out


def ucsc_base(chrom, pos1, genome="hg38"):
    url = (f"https://api.genome.ucsc.edu/getData/sequence?genome={genome}"
           f";chrom={chrom};start={pos1 - 1};end={pos1}")
    d = get_json(url, tries=3)
    return ((d or {}).get("dna") or "").upper() or None


def main():
    x = pd.read_excel(TABLE, sheet_name=SHEET)
    x["NAME"] = x["NAME"].astype(str)
    p = x[x["NAME"].str.startswith("PRAD_")].copy()
    print(f"表 2 PRAD 页 {len(x)} 行，PRAD 专属峰上位点 {len(p)} 行", flush=True)

    rsids = [r for r in dict.fromkeys(p["RSID"].astype(str)) if r.startswith("rs")]
    print(f"唯一 rsID：{len(rsids)}", flush=True)

    d = fetch_ensembl(rsids)
    json.dump(d, open(CACHE, "w"))

    rows, skipped = [], []
    for rs in rsids:
        rec = d.get(rs) or {}
        maps = [m for m in rec.get("mappings", [])
                if m.get("assembly_name") == "GRCh38" and m.get("coord_system") == "chromosome"]
        if not maps:
            skipped.append((rs, "无 GRCh38 定位"))
            continue
        m = maps[0]
        alleles = str(m.get("allele_string", "")).split("/")
        if len(alleles) != 2 or any(len(a) != 1 for a in alleles):
            skipped.append((rs, f"非双等位 SNP {m.get('allele_string')}"))
            continue
        rows.append(dict(rsid=rs, chrom="chr" + str(m["seq_region_name"]),
                         pos=int(m["start"]), ref=alleles[0].upper(),
                         alt=alleles[1].upper(), maf=rec.get("MAF"),
                         minor=rec.get("minor_allele")))
    print(f"Ensembl 注释可用 {len(rows)}，跳过 {len(skipped)}", flush=True)
    for s in skipped[:10]:
        print("  跳过", s, flush=True)

    ann = pd.DataFrame(rows)
    print("校验 GRCh38 参考碱基...", flush=True)

    def check(row):
        return ucsc_base(row["chrom"], row["pos"])
    with ThreadPoolExecutor(max_workers=10) as ex:
        ann["hg38"] = list(ex.map(check, [r for _, r in ann.iterrows()]))
    ann["ref_ok"] = ann["hg38"] == ann["ref"]
    print(f"GRCh38 参考碱基一致：{int(ann['ref_ok'].sum())}/{len(ann)}", flush=True)

    p["abs_dev"] = (p["C1.AF"] - 0.5).abs()
    first = (p.sort_values("abs_dev", ascending=False)
              .drop_duplicates("RSID")[["RSID", "NAME", "C1.AF", "N.HET",
                                        "N.READS", "DIFF.BBINOM.P", "POS", "CHR"]])
    out = ann[ann["ref_ok"]].merge(first, left_on="rsid", right_on="RSID", how="inner")
    out["variant_id"] = out["rsid"]
    out["effect"] = (out["C1.AF"] - 0.5).round(4)
    out["significant"] = 1
    res = out[["variant_id", "chrom", "pos", "ref", "alt", "effect", "significant",
               "rsid", "NAME", "N.HET", "N.READS", "DIFF.BBINOM.P", "maf",
               "minor", "POS", "CHR"]]
    res.to_csv(OUT, index=False)
    print(f"\n写出 {OUT}：{len(res)} 个位点", flush=True)
    print(res.head(8).to_string(index=False), flush=True)
    print("\neffect 分布:\n" + res["effect"].describe().to_string(), flush=True)
    print("正效应（ALT 更开放）比例:", round(float((res["effect"] > 0).mean()), 3), flush=True)


if __name__ == "__main__":
    main()
