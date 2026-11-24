#!/usr/bin/env python3
"""补充对照位点：在 PRAD 里测过但等位不平衡不显著的位点，用来算 AUROC。

数据来自 Zenodo 10.5281/zenodo.6371439 -> stratAS_results/PRAD_06.txt
（无表头，列序按 README：CHR POS RSID P0 P1 NAME CENTER N.HET N.READS
  ALL.AF ALL.BBINOM.P C0.AF C0.BBINOM.P C1.AF C1.BBINOM.P DIFF.BBINOM.P
  IND.C0 IND.C0.COUNT.REF IND.C0.COUNT.ALT IND.C1 IND.C1.COUNT.REF IND.C1.COUNT.ALT）
"""
import json
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor

import pandas as pd

COLS = ["CHR", "POS", "RSID", "P0", "P1", "NAME", "CENTER", "N.HET", "N.READS",
        "ALL.AF", "ALL.BBINOM.P", "C0.AF", "C0.BBINOM.P", "C1.AF",
        "C1.BBINOM.P", "DIFF.BBINOM.P", "IND.C0", "IND.C0.REF", "IND.C0.ALT",
        "IND.C1", "IND.C1.REF", "IND.C1.ALT"]
USE = ["CHR", "POS", "RSID", "NAME", "N.HET", "N.READS", "C1.AF", "C1.BBINOM.P"]
TXT = "stratAS_results/PRAD_06.txt"
HDR = {"Content-Type": "application/json", "User-Agent": "ar-pilot/1.0"}
N_CTRL = 400


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


def get_json(url, tries=3):
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


def main():
    cases = pd.read_csv("variants.csv", dtype={"variant_id": str})
    have = set(cases["rsid"].astype(str))
    print(f"已有病例位点 {len(cases)} 个", flush=True)

    df = pd.read_csv(TXT, sep="\t", names=COLS, usecols=USE,
                     dtype={"RSID": str, "NAME": str}, on_bad_lines="skip")
    print(f"PRAD 全量结果 {len(df)} 行，唯一 rsID {df['RSID'].nunique()}", flush=True)
    for c in ["N.HET", "N.READS", "C1.AF", "C1.BBINOM.P"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    ok = df[(df["RSID"].str.startswith("rs"))
            & (df["N.READS"] >= 200) & (df["N.HET"] >= 20)
            & (df["C1.BBINOM.P"] > 0.2)
            & ((df["C1.AF"] - 0.5).abs() < 0.04)
            & (~df["RSID"].isin(have))]
    ok = ok[ok["NAME"].str.startswith("PRAD_", na=False)]
    print(f"符合对照条件（测过、读段够、不平衡不显著）{len(ok)} 行，"
          f"唯一 rsID {ok['RSID'].nunique()}", flush=True)

    ctrl = (ok.assign(dev=(ok["C1.AF"] - 0.5).abs())
              .sort_values(["N.READS"], ascending=False)
              .drop_duplicates("RSID").head(N_CTRL))
    print(f"取前 {len(ctrl)} 个对照", flush=True)

    rsids = list(ctrl["RSID"])
    ann = {}
    for i in range(0, len(rsids), 150):
        chunk = rsids[i:i + 150]
        print(f"Ensembl 批量 {i + 1}-{i + len(chunk)}", flush=True)
        d = post_json("https://rest.ensembl.org/variation/human?content-type=application/json",
                      {"ids": chunk})
        if d:
            ann.update(d)
        time.sleep(2)

    rows = []
    for _, r in ctrl.iterrows():
        rec = ann.get(r["RSID"]) or {}
        maps = [m for m in rec.get("mappings", [])
                if m.get("assembly_name") == "GRCh38" and m.get("coord_system") == "chromosome"]
        if not maps:
            continue
        m = maps[0]
        al = str(m.get("allele_string", "")).split("/")
        if len(al) != 2 or any(len(a) != 1 for a in al):
            continue
        rows.append(dict(rsid=r["RSID"], chrom="chr" + str(m["seq_region_name"]),
                         pos=int(m["start"]), ref=al[0].upper(), alt=al[1].upper(),
                         C1AF=float(r["C1.AF"]), NAME=r["NAME"],
                         N_HET=int(r["N.HET"]), N_READS=int(r["N.READS"]),
                         P=float(r["C1.BBINOM.P"]), POS19=int(r["POS"]),
                         CHR19=str(r["CHR"])))
    cdf = pd.DataFrame(rows)
    print(f"对照位点拿到注释 {len(cdf)} 个", flush=True)

    with ThreadPoolExecutor(max_workers=10) as ex:
        cdf["hg38"] = list(ex.map(lambda r: ucsc_base(r["chrom"], r["pos"]),
                                  [r for _, r in cdf.iterrows()]))
    cdf = cdf[cdf["hg38"] == cdf["ref"]]
    print(f"参考碱基校验通过 {len(cdf)} 个", flush=True)

    cdf["variant_id"] = cdf["rsid"] + "_ctrl"
    cdf["effect"] = (cdf["C1AF"] - 0.5).round(4)
    cdf["significant"] = 0
    add = cdf.rename(columns={"C1AF": "C1.AF", "N_HET": "N.HET",
                              "N_READS": "N.READS", "P": "DIFF.BBINOM.P",
                              "POS19": "POS", "CHR19": "CHR"})
    add["maf"] = None
    add["minor"] = None
    keep = list(cases.columns) + [c for c in add.columns if c not in cases.columns]
    out = pd.concat([cases, add.reindex(columns=keep)], ignore_index=True)
    out.to_csv("variants.csv", index=False)
    print(f"\nvariants.csv 更新为 {len(out)} 行（病例 {(out['significant'] == 1).sum()}，"
          f"对照 {(out['significant'] == 0).sum()}）", flush=True)


if __name__ == "__main__":
    main()
