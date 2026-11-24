#!/usr/bin/env python3
"""从 Baca et al. Nat Genet 2022 的 LNCaP 等位不平衡数据建表。

数据：Zenodo 10.5281/zenodo.6666796 -> cwas_paper_scripts/LNCaP.ai/data/
  K16h / KVh   : H3K27ac ChIP-seq，加雄激素 16h / 溶剂对照
  AT16h / ATVh : ATAC-seq，同上
  列序 CHR POS RSID P0 P1 NAME CENTER N.HET N.READS AF P COND（无表头）
  AF = ALT 等位基因读段比例（stratAS 同一套代码）=> effect = AF - 0.5

产出 variants_lncap.csv：
  effect       K16h（H3K27ac + 雄激素）的 AF-0.5
  effect_vh    KVh（溶剂对照）同一 rsID 的 AF-0.5，没有就留空
  effect_at    AT16h（ATAC + 雄激素）的 AF-0.5
  significant  1 = K16h 里 q<0.05 且 |effect|>=0.1；0 = 测过但平衡
"""
import json
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd

COLS = ["CHR", "POS", "RSID", "P0", "P1", "NAME", "CENTER", "N.HET",
        "N.READS", "AF", "P", "COND"]
CACHE = "ensembl_cache_lncap.json"
OUT = "variants_lncap.csv"
N_CASE, N_CTRL = 900, 600
MIN_READS = 100
HDR = {"Content-Type": "application/json", "User-Agent": "ar-pilot/1.0"}


def load(fn):
    d = pd.read_csv(f"lncap_ai/{fn}", sep="\t", names=COLS, dtype=str,
                    index_col=False, on_bad_lines="skip")
    for c in ["N.HET", "N.READS", "AF", "P"]:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d = d[d["N.READS"].notna() & d["RSID"].str.startswith("rs", na=False)]
    d["effect"] = d["AF"] - 0.5
    return d


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    q = np.empty(n)
    prev = 1.0
    for rank in range(n - 1, -1, -1):
        i = o[rank]
        prev = min(prev, p[i] * n / (rank + 1))
        q[i] = prev
    return q


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
                return None
            time.sleep(4 * (i + 1))
        except Exception:
            time.sleep(4 * (i + 1))
    return None


def fetch_ids(url, ids, depth=0):
    d = post_json(url, {"ids": ids})
    if d is not None:
        return d
    if len(ids) <= 6 or depth >= 5:
        return {}
    mid = len(ids) // 2
    return {**fetch_ids(url, ids[:mid], depth + 1),
            **fetch_ids(url, ids[mid:], depth + 1)}


def ucsc_base(chrom, pos1):
    url = (f"https://api.genome.ucsc.edu/getData/sequence?genome=hg38"
           f";chrom={chrom};start={pos1 - 1};end={pos1}")
    for _ in range(3):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return (json.load(r).get("dna") or "").upper() or None
        except Exception:
            time.sleep(2)
    return None


def freq_map(rec):
    pops = rec.get("populations") or []
    allf, anyf = {}, {}
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
    maps = [m for m in rec.get("mappings", [])
            if m.get("assembly_name") == "GRCh38" and m.get("coord_system") == "chromosome"]
    if not maps:
        return None
    m = maps[0]
    al = [a.upper() for a in str(m.get("allele_string", "")).split("/")]
    if not al or any(len(a) != 1 for a in al):
        return None
    if len(al) == 2:
        return "chr" + str(m["seq_region_name"]), int(m["start"]), al[0], al[1]
    fm = freq_map(rec)
    cands = sorted([(fm.get(a, 0.0), a) for a in al[1:] if a in fm], reverse=True)
    if not cands:
        return None
    return "chr" + str(m["seq_region_name"]), int(m["start"]), al[0], cands[0][1]


def annotate(rsids):
    cache = {}
    for seed in ["ensembl_cache_v2.json", CACHE]:
        try:
            cache.update({k: v for k, v in json.load(open(seed)).items() if v})
        except Exception:
            pass
    todo = [r for r in rsids if r not in cache]
    print(f"需查询 {len(todo)} / {len(rsids)}", flush=True)
    chunks = [todo[i:i + 40] for i in range(0, len(todo), 40)]
    done = {"n": 0}

    def work(chunk):
        d = fetch_ids("https://rest.ensembl.org/variation/human?content-type=application/json",
                      chunk)
        done["n"] += len(chunk)
        print(f"  {done['n']} / {len(todo)}", flush=True)
        return chunk, d or {}

    with ThreadPoolExecutor(max_workers=6) as ex:
        for chunk, d in ex.map(work, chunks):
            cache.update(d)
            for r in chunk:
                cache.setdefault(r, None)
            json.dump(cache, open(CACHE, "w"))

    def nalt(r):
        rec = cache.get(r) or {}
        maps = [m for m in rec.get("mappings", []) if m.get("assembly_name") == "GRCh38"]
        return len(str(maps[0].get("allele_string", "")).split("/")) if maps else 0

    multi = [r for r in rsids if nalt(r) > 2]
    print(f"多等位 {len(multi)} 个，补取频率", flush=True)
    mc = [multi[i:i + 20] for i in range(0, len(multi), 20)]

    def work2(chunk):
        return chunk, fetch_ids(
            "https://rest.ensembl.org/variation/human?pops=1;content-type=application/json",
            chunk) or {}

    with ThreadPoolExecutor(max_workers=4) as ex:
        for chunk, d in ex.map(work2, mc):
            cache.update(d)
            for r in chunk:
                cache.setdefault(r, None)

    return {k: v for k, v in cache.items() if v}


def main():
    k16 = load("K16h_stratas_results.txt")
    kvh = load("KVh_stratas_results.txt")
    at16 = load("AT16h_stratas_results.txt")

    # 每个 rsID 取 P 最小的一条（对应最强峰的测量）
    k16 = k16.sort_values("P").drop_duplicates("RSID")
    kvh = kvh.sort_values("P").drop_duplicates("RSID")
    at16 = at16.sort_values("P").drop_duplicates("RSID")
    k16["q"] = bh(k16["P"])
    kvh["q"] = bh(kvh["P"])

    print(f"K16h {len(k16)} 个 rsID，q<0.05: {(k16.q < 0.05).sum()}")
    print(f"KVh  {len(kvh)} 个 rsID，q<0.05: {(kvh.q < 0.05).sum()}")

    strong = k16[(k16.q < 0.05) & (k16["N.READS"] >= MIN_READS)
                 & (k16["effect"].abs() >= 0.10)]
    cases = strong.sort_values("P").head(N_CASE)
    ok = k16[(k16["N.READS"] >= MIN_READS) & (k16.q > 0.5)
             & (k16["effect"].abs() < 0.03)
             & (~k16["RSID"].isin(set(cases["RSID"])))]
    ctrl = ok.sort_values("N.READS", ascending=False).head(N_CTRL)
    print(f"病例 {len(cases)}（q<0.05 且 |effect|>=0.1 共 {len(strong)}），对照 {len(ctrl)}")

    sel = pd.concat([cases.assign(sig=1), ctrl.assign(sig=0)], ignore_index=True)
    kvh_map = kvh.set_index("RSID")["effect"]
    at_map = at16.set_index("RSID")["effect"]
    sel["effect_vh"] = sel["RSID"].map(kvh_map)
    sel["effect_at"] = sel["RSID"].map(at_map)
    print(f"其中在 KVh 里也测到的 {sel['effect_vh'].notna().sum()} 个，"
          f"在 AT16h 里也测到的 {sel['effect_at'].notna().sum()} 个")

    sel = sel.rename(columns={"N.READS": "nreads"})
    cache = annotate(list(sel["RSID"]))
    rows, drop = [], 0
    for r in sel.itertuples(index=False):
        got = resolve(cache.get(r.RSID) or {})
        if not got:
            drop += 1
            continue
        chrom, pos, ref, alt = got
        rows.append(dict(variant_id=r.RSID, rsid=r.RSID, chrom=chrom, pos=pos,
                         ref=ref, alt=alt, effect=round(float(r.effect), 4),
                         significant=int(r.sig),
                         effect_vh=(round(float(r.effect_vh), 4)
                                    if pd.notna(r.effect_vh) else None),
                         effect_at=(round(float(r.effect_at), 4)
                                    if pd.notna(r.effect_at) else None),
                         q=float(r.q), N_READS=int(r.nreads),
                         POS19=int(r.POS), CHR19=str(r.CHR)))
    v = pd.DataFrame(rows)
    print(f"注释成功 {len(v)}，失败 {drop}", flush=True)

    with ThreadPoolExecutor(max_workers=10) as ex:
        v["hg38"] = list(ex.map(lambda r: ucsc_base(r["chrom"], r["pos"]),
                                [r for _, r in v.iterrows()]))
    bad = int((v.hg38 != v.ref).sum())
    v = v[v.hg38 == v.ref].drop(columns=["hg38"])
    print(f"参考碱基不一致 {bad} 个，已剔除；剩 {len(v)}", flush=True)

    v.to_csv(OUT, index=False)
    print(f"\n写出 {OUT}：{len(v)} 行（病例 {(v.significant == 1).sum()}，"
          f"对照 {(v.significant == 0).sum()}）")


if __name__ == "__main__":
    main()
