#!/usr/bin/env python3
"""dex 条件的方向检查（MPRA A 路的最后一道前置）。

用 GSE46448（U2OS-GR，vehicle 0h vs dex 24h，各 4 重复）造两个条件，
在 GR 靶基因启动子窗口上打分，看预测的 dex/vehicle 比值是否朝已知方向（上调）。
若方向对，说明 bulk 数据能造出被模型正确理解的 dex 条件。
"""
import gzip
import json
import urllib.request

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from alphagenome.data import genome
from alphagenome.models.v2 import dna_client as c2, dna_output as do, interval_scorers as isc

TARGETS = {"FKBP5": "up", "PER1": "up", "TSC22D3": "up", "DUSP1": "up",
           "KLF15": "up", "ANGPTL4": "up", "PDK4": "up"}
OUT = "b_data/dex_direction.csv"


def parse_matrix(path):
    rows, titles, accs, head = [], None, None, None
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("!Sample_title"):
                titles = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_geo_accession"):
                accs = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            if line.startswith("!") or not line.strip():
                continue
            f2 = [x.strip('"') for x in line.rstrip("\n").split("\t")]
            if f2[0] == "ID_REF":
                head = f2
                continue
            if head is None or len(f2) != len(head):
                continue
            rows.append(f2)
    df = pd.DataFrame(rows, columns=head).set_index("ID_REF")
    df = df.apply(pd.to_numeric, errors="coerce")
    mp = dict(zip(accs, titles))
    df.columns = [mp.get(c, c) for c in df.columns]
    return df


def ensg_map(symbols):
    import re
    clean = [s for s in dict.fromkeys(symbols)
             if isinstance(s, str) and re.fullmatch(r"[A-Za-z0-9\.\-]{1,30}", s)]
    out = {}
    for i in range(0, len(clean), 500):
        chunk = clean[i:i + 500]
        body = json.dumps({"symbols": chunk}).encode()
        req = urllib.request.Request("https://rest.ensembl.org/lookup/symbol/homo_sapiens",
                                     data=body, headers={"Content-Type": "application/json"})
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                d = json.load(r)
            out.update({k: (v or {}).get("id") for k, v in d.items()})
        except Exception as e:
            print(f"  Ensembl 批次 {i} 失败 {type(e).__name__}")
        print(f"  映射 {len(out)}/{len(clean)}", flush=True)
    return out


def main():
    if __import__("os").path.exists("dex/dex_cond_linear.h5ad"):
        adata = ad.read_h5ad("dex/dex_cond_linear.h5ad")
        print(f"读线性尺度条件：{adata.shape}")
        return score(adata)
    if __import__("os").path.exists("dex/dex_cond.h5ad"):
        adata = ad.read_h5ad("dex/dex_cond.h5ad")
        print(f"读缓存条件 AnnData：{adata.shape}")
        return score(adata)
    mat = parse_matrix("dex/GSE46448_series_matrix.txt.gz")
    veh = [c for c in mat.columns if c.startswith("NI_0h")]
    dex = [c for c in mat.columns if c.startswith("NI_24h")]
    print(f"vehicle {len(veh)}：{veh}\ndex {len(dex)}：{dex}")

    with gzip.open("dex/GPL10558.annot.gz", "rt") as f:
        lines = f.readlines()
    hdr = next(i for i, l in enumerate(lines) if l.startswith("ID\t"))
    ann = pd.read_csv("dex/GPL10558.annot.gz", sep="\t", skiprows=hdr,
                      on_bad_lines="skip", low_memory=False)
    print("注释列:", list(ann.columns)[:8])
    probe_col = ann.columns[0]
    sym_col = [c for c in ann.columns if "symbol" in c.lower()][0]
    m = ann[[probe_col, sym_col]].dropna().drop_duplicates(probe_col)
    m.columns = ["probe", "symbol"]
    keep = mat.loc[mat.index.intersection(m.probe)]
    print(f"矩阵探针 {mat.shape[0]}，注释匹配 {len(keep)}")

    genes = sorted(set(m.symbol.astype(str)))
    emap = ensg_map(genes)
    m["ensg"] = m.symbol.astype(str).map(emap)
    m = m.dropna(subset=["ensg"])
    print(f"符号→ENSG 成功 {len(m)}/{len(genes)}")

    k = keep.merge(m, left_index=True, right_on="probe")
    val = k.groupby("ensg")[[*veh, *dex]].mean()      # 同基因多探针取均值
    X = val.T.values.astype(float)                    # 样本 × 基因
    obs = pd.DataFrame({"cond": ["veh"] * len(veh) + ["dex"] * len(dex)},
                       index=[*veh, *dex])
    adata = ad.AnnData(X=sp.csr_matrix(X), obs=obs,
                       var=pd.DataFrame(index=val.index.astype(str)))
    adata.write_h5ad("dex/dex_cond.h5ad")
    print(f"条件 AnnData：{adata.shape}（已缓存 dex/dex_cond.h5ad）")
    return score(adata)


def score(adata):
    cond = {"veh": adata[adata.obs.cond == "veh"].copy(),
            "dex": adata[adata.obs.cond == "dex"].copy()}

    # 靶基因坐标
    body = json.dumps({"symbols": list(TARGETS)}).encode()
    req = urllib.request.Request("https://rest.ensembl.org/lookup/symbol/homo_sapiens",
                                 data=body, headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=120) as r:
        info = json.load(r)
    rows = []
    model = c2.create(__import__("os").environ["ALPHAGENOME_API_KEY"], timeout=600)
    sc = [isc.GeneMaskScorer(requested_output=do.OutputType.SC_GENE_EXPRESSION_DENOISED,
                             width=200001, aggregation_type=isc.IntervalAggregationType.MEAN)]
    for sym, expect in TARGETS.items():
        v = info.get(sym)
        if not v or not v.get("id"):
            print(f"{sym}: 无坐标，跳过")
            continue
        tss = v["start"] if v["strand"] == 1 else v["end"]
        chrom = str(v["seq_region_name"])
        chrom = chrom if chrom.startswith("chr") else "chr" + chrom
        iv = genome.Interval(chrom,
                             max(0, tss - 262144), tss + 262144)
        vals = {}
        for tag, c in cond.items():
            a = model.score_interval(interval=iv, interval_scorers=sc,
                                     cell_conditioning=c, cell_condition_name=tag)[0]
            ids = [str(x).split(".")[0] for x in a.obs.gene_id]
            d = dict(zip(ids, [float(np.nanmean(x)) for x in np.asarray(a.X)]))
            vals[tag] = d.get(v["id"], np.nan)
        lfc = np.log2((vals["dex"] + 1) / (vals["veh"] + 1)) if vals["veh"] == vals["veh"] else np.nan
        rows.append(dict(symbol=sym, ensg=v["id"], veh=vals["veh"], dex=vals["dex"],
                         lfc=lfc, expected=expect,
                         ok=bool(lfc == lfc and ((expect == "up" and lfc > 0) or
                                                 (expect == "down" and lfc < 0)))))
        print(f"{sym:9s} veh={vals['veh']:.3f} dex={vals['dex']:.3f} "
              f"lfc={lfc:+.3f} 预期{expect} {'✔' if rows[-1]['ok'] else '✘'}")
    pd.DataFrame(rows).to_csv(OUT, index=False)
    oks = sum(r["ok"] for r in rows)
    print(f"\n方向一致 {oks}/{len(rows)} → {OUT}")


if __name__ == "__main__":
    main()
