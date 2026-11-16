#!/usr/bin/env python3
"""B：用自定义细胞状态条件，检验模型能否复现 IFN-β 引起的表达变化。

判据（事先写死）：
  - 主指标：预测 log2FC 与实测 log2FC 的相关系数（跨基因）
  - 噪声基线：把对照细胞对半分成两组，同样流程算出「预测 log2FC」，
    它与实测的相关系数应接近 0
  - 通过：主相关系数 > 0.4；不通过：< 0.2

防线：输出带时间戳；读入后断言主键无重复；每步打印行数与空值数。
"""
import glob
import json
import os
import time
import urllib.request

import anndata as ad
import numpy as np
import pandas as pd
from alphagenome.data import genome
from alphagenome.models.v2 import dna_client as c2, dna_output as do, interval_scorers as isc

STAMP = time.strftime("%m%d_%H%M%S")
NBIN = 40          # 用基因最密集的 40 个 1 Mb 窗口
EPS = 0.1


def gene_coords_batch(ids):
    """Ensembl 批量 lookup，一次最多 1000 个 id。"""
    out = {}
    ids = list(ids)
    for i in range(0, len(ids), 900):
        chunk = ids[i:i + 900]
        body = json.dumps({"ids": chunk}).encode()
        for attempt in range(4):
            try:
                req = urllib.request.Request(
                    "https://rest.ensembl.org/lookup/id",
                    data=body,
                    headers={"Content-Type": "application/json",
                             "Accept": "application/json"})
                with urllib.request.urlopen(req, timeout=180) as r:
                    d = json.load(r)
                for k, v in d.items():
                    if isinstance(v, dict) and "seq_region_name" in v:
                        out[k] = ("chr" + str(v["seq_region_name"]),
                                  int(v["start"]), int(v["end"]))
                break
            except Exception as e:
                print(f"  lookup 重试 {attempt+1}: {type(e).__name__}", flush=True)
                time.sleep(3 * (attempt + 1))
    return out


def main():
    ctrl = ad.read_h5ad(sorted(glob.glob("b_data/mono_ctrl_*.h5ad"))[-1])
    stim = ad.read_h5ad(sorted(glob.glob("b_data/mono_stim_*.h5ad"))[-1])
    de = pd.read_csv(sorted(glob.glob("b_data/de_pseudobulk_*.csv"))[-1])
    assert not de.ensg.duplicated().any(), "de 表主键重复"
    print(f"ctrl {ctrl.shape}｜stim {stim.shape}｜表达基因 {(de.ensg.notna()).sum()}")

    # 噪声基线：对照细胞对半分
    idx = np.arange(ctrl.n_obs)
    rng = np.random.default_rng(11)
    rng.shuffle(idx)
    h1, h2 = ctrl[idx[:128]].copy(), ctrl[idx[128:]].copy()

    # 取基因密集的窗口
    expr = de[de.expressed].copy()
    print(f"表达的基因 {len(expr)}")
    coords = gene_coords_batch(expr.ensg.head(1500))
    print(f"拿到坐标 {len(coords)}")
    dfc = pd.DataFrame([(e, c[0], c[1], c[2]) for e, c in coords.items()],
                       columns=["ensg", "chrom", "start", "end"])
    dfc["bin"] = (dfc.start + dfc.end) // 2 // 1048576
    top = dfc.groupby(["chrom", "bin"]).size().sort_values(ascending=False).head(NBIN)
    print(f"选 {len(top)} 个窗口，覆盖基因 {dfc.set_index(['chrom','bin']).loc[top.index].shape[0] if len(top) else 0}")

    m = c2.create(os.environ["ALPHAGENOME_API_KEY"], timeout=600)
    sc = [isc.GeneMaskScorer(requested_output=do.OutputType.SC_GENE_EXPRESSION_DENOISED,
                             width=200001,
                             aggregation_type=isc.IntervalAggregationType.MEAN)]
    conds = {"ctrl": ctrl, "stim": stim, "half1": h1, "half2": h2}

    rows = []
    for i, ((chrom, b), _) in enumerate(top.items(), 1):
        mid = int(b) * 1048576 + 524288
        iv = genome.Interval(chrom, mid - 524288, mid + 524288)
        vals = {}
        for tag, cond in conds.items():
            try:
                a = m.score_interval(interval=iv, interval_scorers=sc,
                                     cell_conditioning=cond, cell_condition_name=tag)[0]
                ids = [str(x).split(".")[0] for x in a.obs.gene_id]
                vals[tag] = dict(zip(ids, [float(np.nanmean(x)) for x in np.asarray(a.X)]))
            except Exception as e:
                print(f"  {chrom}:{mid} {tag} 失败 {type(e).__name__}", flush=True)
                vals[tag] = {}
        for ensg in set(vals["ctrl"]) & set(vals["stim"]):
            rows.append(dict(ensg=ensg,
                             ctrl=vals["ctrl"][ensg], stim=vals["stim"][ensg],
                             half1=vals["half1"].get(ensg, np.nan),
                             half2=vals["half2"].get(ensg, np.nan)))
        if i % 10 == 0:
            print(f"  窗口 {i}/{len(top)}｜基因 {len(rows)}", flush=True)

    p = pd.DataFrame(rows).drop_duplicates("ensg")
    print(f"\n预测结果 {len(p)} 个基因，空值："
          + ", ".join(f"{c}={p[c].isna().sum()}" for c in p.columns))
    p = p.merge(de[["ensg", "log2fc", "cpm_ctrl", "cpm_stim"]], on="ensg", how="left")
    p["pred"] = np.log2((p.stim + EPS) / (p.ctrl + EPS))
    p["noise"] = np.log2((p.half1 + EPS) / (p.half2 + EPS))
    p.to_csv(f"b_data/b_pred_{STAMP}.csv", index=False)
    print(f"写出 b_data/b_pred_{STAMP}.csv")

    ok = p.dropna(subset=["pred", "noise", "log2fc"])
    ok = ok[(ok.ctrl > 1) | (ok.stim > 1)]
    print(f"\n参与相关分析的基因 {len(ok)}（过滤预测值全低的）")
    print(f"实测 log2FC 范围 {ok.log2fc.min():+.1f} ~ {ok.log2fc.max():+.1f}，"
          f"|log2FC|>1 的 {int((ok.log2fc.abs() > 1).sum())}")
    r = np.corrcoef(ok.pred, ok.log2fc)[0, 1]
    rn = np.corrcoef(ok.noise, ok.log2fc)[0, 1]
    from scipy import stats
    sp = stats.spearmanr(ok.pred, ok.log2fc).correlation
    print(f"\n主指标：pred vs 实测 Pearson r = {r:.3f}（Spearman {sp:.3f}）")
    print(f"噪声基线：half1 vs half2 的 r = {rn:.3f}")
    print(f"判据：>0.4 通过，<0.2 不通过 → " +
          ("通过" if r > 0.4 else ("不通过" if r < 0.2 else "灰区")))

    top_g = ok.reindex(ok.pred.abs().sort_values(ascending=False).index).head(12)
    print("\n预测变化最大的基因：")
    print(top_g[["ensg", "ctrl", "stim", "pred", "log2fc"]].round(2).to_string(index=False))


if __name__ == "__main__":
    main()
