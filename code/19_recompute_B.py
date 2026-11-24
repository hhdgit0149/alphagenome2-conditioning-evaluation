#!/usr/bin/env python3
"""B 重算（不调 AlphaGenome）：修实测比值、改 Spearman、去 CXCL 窗口、按窗口重采样。

改动：
  1. 两侧统一 log2((CPM+1)/(CPM+1))，避免分母趋零把比值撑爆
  2. 主指标 = Spearman，Pearson 仅作参考
  3. 去掉含 CXCL9/10/11 的那个窗口后再算一次
  4. 按 1 Mb 窗口做 cluster bootstrap 给 95% CI
  5. 完整打印 310 → 最终 的筛选链
"""
import json
import time
import urllib.request

import numpy as np
import pandas as pd
from scipy import stats

RNG = np.random.default_rng(20261113)
CXCL = {"ENSG00000138755", "ENSG00000169245", "ENSG00000169248"}  # CXCL9/10/11


def coords_batch(ids):
    out = {}
    ids = list(ids)
    for i in range(0, len(ids), 900):
        body = json.dumps({"ids": ids[i:i + 900]}).encode()
        for a in range(4):
            try:
                req = urllib.request.Request("https://rest.ensembl.org/lookup/id", data=body,
                                             headers={"Content-Type": "application/json",
                                                      "Accept": "application/json"})
                with urllib.request.urlopen(req, timeout=180) as r:
                    d = json.load(r)
                for k, v in d.items():
                    if isinstance(v, dict) and "seq_region_name" in v:
                        out[k] = ("chr" + str(v["seq_region_name"]), int(v["start"]))
                break
            except Exception:
                time.sleep(3 * (a + 1))
    return out


def boot_ci(df, B=2000):
    keys = df["win"].unique()
    vals = []
    for _ in range(B):
        pick = RNG.choice(len(keys), len(keys), replace=True)
        sub = pd.concat([df[df.win == keys[i]] for i in pick])
        if len(sub) > 5:
            vals.append(stats.spearmanr(sub["pred"], sub["obs"]).correlation)
    vals = np.array([v for v in vals if np.isfinite(v)])
    return np.percentile(vals, [2.5, 97.5]) if len(vals) > 50 else (np.nan, np.nan)


def main():
    p = pd.read_csv(sorted(__import__("glob").glob("b_data/b_pred_*.csv"))[-1])
    print(f"读入 {len(p)} 个基因（预测表原始行数）")

    # 筛选链
    n0 = len(p)
    p = p.drop_duplicates("ensg")
    print(f"  去重后 {len(p)}（重复 {n0 - len(p)}）")
    has_de = p.log2fc.notna()
    print(f"  有实测 log2FC 的 {int(has_de.sum())}（缺 {int((~has_de).sum())}）")
    p = p[has_de]
    keep = (p.ctrl > 1) | (p.stim > 1)
    print(f"  预测值 >1 CPM 的 {int(keep.sum())}（剔除 {int((~keep).sum())}，"
          f"规则：ctrl>1 或 stim>1）")
    p = p[keep].copy()

    # 统一算法
    p["obs"] = np.log2((p.cpm_stim + 1) / (p.cpm_ctrl + 1))
    p["pred"] = np.log2((p.stim + 1) / (p.ctrl + 1))
    p["noise"] = np.log2((p.half1 + 1) / (p.half2 + 1))

    # 窗口归属
    co = coords_batch(p.ensg)
    p["win"] = [f"{co[e][0]}:{co[e][1] // 1048576}" if e in co else "NA" for e in p.ensg]
    print(f"  能归到窗口的 {int((p.win != 'NA').sum())}／窗口数 {p.win.nunique()}")

    def report(d, tag):
        if len(d) < 10:
            print(f"\n{tag}: 样本 {len(d)}，跳过")
            return
        sp = stats.spearmanr(d["pred"], d["obs"])
        pe = np.corrcoef(d["pred"], d["obs"])[0, 1]
        rn = stats.spearmanr(d["noise"], d["obs"]).correlation
        lo, hi = boot_ci(d)
        print(f"\n{tag}（n={len(d)}，窗口 {d.win.nunique()}）")
        print(f"  Spearman = {sp.correlation:.3f} (p={sp.pvalue:.2g})  ← 主指标")
        print(f"  Pearson  = {pe:.3f}")
        print(f"  按窗口 bootstrap 95% CI: [{lo:.3f}, {hi:.3f}]")
        print(f"  噪声基线 Spearman = {rn:.3f}")

    report(p, "全部")
    report(p[~p.ensg.isin(CXCL)], "去掉 CXCL9/10/11")
    cx_win = p[p.ensg.isin(CXCL)]["win"].unique()
    report(p[~p.win.isin(cx_win)], f"去掉含 CXCL 的整个窗口（{len(cx_win)} 个）")

    print("\n实测 |log2FC| 最大的 8 个（新算法）：")
    t = p.reindex(p.obs.abs().sort_values(ascending=False).index).head(8)
    print(t[["ensg", "cpm_ctrl", "cpm_stim", "obs", "pred"]].round(2).to_string(index=False))
    p.to_csv("b_data/b_recomputed.csv", index=False)
    print("\n写出 b_data/b_recomputed.csv")


if __name__ == "__main__":
    main()
