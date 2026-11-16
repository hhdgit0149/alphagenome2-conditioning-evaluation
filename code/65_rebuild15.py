#!/usr/bin/env python3
"""重出 15 点校准表：单一基因集（base 119）、单一基线口径（各类型全部细胞）。

修两处：
  (1) 未刺激 7 点的实测重叠此前用「排除条件侧细胞」的基线（2785→2273），与刺激后那批不一致；
      现在统一用导出的 CPM（各类型全部细胞）
  (2) 5 个未刺激条件此前只有 110 个基因，已用 15 次调用补齐 → 现在 15 点全在 119 基因上
"""
import glob

import numpy as np
import pandas as pd
from scipy import stats

D = "b_data"
TYPES = ["CD4 T cells", "CD8 T cells", "NK cells", "B cells", "Dendritic cells",
         "Megakaryocytes", "FCGR3A+ Monocytes"]


def main():
    tb = pd.read_csv(f"{D}/null_tb_filled.csv").set_index("ensg")
    near = pd.read_csv(f"{D}/null_near_1115.csv").set_index("ensg")
    m1 = pd.read_csv(f"{D}/calib_mid_pred_1115_121815.csv").set_index("ensg")
    m2 = pd.read_csv(f"{D}/calib_mid_pred_1115_130035.csv").set_index("ensg")
    pred = pd.read_csv(f"{D}/b_pred_disjoint_1115_115447.csv").set_index("ensg")
    ho = pd.read_csv(f"{D}/holdout_pred_1115_121546.csv").set_index("ensg")
    C = pd.read_csv(f"{D}/for_audit_celltype_cpm.csv").set_index("ensg")

    base = pred[(pred[["ctrl", "stim"]].max(axis=1) > 1) & pred.obs.notna()]
    ids = base.index.tolist()
    obs = base.obs.values
    win = base.win
    print(f"base {len(ids)} 基因 / {win.nunique()} 窗口")

    ctrl_c = C.loc[ids, "ctrl|CD14+ Monocytes"].values
    byw = {w: np.where(win.values == w)[0] for w in win.unique()}
    rng = np.random.default_rng(0)

    def meas(kind, t):
        col = f"{'ctrl' if kind == 'unstim' else 'stim'}|{t}"
        v = np.log2((C.loc[ids, col].values + 1) / (ctrl_c + 1))
        r = stats.spearmanr(v, obs).correlation
        bs = []
        for _ in range(2000):
            ii = np.concatenate([byw[w] for w in rng.choice(list(byw), len(byw), True)])
            bs.append(stats.spearmanr(v[ii], obs[ii]).correlation)
        lo, hi = np.percentile(bs, [2.5, 97.5])
        return r, lo, hi

    treat = {}
    for t in TYPES:
        treat[(t, "unstim")] = (tb if t in ("CD4 T cells", "CD8 T cells", "NK cells", "B cells",
                                           "Dendritic cells") else near)
    tbl = {"tcell": "CD4 T cells", "bcell": "B cells", "cd8": "CD8 T cells",
           "nk": "NK cells", "dc": "Dendritic cells"}

    rows = []
    for t in TYPES:
        if t in ("FCGR3A+ Monocytes", "Megakaryocytes"):
            v = near.loc[ids, "fcgr3a" if t.startswith("FCGR") else "mega"].values
            b = pred.loc[ids, "ctrl"].values
        else:
            v = tb.loc[ids, {v_: k for k, v_ in tbl.items()}[t]].values
            b = pred.loc[ids, "ctrl"].values   # 统一基线；已验证与 tb.ctrl 完全相同
        y = stats.spearmanr(np.log2((v + 1) / (b + 1)), obs).correlation
        x, lo, hi = meas("unstim", t)
        rows.append(dict(cond=t, kind="unstim", measured=x, lo=lo, hi=hi, model=y))

    stim_src = {"CD4 T cells": m1, "Dendritic cells": m1, "NK cells": m1,
                "FCGR3A+ Monocytes": m1, "CD8 T cells": m2, "B cells": m2,
                "Megakaryocytes": m2}
    for t in TYPES:
        s = stim_src[t]
        v = s.loc[ids, t].values
        b = s.loc[ids, "cd14ctrl"].values
        y = stats.spearmanr(np.log2((v + 1) / (b + 1)), obs).correlation
        x, lo, hi = meas("stim", t)
        rows.append(dict(cond=t, kind="stim", measured=x, lo=lo, hi=hi, model=y))

    rows.append(dict(cond="IFN-beta (matched)", kind="matched", measured=1.0, lo=np.nan,
                     hi=np.nan, model=stats.spearmanr(
                         np.log2((pred.loc[ids, "stim"] + 1) / (pred.loc[ids, "ctrl"] + 1)),
                         obs).correlation))
    R = pd.DataFrame(rows)
    R["dev"] = R.model - R.measured
    print(R.round(4).to_string(index=False))
    R.to_csv(f"{D}/calibration_final.csv", index=False)

    print()
    for lab, sub in [("全 15 点", R), ("去匹配 14 点", R[R.kind != "matched"]),
                     ("仅未刺激 7", R[R.kind == "unstim"]), ("仅刺激后 7", R[R.kind == "stim"])]:
        x, y = sub.measured.values, sub.model.values
        A = np.column_stack([np.ones(len(x)), x])
        be, *_ = np.linalg.lstsq(A, y, rcond=None)
        pr = stats.pearsonr(x, y); sp = stats.spearmanr(x, y)
        print(f"{lab:12s} n={len(x):2d} 斜率 {be[1]:.3f} 截距 {be[0]:+.3f} "
              f"Pearson {pr.statistic:.3f} (p={pr.pvalue:.1e}) Spearman {sp.statistic:.3f} (p={sp.pvalue:.3f})")

    res = R.model - (np.polyval(np.polyfit(R.measured, R.model, 1), R.measured))
    ok = np.isfinite(res)
    print(f"\n残差 SD {res[ok].std(ddof=1):.3f}，最大 |残差| {np.abs(res[ok]).max():.3f}"
          f"（{R.cond[ok][np.abs(res[ok]).values.argmax()]}）")
    print(f"校正斜率（r_full=0.9832）= {0.709 / 0.9832:.3f} → 本次斜率/0.9832 = "
          f"{np.polyfit(R.measured, R.model, 1)[0] / 0.9832:.3f}")

    # 更新长表
    long = []
    for _, r in R.iterrows():
        long.append(dict(cond=r.cond, kind=r.kind, measured=r.measured,
                         model=r.model, dev=r.dev, lo=r.lo, hi=r.hi))
    pd.DataFrame(long).to_csv(f"{D}/calibration_final.csv", index=False)
    print(f"\n写出 {D}/calibration_final.csv")


if __name__ == "__main__":
    main()
