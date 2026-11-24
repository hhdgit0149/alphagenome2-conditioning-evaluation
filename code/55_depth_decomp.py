#!/usr/bin/env python3
"""深度 vs 循环的分解，以及把已有五个零分布条件重挂到不相交实测侧。

两个问题：
  (1) 0.486 → 0.6155 的升高，是文库变深（5.57e5 → 5.63e6）带来的，
      还是循环本就没有抬分（甚至压分）？
      做法：把不相交 counts 用二项稀释降到原深度，重算 obs 与 rho，重复 20 次。
  (2) 五个零分布条件当初是对旧（循环）实测侧算的，换 obs 后要重挂。
      零分布预测已在 null_tb.csv，重算相关即可，零 API。
"""
import glob

import numpy as np
import pandas as pd
from scipy import stats

D = "b_data"
OLD_LIB = {"ctrl": 557189.0, "stim": 612055.0}
RNG = np.random.default_rng(20261115)


def main():
    de = pd.read_csv(sorted(glob.glob(f"{D}/de_disjoint_*.csv"))[-1])
    pred = pd.read_csv(sorted(glob.glob(f"{D}/b_pred_disjoint_*.csv"))[-1])
    base = pred[(pred[["ctrl", "stim"]].max(axis=1) > 1) & np.isfinite(pred.obs)].copy()
    print(f"不相交 base: {len(base)}；文库 ctrl {de.cnt_ctrl.sum():,.0f} / stim {de.cnt_stim.sum():,.0f}")

    r_obs = stats.spearmanr(base.pred1, base.obs).correlation
    print(f"原深度 rho = {r_obs:+.4f}")

    # --- (1) 二项稀释到旧深度 ---
    sub = de.set_index("ensg").loc[base.ensg]
    fc = OLD_LIB["ctrl"] / de.cnt_ctrl.sum()
    fs = OLD_LIB["stim"] / de.cnt_stim.sum()
    print(f"稀释比例 ctrl {fc:.4f} / stim {fs:.4f}")
    rs = []
    for i in range(20):
        c = RNG.binomial(sub.cnt_ctrl.astype(int).values, fc)
        s = RNG.binomial(sub.cnt_stim.astype(int).values, fs)
        pc = c / c.sum() * 1e6
        ps = s / s.sum() * 1e6
        o = np.log2((ps + 1) / (pc + 1))
        rs.append(stats.spearmanr(base.pred1.values, o).correlation)
    rs = np.array(rs)
    print(f"稀释到旧深度后：rho 均值 {rs.mean():+.4f}（SD {rs.std(ddof=1):.4f}，范围 {rs.min():+.3f}~{rs.max():+.3f}）")
    print(f"  对照：循环实测侧（原数） {0.4860:+.4f}；不相交原深度 {r_obs:+.4f}")

    # --- (2) 五个零分布重挂 ---
    nul = pd.read_csv(f"{D}/null_tb.csv")
    nul = nul[nul.ensg.isin(base.ensg)].copy()
    m = nul.merge(base[["ensg", "obs", "pred1"]], on="ensg", suffixes=("_old", "_new"))
    print(f"\n零分布重挂（{len(m)} 基因）")
    print(f"{'条件':16s} {'旧obs':>8s} {'新obs':>8s} {'中位|pred|':>10s}")
    rows = []
    for c in ["dc", "bcell", "nk", "cd8", "tcell"]:
        col = f"pred_null_{c}"
        old = stats.spearmanr(m[col], m.obs_old).correlation
        new = stats.spearmanr(m[col], m.obs_new).correlation
        rows.append(dict(cond=c, rho_old=old, rho_new=new, mag=m[col].abs().median()))
        print(f"{c:16s} {old:+8.3f} {new:+8.3f} {m[col].abs().median():10.3f}")
    real = stats.spearmanr(m.pred1, m.obs_new).correlation
    print(f"{'真实 IFN-β':16s} {'—':>8s} {real:+8.3f} {m.pred1.abs().median():10.3f}")
    rows.append(dict(cond="real", rho_old=np.nan, rho_new=real, mag=m.pred1.abs().median()))
    pd.DataFrame(rows).to_csv(f"{D}/null_rehang_disjoint.csv", index=False)
    print(f"\n写出 {D}/null_rehang_disjoint.csv")

    # 深度分解的明细也存一份
    pd.DataFrame({"rho_shallow": rs}).to_csv(f"{D}/depth_thinning.csv", index=False)


if __name__ == "__main__":
    main()
