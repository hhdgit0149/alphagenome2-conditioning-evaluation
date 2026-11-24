#!/usr/bin/env python3
"""比较模型对「加雄激素 16h」和「溶剂对照」两种等位不平衡的预测能力。

同一批位点，同一套模型预测，只在实验条件上不同：
  effect      H3K27ac + 雄激素 16h 的 ALT 比例偏移
  effect_vh   H3K27ac + 溶剂对照的 ALT 比例偏移
如果模型真的抓到了雄激素依赖的调控，16h 的一致率应该更高。
"""
import numpy as np
import pandas as pd
from scipy import stats

V, S = "variants_lncap.csv", "scores_lncap.csv"
OUT = "summary_lncap.csv"

v = pd.read_csv(V, dtype={"variant_id": str})
s = pd.read_csv(S, dtype={"input_id": str}, low_memory=False)
df = s.merge(v[["variant_id", "effect", "effect_vh", "effect_at", "significant"]],
             left_on="input_id", right_on="variant_id", suffixes=("_ag", ""))
print(f"位点 {df['input_id'].nunique()}，评分行 {len(df)}")

case = df[df.significant == 1]
rows = []
groups = ["output_type", "biosample_name", "transcription_factor", "histone_mark"]
for keys, g in case.groupby(groups, dropna=False):
    a = g[(g.effect != 0) & (g.raw_score != 0)]
    b = a[a.effect_vh.notna() & (a.effect_vh != 0)]
    if len(a) < 50 or len(b) < 50:
        continue
    ok16 = np.sign(a.raw_score) == np.sign(a.effect)
    okvh = np.sign(b.raw_score) == np.sign(b.effect_vh)
    rec = dict(zip(groups, keys))
    rec.update(n=len(a), n_vh=len(b),
               agree_16h=round(float(ok16.mean()), 3),
               agree_vh=round(float(okvh.mean()), 3),
               rho_16h=round(float(stats.spearmanr(a.raw_score, a.effect).correlation), 3),
               rho_vh=round(float(stats.spearmanr(b.raw_score, b.effect_vh).correlation), 3))
    # 用交集做配对 McNemar
    idx = a.set_index("input_id").index.intersection(b.set_index("input_id").index)
    a2 = a.set_index("input_id").loc[idx]
    b2 = b.set_index("input_id").loc[idx]
    m16 = np.sign(a2.raw_score) == np.sign(a2.effect)
    mvh = np.sign(b2.raw_score) == np.sign(b2.effect_vh)
    rec["only_16h"] = int((m16 & ~mvh).sum())
    rec["only_vh"] = int((~m16 & mvh).sum())
    rec["mcnemar_p"] = round(float(stats.binomtest(rec["only_16h"],
                                                   rec["only_16h"] + rec["only_vh"],
                                                   0.5).pvalue), 4) \
        if rec["only_16h"] + rec["only_vh"] > 0 else np.nan
    # 只保留两条件差别大的位点（雄激素特异的失衡）
    d = a[a.effect_vh.notna()].copy()
    d["delta"] = (d.effect - d.effect_vh).abs()
    disc = d[d.delta >= 0.10]
    if len(disc) >= 30:
        rec["n_disc"] = len(disc)
        rec["agree_16h_disc"] = round(float((np.sign(disc.raw_score) ==
                                             np.sign(disc.effect)).mean()), 3)
        rec["agree_vh_disc"] = round(float((np.sign(disc.raw_score) ==
                                            np.sign(disc.effect_vh)).mean()), 3)
    # ATAC + 雄激素
    at = a[a.effect_at.notna() & (a.effect_at != 0)]
    if len(at) >= 30:
        rec["n_at"] = len(at)
        rec["agree_at"] = round(float((np.sign(at.raw_score) ==
                                       np.sign(at.effect_at)).mean()), 3)
    rows.append(rec)

res = pd.DataFrame(rows).sort_values("agree_16h", ascending=False)
res.to_csv(OUT, index=False)
cols = ["output_type", "biosample_name", "transcription_factor", "histone_mark", "n",
        "agree_16h", "agree_vh", "rho_16h", "rho_vh", "only_16h", "only_vh",
        "mcnemar_p", "n_disc", "agree_16h_disc", "agree_vh_disc", "n_at", "agree_at"]
print(res[[c for c in cols if c in res.columns]].head(18).to_string(index=False))

print("\n=== 全轨道汇总（跨轨道平均） ===")
print(f"平均一致率 16h={res.agree_16h.mean():.3f}  vs  vehicle={res.agree_vh.mean():.3f}")
print(f"16h 占优的轨道 {(res.agree_16h > res.agree_vh).sum()}/{len(res)}")
print(f"\n写出 {OUT}")
