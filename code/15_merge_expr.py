#!/usr/bin/env python3
"""把预测表达量并进长表，产出可直接跑 a_group_check.py 的文件。"""
import numpy as np
import pandas as pd

long = pd.read_csv("tcell_per_variant_long.csv")
expr = pd.read_csv("expr_tcell_clean.csv")
print(f"长表 {len(long)} 行；表达量 {len(expr)} 条")

long = long.drop(columns=["pred_ref_expr"])
e = (expr.groupby(["variant", "gene", "cond_role"])["pred_ref_expr"]
        .mean().reset_index())
m = long.merge(e, on=["variant", "gene", "cond_role"], how="left")
print(f"合并后非空 {m.pred_ref_expr.notna().sum()}/{len(m)}")
m.to_csv("tcell_per_variant_long_expr.csv", index=False)

x = (m.dropna(subset=["pred_ref_expr"])
       .pivot_table(index=["variant", "gene"], columns="cond_role",
                    values="pred_ref_expr", aggfunc="first").dropna())
print(f"两种条件都有表达量的变异-基因：{len(x)}")
if len(x):
    print(f"表达量中位：act {x['act'].median():.4f}｜rest {x['rest'].median():.4f}")
    print(f"两种条件都 >=0.1 的：{((x['act'] >= 0.1) & (x['rest'] >= 0.1)).sum()}")
    print(f"任一条件 >=0.1 的：{((x['act'] >= 0.1) | (x['rest'] >= 0.1)).sum()}")
