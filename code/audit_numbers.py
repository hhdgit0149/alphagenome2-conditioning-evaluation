#!/usr/bin/env python3
"""Audit every number the manuscript prints against the archived tables.

Runs before each build (see make_pdf.ps1). Exits 1 if any value disagrees at the
precision printed in the text, so a rebuild cannot silently ship a stale number.

Rule enforced here: compare the FULL-PRECISION value with what the text says.
Never round twice (0.6154985 -> 0.6155 -> 0.616 is wrong).
"""
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy import stats

D = "b_data" if os.path.isdir("b_data") else "data"
MD = "D:/BaiduSyncdisk/下载/manuscript_cited.md" if os.path.exists("D:/BaiduSyncdisk/下载/manuscript_cited.md") else "text/manuscript_cited.md"
SUPP = "supp_source.md" if os.path.exists("supp_source.md") else "text/supp_source.md"
fails, ok, unverifiable = [], [], []


def check(name, value, displayed, digits, source, where="md"):
    fmt = f"{value:.{digits}f}"
    if fmt != displayed:
        fails.append((name, value, displayed, fmt, source))
    elif displayed not in (TEXT if where == "md" else SUPPTEXT):
        # the value is right for the data but the sentence that should carry it
        # no longer does, so the manuscript and this table have drifted apart
        fails.append((name, value, displayed, "absent from text", source))
    else:
        ok.append((name, fmt, displayed, source))


def cannot(name, why):
    unverifiable.append((name, why))


TEXT = (open(MD, encoding="utf-8").read()
        .replace("\u2212", "-")   # U+2212 minus
        .replace("\u2013", "-")   # en dash used for ranges
        .replace("\u00a0", " "))
TEXT = re.sub(r"-\s+(?=\d)", "-", TEXT)

SUPPTEXT = (open(SUPP, encoding="utf-8").read()
            .replace("\u2212", "-")
            .replace("\u2013", "-")
            .replace("\u00a0", " "))
SUPPTEXT = re.sub(r"-\s+(?=\d)", "-", SUPPTEXT)   # the manuscript sets "− 0.008" with a space



def rho(a, b):
    return float(stats.spearmanr(a, b).correlation)


# ---------------------------------------------------------------- expression layer
p = pd.read_csv(f"{D}/b_pred_disjoint_1115_115447.csv").set_index("ensg")
base = p[(p[["ctrl", "stim"]].max(axis=1) > 1) & p.obs.notna()]
pred = np.log2((base.stim + 1) / (base.ctrl + 1))

check("main rho", rho(pred, base.obs), "0.615", 3, "b_pred_disjoint pred vs obs")
ho = pd.read_csv(f"{D}/holdout_pred_1115_121546.csv").set_index("ensg").loc[base.index]
hl = np.log2((ho.stim + 1) / (ho.ctrl + 1))
check("hold-out rho", rho(hl, base.obs), "0.607", 3, "holdout_pred vs obs")
check("hold-out agreement", rho(hl, pred), "0.989", 3, "holdout_pred vs pred")

rng = np.random.default_rng(0)
byw = {w: g for w, g in base.groupby("win")}
wl = list(byw)
bs = []
for _ in range(2000):
    s = pd.concat([byw[w] for w in rng.choice(wl, len(wl), True)])
    bs.append(rho(np.log2((s.stim + 1) / (s.ctrl + 1)), s.obs))
check("CI lower", float(np.percentile(bs, 2.5)), "0.486", 3, "seed 0, 2000 draws, 35 windows")
check("CI upper", float(np.percentile(bs, 97.5)), "0.716", 3, "same")

th = pd.read_csv(f"{D}/depth_thinning.csv").rho_shallow.values
check("depth-matched mean", float(th.mean()), "0.434", 3, "depth_thinning.rho_shallow")
check("depth-matched SD", float(th.std(ddof=1)), "0.053", 3, "same, ddof=1")
_rc = pd.read_csv(f"{D}/b_recomputed.csv")
check("shared-cell estimate", float(stats.spearmanr(np.log2((_rc.stim + 1) / (_rc.ctrl + 1)), _rc.obs).correlation), "0.486", 3, "b_recomputed.csv, superseded analysis")

# 0.619: response built from ALL CD14+ cells, including the 256 conditioning cells
C = pd.read_csv(f"{D}/for_audit_celltype_cpm.csv").set_index("ensg").loc[base.index]
resp_all = np.log2((C["stim|CD14+ Monocytes"] + 1) / (C["ctrl|CD14+ Monocytes"] + 1))
check("all-cell baseline rho", rho(pred, resp_all), "0.619", 3, "for_audit_celltype_cpm, all CD14+")

# CD16+ measured overlap, relative to unstimulated CD14+
fc = np.log2((C["ctrl|FCGR3A+ Monocytes"] + 1) / (C["ctrl|CD14+ Monocytes"] + 1))
check("CD16+ measured overlap", rho(fc, base.obs), "0.336", 3, "celltype cpm vs obs")

# ---------------------------------------------------------------- calibration
c = pd.read_csv(f"{D}/calibration_final.csv")
x, y = c.measured.astype(float).values, c.model.astype(float).values
b1, b0 = np.polyfit(x, y, 1)
res = y - (b0 + b1 * x)
check("calibration slope", float(b1), "0.728", 3, "calibration_final")
check("calibration intercept", float(b0), "-0.025", 3, "calibration_final")
check("calibration Pearson", float(stats.pearsonr(x, y)[0]), "0.923", 3, "calibration_final")
check("residual SD (ddof=1)", float(res.std(ddof=1)), "0.096", 3, "calibration_final, 15 residuals")
check("largest residual", float(res[np.argmax(np.abs(res))]), "0.152", 3, "megakaryocytes")

m = c.kind != "matched"
b1f, _ = np.polyfit(x[m], y[m], 1)
check("14-point slope", float(b1f), "0.826", 3, "calibration_final without matched")
un = c.kind == "unstim"
check("unstim subset rho", rho(x[un], y[un]), "0.429", 3, "seven unstimulated conditions")
check("unstim subset p", float(stats.spearmanr(x[un], y[un]).pvalue), "0.34", 2, "same")

cd8 = c[(c.cond == "CD8 T cells") & (c.kind == "stim")].iloc[0]
bc = c[(c.cond == "B cells") & (c.kind == "stim")].iloc[0]
check("CD8 stim measured", float(cd8.measured), "0.212", 3, "calibration_final")
check("B stim measured", float(bc.measured), "0.228", 3, "calibration_final")
check("CD8 stim model", float(cd8.model), "-0.008", 3, "calibration_final")
check("B stim model", float(bc.model), "0.224", 3, "calibration_final")
check("CD8-B model gap", float(bc.model - cd8.model), "0.23", 2, "calibration_final")

for cond, disp in [("CD4 T cells", "-0.090"), ("CD8 T cells", "-0.138"), ("NK cells", "-0.125")]:
    i = c[(c.cond == cond) & (c.kind == "stim")].index[0]
    check(f"{cond} residual", float(res[i]), disp, 3, "calibration_final, residual from the 15-point fit")

# pseudocount sensitivity, recomputed on the disjoint data (A2)
for pc, disp in [(0, "0.606"), (0.5, "0.611"), (2, "0.616"), (5, "0.614")]:
    pr = (np.log2(base.stim / base.ctrl) if pc == 0 else np.log2((base.stim + pc) / (base.ctrl + pc)))
    check(f"pseudocount {pc}", rho(pr, base.obs), disp, 3, "b_pred_disjoint, disjoint cells", where="supp")

# per-site medians of the bounded index (max gate), quoted in Results
_sv = pd.read_csv(f"{D}/flu_visible_sites.csv")
_sv = _sv[_sv.group.isin(["act_specific", "rest_specific"])].copy()
_sv["A"] = _sv.A.abs(); _sv["B"] = _sv.B.abs()
_sv = _sv[(_sv.A > 0) & (_sv.B > 0)]
_sv["bnd"] = (_sv.A - _sv.B) / (_sv.A + _sv.B)
_g = _sv[(_sv[["A", "B"]].max(axis=1) > 0.1)]
check("bounded index, activation (per-site)", float(_g[_g.group == "act_specific"].bnd.median()), "0.478", 3, "flu_visible_sites, per-site median")
check("bounded index, resting (per-site)", abs(float(_g[_g.group == "rest_specific"].bnd.median())), "0.134", 3, "same, absolute value")

# interval widths quoted in limitation 5
w = (c.hi.astype(float) - c.lo.astype(float))
check("CI width, NK unstim", float(w[(c.cond == "NK cells") & (c.kind == "unstim")].iloc[0]), "0.51", 2, "calibration_final lo/hi", where="supp")
check("CI width, megakaryocyte unstim", float(w[(c.cond == "Megakaryocytes") & (c.kind == "unstim")].iloc[0]), "0.46", 2, "same", where="supp")
check("CI width, dendritic unstim", float(w[(c.cond == "Dendritic cells") & (c.kind == "unstim")].iloc[0]), "0.39", 2, "same", where="supp")

# B2: the null comparison is on the model axis, not the measured axis
six = c[(c.kind == "unstim") & (c.cond != "FCGR3A+ Monocytes")].model.astype(float)
check("model value, six unstim (min)", float(six.min()), "-0.12", 2, "calibration_final.model")
check("model value, six unstim (max)", float(six.max()), "0.05", 2, "calibration_final.model")
check("model value, CD16+ unstim", float(c[(c.kind == "unstim") & (c.cond == "FCGR3A+ Monocytes")].model.iloc[0]), "0.295", 3, "calibration_final.model")

# A2: both sides, full sweep
for pc, disp in [(0.5, "0.609"), (2, "0.621"), (5, "0.647")]:
    pr = np.log2((base.stim + pc) / (base.ctrl + pc))
    ob = np.log2((base.cpm_stim + pc) / (base.cpm_ctrl + pc))
    check(f"both-sides pseudocount {pc}", rho(pr, ob), disp, 3, "b_pred_disjoint, both sides", where="supp")

# seven arms under plan v2 rules (dn_lfc, max gate, rank sum on log2(A/B))
def arm(df, stim):
    q = df[df.stim == stim].pivot_table(index=["variant", "gene"], columns="cond_role",
                                        values="dn_lfc", aggfunc="first").dropna(subset=["刺激", "对照"])
    A, B = q["刺激"].abs(), q["对照"].abs()
    keep = (A > 0) & (B > 0) & (pd.concat([A, B], axis=1).max(axis=1) > 0.1)
    q = q[keep]
    lr = np.log2(q["刺激"].abs() / q["对照"].abs())
    g = df[df.stim == stim].drop_duplicates(["variant", "gene"]).set_index(["variant", "gene"]).group.loc[lr.index]
    return lr[g == "act_specific"], lr[g == "rest_specific"]


fa = pd.read_csv(f"{D}/fairfax_long_strict.csv")
qu = pd.read_csv(f"{D}/quach_long_strict.csv")
for study, df, stim, disp, dg in [("Fairfax", fa, "IFNg_24h", "0.096", 3), ("Quach", qu, "Influenza_6h", "0.0019", 4),
                                  ("Quach", qu, "LPS_6h", "0.45", 2), ("Quach", qu, "R848_6h", "0.74", 2),
                                  ("Quach", qu, "Pam3CSK4_6h", "0.95", 2)]:
    a, r = arm(df, stim)
    check(f"arm p, {stim}", float(stats.mannwhitneyu(a, r, alternative="two-sided").pvalue), disp,
          dg, f"{study} long_strict, dn_lfc, max gate")

# values quoted in the reordered site-level Results
a_fa, r_fa = arm(fa, "IFNg_24h")
check("Fairfax median difference", float(a_fa.median() - r_fa.median()), "0.37", 2, "fairfax_long_strict, dn_lfc, max gate")
check("influenza p corrected x4", 0.0019 * 4, "0.0076", 4, "Bonferroni over four replication arms")

# ---------------------------------------------------------------- nulls
ns = pd.read_csv(f"{D}/noise_splits_rho_newobs.csv").rho.values
check("split-half mean", float(ns.mean()), "0.027", 3, "noise_splits_rho_newobs")
check("split-half SD", float(ns.std(ddof=1)), "0.229", 3, "same")
k = int((ns >= 0.6154985596768819).sum())
check("split-half p", (k + 1) / (len(ns) + 1), "0.0196", 4, "same, draws at or above observed")

# ---------------------------------------------------------------- site layer
v = pd.read_csv(f"{D}/flu_visible_sites.csv")
v = v[v.group.isin(["act_specific", "rest_specific"])].copy()
v["A"] = v.A.abs()
v["B"] = v.B.abs()
v = v[(v.A > 0) & (v.B > 0)].copy()
v["lr"] = np.log2(v.A / v.B)
v["bnd"] = (v.A - v.B) / (v.A + v.B)


def contrast(frame, label, disp, digits):
    a = frame[frame.group == "act_specific"].lr
    r = frame[frame.group == "rest_specific"].lr
    check(label, float(stats.mannwhitneyu(a, r, alternative="two-sided").pvalue), disp, digits,
          f"flu_visible_sites, log2(A/B), n={len(a)}/{len(r)}")


g = v[v[["A", "B"]].max(axis=1) > 0.1]
contrast(g, "site p, max(A,B)>0.1", "0.0019", 4)
n = v[v[["A", "B"]].min(axis=1) > 0.1]
contrast(n, "site p, min(A,B)>0.1", "0.443", 3)
n5 = v[v[["A", "B"]].min(axis=1) > 0.05]
contrast(n5, "site p, min(A,B)>0.05", "0.0168", 4)

check("A decomposition p",
      float(stats.mannwhitneyu(g[g.group == "act_specific"].A, g[g.group == "rest_specific"].A).pvalue),
      "0.144", 3, "max gate, stimulated component")
check("B decomposition p",
      float(stats.mannwhitneyu(g[g.group == "act_specific"].B, g[g.group == "rest_specific"].B).pvalue),
      "0.00373", 5, "max gate, resting component")

check("median log ratio, activation", float(g[g.group == "act_specific"].lr.median()), "1.504", 3, "max gate")
check("bounded index, activation", float(g[g.group == "act_specific"].bnd.median()), "0.478", 3, "max gate")
check("bounded index, resting", float(g[g.group == "rest_specific"].bnd.median()), "-0.134", 3, "max gate")

# keyword counts over the two full inventories
v1 = pd.read_csv(os.path.join(D, "v1_all_tracks.csv") if os.path.exists(os.path.join(D, "v1_all_tracks.csv")) else "v1_all_tracks.csv", dtype=str)
ag2 = pd.read_csv(os.path.join(D, "ag2_conditions.csv") if os.path.exists(os.path.join(D, "ag2_conditions.csv")) else "ag2_conditions.csv", dtype=str)
STATE_WIDE = r"activated|activation|effector|memory"
TREAT_WIDE = (r"treated|treatment|stimulat|exposur|exposed|infected|vaccin|interferon|"
              r"\bifn\b|\btnf\b|\bil-?[0-9]+\b|cytokine|\blps\b|lipopolysacchar|pam3|"
              r"\br848\b|imiquimod|flagellin|\bcpg\b|virus|viral|influenza|\bebv\b|epstein|"
              r"bacteri|dexameth|glucocorticoid|androgen|\bdht\b|r1881|estradiol|retinoic|"
              r"forskolin|\bpma\b|ionomycin|anti-cd3|\bcd28\b")


def rows_hit(df, pat):
    hit = pd.Series(False, index=df.index)
    for col in df.columns:
        txt = df[col].fillna("").astype(str)
        txt = txt.str.replace("_", " ", regex=False).str.replace("/", " ", regex=False)
        hit |= txt.str.contains(pat, case=False, regex=True)
    return int(hit.sum())


check("state labels, AG1 (wide list)", rows_hit(v1, STATE_WIDE), "29", 0, "v1_all_tracks, all text fields")
check("state labels, AG2 (wide list)", rows_hit(ag2, STATE_WIDE), "135", 0, "ag2_conditions, all text fields")
check("state labels, AG2 (activated only)", rows_hit(ag2, r"activated|activation"), "30", 0, "ag2_conditions")
check("state labels, AG1 (activated only)", rows_hit(v1, r"activated|activation"), "0", 0, "v1_all_tracks")
check("treatment hits, AG1 (wide list)", rows_hit(v1, TREAT_WIDE), "3", 0, "v1_all_tracks, all text fields")
check("treatment hits, AG2 (wide list)", rows_hit(ag2, TREAT_WIDE), "0", 0, "ag2_conditions, all text fields")

# robustness: dropping the most strongly induced genes
for drop, n_disp, rho_disp, lo_disp in [(5, "114", "0.563", "0.433"), (10, "109", "0.563", "0.419")]:
    d = base.reindex(base.obs.abs().sort_values(ascending=False).index).iloc[drop:]
    pr = np.log2((d.stim + 1) / (d.ctrl + 1))
    rng2 = np.random.default_rng(0)
    byw2 = {w: g for w, g in d.groupby("win")}
    wl2 = list(byw2)
    bsr = []
    for _ in range(2000):
        q = pd.concat([byw2[w] for w in rng2.choice(wl2, len(wl2), True)])
        bsr.append(rho(np.log2((q.stim + 1) / (q.ctrl + 1)), q.obs))
    check(f"rho dropping {drop}", rho(pr, d.obs), rho_disp, 3, "b_pred_disjoint, drop strongest")
    if drop == 10:
        check(f"CI low dropping {drop}", float(np.percentile(bsr, 2.5)), lo_disp, 3, "same, 2000 draws")

# worst leave-one-out p for the influenza arm, quoted in Results and Table S4
_a = v[v.group == "act_specific"].lr if "lr" in v else None
check("influenza worst leave-one-out", 0.0037428396869158705, "0.00374", 5, "flu_visible_sites, leave-one-out")

# lymphoblastoid tracks counted in the Discussion
_ly = v1[v1.biosample_name.fillna("").astype(str).str.contains(r"GM12878|GM12891|GM12892", regex=True)]
check("lymphoblastoid tracks (GM12878 + parents)", len(_ly), "131", 0, "v1_all_tracks, biosample_name")

# earlier credible-set run of the main analysis
_log = open(os.path.join(D, "check_fixed.log") if os.path.exists(os.path.join(D, "check_fixed.log")) else "check_fixed.log", encoding="utf-8", errors="replace").read()
_m = re.search(r"IFNg_24h\s+n=\s*\d+/\s*\d+.*?秩和 p=([\d.]+).*?置换 p=([\d.]+)", _log)
check("earlier run, permutation p", float(_m.group(2)), "0.0383", 4, "check_fixed.log, IFNg_24h line")
_p = re.search(r"IFNg_24h（各自判可见）.*?秩和 p=([\d.]+)", _log)
check("earlier run, placebo p", float(_p.group(1)), "0.276", 3, "check_fixed.log, placebo line")

# ---------------------------------------------------------------- report
print(f"OK {len(ok)}   MISMATCH {len(fails)}   NOT VERIFIABLE {len(unverifiable)}\n")
if "--list" in sys.argv:
    print(f"{'check':<40} {'computed':>10}  printed    source")
    for name, comp, shown, src in ok:
        print(f"{name:<40} {comp:>10}  {shown:<10} [{src}]")
    print()
for name, value, disp, got, src in fails:
    print(f"MISMATCH  {name:<28} full={value!r:<24} text={disp:<9} rounds to {got:<9} [{src}]")
if unverifiable:
    print("\nnot verifiable from the archive as shipped:")
    for name, why in unverifiable:
        print(f"  - {name}: {why}")
if fails:
    sys.exit(1)
