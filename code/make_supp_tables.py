#!/usr/bin/env python3
"""Build the three English supplementary tables for submission.

Inputs are the archived working CSVs (Chinese labels in S1/S2); this script holds the
English rendering so the mapping is explicit and reviewable rather than machine-translated.

Outputs (in submission/):
  TableS1_provenance.csv / .xlsx
  TableS2_cells_and_libraries.csv / .xlsx
  TableS3_keyword_counts.csv / .xlsx
  SupplementaryTables.xlsx    (all three as sheets)
"""
import os

import pandas as pd

OUT = "submission"
os.makedirs(OUT, exist_ok=True)

STATUS = {
    "预先规定": "Pre-specified",
    "预先规定（构标签之前）": "Pre-specified (written before labels were built)",
    "审计后新增": "Added after audit",
    "审计后新增（探索性）": "Added after audit (exploratory)",
    "审计后更正": "Corrected after audit",
    "审计后扩展": "Extended after audit",
    "外部复算": "External recomputation",
}

# (item, status, section/hash, time, note) — one row per line of the archived Chinese table
S1 = [
    ("Primary test statistic (log2 ratio A/B; two-sided Mann-Whitney and permutation)",
     "Pre-specified", "PREREG_v2 section 1", "2026-11-14 00:50", ""),
    ("Strict definition of stimulus specificity promoted to primary analysis",
     "Pre-specified", "PREREG_v2 section 3", "2026-11-14 00:50", ""),
    ("Ordering of the four Quach stimuli; Influenza_6h designated the primary arm "
     "(mechanistic reason: influenza strongly induces a type I interferon response)",
     "Pre-specified", "PREREG_v2 section 4", "2026-11-14 00:50",
     "Scoring finished 03:16, after this section; hash 4f91aab9 (01:20) is the next link"),
    ("Operationalisation of the strict definition (three unverifiable classes counted "
     "separately; shared sites require concordant direction)",
     "Pre-specified", "Supplement A, section 7", "2026-11-14 01:20", ""),
    ("Weight of the reversed LPS_2h arm (a reversal across stimuli counts as failure)",
     "Pre-specified", "Supplement B, section 8", "2026-11-14 10:30", ""),
    ("Strict labels matched by molecular trait; validation protocol for the remote client",
     "Pre-specified", "Supplement C, section 9", "2026-11-14 10:50",
     "Supersedes Supplement A section 2 (max |z| across probes)"),
    ("Symmetric application of the match score (Fairfax via E-MTAB-2232); match score "
     "may only downgrade, never upgrade",
     "Pre-specified", "Supplement D, section 10", "2026-11-14 20:15", ""),
    ("Validation criteria changed to sign agreement + Spearman + classification agreement; "
     "|z| < 1.9 sensitivity",
     "Pre-specified (before labels were built)", "Supplement E, section 11", "2026-11-14 22:40", ""),
    ("Visibility gate re-run as min(A,B) rather than max(A,B)",
     "Added after audit", "-", "2026-11-15 morning",
     "Raised by the independent recomputation: the denominator had no lower bound"),
    ("A/B decomposition (denominator-only p = 0.0037; numerator-only p = 0.1442), "
     "primary evidence for methodological result 2",
     "Added after audit", "-", "2026-11-15 morning", ""),
    ("Attrition denominator corrected: 71.5%/75.9% to 57.2%/59.7%",
     "Corrected after audit", "-", "2026-11-15 morning",
     "Numerators were right; denominator was the validation-pair count"),
    ("Bounded index mislabelled: sign(x)|x| - sign(y)|y| is identically x - y, "
     "i.e. the raw effect difference",
     "Corrected after audit", "-", "2026-11-15 morning",
     "The earlier sentence 'p = 0.40 with a bounded index' was reworded"),
    ("Measured response recomputed on cells disjoint from the conditioning input "
     "(256 conditioning cells; 2,529 / 2,334 measured)",
     "Added after audit", "-", "2026-11-15 midday",
     "Previously the measured side was computed from the same cells as the condition "
     "(17_build_B.py, bulk())"),
    ("All-gene hold-out (all 119 target genes removed from the conditioning input)",
     "Added after audit", "-", "2026-11-15 midday", "0.6071 versus 0.6155"),
    ("Twelve-point calibration (including four stimulated conditions)",
     "Added after audit (exploratory)", "-", "2026-11-15 midday",
     "Not pre-registered; flagged as exploratory in the text and figure legend"),
    ("Magnitude-matched null distribution (seven cell-type conditions)",
     "Added after audit", "-", "2026-11-15 early morning",
     "Originally five; FCGR3A+ reclassified from null to a genuine overlap point"),
    ("Split-half null extended from 8 to 50 draws",
     "Extended after audit", "8 draws originally", "2026-11-15 early morning",
     "Not pre-specified; extended at the recomputation's request"),
    ("Withdrawal of the pseudobulk 'no gradient' claim (it is monotonic, 0.891 to 0.928)",
     "Corrected after audit", "-", "2026-11-15 early morning",
     "The original number used the filtered 119 genes, which flattened the gradient. The result is no longer reported in the manuscript (removed on 22 November)"),
    ("Inventory of built-in conditions: 2,799 conditions, none experimentally perturbed",
     "Added after audit", "-", "2026-11-14 night", ""),
    ("Scale check for the dexamethasone condition (near-inert on the log2 scale; "
     "6 of 7 genes correct after conversion to linear scale)",
     "Added after audit", "-", "2026-11-15 early morning", ""),
    ("Attenuation correction of the calibration slope (reliability 0.983; 0.710 to 0.722)",
     "Added after audit", "-", "2026-11-15 midday", ""),
    ("Condition label verification (three layers: annotation, sampled cells, files)",
     "Added after audit", "-", "2026-11-15 afternoon",
     "No API cost; confirms the CD8/B difference is not a labelling error"),
    ("Unstimulated conditions put on the same baseline as the stimulated ones "
     "(all cells; FCGR3A+ 0.268 to 0.336)",
     "Corrected after audit", "-", "2026-11-15 afternoon",
     "The two sets had been computed on different baselines"),
    ("Nine genes missing for five unstimulated conditions filled in "
     "(3 windows x 5 conditions, 15 calls)",
     "Corrected after audit", "-", "2026-11-15 afternoon",
     "52_null_tb.py failed silently; all 15 points now on the same 119 genes"),
    ("Fifteen-point calibration table rebuilt on the common basis "
     "(slope 0.728, corrected 0.741)",
     "Corrected after audit", "-", "2026-11-15 afternoon", "Replaces 0.709 / 0.721"),
    ("Wording corrected for the 2,799 conditions: zero experimental perturbations, "
     "but 30 carry an 'activated' label",
     "Corrected after audit", "-", "2026-11-15 evening",
     "The earlier 'no treated state' was too strong; the full table was exported and searched"),
    ("Full export of the 4,430 reference tracks (zero perturbation matches, activated included)",
     "Added after audit", "-", "2026-11-15 evening",
     "Previously only a 79-row monocyte/prostate subset was available"),
    ("Marker-gene table recomputed on the all-cell baseline (ISG15 0.27 to 0.49, among others)",
     "Corrected after audit", "-", "2026-11-15 evening",
     "The earlier values used the old baseline that excluded conditioning cells"),
    ("Keyword search widened to include the condition and biosample_type fields "
     "(which surfaced the 30 activated labels)",
     "Corrected after audit", "-", "2026-11-15 evening",
     "Previously only biosample_name was searched"),
    ("obs|ctrl and obs|stim columns added to the marker table so log2 IFN response is reproducible",
     "Corrected after audit", "-", "2026-11-15 evening", "Previously max|diff| 0.118. The marker panel is no longer reported in the manuscript (removed on 22 November); the table stays in the deposit"),
    ("Wide keyword list re-checked by the recomputation (five groups: stimulus, "
     "glucocorticoid, hormone, state label, prostate)",
     "External recomputation", "recheck8.py", "2026-11-15 evening",
     "Both the narrow and the wide list appear in Table S3"),
    ("Attrition denominator restored: 57.2%/59.7% to 71.5%/75.9%",
     "Corrected after audit", "Figure 4; section on credible-set attrition",
     "2026-11-22 20:30",
     "The earlier correction had the direction wrong: 57.2%/59.7% divided by all label "
     "rows, including shared ones, and counted untested rows (20/455) as excluded. "
     "Over condition-specific rows the shares are 1,077/1,506 and 2,754/3,628"),
    ("Main correlation printed as 0.615, not 0.616",
     "Corrected after audit", "Abstract, Results, Figure 1 legend, supporting analyses",
     "2026-11-22 20:30",
     "Full-precision value 0.6154985596768819 rounds to 0.615; the earlier 0.616 came "
     "from rounding twice, via the four-decimal 0.6155"),
    ("Robustness analysis added: dropping the ten genes with the largest absolute measured log2 fold change",
     "Added after audit", "Results", "2026-11-23 13:05",
     "Window bootstrap, 35 windows, 2,000 draws, seed 0; rho = 0.563 (n = 109, 95% CI 0.419 to 0.674); computed in audit_numbers.py"),
    ("Pseudocount sensitivity values",
     "Added after audit", "Supporting analyses; Methods", "2026-11-23 13:05",
     "Computed in audit_numbers.py from b_pred_disjoint_1115_115447.csv (stim, ctrl, cpm_stim, cpm_ctrl)"),
    ("Gene count revised: the figure of 310 could not be traced to an archived file",
     "Corrected after audit", "Results; Methods", "2026-11-23 10:40",
     "Replaced by the reproducible funnel: genes returned in the 40 intervals, 246 with a measured "
     "value, 119 passing the activity filter"),
    ("Split-half null recomputed on the disjoint design (0 of 50, p = 0.0196, mean 0.027, SD 0.229)",
     "Corrected after audit", "Supporting analyses; Methods", "2026-11-15 23:35",
     "Recomputed by 53_noise_splits.py (seeds 1000+s, 50 splits) after the observed value moved from "
     "0.486 to 0.615; the v2 preprint still reports the earlier 1 of 50, p = 0.039"),
    ("Lymphoblastoid tracks: the GM prefix does not identify EBV-transformed lines",
     "Corrected after audit", "Discussion; Table S3", "2026-11-22 22:45",
     "GM23338 is an iPSC line derived from the fibroblast line GM23248, and GM23248 is a fibroblast "
     "line; the earlier statement that all 246 GM-prefixed tracks were lymphoblastoid was wrong. "
     "The text now counts only GM12878 with its parental lines GM12891 and GM12892 (131 tracks)"),
    ("Genetically modified tracks: 694 CRISPR insertions plus 40 stable transfections",
     "Added after audit", "Methods, keyword searches", "2026-11-22 22:45",
     "All 734 are epitope-tag constructs; none is a knockdown or knockout"),
    ("Activated-label counts revised to 29 first-generation tracks and 135 second-generation conditions",
     "Corrected after audit", "Table 1; Table S3", "2026-11-22 22:40",
     "The earlier 44 and 136 came from a script that was not archived and could not be reproduced. "
     "The counts now reported come from the full returned tables with the lists given in Table S3"),
    ("Note on the count of 30 corrected: it comes from all text fields with activated|activation, "
     "not from a biosample_name-only search",
     "Corrected after audit", "Table S3", "2026-11-22 22:40", "The earlier note was wrong"),
    ("Supplement E describes itself as written before label construction; the file was hashed at 22:40, "
     "after the labels were built at 22:27",
     "Corrected after audit", "Methods, validation", "2026-11-22 22:40",
     "Recorded, not resolved: the text of Supplement E was not edited"),
    ("Activated-label counts: 44 first-generation tracks and 136 second-generation conditions",
     "Corrected after audit", "Table 1; Table S3", "2026-11-22 20:45",
     "The 30 recorded above is what a biosample_name-only search returned. The final counts "
     "come from all text fields with the wide list; the narrow list returns 0 perturbed "
     "conditions in either generation"),
    ("Stimulated CD8 T measured overlap printed as 0.212, not 0.213",
     "Corrected after audit", "Results; Figure 2 legend", "2026-11-22 20:30",
     "Full-precision value 0.2124651037403786 rounds to 0.212, same double rounding"),
]

s1 = pd.DataFrame(S1, columns=["Item", "Status", "Section or hash", "Date and time", "Note"])


def plan_status(item, sec, when):
    """Where each analysis stood relative to the versions of the plan."""
    low = (item + " " + sec).lower()
    if any(k in low for k in ("calibration", "expression-level", "condition vector",
                              "cell counts", "reliability", "spearman")):
        return "Analysed before plan v1"
    if "prereg_v2" in low:
        return "Specified in plan v2"
    if "supplement a" in low:
        return "Specified in plan v2 (supplement A)"
    if any(k in low for k in ("supplement b", "supplement c", "supplement d", "supplement e")):
        return "After plan v2"
    if "prereg_fairfax" in low or "section 1" in low and "v1" in low:
        return "Specified in plan v1"
    return "After plan v2"


# Statuses assigned per row: a correction inherits the status of the analysis it
# corrects, rather than the date on which the correction was made.
PLAN_STATUS = ['Specified in plan v2', 'Specified in plan v2', 'Specified in plan v2', 'Specified in plan v2 (supplement A)', 'After plan v2', 'After plan v2', 'After plan v2', 'After plan v2', 'After plan v2', 'After plan v2', 'Specified in plan v2', 'After plan v2', 'After plan v2', 'Run before plan v1', 'After plan v2', 'After plan v2', 'Run before plan v1', 'Run before plan v1', 'After plan v2', 'After plan v2', 'After plan v2', 'After plan v2', 'After plan v2', 'After plan v2', 'After plan v2', 'After plan v2', 'After plan v2', 'After plan v2', 'After plan v2', 'Run before plan v1', 'After plan v2', 'After plan v2', 'Run before plan v1', 'After plan v2', 'Run before plan v1', 'After plan v2', 'Run before plan v1', 'After plan v2', 'After plan v2', 'After plan v2', 'After plan v2', 'After plan v2', 'After plan v2', 'After plan v2']

# the version-1 main test, superseded in version 2 but still reported
s1.loc[len(s1)] = ["Main test as specified in version 1 (mean difference of log ratios; "
                   "permutation within study; 1,035 sites)",
                   "Specified in plan v1", "PREREG_fairfax_main, section 1",
                   "2026-11-13 23:19",
                   "Replaced in version 2 by the rank sum; still reported (0.0383 by permutation)"]
PLAN_STATUS = PLAN_STATUS + ["Specified in plan v1"]

s1.insert(1, "Plan status", PLAN_STATUS)

# ---- S2: cells and libraries -------------------------------------------------
ROLE = {"x 轴（CPM）": "x axis (CPM)", "y 轴（obs）": "y axis (observed response)"}
NOTE = {
    "该类型全部 singlet，不排除条件侧": "all singlet cells of this type; conditioning cells not excluded",
    "排除本文库 256 个条件侧细胞": "excluding the 256 conditioning cells from this library",
}
s2 = pd.read_csv("b_data/for_audit_cellcounts.csv")
s2["role"] = s2["role"].map(ROLE).fillna(s2["role"])
s2["note"] = s2["note"].map(NOTE).fillna(s2["note"])
s2.columns = ["Column", "Role", "Cells used", "Library (counts)", "UMI per cell",
              "Cells excluded", "Note"]

# ---- S3: keyword counts ------------------------------------------------------
s3 = pd.read_csv("table1_keyword_counts.csv")
NAME = {"AlphaGenome 1 tracks": "AlphaGenome 1 (reference tracks)",
        "AlphaGenome 2 conditions": "AlphaGenome 2 (built-in conditions)"}
GROUP = {"glucocorticoid": "glucocorticoid",
         "cytokine_or_stimulus": "stimulus / cytokine",
         "stimulus_cytokine": "stimulus / cytokine",
         "activated_state": "activated state label",
         "hormone_treated": "hormone",
         "hormone": "hormone",
         "prostate": "prostate"}
s3["dataset"] = s3["dataset"].map(NAME).fillna(s3["dataset"])
s3["keyword_set"] = s3["keyword_set"].replace({"narrow (ours)": "narrow (this study)",
                                               "broad (auditor)": "broad (independent)"})
s3["group"] = s3["group"].map(GROUP).fillna(s3["group"])
s3["Hits adjudicated as a real perturbation"] = 0
order = {"stimulus / cytokine": 0, "glucocorticoid": 1, "hormone": 2,
         "activated state label": 3, "prostate": 4}
s3 = (s3.assign(_o=s3["group"].map(order))
        .sort_values(["dataset", "keyword_set", "_o"], ascending=[True, True, True])
        .drop(columns="_o"))
s3.columns = ["Dataset", "Items in table", "Keyword list", "Keyword group", "Hits",
              "Hits adjudicated as a real perturbation"]

for name, df in [("TableS1_provenance", s1), ("TableS2_cells_and_libraries", s2),
                 ("TableS3_keyword_counts", s3)]:
    df.to_csv(f"{OUT}/{name}.csv", index=False)
    df.to_excel(f"{OUT}/{name}.xlsx", index=False)

# display copy with thousands separators (used by the PDF; the CSV above stays numeric)
s2d = s2.copy()
s2d["Library (counts)"] = s2d["Library (counts)"].map(lambda v: f"{int(v):,}")
s2d["Cells used"] = s2d["Cells used"].map(lambda v: f"{int(v):,}")
s2d["UMI per cell"] = s2d["UMI per cell"].map(lambda v: f"{v:,.1f}")
s2d["Cells excluded"] = s2d["Cells excluded"].map(lambda v: f"{int(v):,}")
s2d.to_csv(f"{OUT}/TableS2_display.csv", index=False)
print("wrote TableS2_display.csv")

# display copy of S3 for the supplementary PDF (short labels, formatted counts)
s3d = pd.DataFrame({
    "Dataset": s3["Dataset"].replace({"AlphaGenome 1 (reference tracks)": "AlphaGenome 1",
                                      "AlphaGenome 2 (built-in conditions)": "AlphaGenome 2"}),
    "Items": s3["Items in table"].map(lambda v: f"{v:,}"),
    "Keyword list": s3["Keyword list"].replace({"narrow (this study)": "narrow",
                                                "broad (independent)": "broad"}),
    "Keyword group": s3["Keyword group"],
    "Hits": s3["Hits"],
    "Hits adjudicated as a real perturbation":
        s3["Hits adjudicated as a real perturbation"],
})
s3d.to_csv(f"{OUT}/TableS3_display.csv", index=False)
print("wrote TableS3_display.csv")

with pd.ExcelWriter(f"{OUT}/SupplementaryTables.xlsx") as w:
    s1.to_excel(w, sheet_name="S1 provenance", index=False)
    s2.to_excel(w, sheet_name="S2 cells libraries", index=False)
    s3.to_excel(w, sheet_name="S3 keyword counts", index=False)

print(f"S1 {len(s1)} rows, S2 {len(s2)}, S3 {len(s3)}")
print(s2.to_string(index=False, max_colwidth=42))
print()
print(s3.head(6).to_string(index=False))
