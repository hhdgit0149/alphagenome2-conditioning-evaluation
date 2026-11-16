#!/usr/bin/env python3
"""Table 1：两套关键词表 × 两代模型的命中数（窄=我们，宽=复算方 recheck8.py）。"""
import re

import pandas as pd

NARROW = {
    "glucocorticoid": r"glucocort|dexameth|hydrocortis|methylprednis|cortisol|prednisol|steroid",
    "cytokine_or_stimulus": r"ifn\b|interferon|\blps\b|pam3|r848|influenza|cytokine|\btnf\b|"
                            r"\bil-?[0-9]|dex\b|treated|treatment|stimulat",
    "activated_state": r"activated|activation",
    "hormone_treated": r"dihydrotestosterone|\bdht\b|testosteron|estradiol|androgen|tamoxifen|"
                       r"estrogen|progest",
    "prostate": r"prostate|lncap|pc-?3|vcapp?|du145",
}
BROAD = {
    "stimulus_cytokine": ["ifn", "interferon", "ifna", "ifnb", "ifng", "lps", "lipopolysaccharide",
                          "pam3", "pam3csk4", "r848", "resiquimod", "poly(i:c)", "polyic", "poly i:c",
                          "influenza", "flu ", "virus", "viral", "infect", "sars", "cmv", "ebv",
                          "tnf", "il-1", "il-2", "il-4", "il-6", "il-10", "il-13", "interleukin",
                          "tgf", "gm-csf", "m-csf", "cytokine", "chemokine", "treat", "stimul",
                          "challenge", "exposed", "exposure", "induc", "agonist", "antagonist",
                          "inhibitor", "drug", "compound"],
    "glucocorticoid": ["glucocorticoid", "dexamethasone", " dex ", "dex-", "cortisol",
                       "hydrocortisone", "prednis", "steroid", "corticoster"],
    "hormone": ["dht", "dihydrotestosterone", "testosterone", "androgen", "estradiol", "estrogen",
                "oestrogen", "progesterone", "hormone", "r1881", "enzalutamide", "bicalutamide",
                "flutamide"],
    "activated_state": ["activated", "activation", "stimulated", "primed", "polarized", "polarised",
                        "m1 ", "m2 ", "th1", "th2", "th17", "effector", "memory", "exhaust"],
    "prostate": ["prostate", "lncap", "vcap", "22rv1", "c4-2", "du145", "pc-3", "pc3", "rwpe"],
}


def text_of(df):
    cols = [c for c in ["biosample_name", "biosample_type", "condition", "ontology_curie",
                        "transcription_factor", "histone_mark", "Assay title", "assay_title",
                        "name", "data_source", "gtex_tissue"] if c in df.columns]
    sub = df[cols].fillna("").astype(str)
    return sub.apply(lambda r: " | ".join(map(str, r.values)), axis=1)


def counts(df, label):
    txt = text_of(df)
    rows = []
    for name, pat in NARROW.items():
        n = int(txt.str.contains(re.compile(pat, re.I), na=False).sum())
        rows.append(dict(dataset=label, n_total=len(df), keyword_set="narrow (ours)",
                         group=name, hits=n))
    for name, words in BROAD.items():
        pat = re.compile("|".join(re.escape(w) for w in words), re.I)
        rows.append(dict(dataset=label, n_total=len(df), keyword_set="broad (auditor)",
                         group=name, hits=int(txt.str.contains(pat, na=False).sum())))
    return rows


def main():
    v1 = pd.read_csv("v1_all_tracks.csv")
    v2 = pd.read_csv("ag2_conditions.csv")
    rows = counts(v1, "AlphaGenome 1 tracks") + counts(v2, "AlphaGenome 2 conditions")
    t = pd.DataFrame(rows)
    t.to_csv("table1_keyword_counts.csv", index=False)
    for ds in t.dataset.unique():
        print(f"\n== {ds}（n={t[t.dataset==ds].n_total.iloc[0]}）")
        sub = t[t.dataset == ds].pivot(index="group", columns="keyword_set", values="hits")
        print(sub.fillna(0).astype(int).to_string())


if __name__ == "__main__":
    main()
