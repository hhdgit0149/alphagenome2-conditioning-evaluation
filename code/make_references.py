#!/usr/bin/env python3
"""Regenerate references.bib from Europe PMC author lists.

Rationale: expanding initials to first names by hand produced wrong names on five of
sixteen entries (Zhou Jian for Jingtian Zhou; Agarwal Vishal for Vikram; the whole Guo
entry; Bass Joonas/Sohini/Joseph for Joshua/Sach/John; Aksu Elif Dilara for Ekin Deniz).
This script therefore takes every author string from `_authors_epmc.json`, which was
fetched from Europe PMC `authorList` (firstName + lastName fields), and never from an
initial.

The AlphaGenome 2 report is not indexed in Europe PMC; its author list is transcribed
from page 1 of the PDF at the URL below (accessed 2026-11-16).
"""
import json

SNAP = "_authors_epmc.json"

# bibliographic fields, from Europe PMC core records retrieved 2026-11-16
META = {
    "10.1038/s41586-025-10014-0": dict(
        key="alphagenome2026", type="article",
        title="Advancing regulatory variant effect prediction with {AlphaGenome}",
        journal="Nature", volume="649", number="8099", pages="1206--1218",
        year="2026", pmid="41606153",
        note="The first-generation model evaluated here"),
    "10.1038/s41592-021-01252-x": dict(
        key="enformer2021", type="article",
        title="Effective gene expression prediction from sequence by integrating "
              "long-range interactions",
        journal="Nature Methods", volume="18", number="10", pages="1196--1203",
        year="2021"),
    "10.64898/2026.10.19.760026": dict(
        key="bolero2026", type="preprint",
        title="Predicting cell-state-specific gene regulation from {DNA} across species",
        journal="bioRxiv", year="2026", posted="20 October 2026"),
    "10.1038/nbt.4042": dict(
        key="kang2018", type="article",
        title="Multiplexed droplet single-cell {RNA}-sequencing using natural genetic "
              "variation",
        journal="Nature Biotechnology", volume="36", number="1", pages="89--94",
        year="2018", pmid="29227470"),
    "10.1038/s41588-021-00924-w": dict(
        key="kerimov2021", type="article",
        title="A compendium of uniformly processed human gene expression and splicing "
              "quantitative trait loci",
        journal="Nature Genetics", volume="53", number="9", pages="1290--1299",
        year="2021", pmid="34493866"),
    "10.1111/rssb.12388": dict(
        key="susie2020", type="article",
        title="A simple new approach to variable selection in regression, with "
              "application to genetic fine mapping",
        journal="Journal of the Royal Statistical Society: Series B", volume="82",
        number="5", pages="1273--1300", year="2020"),
    "10.1126/science.1246949": dict(
        key="fairfax2014", type="article",
        title="Innate immune activity conditions the effect of regulatory variants upon "
              "monocyte gene expression",
        journal="Science", volume="343", number="6175", pages="1246949",
        year="2014", pmid="24604202"),
    "10.1016/j.cell.2016.09.024": dict(
        key="quach2016", type="article",
        title="Genetic adaptation and {Neandertal} admixture shaped the immune system of "
              "human populations",
        journal="Cell", volume="167", number="3", pages="643--656.e17",
        year="2016", pmid="27768888"),
    "10.64898/2026.11.06.764851": dict(
        key="sreenivasan2026", type="preprint",
        title="Systematic promoter perturbation reveals limitations of {AlphaGenome} gene "
              "expression predictions",
        journal="bioRxiv", year="2026", posted="9 November 2026"),
    "10.64898/2026.10.13.758880": dict(
        key="bass2026", type="preprint",
        title="Genomic foundation models do not yet boost common variant discovery beyond "
              "standard functional annotations",
        journal="bioRxiv", year="2026", posted="19 October 2026"),
    "10.64898/2026.10.22.760090": dict(
        key="guo2026", type="preprint",
        title="Evaluating the added value of sequence models in noncoding {RNA} assays: "
              "a baseline-aware framework",
        journal="bioRxiv", year="2026", posted="25 October 2026"),
    "10.64898/2026.10.28.762238": dict(
        key="scooby2_2026", type="preprint",
        title="Scooby-2: fine-tuning {AlphaGenome} for single-cell genomic track prediction",
        journal="bioRxiv", year="2026", posted="29 October 2026"),
    "10.1038/s41467-021-25960-2": dict(
        key="squair2021", type="article",
        title="Confronting false discoveries in single-cell differential expression",
        journal="Nature Communications", volume="12", pages="5692", year="2021"),
    "10.1073/pnas.2305773120": dict(
        key="pennergoeke2023", type="article",
        title="High-throughput screening of glucocorticoid-induced enhancer activity "
              "reveals mechanisms of stress-related psychiatric disorders",
        journal="Proceedings of the National Academy of Sciences", volume="120",
        number="49", pages="e2305773120", year="2023", pmid="38011552"),
    "10.1093/nar/gkac1198": dict(
        key="morova2023", type="article",
        title="Optimized high-throughput screening of non-coding variants identified from "
              "genome-wide association studies",
        journal="Nucleic Acids Research", volume="51", number="3", pages="e18",
        year="2023", pmid="36546757"),
    "10.1186/s13059-021-02339-6": dict(
        key="huang2021", type="article",
        title="Functional mapping of androgen receptor enhancer activity",
        journal="Genome Biology", volume="22", number="1", pages="149",
        year="2021", pmid="33975627"),
}

# transcribed from page 1 of the report (accessed 2026-11-16); not in Europe PMC
AG2 = dict(
    key="alphagenome2_2026", type="techreport",
    title="{AlphaGenome} 2 predicts regulatory activity and variant effects across cell "
          "states",
    authors=[("Gayoso", "Adam"), ("Mottram", "Anne"), ("Karollus", "Alexander"),
             ("Novati", "Guido"), ("Drotár", "Pavol"), ("Reid", "John"),
             ("Meng-Papaxanthos", "Laetitia"), ("Kosiorek", "Adam"),
             ("Arvaniti", "Eirini"), ("Dutordoir", "Vincent"), ("De", "Soham"),
             ("Ward", "Tom"), ("Wong", "Lai Hong"), ("Perino", "Matteo"),
             ("Taylor", "Kyle"), ("Cheng", "Jun"), ("Latysheva", "Natasha"),
             ("Pan", "Joshua"), ("Hariharan", "Dhavanthi"), ("Bycroft", "Clare"),
             ("Nicolaisen", "Lauren"), ("Hassabis", "Demis"), ("Kohli", "Pushmeet"),
             ("Avsec", "Žiga")],
    institution="Google DeepMind", year="2026",
    url="https://storage.googleapis.com/deepmind-media/DeepMind.com/Blog/"
        "alphagenome-2-deciphering-how-cells-read-dna-to-better-understand-disease/"
        "alphagenome-2.pdf",
    note="Technical report, released 11 November 2026; accessed 16 November 2026")


def fmt_authors(pairs):
    """pairs are (last, first) or (last, first, fullName) from the Europe PMC snapshot."""
    out = []
    for a in pairs:
        last, first = a[0], a[1]
        out.append(f"{last}, {{{first}}}")
    return " and ".join(out)


def main():
    snap = json.load(open(SNAP, encoding="utf-8"))
    lines = [
        "% Generated by make_references.py from Europe PMC authorList records",
        "% (retrieved 2026-11-16). Do not edit author strings by hand: the previous",
        "% hand-expanded initials produced wrong names on five of sixteen entries.",
        "",
    ]
    order = [k for k, v in META.items()]

    def emit(entry):
        pairs = snap.get(entry["doi"]) if entry.get("doi") else None
        if pairs is None:
            pairs = entry.get("authors_pairs") or entry["authors"]
        body = [f"@{entry['type']}{{{entry['key']},",
                f"  author      = {{{fmt_authors(pairs)}}},",
                f"  title       = {{{entry['title']}}},"]
        if entry["type"] == "techreport":
            body.append(f"  institution = {{{entry['institution']}}},")
        else:
            body.append(f"  journal     = {{{entry['journal']}}},")
            for f in ("volume", "number", "pages"):
                if entry.get(f):
                    body.append(f"  {f:11s} = {{{entry[f]}}},")
        body.append(f"  year        = {{{entry['year']}}},")
        if entry.get("doi"):
            body.append(f"  doi         = {{{entry['doi']}}},")
        if entry.get("url"):
            body.append(f"  url         = {{{entry['url']}}},")
        if entry.get("pmid"):
            body.append(f"  pmid        = {{{entry['pmid']}}},")
        if entry.get("note"):
            body.append(f"  note        = {{{entry['note']}}},")
        body[-1] = body[-1].rstrip(",")
        body.append("}")
        return "\n".join(body)

    for doi in order:
        e = dict(META[doi]); e["doi"] = doi
        lines.append(emit(e)); lines.append("")

    lines.append(emit(AG2)); lines.append("")

    open("references.bib", "w", encoding="utf-8").write("\n".join(lines))
    n = len(order) + 1
    print(f"wrote references.bib: {n} entries")
    for doi in order:
        a = snap[doi]
        print(f"  {META[doi]['key']:18s} {len(a):2d} authors | first: {a[0][0]}, {a[0][1]}")
    print(f"  {'alphagenome2_2026':18s} {len(AG2['authors']):2d} authors | first: "
          f"{AG2['authors'][0][0]}, {AG2['authors'][0][1]}")


if __name__ == "__main__":
    main()
