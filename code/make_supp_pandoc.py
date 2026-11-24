#!/usr/bin/env python3
"""Emit a pandoc-compatible source for the supplementary information.

The working source (supp.md) uses [[TABLE:...]] directives for the PDF typesetter.
This script rewrites them as pandoc pipe tables whose column widths follow the same
fractions, so the output can be built with:

  pandoc supp_pandoc.md -o Supplementary.pdf --pdf-engine=xelatex \
    -V mainfont="Liberation Serif" -V fontsize=10pt \
    -V geometry:a4paper -V geometry:margin=2.54cm -V linestretch=1.15 \
    --columns=110

Cells wrap because the separator row carries relative dash widths.
"""
import os
import re

import pandas as pd

SRC = "supp.md"
OUT = "supp_pandoc.md"

# column fractions used by the typesetter, keyed by csv name
FRACS = {
    "TableS1_provenance.csv": [0.31, 0.14, 0.13, 0.16, 0.26],
    "TableS2_display.csv": [0.21, 0.15, 0.08, 0.13, 0.09, 0.09, 0.25],
    "TableS3_display.csv": [0.16, 0.08, 0.13, 0.26, 0.07, 0.30],
}
# width of the text block in characters at --columns=110
TOTAL = 104


def pipe_table(csv_path, fracs):
    df = pd.read_csv(csv_path).fillna("")
    cols = [str(c) for c in df.columns]
    widths = [max(3, int(round(f * TOTAL))) for f in fracs]

    def esc(v):
        return str(v).replace("|", "\\|").strip()

    def row(cells):
        return "| " + " | ".join(esc(c) for c in cells) + " |"

    sep = "|" + "|".join("-" * w for w in widths) + "|"
    lines = [row(cols), sep]
    for _, r in df.iterrows():
        lines.append(row([r[c] for c in df.columns]))
    return "\n".join(lines)


def main():
    s = open(SRC, encoding="utf-8").read()

    def repl(m):
        body = m.group(1)
        path, label, cap, frac = body.split("|", 3)
        path = path.strip()
        table = pipe_table(path, [float(x) for x in frac.split(",")])
        head = cap.strip()
        return f"{head}\n\n{table}\n"

    s = re.sub(r"\[\[TABLE:(.+?)\]\]", repl, s)

    def fig(m):
        parts = m.group(1).split("|")
        path, label, cap = parts[0], parts[1], parts[2]
        width = parts[3].strip() if len(parts) > 3 else "100%"
        return (f"![{label.strip()}. {cap.strip()}]({path.strip()})"
                f"{{width={width}}}\n")

    s = re.sub(r"\[\[FIG:(.+?)\]\]", fig, s)
    s = s.replace("[[PAGEBREAK]]", "\\newpage")
    # the pandoc source carries no page-number footer; LaTeX supplies it
    open(OUT, "w", encoding="utf-8").write(s)
    print(f"wrote {OUT}  ({len(s.splitlines())} lines)")


if __name__ == "__main__":
    main()
