#!/usr/bin/env python3
"""Typeset the submission PDF from a Markdown-subset source, using matplotlib only.

No LaTeX, pandoc or Office needed: text is written as vector glyphs, figures embedded at
300 dpi, tables drawn with booktabs-style rules. Page size A4, 2.54 cm margins.

Markdown subset
---------------
# Title                      -> title block (first one only)
## Section                   -> heading level 1
### Subsection              -> heading level 2
**bold** and *italic*       -> inline runs
- item                      -> bullet
[[FIG:file.png|N|caption]]  -> figure with caption "Figure N | caption"
[[TABLE:file.csv|l|N|caption|0.2,0.5,0.3]] -> table from csv
[[PAGEBREAK]]
Other lines                 -> paragraphs (blank line separates)

Usage:  python build_pdf.py manuscript.md out.pdf
"""
import re
import sys

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42      # embed as CID TrueType, not Type 3
matplotlib.rcParams["ps.fonttype"] = 42
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Rectangle
from PIL import Image

A4 = (595.28, 841.89)          # points
MARGIN = 72.0                  # 2.54 cm
TEXT_W = A4[0] - 2 * MARGIN
TEXT_H = A4[1] - 2 * MARGIN
SERIF = "Times New Roman"

from matplotlib.textpath import TextPath

_FP = {}
_TP = {}
_SPACE = {}


def _fp(size, bold, italic):
    key = (size, bold, italic)
    if key not in _FP:
        _FP[key] = FontProperties(family=SERIF, size=size,
                                  weight="bold" if bold else "normal",
                                  style="italic" if italic else "normal")
    return _FP[key]


def _tp(s, size, bold, italic):
    key = (s, size, bold, italic)
    if key not in _TP:
        # measured with the same font object used for drawing; the Agg renderer falls
        # back to a wider default face and over-measures by ~11%
        _TP[key] = TextPath((0, 0), s, size=size, prop=_fp(size, bold, italic)) \
            .get_extents().width
    return _TP[key]


def _space(size, bold, italic):
    key = (size, bold, italic)
    if key not in _SPACE:
        _SPACE[key] = max(_tp("H H", size, bold, italic)
                          - 2 * _tp("H", size, bold, italic), 1.0)
    return _SPACE[key]


def tw(s, size=10, bold=False, italic=False):
    """Advance width of a string in points (1.02 safety margin)."""
    if not s.strip():
        return _space(size, bold, italic) * len(s)
    return _tp(s, size, bold, italic) * 1.02


TOKEN = re.compile(r"(\*\*.+?\*\*|\*[^*]+?\*|_[^_]+?_)")


def runs(s):
    """Split a paragraph into (text, bold, italic) runs."""
    out = []
    for part in TOKEN.split(s):
        if not part:
            continue
        if part.startswith("**") and part.endswith("**"):
            out.append((part[2:-2], True, False))
        elif part.startswith("*") and part.endswith("*") and len(part) > 2:
            out.append((part[1:-1], False, True))
        elif part.startswith("_") and part.endswith("_") and len(part) > 2:
            out.append((part[1:-1], False, True))
        else:
            out.append((part, False, False))
    return out


def words_of(s, size):
    """Flatten runs into (word, bold, italic, space_before) tokens with widths."""
    toks = []
    for text, b, i in runs(s):
        pieces = text.split(" ")
        for k, p in enumerate(pieces):
            if p == "":
                continue
            toks.append(dict(w=p, b=b, i=i, width=tw(p, size, b, i),
                             space=(k > 0)))
    return toks


class Doc:
    def __init__(self, footnote=None):
        self.pages = []
        self.footnote = footnote
        self._new_page()

    # ---------------------------------------------------------------- pages
    def _new_page(self):
        fig = plt.figure(figsize=(A4[0] / 72.0, A4[1] / 72.0), dpi=72)
        ax = fig.add_axes([0, 0, 1, 1])
        ax.set_xlim(0, A4[0])
        ax.set_ylim(0, A4[1])
        ax.axis("off")
        ax.set_aspect("auto")
        self.pages.append((fig, ax))
        self.y = A4[1] - MARGIN
        self.ax = ax

    def ensure(self, h):
        if self.y - h < MARGIN:
            self._new_page()

    def space(self, h):
        self.y -= h

    # ---------------------------------------------------------------- text
    def _draw_line(self, toks, x, baseline, size, target=None):
        """One text object per style run, spaces inside the run.

        Drawing word by word leaves isolated glyphs in the text layer and lets a
        symbol such as rho be extracted on its own line, so consecutive words of the
        same style are joined into a single string here.
        """
        chunks = []
        for k, t in enumerate(toks):
            txt = (" " if k > 0 else "") + t["w"]
            if chunks and chunks[-1][1] == t["b"] and chunks[-1][2] == t["i"]:
                chunks[-1][0] += txt
            else:
                chunks.append([txt, t["b"], t["i"]])
        cx = x
        for txt, b, i in chunks:
            self.ax.text(cx, baseline, txt, fontsize=size, family=SERIF,
                         fontweight="bold" if b else "normal",
                         fontstyle="italic" if i else "normal",
                         ha="left", va="baseline", color="black")
            cx += tw(txt, size, b, i)

    def paragraph(self, s, size=10.0, leading=1.38, indent=0.0, hang=0.0,
                  space_before=0.0, space_after=6.0, justify=False,
                  align="left", color_note=None):
        s = s.strip()
        if not s:
            return
        self.space(space_before)
        toks = words_of(s, size)
        avail = TEXT_W - indent - hang
        space_w = tw(" ", size)
        lines, cur, curw = [], [], 0.0
        for t in toks:
            add = t["width"] + (space_w if cur else 0.0)
            if cur and curw + add > avail:
                lines.append(cur)
                cur, curw = [t], t["width"]
            else:
                cur.append(t)
                curw += add
        if cur:
            lines.append(cur)
        for idx, line in enumerate(lines):
            lh = size * leading
            self.ensure(lh)
            x = MARGIN + indent + (hang if idx > 0 else 0.0)
            target = None
            if justify and idx < len(lines) - 1 and align == "justify":
                target = TEXT_W - indent - (hang if idx > 0 else 0.0)
            self._draw_line(line, x, self.y - size * 0.82, size, target)
            self.y -= lh
        self.space(space_after)

    def heading(self, s, level=1):
        size = {1: 11.5, 2: 10.2, 3: 10.0}.get(level, 10)
        self.space({1: 12, 2: 9, 3: 7}.get(level, 6))
        self.ensure(size * 2.2)
        self.ax.text(MARGIN, self.y - size * 0.82, s, fontsize=size, family=SERIF,
                     fontweight="bold", ha="left", va="baseline")
        self.y -= size * 1.5
        self.space({1: 4, 2: 3, 3: 2}.get(level, 2))

    def bullets(self, items, size=10.0, leading=1.32):
        for it in items:
            self.paragraph("\u2022  " + it, size=size, leading=leading,
                           indent=14, space_after=2.5, justify=False)

    # ------------------------------------------------------------- figures
    def figure(self, path, label, caption, width=TEXT_W, max_frac=0.62,
               space_before=10, space_after=6):
        im = Image.open(path)
        w, h = im.size
        ww = width
        hh = ww * h / float(w)
        cap = max_frac * TEXT_H
        if hh > cap:
            hh = cap
            ww = hh * w / float(h)
        self.space(space_before)
        if self.y - hh < MARGIN + 40:
            self._new_page()
        x0 = MARGIN + (TEXT_W - ww) / 2.0
        self.ax.imshow(np.asarray(im), extent=(x0, x0 + ww, self.y - hh, self.y),
                       aspect="auto", zorder=2, interpolation="lanczos")
        self.y -= hh
        self.paragraph(f"**{label}**  {caption}", size=8.2, leading=1.30,
                       space_before=5, space_after=space_after, justify=False)

    # -------------------------------------------------------------- tables
    def table(self, df, label, caption, fracs, size=7.6, leading=1.22,
              space_before=8, space_after=8, max_frac=0.75):
        cols = list(df.columns)
        widths = [f * TEXT_W for f in fracs]
        self.space(space_before)
        head = caption if label in ("-", "") else f"**{label}**  {caption}"
        self.paragraph(head, size=8.2, leading=1.30,
                       space_before=0, space_after=5, justify=False)

        def wrap_cell(txt, w):
            toks = words_of(str(txt), size)
            out, cur, curw = [], [], 0.0
            sw = tw(" ", size)
            for t in toks:
                add = t["width"] + (sw if cur else 0)
                if cur and curw + add > w - 4:
                    out.append(cur)
                    cur, curw = [t], t["width"]
                else:
                    cur.append(t)
                    curw += add
            if cur:
                out.append(cur)
            return out or [[dict(w="", b=False, i=False, width=0)]]

        def row_height(cells):
            n = max(len(wrap_cell(c, w)) for c, w in zip(cells, widths))
            return n * size * leading + 3.5

        header = list(cols)
        hh = row_height(header)
        self.ensure(hh + 20)
        top = self.y

        def draw_row(cells, bold=False, ytop=None):
            ytop = self.y if ytop is None else ytop
            x = MARGIN
            rh = row_height(cells)
            for c, w in zip(cells, widths):
                lines = wrap_cell(c, w)
                yy = ytop - size * 0.95
                for ln in lines:
                    xx = x
                    for j, t in enumerate(ln):
                        if bold:
                            t = dict(t, b=True)
                        chunk = t["w"] + (" " if j < len(ln) - 1 else "")
                        self.ax.text(xx, yy, chunk, fontsize=size, family=SERIF,
                                     fontweight="bold" if t["b"] else "normal",
                                     fontstyle="italic" if t["i"] else "normal",
                                     ha="left", va="baseline")
                        xx += t["width"] + tw(" ", size)
                    yy -= size * leading
                x += w
            return rh

        rh = draw_row(header, bold=True, ytop=top)
        self.y = top - rh
        self.ax.plot([MARGIN, MARGIN + TEXT_W], [self.y + 2, self.y + 2],
                     lw=0.8, color="black")
        self.ax.plot([MARGIN, MARGIN + TEXT_W], [top + 4, top + 4],
                     lw=0.8, color="black")

        limit = MARGIN + max(0, 0) + 40
        for _, r in df.iterrows():
            cells = [r[c] for c in cols]
            rh = row_height(cells)
            if self.y - rh < MARGIN + 20:
                self._new_page()
                top2 = self.y
                rh2 = draw_row(header, bold=True, ytop=top2)
                self.ax.plot([MARGIN, MARGIN + TEXT_W], [top2 + 4, top2 + 4],
                             lw=0.8, color="black")
                self.ax.plot([MARGIN, MARGIN + TEXT_W], [top2 - rh2 + 2, top2 - rh2 + 2],
                             lw=0.8, color="black")
                self.y = top2 - rh2
            self.y -= draw_row(cells)
            self.ax.plot([MARGIN, MARGIN + TEXT_W], [self.y + 4, self.y + 4],
                         lw=0.35, color="0.75")
        self.ax.plot([MARGIN, MARGIN + TEXT_W], [self.y + 1, self.y + 1],
                     lw=0.8, color="black")
        self.space(space_after)

    # ---------------------------------------------------------------- save
    def save(self, path):
        n = len(self.pages)
        for i, (fig, ax) in enumerate(self.pages, 1):
            ax.text(A4[0] / 2.0, MARGIN / 2.0, str(i), fontsize=8, family=SERIF,
                    ha="center", va="center", color="0.35")
        with PdfPages(path) as pdf:
            for fig, _ in self.pages:
                pdf.savefig(fig, dpi=300)
        for fig, _ in self.pages:
            plt.close(fig)
        print(f"wrote {path}  ({n} pages)")


# ------------------------------------------------------------------ driver
def build(src, out):
    doc = Doc()
    lines = open(src, encoding="utf-8").read().split("\n")
    buf = []
    i = 0

    def flush():
        nonlocal buf
        if not buf:
            return
        block = " ".join(x.strip() for x in buf).strip()
        buf = []
        if not block:
            return
        if block.startswith("- "):
            doc.bullets([block[2:].strip()])
        else:
            doc.paragraph(block, justify=False)

    while i < len(lines):
        ln = lines[i].rstrip()
        if ln.startswith("# "):
            flush()
            doc.space(18)
            doc.paragraph(ln[2:], size=15, leading=1.25, space_after=10,
                          justify=False)
        elif ln.startswith("## "):
            flush()
            doc.heading(ln[3:], 1)
        elif ln.startswith("### "):
            flush()
            doc.heading(ln[4:], 2)
        elif ln.startswith("[[FIG:"):
            flush()
            body = ln[6:].rstrip("]")
            path, label, cap = body.split("|", 2)
            doc.figure(path.strip(), label.strip(), cap.strip())
        elif ln.startswith("[[TABLE:"):
            flush()
            body = ln[8:].rstrip("]")
            path, label, cap, fracs = body.split("|", 3)
            df = pd.read_csv(path.strip()).fillna("")
            doc.table(df, label.strip(), cap.strip(),
                      [float(x) for x in fracs.split(",")])
        elif ln.startswith("@"):
            flush()
            doc.paragraph(ln[1:], size=8.6, leading=1.30, hang=16,
                          space_after=3.2, justify=False)
        elif ln.startswith("- "):
            flush()
            doc.bullets([ln[2:].strip()])
        elif ln.startswith("[[PAGEBREAK]]"):
            flush()
            doc._new_page()
        elif ln.strip() == "":
            flush()
        else:
            buf.append(ln)
        i += 1
    flush()
    doc.save(out)


if __name__ == "__main__":
    build(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else "submission/manuscript.pdf")
