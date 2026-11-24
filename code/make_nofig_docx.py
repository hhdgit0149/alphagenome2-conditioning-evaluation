#!/usr/bin/env python3
"""Manuscript docx without embedded figures, for journals that want figures uploaded
separately. Captions stay where they are; only the image is removed."""
import re, subprocess, sys, os

src = "manuscript_cited.md"
md = open(src, encoding="utf-8").read()

n = {"n": 0}

def strip(m):
    n["n"] += 1
    return f"**Figure {n['n']}.** " + m.group(1) + "\n"

md = re.sub(r"!\[([^\]]+)\]\([^)]+\)\{width=[^}]*\}", strip, md)
open("manuscript_nofig.md", "w", encoding="utf-8", newline="\n").write(md)

subprocess.run([os.path.join(os.environ["LOCALAPPDATA"], "Pandoc", "pandoc.exe"),
                "manuscript_nofig.md", "-o", "manuscript_submission.docx",
                "--citeproc", "--bibliography=references_final.bib",
                "--csl=biomed-central.csl"], check=True)
print("wrote manuscript_submission.docx")

def bmc_layout(path):
    """Line numbers, page numbers and double spacing, as BMC asks for."""
    from docx import Document
    from docx.shared import Pt
    from docx.oxml.ns import qn
    from docx.oxml import OxmlElement

    doc = Document(path)
    for para in doc.paragraphs:
        para.paragraph_format.line_spacing = 2.0

    for section in doc.sections:
        sectPr = section._sectPr
        for old in sectPr.findall(qn("w:lnNumType")):
            sectPr.remove(old)
        ln = OxmlElement("w:lnNumType")
        ln.set(qn("w:countBy"), "1")
        ln.set(qn("w:restart"), "continuous")
        sectPr.append(ln)

        footer = section.footer
        fp = footer.paragraphs[0] if footer.paragraphs else footer.add_paragraph()
        fp.text = ""
        fp.alignment = 1
        run = fp.add_run()
        fld1 = OxmlElement("w:fldChar"); fld1.set(qn("w:fldCharType"), "begin")
        instr = OxmlElement("w:instrText"); instr.text = " PAGE "
        fld2 = OxmlElement("w:fldChar"); fld2.set(qn("w:fldCharType"), "end")
        run._r.append(fld1); run._r.append(instr); run._r.append(fld2)
    doc.save(path)
    print("layout applied: line numbers, page numbers, double spacing")


bmc_layout("manuscript_submission.docx")
