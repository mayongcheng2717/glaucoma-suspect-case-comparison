"""Manuscript builder: journal-style DOCX with figures, tables and references."""

from __future__ import annotations

import os
from typing import List, Optional, Sequence, Tuple

import pandas as pd
from docx import Document
from docx.enum.section import WD_SECTION
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_BREAK
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Cm, Pt, RGBColor

BODY_FONT = "Times New Roman"
HEAD_FONT = "Arial"


def _set_cell_border(cell, **kwargs):
    tc = cell._tc
    tcPr = tc.get_or_add_tcPr()
    borders = OxmlElement("w:tcBorders")
    for edge in ("top", "left", "bottom", "right"):
        if edge in kwargs:
            el = OxmlElement(f"w:{edge}")
            for k, v in kwargs[edge].items():
                el.set(qn(f"w:{k}"), str(v))
            borders.append(el)
    tcPr.append(borders)


def _add_page_numbers(section):
    footer = section.footer
    p = footer.paragraphs[0]
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = p.add_run()
    fld = OxmlElement("w:fldSimple")
    fld.set(qn("w:instr"), "PAGE")
    run._r.addnext(fld)


class Manuscript:
    """Accumulates a structured manuscript and writes a submission-style DOCX."""

    def __init__(self, title: str, short_title: str = "", authors: str = "",
                 affiliations: str = "", corresponding: str = ""):
        self.doc = Document()
        self.title = title
        self.n_fig = 0
        self.n_tab = 0
        self.references: List[str] = []
        sec = self.doc.sections[0]
        sec.top_margin = Cm(2.2)
        sec.bottom_margin = Cm(2.2)
        sec.left_margin = Cm(2.2)
        sec.right_margin = Cm(2.2)
        _add_page_numbers(sec)
        st = self.doc.styles["Normal"]
        st.font.name = BODY_FONT
        st.font.size = Pt(11)
        st.paragraph_format.line_spacing = 1.5
        st.paragraph_format.space_after = Pt(0)
        self._title(title)
        if authors:
            self._line(authors, size=11, bold=False, italic=False, space_after=4)
        if affiliations:
            for a in affiliations.split(";"):
                if a.strip():
                    self._line(a.strip(), size=9, italic=True, space_after=0)
        if corresponding:
            self._line(corresponding, size=9, italic=False, space_after=10)

    # -- primitives --------------------------------------------------------
    def _line(self, text, size=11, bold=False, italic=False, space_after=6,
              align=WD_ALIGN_PARAGRAPH.LEFT, color=None):
        p = self.doc.add_paragraph()
        p.alignment = align
        p.paragraph_format.space_after = Pt(space_after)
        p.paragraph_format.line_spacing = 1.3
        r = p.add_run(text)
        r.font.size = Pt(size)
        r.bold = bold
        r.italic = italic
        r.font.name = BODY_FONT
        if color:
            r.font.color.rgb = RGBColor(*color)
        return p

    def _title(self, text):
        p = self.doc.add_paragraph()
        p.paragraph_format.space_after = Pt(10)
        p.paragraph_format.line_spacing = 1.2
        r = p.add_run(text)
        r.font.size = Pt(16)
        r.bold = True
        r.font.name = HEAD_FONT

    def heading(self, text: str, level: int = 1):
        p = self.doc.add_paragraph()
        p.paragraph_format.space_before = Pt(12 if level == 1 else 9)
        p.paragraph_format.space_after = Pt(4)
        p.paragraph_format.line_spacing = 1.15
        r = p.add_run(text)
        r.bold = True
        r.font.name = HEAD_FONT
        r.font.size = Pt(12 if level == 1 else 10.5)
        return p

    def paragraph(self, text: str, first_indent: bool = False, size: int = 11,
                  italic: bool = False):
        p = self.doc.add_paragraph()
        p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
        p.paragraph_format.space_after = Pt(6)
        p.paragraph_format.line_spacing = 1.5
        if first_indent:
            p.paragraph_format.first_line_indent = Cm(0.6)
        # support *italic* segments and **bold** segments minimally
        self._rich(p, text, size, italic)
        return p

    @staticmethod
    def _rich(p, text, size, italic=False):
        import re
        tokens = re.split(r"(\*\*[^*]+\*\*|\*[^*]+\*)", text)
        for tok in tokens:
            if not tok:
                continue
            bold = tok.startswith("**") and tok.endswith("**") and len(tok) > 4
            ital = tok.startswith("*") and tok.endswith("*") and len(tok) > 2 and not bold
            body = tok[2:-2] if bold else (tok[1:-1] if ital else tok)
            r = p.add_run(body)
            r.font.size = Pt(size)
            r.font.name = BODY_FONT
            r.bold = bold
            r.italic = italic or ital

    def bullets(self, items: Sequence[str], size: int = 11):
        for it in items:
            p = self.doc.add_paragraph(style="List Bullet")
            p.paragraph_format.space_after = Pt(2)
            p.paragraph_format.line_spacing = 1.4
            self._rich(p, it, size)

    def keypoints(self, items: Sequence[str]):
        for it in items:
            p = self.doc.add_paragraph()
            p.paragraph_format.left_indent = Cm(0.4)
            p.paragraph_format.first_line_indent = Cm(-0.4)
            p.paragraph_format.space_after = Pt(3)
            p.paragraph_format.line_spacing = 1.35
            self._rich(p, "\u2022  " + it, 10.5)

    # -- tables ------------------------------------------------------------
    def table(self, df: pd.DataFrame, caption: str, font_size: int = 8,
              col_widths: Optional[Sequence[float]] = None, note: str = ""):
        """Three-line table: top rule, rule under the header, bottom rule only."""
        self.n_tab += 1
        cap = self.doc.add_paragraph()
        cap.paragraph_format.space_before = Pt(10)
        cap.paragraph_format.space_after = Pt(3)
        cap.paragraph_format.line_spacing = 1.15
        r = cap.add_run(f"Table {self.n_tab}. ")
        r.bold = True
        r.font.size = Pt(9)
        r.font.name = HEAD_FONT
        r2 = cap.add_run(caption)
        r2.font.size = Pt(9)
        r2.font.name = BODY_FONT

        n_cols = len(df.columns)
        t = self.doc.add_table(rows=1, cols=n_cols)
        t.alignment = WD_TABLE_ALIGNMENT.CENTER
        # no table-level borders; the three rules are applied per cell
        hdr = t.rows[0].cells
        for j, c in enumerate(df.columns):
            hdr[j].text = ""
            p = hdr[j].paragraphs[0]
            p.paragraph_format.space_after = Pt(1)
            run = p.add_run(str(c))
            run.bold = True
            run.font.size = Pt(font_size)
            run.font.name = HEAD_FONT
            # top rule (thick) and rule under the header (thin)
            _set_cell_border(hdr[j], top={"val": "single", "sz": 12, "space": 0,
                                          "color": "000000"},
                             bottom={"val": "single", "sz": 8, "space": 0,
                                     "color": "000000"})
        last = None
        for _, row in df.iterrows():
            cells = t.add_row().cells
            last = cells
            for j, c in enumerate(df.columns):
                cells[j].text = ""
                p = cells[j].paragraphs[0]
                p.paragraph_format.space_after = Pt(1)
                p.paragraph_format.line_spacing = 1.0
                run = p.add_run(str(row[c]))
                run.font.size = Pt(font_size)
                run.font.name = BODY_FONT
        if last is not None:
            for c in last:
                _set_cell_border(c, bottom={"val": "single", "sz": 12, "space": 0,
                                            "color": "000000"})
        if note:
            n = self.doc.add_paragraph()
            n.paragraph_format.space_before = Pt(3)
            n.paragraph_format.line_spacing = 1.15
            rr = n.add_run(note)
            rr.font.size = Pt(7.5)
            rr.italic = True
            rr.font.name = BODY_FONT
        return t

    @staticmethod
    def _cell_rule(cell, top=False, bottom=False, thin=False):
        kw = {}
        w = 2 if thin else 8
        if top:
            kw["top"] = {"val": "single", "sz": w, "space": 0, "color": "000000"}
        if bottom:
            kw["bottom"] = {"val": "single", "sz": w, "space": 0, "color": "000000"}
        if kw:
            _set_cell_border(cell, **kw)

    # -- figures -----------------------------------------------------------
    def figure(self, path_png: str, caption: str, width_cm: float = 16.5):
        self.n_fig += 1
        p = self.doc.add_paragraph()
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        p.paragraph_format.space_before = Pt(10)
        p.paragraph_format.space_after = Pt(3)
        p.add_run().add_picture(path_png, width=Cm(width_cm))
        cap = self.doc.add_paragraph()
        cap.paragraph_format.space_after = Pt(10)
        cap.paragraph_format.line_spacing = 1.2
        r = cap.add_run(f"Figure {self.n_fig}. ")
        r.bold = True
        r.font.size = Pt(9)
        r.font.name = HEAD_FONT
        self._rich(cap, caption, 9)

    # -- references --------------------------------------------------------
    def add_references(self, refs: Sequence[str]):
        self.references = list(refs)
        self.heading("References")
        for i, r in enumerate(self.references, 1):
            p = self.doc.add_paragraph()
            p.paragraph_format.left_indent = Cm(0.6)
            p.paragraph_format.first_line_indent = Cm(-0.6)
            p.paragraph_format.space_after = Pt(3)
            p.paragraph_format.line_spacing = 1.25
            run = p.add_run(f"{i}.  {r}")
            run.font.size = Pt(9)
            run.font.name = BODY_FONT

    def page_break(self):
        self.doc.add_paragraph().add_run().add_break(WD_BREAK.PAGE)

    def save(self, path: str):
        os.makedirs(os.path.dirname(path), exist_ok=True)
        self.doc.save(path)
        return path
