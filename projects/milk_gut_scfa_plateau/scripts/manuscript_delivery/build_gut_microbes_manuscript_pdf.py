import os
from pathlib import Path
import re

from pypdf import PdfReader, PdfWriter
from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import cm
from reportlab.platypus import (
    BaseDocTemplate,
    Frame,
    PageTemplate,
    Paragraph,
    Preformatted,
    Spacer,
    Table,
    TableStyle,
)


ROOT = Path(os.environ.get("FANLAB_RESEARCH2_ROOT", Path.cwd())).resolve()
PKG = ROOT / "final_manuscript_planning_20260605"
REV = PKG / "gut_microbes_revision"
MD_IN = REV / "gut_microbes_manuscript_expert_review.md"
TEXT_PDF = REV / "gut_microbes_manuscript_expert_review_text.pdf"
FINAL_PDF = REV / "gut_microbes_manuscript_expert_review_with_figures.pdf"
FIGURES_PDF = PKG / "compiled_figures" / "all_figures_one_pdf_each_page.pdf"


styles = getSampleStyleSheet()
styles.add(ParagraphStyle(name="GMTitle", fontName="Helvetica-Bold", fontSize=15, leading=20, spaceAfter=10, textColor=colors.HexColor("#102A43")))
styles.add(ParagraphStyle(name="GMH1", fontName="Helvetica-Bold", fontSize=12.5, leading=16, spaceBefore=10, spaceAfter=5, textColor=colors.HexColor("#0B3D4A")))
styles.add(ParagraphStyle(name="GMH2", fontName="Helvetica-Bold", fontSize=10.2, leading=13.5, spaceBefore=7, spaceAfter=4, textColor=colors.HexColor("#174A7C")))
styles.add(ParagraphStyle(name="GMBody", fontName="Helvetica", fontSize=8.4, leading=11.7, alignment=TA_LEFT, spaceAfter=4))
styles.add(ParagraphStyle(name="GMSmall", fontName="Helvetica", fontSize=7.2, leading=9.6, textColor=colors.HexColor("#444444"), spaceAfter=3))
styles.add(ParagraphStyle(name="GMBullet", fontName="Helvetica", fontSize=8.1, leading=11.2, leftIndent=12, firstLineIndent=-8, spaceAfter=2))


def inline_md(text: str) -> str:
    text = text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    text = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", text)
    text = re.sub(r"`(.+?)`", r"<font face='Courier'>\1</font>", text)
    return text


def flush_table(story, rows):
    if not rows:
        return
    data = []
    for row in rows:
        cells = [c.strip() for c in row.strip().strip("|").split("|")]
        if all(set(c) <= {"-", ":", " "} for c in cells):
            continue
        data.append([Paragraph(inline_md(c), styles["GMSmall"]) for c in cells])
    if not data:
        return
    ncols = max(len(r) for r in data)
    widths = [16.3 * cm / ncols] * ncols
    table = Table(data, colWidths=widths, repeatRows=1)
    table.setStyle(
        TableStyle(
            [
                ("FONTNAME", (0, 0), (-1, -1), "Helvetica"),
                ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#E7F0F2")),
                ("GRID", (0, 0), (-1, -1), 0.25, colors.HexColor("#C8D7DA")),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("LEFTPADDING", (0, 0), (-1, -1), 4),
                ("RIGHTPADDING", (0, 0), (-1, -1), 4),
                ("TOPPADDING", (0, 0), (-1, -1), 3),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
            ]
        )
    )
    story.append(table)
    story.append(Spacer(1, 0.08 * cm))


def build_story(markdown: str):
    story = []
    table_rows = []
    quote_lines = []

    def flush_quote():
        nonlocal quote_lines
        if quote_lines:
            story.append(Preformatted("\n".join(quote_lines), styles["GMSmall"], maxLineLength=110))
            story.append(Spacer(1, 0.08 * cm))
        quote_lines = []

    for raw in markdown.splitlines():
        line = raw.rstrip()
        if line.startswith("|"):
            flush_quote()
            table_rows.append(line)
            continue
        flush_table(story, table_rows)
        table_rows = []

        if line.startswith(">"):
            quote_lines.append(line[1:].strip())
            continue
        flush_quote()

        if not line.strip():
            story.append(Spacer(1, 0.07 * cm))
        elif line.startswith("# "):
            story.append(Paragraph(inline_md(line[2:]), styles["GMTitle"]))
        elif line.startswith("## "):
            story.append(Paragraph(inline_md(line[3:]), styles["GMH1"]))
        elif line.startswith("### "):
            story.append(Paragraph(inline_md(line[4:]), styles["GMH2"]))
        elif line.startswith("- "):
            story.append(Paragraph("- " + inline_md(line[2:]), styles["GMBullet"]))
        elif re.match(r"^\d+\.\s+", line):
            story.append(Paragraph(inline_md(line), styles["GMBullet"]))
        else:
            story.append(Paragraph(inline_md(line), styles["GMBody"]))

    flush_table(story, table_rows)
    flush_quote()
    return story


def footer(canvas, doc):
    canvas.saveState()
    canvas.setFont("Helvetica", 7)
    canvas.setFillColor(colors.HexColor("#666666"))
    canvas.drawString(1.7 * cm, 1.05 * cm, "Gut Microbes expert-review manuscript draft")
    canvas.drawRightString(A4[0] - 1.7 * cm, 1.05 * cm, str(doc.page))
    canvas.restoreState()


def build_text_pdf():
    story = build_story(MD_IN.read_text(encoding="utf-8"))
    doc = BaseDocTemplate(
        str(TEXT_PDF),
        pagesize=A4,
        leftMargin=1.65 * cm,
        rightMargin=1.65 * cm,
        topMargin=1.4 * cm,
        bottomMargin=1.55 * cm,
    )
    frame = Frame(doc.leftMargin, doc.bottomMargin, doc.width, doc.height, id="normal")
    doc.addPageTemplates([PageTemplate(id="normal", frames=[frame], onPage=footer)])
    doc.build(story)


def merge_with_figures():
    writer = PdfWriter()
    for pdf in [TEXT_PDF, FIGURES_PDF]:
        reader = PdfReader(str(pdf))
        for page in reader.pages:
            writer.add_page(page)
    with FINAL_PDF.open("wb") as f:
        writer.write(f)


def main():
    build_text_pdf()
    merge_with_figures()
    print(FINAL_PDF)


if __name__ == "__main__":
    main()
