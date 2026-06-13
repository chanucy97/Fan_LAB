import os
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import cm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
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
REV = ROOT / "final_manuscript_planning_20260605" / "gut_microbes_revision"
MD_IN = REV / "gut_microbes_expert_review_brief.md"
PDF_OUT = REV / "gut_microbes_expert_review_brief.pdf"

pdfmetrics.registerFont(TTFont("MSYH", r"C:\Windows\Fonts\msyh.ttc"))
pdfmetrics.registerFont(TTFont("MSYH-Bold", r"C:\Windows\Fonts\msyhbd.ttc"))

styles = getSampleStyleSheet()
styles.add(ParagraphStyle(name="TitleCN", fontName="MSYH-Bold", fontSize=16, leading=22, spaceAfter=10, textColor=colors.HexColor("#102A43")))
styles.add(ParagraphStyle(name="H1CN", fontName="MSYH-Bold", fontSize=13, leading=18, spaceBefore=10, spaceAfter=6, textColor=colors.HexColor("#0B3D4A")))
styles.add(ParagraphStyle(name="H2CN", fontName="MSYH-Bold", fontSize=10.5, leading=15, spaceBefore=8, spaceAfter=4, textColor=colors.HexColor("#174A7C")))
styles.add(ParagraphStyle(name="BodyCN", fontName="MSYH", fontSize=8.4, leading=12.2, alignment=TA_LEFT, spaceAfter=4))
styles.add(ParagraphStyle(name="SmallCN", fontName="MSYH", fontSize=7.4, leading=10.2, textColor=colors.HexColor("#555555"), spaceAfter=3))
styles.add(ParagraphStyle(name="BulletCN", fontName="MSYH", fontSize=8.2, leading=11.6, leftIndent=12, firstLineIndent=-8, spaceAfter=2))


def esc(text: str) -> str:
    return (
        text.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace("**", "")
    )


def inline_md(text: str) -> str:
    text = esc(text)
    # Minimal bold support after escaping.
    while "**" in text:
        text = text.replace("**", "", 1)
    return text


def flush_table(story, rows):
    if not rows:
        return
    data = []
    for row in rows:
        cells = [c.strip() for c in row.strip().strip("|").split("|")]
        if all(set(c) <= {"-", ":", " "} for c in cells):
            continue
        data.append([Paragraph(inline_md(c), styles["SmallCN"]) for c in cells])
    if not data:
        return
    ncols = max(len(r) for r in data)
    widths = [16.3 * cm / ncols] * ncols
    table = Table(data, colWidths=widths, repeatRows=1)
    table.setStyle(
        TableStyle(
            [
                ("FONTNAME", (0, 0), (-1, -1), "MSYH"),
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
    story.append(Spacer(1, 0.12 * cm))


def build_story(markdown: str):
    story = []
    table_rows = []
    in_quote = False
    quote_lines = []

    def flush_quote():
        nonlocal quote_lines, in_quote
        if quote_lines:
            story.append(Preformatted("\n".join(quote_lines), styles["SmallCN"], maxLineLength=110))
            story.append(Spacer(1, 0.08 * cm))
        quote_lines = []
        in_quote = False

    for raw in markdown.splitlines():
        line = raw.rstrip()
        if line.startswith("|"):
            flush_quote()
            table_rows.append(line)
            continue
        flush_table(story, table_rows)
        table_rows = []

        if line.startswith(">"):
            in_quote = True
            quote_lines.append(line[1:].strip())
            continue
        if in_quote:
            flush_quote()

        if not line.strip():
            story.append(Spacer(1, 0.08 * cm))
        elif line.startswith("# "):
            story.append(Paragraph(inline_md(line[2:]), styles["TitleCN"]))
        elif line.startswith("## "):
            story.append(Paragraph(inline_md(line[3:]), styles["H1CN"]))
        elif line.startswith("### "):
            story.append(Paragraph(inline_md(line[4:]), styles["H2CN"]))
        elif line.startswith("- "):
            story.append(Paragraph("• " + inline_md(line[2:]), styles["BulletCN"]))
        elif len(line) > 2 and line[0].isdigit() and line[1:3] in {". ", ") "}:
            story.append(Paragraph(inline_md(line), styles["BulletCN"]))
        else:
            story.append(Paragraph(inline_md(line), styles["BodyCN"]))

    flush_table(story, table_rows)
    if in_quote:
        flush_quote()
    return story


def footer(canvas, doc):
    canvas.saveState()
    canvas.setFont("MSYH", 7)
    canvas.setFillColor(colors.HexColor("#666666"))
    canvas.drawString(1.7 * cm, 1.05 * cm, "Gut Microbes expert-review strengthening brief")
    canvas.drawRightString(A4[0] - 1.7 * cm, 1.05 * cm, str(doc.page))
    canvas.restoreState()


def main():
    story = build_story(MD_IN.read_text(encoding="utf-8"))
    doc = BaseDocTemplate(
        str(PDF_OUT),
        pagesize=A4,
        leftMargin=1.65 * cm,
        rightMargin=1.65 * cm,
        topMargin=1.4 * cm,
        bottomMargin=1.55 * cm,
    )
    frame = Frame(doc.leftMargin, doc.bottomMargin, doc.width, doc.height, id="normal")
    doc.addPageTemplates([PageTemplate(id="normal", frames=[frame], onPage=footer)])
    doc.build(story)
    print(PDF_OUT)


if __name__ == "__main__":
    main()
