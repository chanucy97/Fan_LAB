import os
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import cm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import BaseDocTemplate, Frame, PageTemplate, Paragraph, Spacer


ROOT = Path(os.environ.get("FANLAB_RESEARCH2_ROOT", Path.cwd())).resolve()
REV = ROOT / "final_manuscript_planning_20260605" / "gut_microbes_revision"
MD_IN = REV / "README_expert_review_CN.md"
PDF_OUT = REV / "README_expert_review_CN.pdf"

pdfmetrics.registerFont(TTFont("MSYH", r"C:\Windows\Fonts\msyh.ttc"))
pdfmetrics.registerFont(TTFont("MSYH-Bold", r"C:\Windows\Fonts\msyhbd.ttc"))

styles = getSampleStyleSheet()
styles.add(ParagraphStyle(name="CNTitle", fontName="MSYH-Bold", fontSize=16, leading=22, spaceAfter=10, textColor=colors.HexColor("#102A43")))
styles.add(ParagraphStyle(name="CNH1", fontName="MSYH-Bold", fontSize=12.5, leading=17, spaceBefore=9, spaceAfter=5, textColor=colors.HexColor("#0B3D4A")))
styles.add(ParagraphStyle(name="CNH2", fontName="MSYH-Bold", fontSize=10.2, leading=14, spaceBefore=7, spaceAfter=4, textColor=colors.HexColor("#174A7C")))
styles.add(ParagraphStyle(name="CNBody", fontName="MSYH", fontSize=8.7, leading=12.6, spaceAfter=4))
styles.add(ParagraphStyle(name="CNBullet", fontName="MSYH", fontSize=8.5, leading=12.2, leftIndent=13, firstLineIndent=-9, spaceAfter=2))


def esc(s: str) -> str:
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def inline_md(s: str) -> str:
    s = esc(s)
    s = s.replace("**", "")
    return s


def build_story(text: str):
    story = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            story.append(Spacer(1, 0.06 * cm))
        elif line.startswith("# "):
            story.append(Paragraph(inline_md(line[2:]), styles["CNTitle"]))
        elif line.startswith("## "):
            story.append(Paragraph(inline_md(line[3:]), styles["CNH1"]))
        elif line.startswith("### "):
            story.append(Paragraph(inline_md(line[4:]), styles["CNH2"]))
        elif line.startswith("- "):
            story.append(Paragraph("• " + inline_md(line[2:]), styles["CNBullet"]))
        elif len(line) > 2 and line[0].isdigit() and line[1] == ".":
            story.append(Paragraph(inline_md(line), styles["CNBullet"]))
        elif line.startswith(">"):
            story.append(Paragraph(inline_md(line[1:].strip()), styles["CNBullet"]))
        else:
            story.append(Paragraph(inline_md(line), styles["CNBody"]))
    return story


def footer(canvas, doc):
    canvas.saveState()
    canvas.setFont("MSYH", 7)
    canvas.setFillColor(colors.HexColor("#666666"))
    canvas.drawString(1.7 * cm, 1.05 * cm, "Gut Microbes 投稿前专家审阅说明")
    canvas.drawRightString(A4[0] - 1.7 * cm, 1.05 * cm, str(doc.page))
    canvas.restoreState()


def main():
    doc = BaseDocTemplate(
        str(PDF_OUT),
        pagesize=A4,
        leftMargin=1.7 * cm,
        rightMargin=1.7 * cm,
        topMargin=1.45 * cm,
        bottomMargin=1.55 * cm,
    )
    frame = Frame(doc.leftMargin, doc.bottomMargin, doc.width, doc.height, id="normal")
    doc.addPageTemplates([PageTemplate(id="normal", frames=[frame], onPage=footer)])
    doc.build(build_story(MD_IN.read_text(encoding="utf-8")))
    print(PDF_OUT)


if __name__ == "__main__":
    main()
