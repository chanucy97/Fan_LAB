import os
from pathlib import Path

import fitz
from PIL import Image
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import cm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (
    BaseDocTemplate,
    Frame,
    Image as RLImage,
    KeepTogether,
    NextPageTemplate,
    PageBreak,
    PageTemplate,
    Paragraph,
    Spacer,
    Table,
    TableStyle,
)


ROOT = Path(os.environ.get("FANLAB_RESEARCH2_ROOT", Path.cwd())).resolve()
PKG = ROOT / "final_manuscript_planning_20260605"
OUT_DIR = PKG / "leadership_manuscript"
FIG_DIR = OUT_DIR / "rendered_figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

PDF_OUT = OUT_DIR / "milk_gut_scfa_leadership_manuscript.pdf"

FONT_REG = r"C:\Windows\Fonts\msyh.ttc"
FONT_BOLD = r"C:\Windows\Fonts\msyhbd.ttc"
pdfmetrics.registerFont(TTFont("MSYH", FONT_REG))
pdfmetrics.registerFont(TTFont("MSYH-Bold", FONT_BOLD))


FIGURE_SOURCES = [
    (
        "Graphical abstract",
        PKG / "compiled_figures" / "00_graphical_abstract.pdf",
    ),
    (
        "Fig. 1",
        PKG
        / "figures_redrawn"
        / "fig1_cohort_design_overview"
        / "Fig1_complete_cohort_coverage.pdf",
    ),
    (
        "Fig. 2",
        PKG
        / "figures_redrawn"
        / "fig2_milk_bioactive_remodeling"
        / "Fig2_milk_bioactive_remodeling_redrawn.pdf",
    ),
    (
        "Fig. 3",
        PKG
        / "figures_redrawn"
        / "fig3_infant_microbiome_maturation"
        / "Fig3_infant_microbiome_maturation_redrawn.pdf",
    ),
    (
        "Fig. 4",
        PKG
        / "figures_redrawn"
        / "fig4_scfa_maturation"
        / "Fig4_scfa_maturation_redrawn.pdf",
    ),
    (
        "Fig. 5",
        PKG
        / "figures_redrawn"
        / "fig5_milk_scfa_lagged_candidates"
        / "Fig5_milk_scfa_lagged_candidates_redrawn.pdf",
    ),
    (
        "Fig. S1",
        PKG
        / "figures_redrawn"
        / "fig6_maternal_fecal_context"
        / "Fig6_maternal_fecal_context_redrawn.pdf",
    ),
    (
        "Fig. S2",
        PKG
        / "figures_redrawn"
        / "fig7_clinical_modifier_layer"
        / "Fig7_clinical_modifier_layer_redrawn.pdf",
    ),
]


def render_figures() -> dict[str, Path]:
    rendered = {}
    for label, pdf_path in FIGURE_SOURCES:
        if not pdf_path.exists():
            raise FileNotFoundError(pdf_path)
        out = FIG_DIR / f"{label.replace('.', '').replace(' ', '_')}.png"
        doc = fitz.open(str(pdf_path))
        page = doc[0]
        pix = page.get_pixmap(matrix=fitz.Matrix(2.0, 2.0), alpha=False)
        pix.save(str(out))
        rendered[label] = out
    return rendered


def fit_image(path: Path, max_w: float, max_h: float) -> RLImage:
    with Image.open(path) as im:
        w, h = im.size
    scale = min(max_w / w, max_h / h)
    return RLImage(str(path), width=w * scale, height=h * scale)


def header_footer(canvas, doc):
    canvas.saveState()
    canvas.setFont("MSYH", 7)
    canvas.setFillColor(colors.HexColor("#666666"))
    canvas.drawString(1.8 * cm, 1.05 * cm, "领导审阅稿 | 内部论文方向评估版")
    canvas.drawRightString(A4[0] - 1.8 * cm, 1.05 * cm, f"{doc.page}")
    canvas.restoreState()


def cover_footer(canvas, doc):
    canvas.saveState()
    canvas.setFont("MSYH", 7)
    canvas.setFillColor(colors.HexColor("#777777"))
    canvas.drawRightString(A4[0] - 1.8 * cm, 1.05 * cm, "Confidential review draft")
    canvas.restoreState()


styles = getSampleStyleSheet()
styles.add(
    ParagraphStyle(
        name="TitleCN",
        fontName="MSYH-Bold",
        fontSize=20,
        leading=27,
        alignment=TA_CENTER,
        textColor=colors.HexColor("#102A43"),
        spaceAfter=12,
    )
)
styles.add(
    ParagraphStyle(
        name="Deck",
        fontName="MSYH",
        fontSize=10.5,
        leading=16,
        alignment=TA_CENTER,
        textColor=colors.HexColor("#52616B"),
        spaceAfter=18,
    )
)
styles.add(
    ParagraphStyle(
        name="H1CN",
        fontName="MSYH-Bold",
        fontSize=14,
        leading=19,
        textColor=colors.HexColor("#0B3D4A"),
        spaceBefore=12,
        spaceAfter=7,
    )
)
styles.add(
    ParagraphStyle(
        name="H2CN",
        fontName="MSYH-Bold",
        fontSize=11.5,
        leading=16,
        textColor=colors.HexColor("#174A7C"),
        spaceBefore=8,
        spaceAfter=5,
    )
)
styles.add(
    ParagraphStyle(
        name="BodyCN",
        fontName="MSYH",
        fontSize=9.1,
        leading=13.7,
        alignment=TA_LEFT,
        firstLineIndent=0,
        textColor=colors.HexColor("#202124"),
        spaceAfter=4,
    )
)
styles.add(
    ParagraphStyle(
        name="BodyNoIndent",
        parent=styles["BodyCN"],
        firstLineIndent=0,
    )
)
styles.add(
    ParagraphStyle(
        name="CaptionCN",
        fontName="MSYH",
        fontSize=8.1,
        leading=12,
        alignment=TA_LEFT,
        textColor=colors.HexColor("#333333"),
        spaceBefore=5,
        spaceAfter=9,
    )
)
styles.add(
    ParagraphStyle(
        name="Callout",
        fontName="MSYH",
        fontSize=9.2,
        leading=14,
        textColor=colors.HexColor("#1B3A4B"),
        backColor=colors.HexColor("#EEF7F6"),
        borderColor=colors.HexColor("#A7D8D3"),
        borderWidth=0.6,
        borderPadding=7,
        spaceBefore=6,
        spaceAfter=9,
    )
)
styles.add(
    ParagraphStyle(
        name="Takeaway",
        fontName="MSYH",
        fontSize=8.9,
        leading=13,
        textColor=colors.HexColor("#12323F"),
        backColor=colors.HexColor("#F1F8F6"),
        borderColor=colors.HexColor("#70B6AD"),
        borderWidth=0.7,
        borderPadding=7,
        spaceBefore=6,
        spaceAfter=8,
    )
)
styles.add(
    ParagraphStyle(
        name="Small",
        fontName="MSYH",
        fontSize=7.6,
        leading=10,
        textColor=colors.HexColor("#555555"),
    )
)


def P(text, style="BodyCN"):
    return Paragraph(text, styles[style])


def figure_block(label: str, image_path: Path, caption: str, max_h=14.2 * cm):
    return KeepTogether(
        [
            fit_image(image_path, 17.2 * cm, max_h),
            P(f"<b>{label}.</b> {caption}", "CaptionCN"),
        ]
    )


TAKEAWAYS = {
    "framework": (
        "本节结论：这不是一篇单纯描述母乳成分的文章，而是一篇关于“早期乳汁暴露如何与婴儿肠道生态成熟同步变化”的纵向多组学文章。",
        "意义：顶刊叙事的核心应从单一营养指标升级为母乳-肠道微生物-SCFA成熟轴，同时主动声明目前是关联框架而不是因果机制证明。",
    ),
    "fig1": (
        "本节结论：样本结构不支持把所有组学硬合成一个 complete-case 模型；最合理、也最诚实的策略是模块化整合。",
        "意义：这把一个潜在弱点转化为方法学合理性。文章可以解释为什么分别分析乳汁、宏基因组和SCFA，再用重叠样本做候选桥接。",
    ),
    "fig2": (
        "本节结论：早期母乳暴露层确实在变化，最可靠的两个入口是覆盖度好的乳铁蛋白和动态最强的HMO特征。",
        "意义：这为后文寻找“乳汁到婴儿代谢”的时间联系提供了上游依据；LCFA可作为背景层，但不宜承担核心结论。",
    ),
    "fig3": (
        "本节结论：婴儿肠道成熟的主要微生物信号集中在Bifidobacterium，尤其是B. longum在D05-D30呈稳定上升和更高检出。",
        "意义：这说明婴儿肠道生态不是随机波动，而有可识别的优势菌群成熟轨迹；但物种标签需复核，不能写成菌株传播。",
    ),
    "fig4": (
        "本节结论：婴儿粪便SCFA存在随日龄上升的代谢成熟轴，这是本文最稳的婴儿侧结局。",
        "意义：SCFA成熟分数可作为整篇文章的下游锚点，把乳汁动态和肠道微生物成熟连接到一个可量化的代谢表型。",
    ),
    "fig5": (
        "本节结论：最值得作为主文亮点的是早期LTF/HMO与后续SCFA的候选滞后联系，而不是同时间点相关。",
        "意义：这让文章从“多组学并列描述”进入“有时间顺序的整合发现”。但这些信号仍是候选关联，不能写成因果调控。",
    ),
    "figs1": (
        "本节结论：母亲粪便数据提示部分母婴配对存在生态相似性，但总体分布重叠，不足以证明垂直传播。",
        "意义：这层数据的正确用途是提供母体肠道生态背景，增强母婴系统完整性，而不是作为传播机制的主证据。",
    ),
    "figs2": (
        "本节结论：临床层面最清楚的信号是母亲BMI可能修饰婴儿SCFA成熟随年龄变化的轨迹，喂养、分娩和出生体重信号不稳定。",
        "意义：它提升临床可读性，但应放在补充或讨论中。领导可把它理解为下一步临床扩展方向，而不是当前文章主结论。",
    ),
    "discussion": (
        "总体结论：本文最有价值的发现是建立了一个可被数据支撑的“乳汁生物活性动态-婴儿Bifidobacterium结构-SCFA成熟”的纵向关联框架。",
        "意义：冲击高水平期刊时，卖点不是单个P值，而是高原母婴队列中多层数据按时间顺序拼成了一个清晰、克制、可验证的生物学故事。",
    ),
}


def takeaway_block(story, key: str):
    finding, meaning = TAKEAWAYS[key]
    story.append(P(f"<b>{finding}</b><br/>{meaning}", "Takeaway"))


def add_cover(story, figs):
    story.append(Spacer(1, 0.7 * cm))
    story.append(
        Paragraph(
            "青藏高原母婴队列中乳汁生物活性成分、婴儿肠道微生物与粪便短链脂肪酸成熟的纵向多组学关联",
            styles["TitleCN"],
        )
    )
    story.append(
        P(
            "领导审阅版 | 仅用于内部汇报和论文方向判断 | 不含作者、单位及正式投稿信息",
            "Deck",
        )
    )
    story.append(fit_image(figs["Graphical abstract"], 16.5 * cm, 7.7 * cm))
    story.append(Spacer(1, 0.35 * cm))
    story.append(
        P(
            "<b>一句话结论：</b>本研究支持一个以“早期乳汁 LTF/HMO 动态 - 婴儿 Bifidobacterium 优势成熟 - 粪便 SCFA 成熟 - 母体肠道生态背景”为主线的纵向关联模型；当前证据强度适合定位为候选时间联系和多组学整合发现，而不应表述为因果调控或垂直传播证明。",
            "Callout",
        )
    )
    story.append(
        P(
            "版本定位：本文按高水平期刊研究论文的叙事结构组织，但保留内部审阅稿的表达方式；背景文献、正式方法细节、伦理编号、作者贡献和参考文献待投稿前补齐。",
            "Small",
        )
    )
    story.append(PageBreak())


def add_summary(story):
    story.append(P("摘要", "H1CN"))
    story.append(
        P(
            "<b>背景：</b>生命早期母乳暴露、婴儿肠道微生物定植和粪便代谢成熟共同塑造婴儿肠道生态，但在高海拔地区母婴人群中，这些层面如何以纵向方式相互衔接仍缺少系统证据。"
            "<b>方法：</b>本研究整合青藏高原母婴队列中的基线临床信息、D05/D14/D30/D90 人乳生物活性成分、婴儿粪便宏基因组、婴儿粪便短链脂肪酸以及母亲粪便宏基因组资料。考虑到严格 complete-case 多组学重叠有限，分析采用模块化纵向框架，分别评估乳汁暴露动态、婴儿微生物成熟、SCFA 代谢成熟、候选滞后关联、母体肠道生态背景和临床修饰层。"
            "<b>结果：</b>基线临床信息覆盖 52 对母婴；乳铁蛋白覆盖最广，HMO 主要集中在早期乳汁窗口，LCFA 覆盖相对有限。婴儿侧结果显示，以 Bifidobacterium 为中心的肠道微生物成熟结构和随年龄上升的 SCFA 成熟轴。乳汁侧动态中，乳铁蛋白和多个 HMO 特征呈现明确的早期重塑；候选滞后模型显示，D05 乳铁蛋白与 D30 异戊酸、D14 LNT 与 D90 乙酸等联系具有较强名义统计信号。母亲粪便宏基因组结果为婴儿肠道装配提供生态背景，但不足以支持菌株传播结论。"
            "<b>结论：</b>本研究形成了一个可投稿打磨的母乳-婴儿肠道-SCFA 成熟整合框架。当前最稳健的主文证据是婴儿 SCFA 成熟轴、Bifidobacterium 结构和乳汁 LTF/HMO 到后续 SCFA 的候选时间联系；临床修饰、母体粪便背景和物种-代谢物桥接适合补充或谨慎主文呈现。",
            "BodyNoIndent",
        )
    )
    story.append(P("核心发现一览", "H1CN"))
    core_cards = [
        ("1. 本文真正发现了什么", "高原母婴队列中，早期乳汁LTF/HMO动态、婴儿Bifidobacterium优势结构和粪便SCFA成熟之间形成了有时间顺序的纵向关联框架。", "意义：文章的主线不是‘测了很多组学’，而是这些组学共同指向婴儿肠道生态成熟。"),
        ("2. 为什么不能硬做全组学模型", "同一时间点完整覆盖LTF、HMO、LCFA、SCFA和宏基因组的dyad很少，D90甚至没有完整重叠。", "意义：模块化分析不是退而求其次，而是最符合数据结构的稳健策略。"),
        ("3. 乳汁端最有价值的发现", "LTF覆盖最好，HMO早期变化最强；6-SL、3-SL和LNnT等HMO呈清晰早期下降轨迹。", "意义：LTF/HMO是后续连接婴儿SCFA的最佳上游候选层，LCFA因覆盖少应降权。"),
        ("4. 微生物端最有价值的发现", "Bifidobacterium，尤其B. longum，在D05-D30表现出稳定上升和更高检出。", "意义：婴儿肠道成熟有可识别的优势菌群轨迹，而不是随机波动。"),
        ("5. 代谢端最稳的发现", "SCFA成熟分数随日龄显著上升，并由多个SCFA特征的时间变化共同支撑。", "意义：这是全文最稳的婴儿侧结局，可以作为整合故事的下游锚点。"),
        ("6. 最适合作为主文亮点的整合发现", "早期LTF/HMO与后续SCFA的候选滞后模型强于同时间点模型。", "意义：这让文章从描述性多组学进入‘有时间顺序的候选机制’，但仍不能写成因果调控。"),
        ("7. 母亲粪便数据说明了什么", "部分母婴对有生态相似性，但总体own-pair与背景分布重叠明显。", "意义：它提供母体肠道生态背景，不是垂直传播证明。"),
        ("8. 临床修饰层说明了什么", "母亲BMI与婴儿日龄的交互是最清楚的探索性临床信号，喂养、分娩和出生体重信号不稳定。", "意义：它适合提示下一步临床扩展方向，不宜作为当前文章主结论。"),
    ]
    for title, finding, meaning in core_cards:
        story.append(P(f"<b>{title}</b><br/>{finding}<br/><font color='#315A67'>{meaning}</font>", "Takeaway"))
    story.append(PageBreak())
    story.append(P("战略价值", "H1CN"))
    strategic = [
        ["层面", "审阅判断"],
        ["科学问题", "从“单一母乳成分”升级为“母乳生物活性组分-婴儿肠道微生物-代谢成熟”的纵向生态系统问题。"],
        ["队列特色", "高原母婴人群和多时间点采样构成差异化场景，可作为文章辨识度来源。"],
        ["最强证据", "SCFA 成熟分数随年龄上升；Bifidobacterium 结构成熟；乳铁蛋白/HMO 与后续 SCFA 的候选滞后联系。"],
        ["最大风险", "样本重叠不支持严格完整多组学模型；滞后关联未达到确认性因果证据；物种标签和母婴传播需谨慎。"],
        ["建议定位", "以纵向关联和候选机制框架冲击高水平综合或营养/微生物方向期刊；投稿前补足文献、方法和外部验证叙事。"],
    ]
    table = Table(strategic, colWidths=[3.0 * cm, 13.4 * cm])
    table.setStyle(
        TableStyle(
            [
                ("FONTNAME", (0, 0), (-1, -1), "MSYH"),
                ("FONTNAME", (0, 0), (-1, 0), "MSYH-Bold"),
                ("FONTSIZE", (0, 0), (-1, -1), 8.2),
                ("LEADING", (0, 0), (-1, -1), 11),
                ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#E9F2F6")),
                ("TEXTCOLOR", (0, 0), (-1, 0), colors.HexColor("#102A43")),
                ("GRID", (0, 0), (-1, -1), 0.35, colors.HexColor("#C9D6DF")),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("LEFTPADDING", (0, 0), (-1, -1), 6),
                ("RIGHTPADDING", (0, 0), (-1, -1), 6),
                ("TOPPADDING", (0, 0), (-1, -1), 5),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
            ]
        )
    )
    story.append(table)
    story.append(PageBreak())


def add_main_text(story, figs):
    story.append(P("研究背景与核心问题", "H1CN"))
    takeaway_block(story, "framework")
    story.append(
        P(
            "生命早期肠道生态成熟由母乳暴露、微生物定植、代谢产物积累和宿主临床背景共同驱动。传统研究往往分别讨论母乳营养、婴儿肠道菌群或粪便代谢物，因而难以解释这些层面在同一母婴队列内如何按时间顺序衔接。本研究的价值在于把乳汁生物活性组分、婴儿粪便宏基因组、粪便 SCFA 和母亲粪便宏基因组放入同一个纵向分析框架中，并明确区分强主文证据、候选关联和探索性补充证据。",
        )
    )
    story.append(
        P(
            "在顶刊叙事上，本文不宜被包装成“单一成分导致单一结局”的线性因果故事，而应定位为一个有边界的纵向生态模型：母乳 LTF/HMO 等生物活性成分在早期快速重塑，婴儿肠道 Bifidobacterium 结构与粪便 SCFA 成熟同步推进，母体肠道则提供可比较的生态背景。这个框架的优势是解释力强、图形结构清晰；限制是多组学 complete-case 重叠较小，不能过度声明因果、调控或菌株传播。",
        )
    )
    story.append(figure_block("Graphical abstract", figs["Graphical abstract"], "研究工作模型。图中箭头表示基于纵向采样和候选模型形成的关联主线，而非因果通路证明。"))

    story.append(P("结果一：模块化覆盖决定了合理的分析策略", "H1CN"))
    takeaway_block(story, "fig1")
    story.append(
        P(
            "队列层面，基线临床信息覆盖 52 对母婴。人乳模块中，乳铁蛋白在 D05、D14、D30 和 D90 分别覆盖 51、46、41 和 24 对母婴，是覆盖最广的乳汁特征；HMO 覆盖主要集中在 D05-D30，分别为 37、33 和 37 对母婴；LCFA 覆盖较小，在四个窗口分别为 22、17、6 和 10 对母婴。婴儿粪便 SCFA 覆盖 21、26、26 和 18 对母婴，婴儿粪便宏基因组覆盖 29、25、16 和 9 对母婴。",
        )
    )
    story.append(
        P(
            "严格要求同一时间点同时具备 LTF、HMO、LCFA、SCFA 和宏基因组数据时，可用样本迅速下降：D05 仅 5 对、D14 为 4 对、D30 为 2 对、D90 为 0 对。因此，本文的设计不是把所有数据压缩为一个完整病例模型，而是以模块为单位建立互相呼应的纵向证据链。这个判断非常关键，它让文章的方法选择与数据实际一致，也避免了为追求“全组学整合”而牺牲统计稳定性。",
        )
    )
    story.append(figure_block("Figure 1", figs["Fig. 1"], "队列设计和模块化多组学覆盖。A 展示采样时间线、乳汁/婴儿粪便/母亲粪便模块和分析层；B 以 dyad 为列、模块-时间窗为行展示覆盖度，支持模块化纵向整合策略。", max_h=18.2 * cm))
    story.append(PageBreak())

    story.append(P("结果二：乳汁生物活性模块呈现早期重塑", "H1CN"))
    takeaway_block(story, "fig2")
    story.append(
        P(
            "乳汁侧结果显示，早期母乳不是稳定背景，而是快速变化的暴露层。乳铁蛋白从 D05 到 D30 中位水平轻度下降，D05、D14 和 D30 分别约为 2492.1、2344.2 和 2275.4；动态筛选中乳铁蛋白时间趋势达到 FDR 支持。相比之下，多个 HMO 特征显示更显著的早期重塑，其中 6-SL、3-SL 和 LNnT 为最清楚的下降轨迹，FDR 分别约为 6.55e-10、1.13e-9 和 8.70e-9。",
        )
    )
    story.append(
        P(
            "这一结果为后续整合提供了暴露端基础：乳铁蛋白提供较好的覆盖度和解释稳定性，HMO 提供更强的时间动态。LCFA 也呈现时间变化，但由于与婴儿结局的重叠较少，更适合在本文中作为补充或背景模块，而不是核心推断层。",
        )
    )
    story.append(figure_block("Figure 2", figs["Fig. 2"], "人乳生物活性组分动态。图中整合乳铁蛋白轨迹、HMO/LCFA 模块变化和候选暴露特征，强调早期乳汁暴露层的快速重塑。"))

    story.append(P("结果三：婴儿肠道微生物成熟以 Bifidobacterium 结构为中心", "H1CN"))
    takeaway_block(story, "fig3")
    story.append(
        P(
            "婴儿粪便宏基因组显示出以 Bifidobacterium 为中心的成熟结构。物种层面，Bifidobacterium longum 是最稳定的高丰度标签，其平均丰度从 D05 的 10.9% 上升到 D14 的 16.1% 和 D30 的 30.7%，检出比例也从 D05 的 55.2% 上升到 D30 的 75.0%。D90 时 B. longum 平均丰度较低，但该时间点仅 9 个宏基因组样本，因此不应过度解释。",
        )
    )
    story.append(
        P(
            "B. breve 和 B. pseudocatenulatum 等晚期信号具有潜在生物学意义，但当前 D90 样本数较少，适合描述为探索性现象。顶刊写法中，这一部分最稳妥的贡献不是宣称某个物种的确定性功能，而是说明婴儿肠道成熟并非单一属水平变化，而存在可进一步精修的物种层结构。",
        )
    )
    story.append(figure_block("Figure 3", figs["Fig. 3"], "婴儿肠道微生物成熟。图中展示整体结构变化、主要属/物种轨迹和 Bifidobacterium 物种层差异；物种标签需在投稿前做分类学复核。"))
    story.append(PageBreak())

    story.append(P("结果四：粪便 SCFA 形成可测量的代谢成熟轴", "H1CN"))
    takeaway_block(story, "fig4")
    story.append(
        P(
            "婴儿粪便 SCFA 是本文最稳健的婴儿侧代谢结局。基于 91 个样本、32 对母婴构建的 SCFA 成熟分数随日龄增加而上升，线性模型斜率为 0.0081/天，P = 2.11e-4，R2 = 0.144。这一结果为文章提供了一个可被图形和模型共同支撑的代谢成熟轴。",
        )
    )
    story.append(
        P(
            "单一 SCFA 层面，琥珀酸、异丁酸、丁酸、2-甲基丁酸、戊酸和异戊酸均表现出 FDR 支持的时间相关变化；乙酸呈名义上升但未达到相同 FDR 阈值，丙酸没有清晰时间趋势。这里需要保持术语克制：该分数是数据驱动的成熟指标，不是已验证的临床成熟指数。",
        )
    )
    story.append(figure_block("Figure 4", figs["Fig. 4"], "婴儿粪便 SCFA 成熟。图中展示 SCFA 特征轨迹、复合成熟分数和时间趋势，构成本文婴儿代谢端的核心证据。"))

    story.append(P("结果五：早期乳汁 LTF/HMO 与后续 SCFA 存在候选时间联系", "H1CN"))
    takeaway_block(story, "fig5")
    story.append(
        P(
            "在整合层面，本文最有潜力的主结果是乳汁特征与后续婴儿粪便 SCFA 的候选滞后联系。滞后模型更符合“早期乳汁暴露 - 后续婴儿代谢表型”的时间顺序。40 个候选滞后基础模型中，14 个达到名义 P < 0.05，20 个达到名义 P < 0.10；在具备临床协变量的调整模型中，6 个仍达到名义 P < 0.05，11 个达到名义 P < 0.10。相比之下，同时间点模型在临床调整后没有候选项保留名义 P < 0.05。",
        )
    )
    story.append(
        P(
            "代表性信号包括 D05 乳铁蛋白与 D30 异戊酸的正相关候选联系（N = 26，beta = 0.538，P = 0.00457），以及 D14 LNT 与 D90 乙酸的正相关候选联系（N = 15，beta = 0.684，P = 0.00489）。此外，D14 LNFP 相关特征、3-SL 和 6-SL 与 D30 琥珀酸呈候选负向联系，D14 LDFT 与 D90 2-甲基丁酸和异丁酸呈候选负向联系。该部分应作为“候选时间联系”来写，不应声称因果调控或中介机制。",
        )
    )
    story.append(figure_block("Figure 5", figs["Fig. 5"], "乳汁 LTF/HMO 与后续婴儿 SCFA 的候选滞后联系。图中突出基础模型、临床调整模型和候选证据计数；该结果是整合叙事核心，但仍属候选关联。"))
    story.append(PageBreak())

    story.append(P("结果六：母亲粪便宏基因组提供生态背景，而非传播证明", "H1CN"))
    takeaway_block(story, "figs1")
    story.append(
        P(
            "母亲粪便宏基因组用于评估婴儿肠道装配的母体生态背景。73 个真实母婴粪便比较被置于非配对母体背景分布中评估，部分 dyad 显示高于背景的 own-pair 相似性，但总体 own-dyad 与 other-dyad 分布重叠明显。当前汇总中，24/73 个比较具有正向 own-versus-other Z 分数，49/73 为负向。",
        )
    )
    story.append(
        P(
            "因此，这一层的价值不是证明垂直传播，而是证明研究已经把母体肠道生态作为可比较背景纳入框架。这种表述既保留了母亲粪便数据的工作量和生态意义，也避免了宏基因组分辨率不足时最容易被审稿人质疑的传播叙事。",
        )
    )
    story.append(figure_block("Figure S1", figs["Fig. S1"], "母亲粪便宏基因组背景。图中展示 own-vs-other 相似性、dyad 层差异和覆盖度；应作为生态背景证据，而非菌株传播证明。"))

    story.append(P("结果七：临床修饰层增强应用意义，但应保持探索性", "H1CN"))
    takeaway_block(story, "figs2")
    story.append(
        P(
            "临床修饰模型显示，婴儿日龄仍是 SCFA 成熟分数最稳定的解释变量。母亲 BMI 与日龄的交互是最清楚的探索性临床信号（DayScaled × maternal BMI z：beta = 0.622，P = 0.0137；N = 46，dyads = 17）。相比之下，喂养方式、分娩方式和出生体重在当前模型中没有形成稳定主效应。",
        )
    )
    story.append(
        P(
            "这一结果可以提升文章的临床可读性，但不宜推到主结论层面。建议在正文讨论中简要提及，在补充图中呈现完整模型森林图、时间点对比和覆盖热图。这样既回应领导和审稿人对临床意义的期待，又不让小样本修饰模型承担过重结论。",
        )
    )
    story.append(figure_block("Figure S2", figs["Fig. S2"], "临床修饰层。图中展示 BMI-age 模型、混合模型项、时间点对比和临床覆盖；定位为探索性支持。"))

    story.append(P("讨论：顶刊叙事应如何收束", "H1CN"))
    takeaway_block(story, "discussion")
    story.append(
        P(
            "本文最强的顶刊潜力来自三点组合。第一，研究对象不是普通单点营养队列，而是具有高原母婴人群特色的纵向多组学框架。第二，结果不是零散相关，而是形成了从乳汁生物活性重塑到婴儿微生物和代谢成熟的时间顺序。第三，文章主动承认证据边界，把 complete-case 限制转化为模块化整合策略，而不是掩盖样本结构问题。",
        )
    )
    story.append(
        P(
            "需要避免的写法同样明确：不要写“母乳成分调控 SCFA 成熟”，不要写“母体菌群垂直传播驱动婴儿肠道”，也不要把 SCFA 成熟分数包装成临床诊断工具。更好的表述是：在高原母婴队列中，早期乳汁 LTF/HMO 动态、婴儿 Bifidobacterium 优势结构和粪便 SCFA 成熟之间存在可复核的纵向关联网络，其中若干候选滞后联系值得在更大样本或机制研究中验证。",
        )
    )
    story.append(
        P(
            "如果后续目标是真正冲击顶刊，建议补三类工作：一是正式文献支持和 novelty 定位，尤其是高海拔/民族地区母婴营养、HMO/LTF、Bifidobacterium 和 SCFA 成熟的交叉空白；二是方法细节和质控透明化，包括宏基因组物种标签复核、SCFA 成熟分数构建、滞后模型多重检验策略；三是把补充分析组织成审稿人容易检查的 source-data package，而不是只靠主文叙事。",
        )
    )

    story.append(P("方法概述", "H1CN"))
    story.append(
        P(
            "本审阅稿基于项目内已完成分析输出整理。研究对象为青藏高原母婴队列，核心时间窗为 D05、D14、D30 和 D90。数据层包括基线/围产期临床信息、人乳生物活性成分（LTF、HMO、LCFA）、婴儿粪便宏基因组、婴儿粪便 SCFA、母亲粪便宏基因组及临床修饰变量。由于各模块覆盖度不完全一致，分析采用模块化策略：先分别建立乳汁动态、婴儿微生物成熟和 SCFA 成熟结果，再用重叠样本评估候选滞后关联和补充桥接模型。",
        )
    )
    story.append(
        P(
            "统计解释遵循预先设定的证据边界：FDR 支持的时间趋势可作为主文发现；名义显著但未获 FDR 支持的跨层模型仅作为候选关联；临床修饰和物种-代谢物桥接因样本数较小，定位为探索性支持。正式投稿前仍需补充完整伦理信息、实验方法、软件版本、外部文献引用和参考文献列表。",
        )
    )


def build_pdf():
    figs = render_figures()
    doc = BaseDocTemplate(
        str(PDF_OUT),
        pagesize=A4,
        leftMargin=1.7 * cm,
        rightMargin=1.7 * cm,
        topMargin=1.45 * cm,
        bottomMargin=1.55 * cm,
        title="milk-gut-SCFA leadership manuscript",
        author="",
    )
    frame = Frame(
        doc.leftMargin,
        doc.bottomMargin,
        doc.width,
        doc.height,
        id="normal",
    )
    doc.addPageTemplates(
        [
            PageTemplate(id="cover", frames=[frame], onPage=cover_footer),
            PageTemplate(id="body", frames=[frame], onPage=header_footer),
        ]
    )
    story = [NextPageTemplate("body")]
    add_cover(story, figs)
    add_summary(story)
    add_main_text(story, figs)
    doc.build(story)
    print(PDF_OUT)


if __name__ == "__main__":
    build_pdf()
