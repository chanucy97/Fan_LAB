import os
from pathlib import Path

from PIL import Image
from pypdf import PdfReader, PdfWriter
from reportlab.lib.utils import ImageReader
from reportlab.pdfgen import canvas


ROOT = Path(os.environ.get("FANLAB_RESEARCH2_ROOT", Path.cwd())).resolve()
PKG = ROOT / "final_manuscript_planning_20260605"
OUT_DIR = PKG / "compiled_figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)


FIGURES = [
    (
        "Graphical abstract",
        PKG / "graphical_abstract" / "graphical_abstract_model_summary_user_final.png",
        "image",
    ),
    (
        "Fig. 1",
        PKG
        / "figures_redrawn"
        / "fig1_cohort_design_overview"
        / "Fig1_complete_cohort_coverage.pdf",
        "pdf",
    ),
    (
        "Fig. 2",
        PKG
        / "figures_redrawn"
        / "fig2_milk_bioactive_remodeling"
        / "Fig2_milk_bioactive_remodeling_redrawn.pdf",
        "pdf",
    ),
    (
        "Fig. 3",
        PKG
        / "figures_redrawn"
        / "fig3_infant_microbiome_maturation"
        / "Fig3_infant_microbiome_maturation_redrawn.pdf",
        "pdf",
    ),
    (
        "Fig. 4",
        PKG
        / "figures_redrawn"
        / "fig4_scfa_maturation"
        / "Fig4_scfa_maturation_redrawn.pdf",
        "pdf",
    ),
    (
        "Fig. 5",
        PKG
        / "figures_redrawn"
        / "fig5_milk_scfa_lagged_candidates"
        / "Fig5_milk_scfa_lagged_candidates_redrawn.pdf",
        "pdf",
    ),
    (
        "Fig. S1",
        PKG
        / "figures_redrawn"
        / "fig6_maternal_fecal_context"
        / "Fig6_maternal_fecal_context_redrawn.pdf",
        "pdf",
    ),
    (
        "Fig. S2",
        PKG
        / "figures_redrawn"
        / "fig7_clinical_modifier_layer"
        / "Fig7_clinical_modifier_layer_redrawn.pdf",
        "pdf",
    ),
]


def image_to_pdf(image_path: Path, pdf_path: Path) -> None:
    with Image.open(image_path) as img:
        img = img.convert("RGB")
        width_px, height_px = img.size

    width_pt = width_px * 72 / 300
    height_pt = height_px * 72 / 300
    c = canvas.Canvas(str(pdf_path), pagesize=(width_pt, height_pt))
    c.drawImage(
        ImageReader(str(image_path)),
        0,
        0,
        width=width_pt,
        height=height_pt,
        preserveAspectRatio=True,
        mask="auto",
    )
    c.showPage()
    c.save()


def main() -> None:
    missing = [str(path) for _, path, _ in FIGURES if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing figure input(s):\n" + "\n".join(missing))

    ga_pdf = OUT_DIR / "00_graphical_abstract.pdf"
    image_to_pdf(FIGURES[0][1], ga_pdf)

    writer = PdfWriter()
    manifest_lines = ["Page,Figure,Source"]
    for index, (label, path, kind) in enumerate(FIGURES, start=1):
        input_pdf = ga_pdf if kind == "image" else path
        reader = PdfReader(str(input_pdf))
        if len(reader.pages) < 1:
            raise ValueError(f"No pages found in {input_pdf}")
        writer.add_page(reader.pages[0])
        manifest_lines.append(f'{index},"{label}","{path}"')

    out_pdf = OUT_DIR / "all_figures_one_pdf_each_page.pdf"
    with out_pdf.open("wb") as f:
        writer.write(f)

    manifest = OUT_DIR / "all_figures_one_pdf_manifest.csv"
    manifest.write_text("\n".join(manifest_lines) + "\n", encoding="utf-8")

    reader = PdfReader(str(out_pdf))
    print(f"Wrote: {out_pdf}")
    print(f"Pages: {len(reader.pages)}")
    print(f"Manifest: {manifest}")


if __name__ == "__main__":
    main()
