from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

from PIL import Image, ImageDraw


ROOT = Path(os.environ.get("HCC_AE_ROOT", Path.cwd())).resolve()
DEFAULT_DOCX = ROOT / "hcc_ae_submission_package" / "HCC_AE_submission.docx"
DEFAULT_OUT = ROOT / "hcc_ae_submission_package" / "qa_pdf_pages"
SOFFICE = Path(
    os.environ.get(
        "LIBREOFFICE_SOFFICE",
        r"C:\Program Files\LibreOffice\program\soffice.com",
    )
)
LOCAL_PDFIUM = ROOT / ".pdf_render_lib"


def ensure_pdfium(install: bool):
    if LOCAL_PDFIUM.exists():
        sys.path.insert(0, str(LOCAL_PDFIUM))
    try:
        import pypdfium2 as pdfium  # type: ignore

        return pdfium
    except ImportError:
        if not install:
            raise RuntimeError(
                "pypdfium2 is not installed. Install requirements.txt or pass --install-pdfium."
            )
        LOCAL_PDFIUM.mkdir(parents=True, exist_ok=True)
        subprocess.run(
            [
                sys.executable,
                "-m",
                "pip",
                "install",
                "--target",
                str(LOCAL_PDFIUM),
                "pypdfium2",
            ],
            check=True,
        )
        sys.path.insert(0, str(LOCAL_PDFIUM))
        import pypdfium2 as pdfium  # type: ignore

        return pdfium


def prepare_output_dir(out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    for pattern in ("*.pdf", "page_*.png", "contact_sheet.jpg", "render_manifest.txt"):
        for path in out_dir.glob(pattern):
            if path.is_file():
                path.unlink()
            elif path.is_dir():
                shutil.rmtree(path)


def convert_docx_to_pdf(docx: Path, out_dir: Path) -> Path:
    if not SOFFICE.exists():
        raise FileNotFoundError(f"LibreOffice was not found: {SOFFICE}")
    cmd = [
        str(SOFFICE),
        "--headless",
        "--norestore",
        "--convert-to",
        "pdf",
        "--outdir",
        str(out_dir),
        str(docx),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.stdout.strip():
        print(result.stdout.strip())
    if result.stderr.strip():
        print(result.stderr.strip())
    if result.returncode != 0:
        raise SystemExit(result.returncode)

    pdf = out_dir / f"{docx.stem}.pdf"
    if not pdf.exists():
        raise FileNotFoundError(pdf)
    return pdf


def render_pdf_pages(pdfium, pdf: Path, out_dir: Path, scale: float) -> list[Path]:
    doc = pdfium.PdfDocument(str(pdf))
    rendered: list[Path] = []
    for idx in range(len(doc)):
        page = doc[idx]
        bitmap = page.render(scale=scale)
        image = bitmap.to_pil()
        page_path = out_dir / f"page_{idx + 1:02d}.png"
        image.save(page_path, optimize=True)
        rendered.append(page_path)
    return rendered


def make_contact_sheet(page_paths: list[Path], out_dir: Path, thumb_width: int = 360) -> Path:
    thumbs = []
    for path in page_paths:
        image = Image.open(path).convert("RGB")
        ratio = thumb_width / image.width
        thumb = image.resize((thumb_width, max(1, int(image.height * ratio))))
        thumbs.append((path, thumb))

    if not thumbs:
        raise ValueError("No pages were rendered.")

    label_height = 28
    gap = 18
    columns = 2 if len(thumbs) > 1 else 1
    rows = (len(thumbs) + columns - 1) // columns
    cell_width = thumb_width
    cell_height = max(thumb.height for _, thumb in thumbs) + label_height
    sheet = Image.new(
        "RGB",
        (
            columns * cell_width + (columns + 1) * gap,
            rows * cell_height + (rows + 1) * gap,
        ),
        "white",
    )
    draw = ImageDraw.Draw(sheet)
    for idx, (path, thumb) in enumerate(thumbs):
        row = idx // columns
        col = idx % columns
        x = gap + col * (cell_width + gap)
        y = gap + row * (cell_height + gap)
        draw.text((x, y), path.stem, fill=(40, 40, 40))
        sheet.paste(thumb, (x, y + label_height))

    contact = out_dir / "contact_sheet.jpg"
    sheet.save(contact, quality=82, optimize=True)
    return contact


def write_manifest(out_dir: Path, docx: Path, pdf: Path, page_paths: list[Path], contact: Path) -> Path:
    manifest = out_dir / "render_manifest.txt"
    lines = [
        "DOCX/PDF QA render completed.",
        "",
        f"DOCX: {docx}",
        f"PDF: {pdf}",
        f"CONTACT_SHEET: {contact}",
        f"PAGE_COUNT: {len(page_paths)}",
        "",
        "Rendered pages:",
    ]
    lines.extend(str(path) for path in page_paths)
    lines.extend(
        [
            "",
            "Safety note:",
            "Rendered pages may contain controlled manuscript text.",
            "Keep this QA directory outside public Git commits.",
        ]
    )
    manifest.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Render a controlled DOCX to PDF and lightweight QA images on Windows."
    )
    parser.add_argument("--docx", default=str(DEFAULT_DOCX), help="DOCX file to render")
    parser.add_argument("--out", default=str(DEFAULT_OUT), help="Output directory")
    parser.add_argument("--scale", type=float, default=1.15, help="PDF render scale for page PNGs")
    parser.add_argument(
        "--install-pdfium",
        action="store_true",
        help="Install pypdfium2 into a local cache under HCC_AE_ROOT if missing.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    docx = Path(args.docx).expanduser().resolve()
    out_dir = Path(args.out).expanduser().resolve()
    if not docx.exists():
        raise FileNotFoundError(docx)

    prepare_output_dir(out_dir)
    pdfium = ensure_pdfium(args.install_pdfium)
    pdf = convert_docx_to_pdf(docx, out_dir)
    page_paths = render_pdf_pages(pdfium, pdf, out_dir, args.scale)
    contact = make_contact_sheet(page_paths, out_dir)
    manifest = write_manifest(out_dir, docx, pdf, page_paths, contact)

    print(f"PDF: {pdf}")
    print(f"Pages rendered: {len(page_paths)}")
    print(f"Contact sheet: {contact}")
    print(f"Manifest: {manifest}")
    print("Public-safety note: QA images may contain controlled manuscript text.")


if __name__ == "__main__":
    main()
