"""
build the submission-ready pdf deliverables for the jme package.

writes three files into jme/:

    figures_JME.pdf                   figures 1-4, one per page, caption beneath
    supplementary_materials_JME.pdf   tables s1/s2, figure s1, file manifest
    figures_and_supplementary_JME.pdf the two above concatenated

two deliberate toolchain choices:

1. figure artwork is taken from figures/figN.pdf (vector), never from the
   300 dpi png. text inside the figures therefore stays selectable and
   scales without pixelation, which is what jme wants for line art.

2. prose and tables go through pandoc -> docx -> libreoffice. this machine
   has no xelatex or lualatex, and pandoc's pdflatex path dies inside the
   bookmark package ("unsupported driver pdftex"). libreoffice is the only
   local renderer that handles the unicode this manuscript is full of:
   superscripts, the norm bars, greek, arrows, and u+2212 minus.

captions are composed with pymupdf using dejavu sans, which was verified to
carry every non-ascii glyph the captions use.

run from the repo root:

    python scripts/build_submission_pdfs.py
"""

from __future__ import annotations

import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import fitz  # pymupdf
from pypdf import PdfWriter

REPO = Path(__file__).resolve().parent.parent
OUT_DIR = REPO / "jme"

# us letter, in points
PAGE_W, PAGE_H = 612.0, 792.0
MARGIN = 54.0
CAPTION_SIZE = 8.6
CAPTION_LEAD = 11.2
FIG_CAPTION_GAP = 20.0

FONT_REGULAR = Path("/usr/share/fonts/dejavu-sans-fonts/DejaVuSans.ttf")
FONT_BOLD = Path("/usr/share/fonts/dejavu-sans-fonts/DejaVuSans-Bold.ttf")

# main-text figures, in manuscript order
MAIN_FIGURES = [
    ("figures/fig1.pdf", "figures/fig1_caption.txt"),
    ("figures/fig2.pdf", "figures/fig2_caption.txt"),
    ("figures/fig3.pdf", "figures/fig3_caption.txt"),
    ("figures/fig4.pdf", "figures/fig4_caption.txt"),
]

SUPP_FIGURE = (
    "figures/figS1_raw_distributions.pdf",
    "figures/figS1_caption.txt",
)


def readCaption(relPath: str) -> str:
    """load a caption file and unwrap its hard line breaks into one paragraph."""
    raw = (REPO / relPath).read_text(encoding="utf-8")
    return re.sub(r"\s+", " ", raw).strip()


def splitLabel(caption: str) -> tuple[str, str]:
    """separate a leading 'Figure 3.' / 'Figure S1.' label from the caption body.

    the label is drawn in bold, so it has to be measured with the bold font
    while the body is measured with the regular one.
    """
    match = re.match(r"^(Figure\s+S?\d+\.)\s*(.*)$", caption, flags=re.S)
    if not match:
        return "", caption
    return match.group(1), match.group(2)


def wrapText(
    text: str,
    font: fitz.Font,
    size: float,
    fullWidth: float,
    firstWidth: float | None = None,
) -> list[str]:
    """greedy word wrap using real glyph metrics.

    firstWidth lets the opening line be narrower than the rest, which is how
    room is made for the bold 'Figure N.' label sitting on that same line.
    """
    limitFirst = fullWidth if firstWidth is None else firstWidth
    lines: list[str] = []
    current = ""
    for word in text.split():
        limit = limitFirst if not lines else fullWidth
        candidate = f"{current} {word}".strip()
        if current and font.text_length(candidate, fontsize=size) > limit:
            lines.append(current)
            current = word
        else:
            current = candidate
    if current:
        lines.append(current)
    return lines


def composeFigurePage(outDoc: fitz.Document, figRel: str, capRel: str) -> None:
    """append one page holding a vector figure with its caption underneath.

    the caption is laid out first, because its height determines how much
    vertical space is left for the artwork. the figure is then scaled to fit
    that remaining box while preserving its aspect ratio -- show_pdf_page
    stretches to whatever rect it is given, so the rect must already have the
    correct proportions or the figure comes out distorted.
    """
    regular = fitz.Font(fontfile=str(FONT_REGULAR))
    bold = fitz.Font(fontfile=str(FONT_BOLD))

    caption = readCaption(capRel)
    label, body = splitLabel(caption)

    usableW = PAGE_W - 2 * MARGIN
    labelW = bold.text_length(label, fontsize=CAPTION_SIZE) if label else 0.0
    spaceW = regular.text_length(" ", fontsize=CAPTION_SIZE) if label else 0.0

    lines = wrapText(
        body,
        regular,
        CAPTION_SIZE,
        usableW,
        firstWidth=usableW - labelW - spaceW,
    )
    captionH = len(lines) * CAPTION_LEAD

    # whatever vertical space the caption does not claim belongs to the figure
    figMaxH = PAGE_H - 2 * MARGIN - captionH - FIG_CAPTION_GAP
    if figMaxH < 120:
        raise RuntimeError(f"caption for {figRel} leaves too little room for art")

    src = fitz.open(REPO / figRel)
    srcRect = src[0].rect
    scale = min(usableW / srcRect.width, figMaxH / srcRect.height)
    figW, figH = srcRect.width * scale, srcRect.height * scale

    page = outDoc.new_page(width=PAGE_W, height=PAGE_H)
    figX = MARGIN + (usableW - figW) / 2  # centre horizontally
    figRect = fitz.Rect(figX, MARGIN, figX + figW, MARGIN + figH)
    page.show_pdf_page(figRect, src, 0)
    src.close()

    # baseline of the first caption line
    y = MARGIN + figH + FIG_CAPTION_GAP + CAPTION_SIZE
    for index, line in enumerate(lines):
        x = MARGIN
        if index == 0 and label:
            page.insert_text(
                (x, y),
                label,
                fontname="djb",
                fontfile=str(FONT_BOLD),
                fontsize=CAPTION_SIZE,
            )
            x += labelW + spaceW
        page.insert_text(
            (x, y),
            line,
            fontname="djr",
            fontfile=str(FONT_REGULAR),
            fontsize=CAPTION_SIZE,
        )
        y += CAPTION_LEAD


def buildReferenceDocx(outPath: Path) -> Path:
    """derive a pandoc reference document with manuscript-appropriate headings.

    pandoc's stock reference.docx styles headings with the theme accent colour
    (0F4761, a dark teal) at 20 pt and 16 pt. that reads as a slide deck rather
    than a journal supplement, so the colour is forced to black and the sizes
    are brought down. sizes are half-points, hence 28 = 14 pt.

    generated at build time rather than committed, so no opaque binary blob
    lives in the repository.
    """
    raw = subprocess.run(
        ["pandoc", "--print-default-data-file", "reference.docx"],
        check=True,
        capture_output=True,
    ).stdout

    with tempfile.TemporaryDirectory() as tmp:
        tmpDir = Path(tmp)
        stock = tmpDir / "stock.docx"
        stock.write_bytes(raw)
        unpacked = tmpDir / "unpacked"
        unpacked.mkdir()
        subprocess.run(
            ["unzip", "-q", str(stock), "-d", str(unpacked)],
            check=True,
            capture_output=True,
        )

        stylesPath = unpacked / "word" / "styles.xml"
        styles = stylesPath.read_text(encoding="utf-8")

        # the colour tag spans newlines in the stock file, so match loosely.
        # bold is added at the same time: pandoc's heading styles carry no
        # <w:b/>, relying on size and the accent colour alone to signal
        # hierarchy, so once the colour goes black the headings would
        # otherwise be indistinguishable from body text.
        styles = re.sub(
            r'<w:color\s+w:val="0F4761"[^/]*/>',
            '<w:color w:val="000000"/><w:b/>',
            styles,
        )
        for oldSize, newSize in (("40", "28"), ("32", "24")):
            styles = styles.replace(
                f'<w:sz w:val="{oldSize}" />\n      <w:szCs w:val="{oldSize}" />',
                f'<w:sz w:val="{newSize}" /><w:szCs w:val="{newSize}" />',
            )
        stylesPath.write_text(styles, encoding="utf-8")

        outPath.parent.mkdir(parents=True, exist_ok=True)
        if outPath.exists():
            outPath.unlink()
        subprocess.run(
            ["zip", "-q", "-r", "-X", str(outPath), "."],
            cwd=unpacked,
            check=True,
            capture_output=True,
        )
    return outPath


def renderMarkdown(markdown: str, outPdf: Path, referenceDocx: Path | None = None) -> Path:
    """render markdown to pdf via pandoc -> docx -> libreoffice.

    libreoffice is handed a private profile directory; without it a second
    invocation can collide with an existing user profile lock and silently
    produce nothing.
    """
    with tempfile.TemporaryDirectory() as tmp:
        tmpDir = Path(tmp)
        mdPath = tmpDir / "part.md"
        docxPath = tmpDir / "part.docx"
        mdPath.write_text(markdown, encoding="utf-8")

        pandocCmd = ["pandoc", str(mdPath), "-o", str(docxPath)]
        if referenceDocx is not None:
            pandocCmd += ["--reference-doc", str(referenceDocx)]
        subprocess.run(pandocCmd, check=True, capture_output=True)
        subprocess.run(
            [
                "soffice",
                "-env:UserInstallation=file://" + str(tmpDir / "loprofile"),
                "--headless",
                "--convert-to",
                "pdf",
                "--outdir",
                str(tmpDir),
                str(docxPath),
            ],
            check=True,
            capture_output=True,
            timeout=300,
        )

        produced = tmpDir / "part.pdf"
        if not produced.exists():
            raise RuntimeError(f"libreoffice produced no pdf for {outPdf.name}")
        outPdf.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(produced, outPdf)
    return outPdf


def mergePdfs(parts: list[Path], outPath: Path) -> None:
    writer = PdfWriter()
    for part in parts:
        writer.append(str(part))
    with open(outPath, "wb") as handle:
        writer.write(handle)
    writer.close()


def buildFigures(workDir: Path) -> Path:
    """figures 1-4, one per page."""
    doc = fitz.open()
    for figRel, capRel in MAIN_FIGURES:
        composeFigurePage(doc, figRel, capRel)
    outPath = OUT_DIR / "figures_JME.pdf"
    doc.save(str(outPath))
    doc.close()
    return outPath


def buildSupplementary(workDir: Path) -> Path:
    """tables s1/s2, then figure s1 as a composed page, then the manifest.

    the markdown is cut at its own section headings so the composed vector
    figure page lands in document order rather than being bolted onto the end.
    """
    source = (REPO / "supplementary_materials_JME.md").read_text(encoding="utf-8")

    figHeading = "## Figure S1"
    manifestHeading = "## File manifest"
    if figHeading not in source or manifestHeading not in source:
        raise RuntimeError("supplementary_materials_JME.md headings moved; update the split")

    head = source.split(figHeading)[0].rstrip() + "\n"
    tail = manifestHeading + source.split(manifestHeading, 1)[1]

    reference = buildReferenceDocx(workDir / "reference.docx")
    headPdf = renderMarkdown(head, workDir / "supp_head.pdf", reference)

    figDoc = fitz.open()
    composeFigurePage(figDoc, *SUPP_FIGURE)
    figPdf = workDir / "supp_figS1.pdf"
    figDoc.save(str(figPdf))
    figDoc.close()

    tailPdf = renderMarkdown(tail, workDir / "supp_tail.pdf", reference)

    outPath = OUT_DIR / "supplementary_materials_JME.pdf"
    mergePdfs([headPdf, figPdf, tailPdf], outPath)
    return outPath


def main() -> int:
    for font in (FONT_REGULAR, FONT_BOLD):
        if not font.exists():
            print(f"missing font: {font}", file=sys.stderr)
            return 1

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        workDir = Path(tmp)
        figuresPdf = buildFigures(workDir)
        supplementaryPdf = buildSupplementary(workDir)

        combinedPdf = OUT_DIR / "figures_and_supplementary_JME.pdf"
        mergePdfs([figuresPdf, supplementaryPdf], combinedPdf)

    for path in (figuresPdf, supplementaryPdf, combinedPdf):
        pages = fitz.open(path).page_count
        size = path.stat().st_size / 1024
        print(f"{path.relative_to(REPO)}  {pages} pages  {size:.0f} KB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
