#!/usr/bin/env python
"""
parse Katoh & Suga 2023 NAR 51:10642, Supplementary Figure 8 (pdf pages 12-15)
into one row per codon x non-cognate amino acid (61 x 19 = 1159).

the figure is vector text, so each cell is read as the literal glyph the
authors typeset (no pixel digitising). each codon block has five rows:
  intensity : 1..5 printed, blank = 0 (no peak detected)
  1st/2nd/3rd : codon position match '-', mismatch '+', '±' = isoacceptors disagree
  total     : number of mismatches, e.g. '1', '2', '1/2', '1/2/3'
cognate columns are grey with no glyphs and are dropped.

output: ../katoh2023_s8_parsed.csv
"""
import re
import fitz
import pandas as pd

PDF = "raw/katoh2023/suppl/Supplementary Figures.pdf"
PAGES = [11, 12, 13, 14]            # 0-based index of pdf pages 12-15
AAS = list("FLIMVSPTAYHQNKDECWRG")
ROWS = ["intensity", "1st", "2nd", "3rd", "total"]
CODON = re.compile(r"^[UCAG]{3}$")


def parsePage(page):
    words = page.get_text("words")
    # column centres from the single-letter aa header row (topmost y)
    # (the U/C/A/G group banner and panel letter also match, so take the y
    # line holding the most aa letters)
    heads = [w for w in words if w[4] in AAS and len(w[4]) == 1]
    ys = [round(w[1]) for w in heads]
    ytop = max(set(ys), key=ys.count)
    cols = {w[4]: (w[0] + w[2]) / 2 for w in heads if round(w[1]) == ytop}
    assert len(cols) == 20, cols
    # row anchors: each row label word gives a y; codon label sits on the intensity row
    labels = [w for w in words if w[4] in ROWS]
    codons = sorted([w for w in words if CODON.match(w[4])], key=lambda w: w[1])
    cells = {}
    for lab in labels:
        yc = (lab[1] + lab[3]) / 2
        # owning codon = last codon label at or above this row
        own = [c for c in codons if (c[1] + c[3]) / 2 <= yc + 1]
        if not own:
            continue
        codon = own[-1][4]
        for w in words:
            if w[0] <= lab[2] or abs((w[1] + w[3]) / 2 - yc) > 1.5:
                continue
            xc = (w[0] + w[2]) / 2
            aa = min(cols, key=lambda a: abs(cols[a] - xc))
            assert abs(cols[aa] - xc) < 5, (codon, lab[4], w)
            key = (codon, aa, lab[4])
            assert key not in cells, key
            cells[key] = w[4]
    return cells


def main():
    doc = fitz.open(PDF)
    cells = {}
    for p in PAGES:
        cells.update(parsePage(doc[p]))
    codons = sorted({k[0] for k in cells})
    recs = []
    for c in codons:
        for aa in AAS:
            v = {r: cells.get((c, aa, r)) for r in ROWS}
            # cognate column carries no position glyphs at all
            if v["1st"] is None and v["2nd"] is None and v["3rd"] is None:
                assert v["intensity"] is None and v["total"] is None, (c, aa, v)
                continue
            assert all(v[r] is not None for r in ROWS[1:]), (c, aa, v)
            recs.append(dict(codon=c, aa=aa,
                             intensity=int(v["intensity"]) if v["intensity"] else 0,
                             p1=v["1st"], p2=v["2nd"], p3=v["3rd"], total=v["total"]))
    df = pd.DataFrame(recs)
    df.to_csv("katoh2023_s8_parsed.csv", index=False)
    print(f"codons={len(codons)} combos={len(df)}")


if __name__ == "__main__":
    main()
