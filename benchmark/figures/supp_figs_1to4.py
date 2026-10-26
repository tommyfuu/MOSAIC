"""Assemble Supplemental Figures 1-4 by stacking the Fig 2-style pages (fig2_style_plots.py) for OR 1, 1.25 and
1.5 vertically on one page, as in the bioRxiv v1 supplement. Vector content and editable text are kept.
  Supp Fig 1: count, batch-library size confounding      Supp Fig 2: count, no confounding
  Supp Fig 3: relative abundance, confounding            Supp Fig 4: relative abundance, no confounding
Needs PyMuPDF (fitz). Usage: python3 supp_figs_1to4.py <fig2_style_dir> <out_dir>
"""
import os
import sys

import fitz

SRC, OUT = sys.argv[1:3]
os.makedirs(OUT, exist_ok=True)
FIGS = {1: ('count', 'libsize_confounded'), 2: ('count', 'no_libsize_confounding'),
        3: ('relab', 'libsize_confounded'), 4: ('relab', 'no_libsize_confounding')}
for n, (dt, conf) in FIGS.items():
    pages = [fitz.open(f'{SRC}/Fig2style_{dt}_OR{OR}_{conf}.pdf') for OR in ('1', '1.25', '1.5')]
    w, h = pages[0][0].rect.width, pages[0][0].rect.height
    gap = 18
    out = fitz.open()
    page = out.new_page(width=w, height=3 * h + 2 * gap)
    for i, src in enumerate(pages):
        y0 = i * (h + gap)
        page.show_pdf_page(fitz.Rect(0, y0, w, y0 + h), src, 0)
        if i:
            page.draw_line((40, y0 - gap / 2), (w - 40, y0 - gap / 2), color=(0.75, 0.75, 0.75), width=0.6)
    name = f'SuppFig{n}_{dt}_{conf}_OR1_1.25_1.5.pdf'
    out.save(f'{OUT}/{name}', garbage=3, deflate=True)
    print(name)
