# Fig S6 (docking poses + interaction maps) — combine the author's pre-rendered 2x5 panels
# (rows = quercetin/kaempferol; cols = AKT1/IL1B/IL6/MMP9/TP53) into one supplementary figure:
#   A = 3D binding poses (PyMOL), B = 2D ligand-residue interaction diagrams (LigPlot style).
import os
from PIL import Image, ImageDraw, ImageFont

ROOT = r"H:\毕业设计\网药部分\JingTong_CSR_paper"
DOCK = os.path.join(ROOT, "JINGTONG_RESULTS", "docking_figures")
OUT = os.path.join(ROOT, "figures", "composite", "FigS6_docking.png")
try:
    FONT = ImageFont.truetype("C:/Windows/Fonts/arialbd.ttf", 60)
except Exception:
    FONT = ImageFont.load_default()

PAD, LAB, GAP, BG = 46, 78, 54, (255, 255, 255)

p3d = Image.open(os.path.join(DOCK, "panel_3D_all.png")).convert("RGB")
p2d = Image.open(os.path.join(DOCK, "panel_2D_all.png")).convert("RGB")
W = max(p3d.width, p2d.width)                      # both 3172
H = PAD * 2 + LAB + p3d.height + GAP + LAB + p2d.height
canvas = Image.new("RGB", (W + 2 * PAD, H), BG)
d = ImageDraw.Draw(canvas)

y = PAD
for letter, im in [("A", p3d), ("B", p2d)]:
    d.text((PAD, y), letter, font=FONT, fill=(0, 0, 0))
    canvas.paste(im, (PAD + (W - im.width) // 2, y + LAB))
    y += LAB + im.height + GAP

canvas.save(OUT, dpi=(200, 200))
print("wrote FigS6_docking.png", canvas.size)
