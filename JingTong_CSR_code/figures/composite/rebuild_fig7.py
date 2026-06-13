# Fig 7 (single-cell) — UMAP coords/adata are server-only, so panels are composited from the
# pre-rendered PNGs. Panel A: the source 3-in-1 UMAP strip (umap_clusters.png) has its legends
# overlapping across panel boundaries (a ghost marker column bleeds into cell_type, and the long
# cell_type labels overflow into the tissue panel). We therefore crop each panel to its SCATTER
# FRAME only (frame x-bounds measured from the strip) and redraw crisp, non-overlapping legends
# with colours sampled from the source markers. Then two rows of 4 equal squares (the 4 TF KO/OE).
import os
from PIL import Image, ImageDraw, ImageFont

F = r"H:\毕业设计\网药部分\JingTong_CSR_paper\figures"
SC = os.path.join(F, "scRNA_dualvirtual")
OUT = os.path.join(F, "composite", "Fig7_singlecell.png")
try:
    FONT = ImageFont.truetype("C:/Windows/Fonts/arialbd.ttf", 48)
except Exception:
    FONT = ImageFont.load_default()
try:
    LFONT = ImageFont.truetype("C:/Windows/Fonts/arial.ttf", 20)
except Exception:
    LFONT = ImageFont.load_default()

W = 2200
GAP = 44
PAD = 52
LAB = 64
BG = (255, 255, 255)

# Scatter-frame x-bounds for the 3 UMAP panels in umap_clusters.png (measured from the plot
# borders); these exclude every source legend, the ghost marker column, and all label bleed.
PANEL_CROP = {
    "leiden":    (30, 524),
    "cell_type": (637, 1126),
    "tissue":    (1243, 1731),   # from the frame-left line (1243) so the tissue panel keeps its left border
}
# Redrawn legends: (label, RGB) with colours sampled from the source legend markers.
PANEL_LEGEND = {
    "leiden": [(str(i), c) for i, c in enumerate([
        (57, 135, 188), (255, 142, 42), (64, 169, 121), (218, 64, 65), (180, 86, 252),
        (140, 86, 75), (227, 119, 194), (181, 189, 97), (50, 197, 212)])],
    "cell_type": [("Chondrocyte", (57, 135, 188)), ("Endothelial", (255, 142, 42)),
                  ("Fibrochondrocyte", (44, 160, 44)), ("Immune/Macrophage", (218, 64, 65))],
    "tissue": [("AF", (57, 135, 188)), ("CEP", (255, 142, 42)), ("NP", (44, 160, 44))],
}


def draw_legend(entries, font, dr=7, rh=32, padr=16):
    """Render a clean single-column legend (coloured dot + label) on white."""
    meas = ImageDraw.Draw(Image.new("RGB", (4, 4)))
    tw = max(meas.textlength(lab, font=font) for lab, _ in entries)
    w = int(2 * dr + 16 + tw + padr)
    h = rh * len(entries)
    im = Image.new("RGB", (w, h), BG)
    d = ImageDraw.Draw(im)
    for i, (lab, col) in enumerate(entries):
        cy = i * rh + rh // 2
        d.ellipse((0, cy - dr, 2 * dr, cy + dr), fill=tuple(col))
        d.text((2 * dr + 16, cy), lab, font=font, fill=(35, 35, 35), anchor="lm")
    return im


def build_panelA(target_w, inner_gap=16, scat_h=440):
    """Panel A: 3 clean scatter frames + redrawn legends, laid out so the row spans the FULL
    width [0, target_w] (matching the B-I rows below) while keeping the three scatters EQUAL-sized
    and EVENLY spaced. The leftmost scatter starts at x=0 and the rightmost legend ends at
    target_w, so panel A's left/right edges align with the panels below and the group is centred."""
    strip = Image.open(os.path.join(SC, "umap_clusters.png")).convert("RGB")
    hs = strip.height
    keys = ("leiden", "cell_type", "tissue")
    crops = {k: strip.crop((PANEL_CROP[k][0], 0, PANEL_CROP[k][1], hs)) for k in keys}
    legs = {k: draw_legend(PANEL_LEGEND[k], LFONT) for k in keys}
    leg_t, leg_c = legs["tissue"].width, legs["cell_type"].width
    # even stride with leiden scatter at x=0 and tissue legend ending at target_w:
    #   stride = (W - slot_w - inner_gap - leg_tissue) / 2,  feasible when stride >= slot_w + inner_gap + leg_cell
    def geom(h):
        sws = {k: round(crops[k].width * h / hs) for k in keys}
        slot = max(sws.values())
        stride = (target_w - slot - inner_gap - leg_t) / 2.0
        return sws, slot, stride
    sws, slot_w, stride = geom(scat_h)
    while stride < slot_w + inner_gap + leg_c and scat_h > 220:   # shrink scatters until the middle legend clears
        scat_h -= 8
        sws, slot_w, stride = geom(scat_h)
    col_h = max(scat_h, max(l.height for l in legs.values()))
    banner = Image.new("RGB", (target_w, col_h), BG)
    for i, k in enumerate(keys):
        sc = crops[k].resize((sws[k], scat_h), Image.LANCZOS)
        slot_x = round(i * stride)
        banner.paste(sc, (slot_x + (slot_w - sws[k]) // 2, (col_h - scat_h) // 2))
        banner.paste(legs[k], (slot_x + slot_w + inner_gap, (col_h - legs[k].height) // 2))
    return banner


def rwid(name, w):
    im = Image.open(os.path.join(SC, name)).convert("RGB")
    return im.resize((w, max(1, round(im.height * w / im.width))), Image.LANCZOS)


rows_spec = [
    [("A", None)],   # banner (3 clean UMAP frames + redrawn legends)
    [("B", "MYC_KO.png"), ("C", "MYC_OE.png"), ("D", "JUN_KO.png"), ("E", "JUN_OE.png")],
    [("F", "ESR1_KO.png"), ("G", "ESR1_OE.png"), ("H", "HIF1A_KO.png"), ("I", "HIF1A_OE.png")],
]

grid = []
for row in rows_spec:
    n = len(row); pw = (W - (n - 1) * GAP) // n
    built = []
    for letter, fn in row:
        built.append((letter, build_panelA(W) if fn is None else rwid(fn, pw)))
    grid.append(built)

rowh = [max(im.height for _, im in r) for r in grid]
H = PAD * 2 + sum(rowh) + len(grid) * LAB + (len(grid) - 1) * GAP
canvas = Image.new("RGB", (PAD * 2 + W, H), BG)
d = ImageDraw.Draw(canvas)

y = PAD
for ri, row in enumerate(grid):
    x = PAD
    for letter, im in row:
        d.text((x, y), letter, font=FONT, fill=(0, 0, 0))
        canvas.paste(im, (x, y + LAB + (rowh[ri] - im.height) // 2))
        x += im.width + GAP
    y += LAB + rowh[ri] + GAP

canvas.save(OUT, dpi=(200, 200))
print("wrote Fig7_singlecell.png", canvas.size)
