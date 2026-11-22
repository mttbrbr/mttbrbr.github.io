"""Charts for the expansion-ratio article. Usage: python3 plot_ex.py datadir outdir"""
import sys, os, csv
_here = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [_here, os.path.join(_here, "..", "common")]
from svgchart import *

src, out = sys.argv[1], sys.argv[2]
set_out(out)
MESH = [("uniform", "uniform", SER[0], None), ("jump 2:1", "jump 2:1", SER[1], None), ("graded r=1.2 to 2x", "graded 1.2 to 2x", SER[2], None),
        ("graded r=1.1 to 2x", "graded 1.1 to 2x", SER[3], None), ("jump 4:1", "jump 4:1", SER[4], "6 4"), ("graded r=1.2 to 4x", "graded 1.2 to 4x", SER[6], "6 4")]
d = list(csv.DictReader(open(os.path.join(src, "diffusion.csv"))))

# 1. diffusion: truncation error and global error vs fine cell size
parts = []
for q, (key, title, ylab, y0, y1, yt) in enumerate([
        ("max_truncation", "Local truncation error (max)", "truncation error (log)", 1e-4, 100, [1e-4, 1e-2, 1, 100]),
        ("l2_error", "Error of the solution (L2)", "error (log)", 1e-6, 0.05, [1e-6, 1e-4, 1e-2])]):
    c = Chart(410, 420, title, 1 / 1400, 1 / 30, y0, y1, L=64, R=20, T=110, xlog=True, ylog=True, ox=q * 410)
    c.axes([1 / 1280, 1 / 320, 1 / 80], yt, "fine cell size h (log)", ylab, fx=lambda v: f"1/{round(1/v)}", fy=LOGF)
    for name, lab, col, dash in MESH:
        p = [(float(r["h_fine"]), float(r[key])) for r in d if r["mesh"] == name]
        c.line([x for x, _ in p], [y for _, y in p], col, w=1.6, dash=dash); c.dots([x for x, _ in p], [y for _, y in p], col)
    parts.append(c)
lx, ly = 74, 50
for i, (name, lab, col, dash) in enumerate(MESH):
    if i == 3: lx, ly = 74, 72
    parts[0].label(lx, ly, lab, col); lx += 30 + 7.4 * len(lab)
save("diffusion.svg", 820, 420, parts)

# 2. waves: reflected energy (log bars), two resolutions
w = list(csv.DictReader(open(os.path.join(src, "waves.csv"))))
W, H, L, T = 820, 400, 70, 84
o = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" %s>' % (W, H, FONT), f'<rect width="{W}" height="{H}" fill="{BG}"/>',
     f'<text x="{L}" y="24" fill="{FG}" font-size="14">Fraction of the wave energy reflected by the transition (log)</text>']
import math
y0, y1 = -9, -1
Y = lambda v: T + (H - T - 60) * (1 - (math.log10(max(v, 1e-9)) - y0) / (y1 - y0))
for e in range(y0, y1 + 1, 2):
    o.append(f'<line x1="{L}" x2="{W-20}" y1="{Y(10**e):.1f}" y2="{Y(10**e):.1f}" stroke="{GRID}"/>')
    o.append(f'<text x="{L-8}" y="{Y(10**e)+4:.1f}" fill="{MUTED}" text-anchor="end">1e{e}</text>')
groups = [m for m in MESH if m[0] != "uniform"]
gw = (W - 20 - L) / len(groups)
for g, (name, lab, col, dash) in enumerate(groups):
    for k, ppw in enumerate(["16", "32"]):
        r = next(x for x in w if x["mesh"] == name and x["cells_per_wavelength_fine"] == ppw)
        v = float(r["reflected_fraction"]); destroyed = float(r["transmitted_fraction"]) < 0.05
        x = L + g * gw + gw / 2 - 34 + k * 36
        o.append(f'<rect x="{x:.1f}" y="{Y(v):.1f}" width="32" height="{Y(1e-9)-Y(v):.1f}" rx="3" fill="{SER[0] if k == 0 else SER[2]}"/>')
        if destroyed: o.append(f'<text x="{x+16:.1f}" y="{Y(v)-6:.1f}" fill="{FG}" text-anchor="middle">✕</text>')
    o.append(f'<text x="{L+g*gw+gw/2:.1f}" y="{Y(1e-9)+20:.1f}" fill="{MUTED}" text-anchor="middle">{lab}</text>')
o.append(f'<rect x="{L}" y="40" width="10" height="10" rx="2" fill="{SER[0]}"/><text x="{L+15}" y="49" fill="{FG}">16 fine cells per wavelength</text>')
o.append(f'<rect x="{L+250}" y="40" width="10" height="10" rx="2" fill="{SER[2]}"/><text x="{L+265}" y="49" fill="{FG}">32 fine cells per wavelength</text>')
o.append(f'<text x="{L+500}" y="49" fill="{MUTED}">✕ = wave destroyed</text>')
o.append("</svg>")
open(os.path.join(out, "waves.svg"), "w").write("\n".join(o))

# 3. profiles: reflected waves in front of the transition, 16 cells per wavelength
pr = list(csv.DictReader(open(os.path.join(src, "wave_profiles.csv"))))
c = Chart(820, 380, "Waves left behind the transition (16 fine cells per wavelength), zoomed", 1.55, 2.0, -0.06, 0.06, R=40, T=84)
c.axes([1.6, 1.7, 1.8, 1.9, 2.0], [-0.06, -0.03, 0, 0.03, 0.06], "x (the transition starts at x = 2)", "u (packet amplitude = 1)")
for name, lab, col in [("jump 2:1", "jump 2:1", SER[1]), ("graded r=1.1 to 2x", "graded 1.1 to 2x", SER[3])]:
    p = [(float(r["x"]), float(r["u"])) for r in pr if r["mesh"] == name and r["ppw"] == "16" and 1.55 <= float(r["x"]) <= 2.0]
    c.line([x for x, _ in p], [y for _, y in p], col, w=1.4)
c.label(c.X(1.55), 50, "jump 2:1", SER[1]); c.label(c.X(1.55) + 120, 50, "graded 1.1 to 2x", SER[3])
save("reflections.svg", 820, 380, [c])
print("ok")
