"""Dark-theme SVG/PNG figures for the low-Mach article (no dependencies
except ImageMagick for the 2D fields). Usage: python3 plot_lm.py datadir outdir"""
import csv, math, os, subprocess, sys
from collections import defaultdict

src, out = sys.argv[1], sys.argv[2]
FG, MUTED, GRID, BG = "#e6e6e6", "#8a8f94", "#23272b", "#08090a"
SER = ["#3987e5", "#d95926", "#199e70", "#c98500"]          # validated on #08090a
FONT = 'font-family="ui-monospace, SFMono-Regular, Menlo, monospace" font-size="12"'

class Chart:
    def __init__(s, W, H, title, x0, x1, y0, y1, L=70, R=150, T=48, B=56, xlog=False, ylog=False, ox=0, oy=0):
        s.W, s.H, s.L, s.R, s.T, s.B, s.ox, s.oy = W, H, L, R, T, B, ox, oy
        s.x0, s.x1, s.y0, s.y1, s.xlog, s.ylog = x0, x1, y0, y1, xlog, ylog
        s.o = [f'<text x="{ox+L}" y="{oy+24}" fill="{FG}" font-size="14">{title}</text>']
    def tx(s, v): return math.log10(v) if s.xlog else v
    def ty(s, v): return math.log10(v) if s.ylog else v
    def X(s, v): return s.ox + s.L + (s.tx(v) - s.tx(s.x0)) / (s.tx(s.x1) - s.tx(s.x0)) * (s.W - s.L - s.R)
    def Y(s, v): return s.oy + s.H - s.B - (s.ty(v) - s.ty(s.y0)) / (s.ty(s.y1) - s.ty(s.y0)) * (s.H - s.T - s.B)
    def axes(s, xt, yt, xl, yl, fx=lambda v: f"{v:g}", fy=lambda v: f"{v:g}"):
        for v in yt:
            s.o.append(f'<line x1="{s.X(s.x0):.1f}" x2="{s.X(s.x1):.1f}" y1="{s.Y(v):.1f}" y2="{s.Y(v):.1f}" stroke="{GRID}"/>')
            s.o.append(f'<text x="{s.X(s.x0)-8:.1f}" y="{s.Y(v)+4:.1f}" fill="{MUTED}" text-anchor="end">{fy(v)}</text>')
        for v in xt:
            s.o.append(f'<text x="{s.X(v):.1f}" y="{s.oy+s.H-s.B+18}" fill="{MUTED}" text-anchor="middle">{fx(v)}</text>')
        cy = s.oy + (s.T + s.H - s.B) / 2
        s.o.append(f'<text x="{s.ox+16}" y="{cy:.0f}" fill="{MUTED}" transform="rotate(-90 {s.ox+16} {cy:.0f})" text-anchor="middle">{yl}</text>')
        s.o.append(f'<text x="{(s.X(s.x0)+s.X(s.x1))/2:.0f}" y="{s.oy+s.H-12}" fill="{MUTED}" text-anchor="middle">{xl}</text>')
    def line(s, xs, ys, col, dash=None, w=2):
        pts = " ".join(f"{s.X(x):.1f},{s.Y(y):.1f}" for x, y in zip(xs, ys))
        d = f' stroke-dasharray="{dash}"' if dash else ""
        s.o.append(f'<polyline points="{pts}" fill="none" stroke="{col}" stroke-width="{w}" stroke-linejoin="round"{d}/>')
    def dots(s, xs, ys, col):
        for x, y in zip(xs, ys):
            s.o.append(f'<circle cx="{s.X(x):.1f}" cy="{s.Y(y):.1f}" r="4" fill="{col}" stroke="{BG}" stroke-width="2"/>')
    def label(s, x, y, text, col):
        s.o.append(f'<rect x="{x:.1f}" y="{y-5:.1f}" width="10" height="10" rx="2" fill="{col}"/>')
        s.o.append(f'<text x="{x+15:.1f}" y="{y+4:.1f}" fill="{FG}">{text}</text>')

def save(name, W, H, parts):
    body = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" {FONT}>', f'<rect width="{W}" height="{H}" fill="{BG}"/>']
    for p in parts: body += p.o
    body.append("</svg>")
    open(os.path.join(out, name), "w").write("\n".join(body))

def legend_right(c, items):
    """direct labels at line ends, de-collided"""
    items = sorted(items)
    for i in range(1, len(items)):
        if items[i][0] - items[i-1][0] < 16: items[i] = (items[i-1][0] + 16,) + items[i][1:]
    for y, text, col in items: c.label(c.X(c.x1) + 10, y, text, col)


import math
def rd(name): return list(csv.DictReader(open(os.path.join(src, name))))
LOGF = lambda v: (f"{v:g}" if 1e-3 <= v < 1e4 else f"1e{int(round(math.log10(v)))}")

# 1. steps vs Mach
r = rd("crossing.csv")
c = Chart(820, 400, "Time steps to carry the flow once across the domain", 0.0008, 0.6, 100, 1e6, R=170, xlog=True, ylog=True)
c.axes([0.001, 0.01, 0.1, 0.5], [100, 1e3, 1e4, 1e5, 1e6], "Mach number (log)", "time steps (log)", fy=LOGF)
items = []
for k, (sol, lab) in enumerate([("density", "density-based"), ("incompressible", "incompressible")]):
    d = [(float(x["M"]), float(x["steps"])) for x in r if x["solver"] == sol]
    c.line([a for a, _ in d], [b for _, b in d], SER[k]); c.dots([a for a, _ in d], [b for _, b in d], SER[k])
    items.append((c.Y(d[-1][1]), lab, SER[k]))
c.x1 = 0.6
for y, t, col in items: c.label(c.X(0.0008) + 8 if False else c.X(0.6) + 10, y, t, col)
save("steps.svg", 820, 400, [c])

# 2. profiles: three panels
pr = rd("profiles.csv")
parts = []
for p, M in enumerate(["0.5", "0.1", "0.01"]):
    d = [x for x in pr if x["M"] == M]
    c = Chart(275, 330, f"Mach {M}", 0, 1, 0.98, 1.11, L=48, R=12, ox=p * 273)
    c.axes([0, 0.5, 1], [1.0, 1.05, 1.1], "x", "density" if p == 0 else "")
    xs = [float(x["x"]) for x in d]
    c.line(xs, [float(x["initial"]) for x in d], FG, dash="5 4", w=1.5)
    c.line(xs, [float(x["incompressible"]) for x in d], SER[1])
    c.line(xs, [float(x["density_based"]) for x in d], SER[0])
    parts.append(c)
lx = 60
parts[0].o.append(f'<line x1="{lx}" x2="{lx+18}" y1="48" y2="48" stroke="{FG}" stroke-dasharray="5 4"/>')
parts[0].o.append(f'<text x="{lx+24}" y="52" fill="{FG}">exact</text>')
parts[0].label(lx + 90, 48, "density-based", SER[0]); parts[0].label(lx + 230, 48, "incompressible", SER[1])
for c in parts: c.T = 66
save("profiles.svg", 820, 330, parts)

# 3. viscosity ratio
v = rd("viscosity.csv")
c = Chart(820, 380, "How much more viscous: density-based (Rusanov) vs incompressible upwind", 0.0008, 0.6, 1, 2000, R=40, xlog=True, ylog=True)
c.axes([0.001, 0.01, 0.1, 0.5], [1, 10, 100, 1000], "Mach number (log)", "ratio of numerical viscosities (log)")
ms = [10 ** (e / 20) for e in range(-62, -5)]
c.line(ms, [1 + 1 / m for m in ms], MUTED, dash="5 4", w=1.5)
d = [(float(x["M"]), float(x["nu_measured"]) / float(x["nu_upwind_incomp"])) for x in v]
c.dots([a for a, _ in d], [b for _, b in d], SER[0])
c.o.append(f'<text x="{c.X(0.02):.1f}" y="{c.Y(300):.1f}" fill="{MUTED}">dashed: 1 + 1/M · dots: measured</text>')
save("viscosity_ratio.svg", 820, 380, [c])

# 4. round-off
ro = rd("roundoff.csv")
c = Chart(820, 420, "Error on the pressure fluctuation recovered from (ρ, ρu, E)", 2e-9, 0.2, 1e-14, 1e10, R=150, xlog=True, ylog=True)
c.axes([1e-8, 1e-6, 1e-4, 1e-2], [1e-14, 1e-10, 1e-6, 1e-2, 1e2, 1e6, 1e10], "Mach number (log)", "relative error, 1 = 100% (log)", fx=LOGF, fy=LOGF)
c.o.append(f'<line x1="{c.X(2e-9):.1f}" x2="{c.X(0.2):.1f}" y1="{c.Y(0.01):.1f}" y2="{c.Y(0.01):.1f}" stroke="{MUTED}" stroke-dasharray="4 4"/>')
c.o.append(f'<text x="{c.X(0.03):.1f}" y="{c.Y(0.01)-6:.1f}" fill="{MUTED}">1% error</text>')
items = []
for k, (col_, lab) in enumerate([("rel_err_float", "single (float)"), ("rel_err_double", "double")]):
    d = [(float(x["M"]), max(float(x[col_]), 1e-14)) for x in ro]
    c.line([a for a, _ in d], [b for _, b in d], SER[k]); c.dots([a for a, _ in d], [b for _, b in d], SER[k])
    items.append((c.Y(d[0][1]), lab, SER[k]))
legend_right(c, items)
save("roundoff.svg", 820, 420, [c])

# 5. eigenvalue ratio
e = rd("eigen.csv")
c = Chart(820, 380, "Fastest / slowest wave speed (stiffness)", 0.0008, 0.6, 1, 2000, R=170, xlog=True, ylog=True)
c.axes([0.001, 0.01, 0.1, 0.5], [1, 10, 100, 1000], "Mach number (log)", "ratio (log)")
items = []
for k, (col_, lab) in enumerate([("kappa", "no preconditioning"), ("kappa_precond", "preconditioned")]):
    d = [(float(x["M"]), float(x[col_])) for x in e]
    c.line([a for a, _ in d], [b for _, b in d], SER[k]); c.dots([a for a, _ in d], [b for _, b in d], SER[k])
    items.append((c.Y(d[0][1]), lab, SER[k]))
legend_right(c, items)
save("eigen.svg", 820, 380, [c])
print("ok")
