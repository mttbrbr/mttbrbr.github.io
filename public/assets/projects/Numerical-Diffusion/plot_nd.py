"""Dark-theme SVG/PNG figures for the numerical diffusion article (no dependencies
except ImageMagick for the 2D fields). Usage: python3 plot_nd.py datadir outdir"""
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

# ---------- 1. profiles ----------
rows = list(csv.DictReader(open(os.path.join(src, "profiles.csv"))))
x = [float(r["x"]) for r in rows]
c = Chart(820, 420, "After one revolution: N = 200 cells, CFL = 0.5", 0, 1, -0.3, 1.3, R=40)
c.axes([0, 0.2, 0.4, 0.6, 0.8, 1], [-0.25, 0, 0.25, 0.5, 0.75, 1, 1.25], "x", "u")
c.line(x, [float(r["exact"]) for r in rows], FG, dash="5 4", w=1.5)
names = [("upwind", "upwind (1st order)"), ("lax-wendroff", "Lax-Wendroff (2nd order)"), ("van-leer", "van Leer (TVD)")]
for k, (key, lab) in enumerate(names):
    c.line(x, [float(r[key]) for r in rows], SER[k])
lx = c.X(0) + 4
c.o.append(f'<line x1="{lx}" x2="{lx+18}" y1="{c.oy+40}" y2="{c.oy+40}" stroke="{FG}" stroke-dasharray="5 4"/>')
c.o.append(f'<text x="{lx+24}" y="{c.oy+44}" fill="{FG}">exact</text>')
lx += 80
for k, (key, lab) in enumerate(names):
    c.label(lx, c.oy + 40, lab, SER[k]); lx += 30 + 7.4 * len(lab)
c.T = 48
save("profiles.svg", 820, 420, [c])

# ---------- 2. convergence (two panels) ----------
conv = defaultdict(list)
for r in csv.DictReader(open(os.path.join(src, "convergence.csv"))):
    conv[(r["profile"], r["scheme"])].append((int(r["N"]), float(r["L1"])))
parts = []
for p, (prof, title) in enumerate([("gaussian", "Smooth profile (gaussian)"), ("square", "Discontinuity (square wave)")]):
    c = Chart(410, 380, title, 40, 2000, 1e-5, 0.3, L=64, R=20, xlog=True, ylog=True, ox=p * 410)
    c.axes([50, 100, 200, 400, 800, 1600], [1e-5, 1e-4, 1e-3, 1e-2, 1e-1], "cells N (log)", "L1 error (log)",
           fy=lambda v: f"1e{int(round(math.log10(v)))}")
    for k, (key, lab) in enumerate(names):
        d = conv[(prof, key)]
        c.line([a for a, _ in d], [b for _, b in d], SER[k]); c.dots([a for a, _ in d], [b for _, b in d], SER[k])
    parts.append(c)
lx = 74
for k, (key, lab) in enumerate(names):
    parts[0].o.append(""); parts[0].label(lx, 44, lab, SER[k]); lx += 30 + 7.4 * len(lab)
for c in parts: c.T = 60
save("convergence.svg", 820, 380, parts)

# ---------- 3. viscosity vs CFL ----------
vis = defaultdict(list)
for r in csv.DictReader(open(os.path.join(src, "viscosity.csv"))):
    vis[int(r["N"])].append((float(r["C"]), float(r["nu_predicted"]), float(r["nu_measured"])))
c = Chart(820, 400, "Upwind: viscosity measured from the spreading of a gaussian", 0, 1.05, 0, 5e-3, R=120)
c.axes([0, 0.2, 0.4, 0.6, 0.8, 1.0], [0, 1e-3, 2e-3, 3e-3, 4e-3, 5e-3], "CFL number C = a dt / dx",
       "numerical viscosity ν (a = 1, L = 1)", fy=lambda v: f"{v*1e3:g}e-3" if v else "0")
items = []
for k, n in enumerate(sorted(vis)):
    d = vis[n]
    xs = [0, 1]; c.line(xs, [ (1.0 / n) * (1 - xx) / 2 for xx in xs], SER[k], w=1.5)
    c.dots([a for a, _, _ in d], [m for _, _, m in d], SER[k])
    items.append((c.Y((1.0 / n) * (1 - 0) / 2), f"N = {n}", SER[k]))
for y, text, col in items: c.label(c.X(0) + 8, y - 12, text, col)
c.o.append(f'<text x="{c.X(0.55):.1f}" y="{c.Y(4.4e-3):.1f}" fill="{MUTED}">lines: theory dx (1 − C) / 2 · dots: measured</text>')
save("viscosity.svg", 820, 400, [c])

# ---------- 4. false diffusion width ----------
fd = defaultdict(list)
for r in csv.DictReader(open(os.path.join(src, "falsediff.csv"))):
    fd[int(r["N"])].append((float(r["s"]), float(r["width_measured"]), float(r["width_predicted"])))
c = Chart(820, 400, "45 degrees: width of the smeared step vs distance travelled", 0, 1.05, 0, 0.25, R=120)
c.axes([0, 0.2, 0.4, 0.6, 0.8, 1.0], [0, 0.05, 0.1, 0.15, 0.2, 0.25], "distance travelled from the corner s / L",
       "10–90% width / L")
items = []
for k, n in enumerate(sorted(fd)):
    d = fd[n]
    ss = [i / 100 for i in range(1, 101)]
    c.line(ss, [2.563 * math.sqrt(s * (1.0 / n) / math.sqrt(2)) for s in ss], SER[k], w=1.5)
    c.dots([a for a, _, _ in d], [m for _, m, _ in d], SER[k])
    items.append((c.Y(d[-1][2]), f"N = {n}", SER[k]))
legend_right(c, items)
c.o.append(f'<text x="{c.X(0.05):.1f}" y="{c.Y(0.235):.1f}" fill="{MUTED}">lines: modified equation · dots: measured</text>')
save("falsediff.svg", 820, 400, [c])

# ---------- 5. 2D fields -> PNG with a sequential blue ramp ----------
ramp = ["#184f95", "#1c5cab", "#256abf", "#2a78d6", "#3987e5", "#5598e7",
        "#6da7ec", "#86b6ef", "#9ec5f4", "#b7d3f6", "#cde2fb"]
rgb = [tuple(int(h[i:i+2], 16) for i in (1, 3, 5)) for h in ramp]
def color(v):
    v = min(max(v, 0.0), 1.0) * (len(rgb) - 1); i = min(int(v), len(rgb) - 2); t = v - i
    return bytes(int(rgb[i][k] + (rgb[i+1][k] - rgb[i][k]) * t) for k in range(3))
for th in ("00", "22", "45"):
    f = os.path.join(src, f"field2d_{th}.csv")
    if not os.path.exists(f): continue
    rows = [list(map(float, l.split(","))) for l in open(f)]
    n = len(rows)
    ppm = os.path.join(out, f"field_{th}.ppm")
    with open(ppm, "wb") as g:
        g.write(f"P6 {n} {n} 255\n".encode())
        for j in range(n - 1, -1, -1):
            g.write(b"".join(color(v) for v in rows[j]))
    subprocess.run(["magick", ppm, "-filter", "point", "-resize", "300x300", os.path.join(out, f"field_{th}.png")], check=True)
    os.remove(ppm)
print("ok")
