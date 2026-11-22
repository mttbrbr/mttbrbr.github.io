"""Dark-theme SVG/PNG figures for the compiler/precision article (no dependencies
except ImageMagick for the 2D fields). Usage: python3 plot_cp.py datadir outdir"""
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


import statistics
LOGF = lambda v: (f"{v:g}" if 1e-3 <= v < 1e5 else f"1e{int(round(math.log10(v)))}")
rows = list(csv.DictReader(open(os.path.join(src, "results.csv"))))
med = defaultdict(list)
for r in rows: med[(r["kind"], r["precision"], r["flags"], int(r["n"]), r["procs"])].append(float(r["ns_per_cell_iter"]))
M = {k: statistics.median(v) for k, v in med.items()}

# 1. compiler flags: two panels of horizontal bars
FL = [("gcc -O0", "gcc -O0"), ("gcc -O2 -march=x86-64 -mtune=generic", "gcc -O2, generic x86-64"), ("gcc -O3", "gcc -O3"),
      ("gcc -O3 -march=native", "gcc -O3 -march=native"), ("gcc -O3 -march=native -mprefer-vector-width=512", "  + 512-bit vectors"),
      ("gcc -O3 -march=native -ffast-math", "  + -ffast-math"), ("clang -O3", "clang -O3"),
      ("clang -O3 -march=native", "clang -O3 -march=native"), ("clang -O3 -march=native -ffast-math", "  + -ffast-math")]
W, rowh, T = 820, 30, 92
H = T + rowh * len(FL) + 50
o = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" %s>' % (W, H, FONT), f'<rect width="{W}" height="{H}" fill="{BG}"/>',
     f'<text x="20" y="24" fill="{FG}" font-size="14">Solver kernel (PCG, LDU matrix): time per cell per iteration, double precision</text>']
L0, pw = 250, 235
for p, (n, title) in enumerate([(24, "24³ cells: fits in cache"), (96, "96³ cells: in RAM")]):
    x0 = L0 + p * (pw + 25)
    vals = [M[("time", "DP", f, n, "1")] for f, _ in FL]
    xmax = 55
    o.append(f'<text x="{x0}" y="{T-30}" fill="{FG}">{title}</text>')
    for t in (0, 10, 20, 30, 40, 50):
        xx = x0 + t / xmax * pw
        o.append(f'<line x1="{xx:.1f}" x2="{xx:.1f}" y1="{T-14}" y2="{T+rowh*len(FL)}" stroke="{GRID}"/>')
        o.append(f'<text x="{xx:.1f}" y="{T+rowh*len(FL)+16}" fill="{MUTED}" text-anchor="middle">{t}</text>')
    for i, v in enumerate(vals):
        y = T + i * rowh
        if p == 0: o.append(f'<text x="{L0-12}" y="{y+18}" fill="{FG}" text-anchor="end">{FL[i][1]}</text>')
        o.append(f'<rect x="{x0}" y="{y+5}" width="{max(v/xmax*pw,1):.1f}" height="18" rx="2" fill="{SER[0]}"/>')
        o.append(f'<text x="{x0+v/xmax*pw+6:.1f}" y="{y+18}" fill="{FG}">{v:.1f}</text>')
o.append(f'<text x="{L0+pw+12}" y="{H-10}" fill="{MUTED}" text-anchor="middle">ns per cell per iteration (lower is better)</text>')
o.append("</svg>")
open(os.path.join(out, "flags.svg"), "w").write("\n".join(o))

# 2. precision vs mesh size: two panels (solver, fields)
Ns = [16, 24, 32, 48, 64, 96, 128, 160]
parts = []
for p, (kind, title) in enumerate([("time", "Solver kernel (PCG)"), ("fields", "Field operations")]):
    c = Chart(410, 400, title, 3000, 6e6, 0, 35, L=56, R=20, T=84, xlog=True, ox=p * 410)
    yt = [0, 5, 10, 15, 20, 25, 30, 35]
    c.axes([4096, 32768, 262144, 2097152], yt, "cells (log)", "ns per cell per iteration", fx=lambda v: {4096: "4k", 32768: "33k", 262144: "262k", 2097152: "2.1M"}[v])
    for k, prec in enumerate(["DP", "SP", "SPDP"]):
        fl = f"gcc -O3 -march=native -DPREC_{prec}"
        d = [(n ** 3, M[(kind, prec, fl, n, "1")]) for n in Ns if (kind, prec, fl, n, "1") in M]
        c.line([a for a, _ in d], [b for _, b in d], SER[k], dash="6 4" if prec == "SPDP" else None)
        c.dots([a for a, _ in d], [b for _, b in d], SER[k])
    if kind == "time":
        c.o.append(f'<text x="{c.X(128**3)-8:.1f}" y="{c.Y(33.5):.1f}" fill="{MUTED}" text-anchor="end">128³: power-of-two trap</text>')
    parts.append(c)
lx = 66
for k, lab in enumerate(["double (DP)", "single (SP)", "mixed (SPDP)"]):
    parts[0].label(lx, 52, lab, SER[k]); lx += 30 + 7.4 * len(lab)
save("precision_size.svg", 820, 400, parts)

# 3. speedup over DP: 1 vs 4 processes
groups = [("time", "1", "solver, 1 process"), ("time", "4", "solver, 4 processes"), ("fields", "1", "fields, 1 process"), ("fields", "4", "fields, 4 processes")]
W, H, L, T = 820, 360, 70, 80
o = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" %s>' % (W, H, FONT), f'<rect width="{W}" height="{H}" fill="{BG}"/>',
     f'<text x="{L}" y="24" fill="{FG}" font-size="14">Speedup over double precision, 96³ cells per process</text>']
o.append(f'<rect x="{L}" y="40" width="10" height="10" rx="2" fill="{SER[1]}"/><text x="{L+15}" y="49" fill="{FG}">single (SP)</text>')
o.append(f'<rect x="{L+130}" y="40" width="10" height="10" rx="2" fill="{SER[2]}"/><text x="{L+145}" y="49" fill="{FG}">mixed (SPDP)</text>')
ymax, ph = 2.2, H - T - 60
Y = lambda v: T + ph - v / ymax * ph
for t in (0, 0.5, 1, 1.5, 2):
    o.append(f'<line x1="{L}" x2="{W-20}" y1="{Y(t):.1f}" y2="{Y(t):.1f}" stroke="{GRID if t != 1 else MUTED}"' + (' stroke-dasharray="4 4"' if t == 1 else '') + '/>')
    o.append(f'<text x="{L-8}" y="{Y(t)+4:.1f}" fill="{MUTED}" text-anchor="end">{t:g}×</text>')
gw = (W - 20 - L) / len(groups)
def med_all(kind, prec, procs):
    fl = f"gcc -O3 -march=native -DPREC_{prec}"
    return statistics.median(med[(kind, prec, fl, 96, procs)])
for g, (kind, procs, lab) in enumerate(groups):
    dp = med_all(kind, "DP", procs)
    for k, prec in enumerate(["SP", "SPDP"]):
        sp = dp / med_all(kind, prec, procs)
        x = L + g * gw + gw / 2 - 46 + k * 48
        o.append(f'<rect x="{x:.1f}" y="{Y(sp):.1f}" width="44" height="{Y(0)-Y(sp):.1f}" rx="3" fill="{SER[1+k]}"/>')
        o.append(f'<text x="{x+22:.1f}" y="{Y(sp)-6:.1f}" fill="{FG}" text-anchor="middle">{sp:.2f}×</text>')
    o.append(f'<text x="{L+g*gw+gw/2:.1f}" y="{Y(0)+20:.1f}" fill="{MUTED}" text-anchor="middle">{lab}</text>')
o.append("</svg>")
open(os.path.join(out, "precision_speedup.svg"), "w").write("\n".join(o))

# 4. convergence at 96^3: true residual, plus the SP internal residual
c = Chart(820, 420, "PCG convergence, 96³ cells: true residual (computed in double)", 1, 900, 1e-13, 10, R=200, ylog=True)
c.axes([0, 200, 400, 600, 800], [1e-12, 1e-9, 1e-6, 1e-3, 1], "CG iteration", "relative residual (log)", fy=LOGF)
items = []
for k, prec in enumerate(["dp", "sp", "spdp"]):
    h = [l.strip().split(",") for l in open(os.path.join(src, f"solve_{prec}_96.csv")) if l.startswith("hist")]
    xs = [int(r[2]) for r in h]
    c.line(xs, [max(float(r[4]), 1e-13) for r in h], SER[k])
    items.append((c.Y(max(float(h[-1][4]), 1e-13)), {"dp": "DP", "sp": "SP", "spdp": "SPDP"}[prec] + " (true)", SER[k]))
    if prec == "sp":
        c.line(xs, [max(float(r[3]), 1e-13) for r in h], SER[k], dash="5 4", w=1.5)
        items.append((c.Y(max(float(h[-1][3]), 1e-13)), "SP (what CG sees)", SER[k]))
legend_right(c, items)
save("convergence.svg", 820, 420, [c])
print("ok")
