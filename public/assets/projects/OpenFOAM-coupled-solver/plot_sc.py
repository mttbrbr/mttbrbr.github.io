"""Dark-theme SVG/PNG figures for the segregated vs coupled article (no dependencies
except ImageMagick for the 2D fields). Usage: python3 plot_sc.py datadir outdir"""
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


def rd(name): return list(csv.DictReader(open(os.path.join(src, name))))
LOGF = lambda v: (f"{v:g}" if 1e-2 <= v < 1e5 else f"1e{int(round(math.log10(v)))}")
it = rd("iterative.csv")
M3 = [("simple_iter_07_03", "SIMPLE 0.7 / 0.3"), ("simple_iter_09_02", "SIMPLE 0.9 / 0.2 (best)"), ("coupled_iter", "coupled")]

# 1. outer iterations vs grid
c = Chart(820, 400, "Outer (non-linear) iterations to converge, Re = 100", 25, 160, 5, 20000, R=200, xlog=True, ylog=True)
c.axes([32, 64, 128], [10, 100, 1000, 10000], "cells per side (log)", "outer iterations (log)", fy=LOGF)
items = []
for k, (key, lab) in enumerate(M3):
    d = [(int(x["N"]), int(x["iters"])) for x in it if x["method"] == key]
    c.line([a for a, _ in d], [b for _, b in d], SER[k]); c.dots([a for a, _ in d], [b for _, b in d], SER[k])
    items.append((c.Y(d[-1][1]), lab, SER[k]))
legend_right(c, items)
save("outer_iterations.svg", 820, 400, [c])

# 2. time and memory vs grid (two panels)
parts = []
for p, (col_, title, ylab, y0, y1, ticks) in enumerate([
        ("seconds", "Wall time (one core)", "seconds (log)", 0.05, 1000, [0.1, 1, 10, 100, 1000]),
        ("MB", "Memory of the linear solvers", "MB (log)", 0.1, 100, [0.1, 1, 10, 100])]):
    c = Chart(410, 410, title, 25, 160, y0, y1, L=64, R=20, T=84, xlog=True, ylog=True, ox=p * 410)
    c.axes([32, 64, 128], ticks, "cells per side (log)", ylab, fy=LOGF)
    for k, (key, lab) in enumerate(M3):
        if col_ == "MB" and k == 0: continue          # same memory as the other SIMPLE
        d = [(int(x["N"]), float(x[col_])) for x in it if x["method"] == key]
        c.line([a for a, _ in d], [b for _, b in d], SER[k]); c.dots([a for a, _ in d], [b for _, b in d], SER[k])
    parts.append(c)
lx = 74
for k, (key, lab) in enumerate(M3):
    parts[0].label(lx, 52, lab, SER[k]); lx += 30 + 7.4 * len(lab)
save("time_memory.svg", 820, 410, parts)

# 3. convergence history N = 64 (iterative)
hist = defaultdict(list)
for x in rd("history_iter.csv" if os.path.exists(os.path.join(src, "history_iter.csv")) else "history.csv"): hist[x["case"]].append((int(x["iter"]), float(x["residual"])))
c = Chart(820, 400, "Convergence history, 64 x 64 cells, Re = 100", 1, 3000, 1e-9, 2, R=40, xlog=True, ylog=True)
c.axes([1, 10, 100, 1000], [1e-8, 1e-6, 1e-4, 1e-2, 1], "outer iteration (log)", "residual / initial residual (log)", fy=LOGF)
for k, (key, lab) in enumerate([("simple_iter_N64", "SIMPLE 0.7 / 0.3"), ("coupled_iter_N64", "coupled")]):
    d = hist[key]; col = SER[0 if k == 0 else 2]
    c.line([a for a, _ in d], [max(b, 1e-9) for _, b in d], col)
    ex, ey = c.X(d[-1][0]), c.Y(max(d[-1][1], 1e-9))
    if k == 0: c.label(ex - 250, ey + 16, f"{lab}: {d[-1][0]} iterations", col)
    else:      c.label(ex + 10, ey, f"{lab}: {d[-1][0]} iterations", col)
save("history.svg", 820, 400, [c])

# 4. relaxation map (heatmap, sequential blue: darker = more iterations)
rl = rd("relaxation.csv")
aus = sorted(set(float(x["alpha_u"]) for x in rl)); aps = sorted(set(float(x["alpha_p"]) for x in rl))
W, H, L, T, cw, ch = 820, 360, 120, 70, 110, 48
o = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" %s>' % (W, H, FONT), f'<rect width="{W}" height="{H}" fill="{BG}"/>',
     f'<text x="{L}" y="24" fill="{FG}" font-size="14">SIMPLE: outer iterations vs under-relaxation (32 x 32, Re = 100)</text>',
     f'<text x="{L}" y="46" fill="{MUTED}">darker = more iterations · empty = diverged · the coupled solver has no such parameter</text>']
ramp = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95"]
conv = [int(x["iters"]) for x in rl if x["converged"] == "1"]
lo, hi = math.log(min(conv)), math.log(max(conv))
for a, au in enumerate(aus):
    y = T + a * ch
    o.append(f'<text x="{L-12}" y="{y+ch/2+4:.0f}" fill="{MUTED}" text-anchor="end">α_u = {au:g}</text>')
    for b, ap in enumerate(aps):
        x = L + b * cw
        r = next(z for z in rl if float(z["alpha_u"]) == au and float(z["alpha_p"]) == ap)
        if r["converged"] == "1":
            n = int(r["iters"]); t = (math.log(n) - lo) / (hi - lo); col = ramp[min(int(t * (len(ramp) - 1) + 0.5), len(ramp) - 1)]
            ink = BG if t < 0.5 else FG
            o.append(f'<rect x="{x+1}" y="{y+1}" width="{cw-2}" height="{ch-2}" rx="3" fill="{col}"/>')
            o.append(f'<text x="{x+cw/2:.0f}" y="{y+ch/2+4:.0f}" fill="{ink}" text-anchor="middle">{n}</text>')
        else:
            o.append(f'<rect x="{x+1}" y="{y+1}" width="{cw-2}" height="{ch-2}" rx="3" fill="none" stroke="{GRID}"/>')
            o.append(f'<text x="{x+cw/2:.0f}" y="{y+ch/2+4:.0f}" fill="{MUTED}" text-anchor="middle">diverged</text>')
for b, ap in enumerate(aps):
    o.append(f'<text x="{L+b*cw+cw/2:.0f}" y="{T+len(aus)*ch+20}" fill="{MUTED}" text-anchor="middle">α_p = {ap:g}</text>')
o.append("</svg>")
open(os.path.join(out, "relaxation_map.svg"), "w").write("\n".join(o))
print("ok")
