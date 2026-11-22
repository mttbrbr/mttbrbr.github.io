"""Shared dark-theme SVG chart helpers for the articles (no dependencies).
Use: sys.path.insert(0, "../common"); from svgchart import *; set_out(outdir)"""
import csv, math, os, subprocess, sys
from collections import defaultdict

FG, MUTED, GRID, BG = "#e6e6e6", "#8a8f94", "#23272b", "#08090a"
SER = ["#3987e5", "#d95926", "#199e70", "#c98500", "#d55181", "#0ca30c", "#9085e9", "#e66767"]   # reference categorical order, dark mode
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

OUT = "."
def set_out(d):
    global OUT
    OUT = d

def save(name, W, H, parts):
    body = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" {FONT}>', f'<rect width="{W}" height="{H}" fill="{BG}"/>']
    for p in parts: body += p.o
    body.append("</svg>")
    open(os.path.join(OUT, name), "w").write("\n".join(body))

def legend_right(c, items):
    """direct labels at line ends, de-collided"""
    items = sorted(items)
    for i in range(1, len(items)):
        if items[i][0] - items[i-1][0] < 16: items[i] = (items[i-1][0] + 16,) + items[i][1:]
    for y, text, col in items: c.label(c.X(c.x1) + 10, y, text, col)


LOGF = lambda v: (f"{v:g}" if 1e-3 <= v < 1e5 else f"1e{int(round(math.log10(v)))}")
RAMP = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7", "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95"]
DIV_NEG, DIV_POS = "#3987e5", "#d95926"

def field_png(values, nx, ny, path, scale=4, vmin=None, vmax=None, diverging=False):
    """Write a 2D field (row-major, row 0 at the bottom) as a PNG via ImageMagick."""
    import subprocess
    vmin = min(values) if vmin is None else vmin
    vmax = max(values) if vmax is None else vmax
    def lerp(a, b, t): return tuple(int(a[k] + (b[k] - a[k]) * t) for k in range(3))
    hexrgb = lambda h: tuple(int(h[i:i+2], 16) for i in (1, 3, 5))
    ramp = [hexrgb(h) for h in RAMP]
    neg, pos, mid = hexrgb(DIV_NEG), hexrgb(DIV_POS), hexrgb("#e6e6e6")
    rows = []
    for j in range(ny - 1, -1, -1):
        row = bytearray()
        for i in range(nx):
            v = values[j * nx + i]
            if diverging:
                m = max(abs(vmin), abs(vmax)) or 1
                t = max(-1, min(1, v / m))
                c = lerp(mid, pos, t) if t >= 0 else lerp(mid, neg, -t)
            else:
                t = 0 if vmax == vmin else max(0, min(1, (v - vmin) / (vmax - vmin)))
                x = t * (len(ramp) - 1); k = min(int(x), len(ramp) - 2)
                c = lerp(ramp[k], ramp[k + 1], x - k)
            row += bytes(c)
        rows.append(bytes(row))
    ppm = path + ".ppm"
    with open(ppm, "wb") as f:
        f.write(f"P6 {nx} {ny} 255\n".encode()); f.write(b"".join(rows))
    subprocess.run(["magick", ppm, "-filter", "point", "-resize", f"{nx*scale}x{ny*scale}", path], check=True)
    os.remove(ppm)
