"""Charts for the checkerboard / Rhie-Chow article. Usage: python3 plot_cb.py datadir outdir"""
import sys, os, csv
_here = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [_here, os.path.join(_here, "..", "common")]
from svgchart import *

src, out = sys.argv[1], sys.argv[2]
set_out(out)
N, g = 32, 10.0
h = 1.0 / N
V = [("linear", "linear interpolation"), ("rhie-chow", "Rhie-Chow"), ("balanced", "Rhie-Chow + balanced force")]

def field(name, tag="a0_g10"):
    return [list(map(float, l.split(","))) for l in open(os.path.join(src, f"p_{name}_N{N}_{tag}.csv"))]

ys = [(j + 0.5) * h for j in range(N)]
pex = [g * y if y < 0.5 else g * (1 - y) for y in ys]
m = sum(pex) / N
pex = [p - m for p in pex]

# 1. pressure profile on the vertical centreline + deviation from the exact one
parts = []
c = Chart(410, 420, "Pressure on the vertical centreline", 0, 1, -3, 3, L=56, R=20, T=84)
c.axes([0, 0.25, 0.5, 0.75, 1], [-3, -2, -1, 0, 1, 2, 3], "height y", "pressure")
c.line(ys, pex, FG, dash="5 4", w=1.5)
for k, (name, lab) in enumerate(V):
    P = field(name); col = [0.5 * (P[j][N // 2 - 1] + P[j][N // 2]) for j in range(N)]
    c.line(ys, col, SER[k], w=1.5)
parts.append(c)
c = Chart(410, 420, "Deviation from the exact pressure", 0, 1, -0.2, 0.2, L=64, R=20, T=84, ox=410)
c.axes([0, 0.25, 0.5, 0.75, 1], [-0.2, -0.1, 0, 0.1, 0.2], "height y", "p - p_exact")
for k, (name, lab) in enumerate(V):
    P = field(name); col = [0.5 * (P[j][N // 2 - 1] + P[j][N // 2]) - pex[j] for j in range(N)]
    c.line(ys, col, SER[k], w=1.5); c.dots(ys, col, SER[k])
parts.append(c)
lx = 66
parts[0].o.append(f'<line x1="{lx}" x2="{lx+18}" y1="52" y2="52" stroke="{FG}" stroke-dasharray="5 4"/>')
parts[0].o.append(f'<text x="{lx+24}" y="56" fill="{FG}">exact</text>'); lx += 90
for k, (name, lab) in enumerate(V):
    parts[0].label(lx, 52, lab, SER[k]); lx += 30 + 7.4 * len(lab)
save("pressure_profile.svg", 820, 420, parts)

# 2. pressure maps (hydrostatic test): linear vs Rhie-Chow vs balanced, deviation from exact
for name, _ in V:
    P = field(name)
    vals = [P[j][i] - pex[j] for j in range(N) for i in range(N)]
    field_png(vals, N, N, os.path.join(out, f"perr_{name}.png"), scale=8, vmin=-0.16, vmax=0.16, diverging=True)
print("ok")
