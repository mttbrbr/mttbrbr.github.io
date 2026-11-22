"""Charts for the Taylor-Green article. Usage: python3 plot_tgv.py datadir outdir"""
import sys, os, csv
_here = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [_here, os.path.join(_here, "..", "common")]
from svgchart import *

src, out = sys.argv[1], sys.argv[2]
set_out(out)
SCH = [("central", "central", SER[0]), ("linearUpwind", "linearUpwind", SER[2]), ("upwind", "upwind", SER[1])]

def load(N, s):
    p = os.path.join(src, f"tgv_N{N}_{s}_Re1600.csv")
    if not os.path.exists(p): return None
    r = [(float(x["t"]), float(x["energy"]), float(x["physical_dissipation"])) for x in csv.DictReader(open(p))]
    d = [((r[i + 1][0] + r[i - 1][0]) / 2, -(r[i + 1][1] - r[i - 1][1]) / (r[i + 1][0] - r[i - 1][0]), r[i][2]) for i in range(1, len(r) - 1)]
    return r, d

# 1. total dissipation vs time for each scheme, 64^3 and 128^3 (two panels)
parts = []
for q, N in enumerate(["64", "128"]):
    c = Chart(410, 420, f"{N}³: total dissipation −dE/dt", 0, 20, 0, 0.016, L=64, R=20, T=96, ox=q * 410)
    c.axes([0, 5, 10, 15, 20], [0, 0.004, 0.008, 0.012, 0.016], "time", "−dE/dt")
    c.o.append(f'<circle cx="{c.X(9):.1f}" cy="{c.Y(0.0127):.1f}" r="5" fill="none" stroke="{FG}" stroke-width="1.5"/>')
    c.o.append(f'<text x="{c.X(9)+8:.1f}" y="{c.Y(0.0127)-8:.1f}" fill="{FG}">reference DNS</text>')
    for s, lab, col in SCH:
        L = load(N, s)
        if not L: continue
        c.line([x[0] for x in L[1]], [x[1] for x in L[1]], col, w=1.8)
    parts.append(c)
lx = 74
for s, lab, col in SCH:
    parts[0].label(lx, 52, lab, col); lx += 30 + 7.4 * len(lab)
save("dissipation_total.svg", 820, 420, parts)

# 2. physical vs numerical, 64^3, one panel per scheme
parts = []
for q, (s, lab, col) in enumerate(SCH):
    c = Chart(273, 360, lab + ", 64³", 0, 20, 0, 0.016, L=56, R=12, T=84, ox=q * 273)
    c.axes([0, 10, 20], [0, 0.005, 0.01, 0.015], "time", "dissipation" if q == 0 else "")
    L = load("64", s)
    t = [x[0] for x in L[1]]; tot = [x[1] for x in L[1]]; ph = [x[2] for x in L[1]]
    c.line(t, tot, MUTED, w=1.5, dash="5 4")
    c.line(t, ph, col, w=2)
    parts.append(c)
parts[0].o.append(f'<line x1="66" x2="84" y1="50" y2="50" stroke="{MUTED}" stroke-dasharray="5 4"/><text x="90" y="54" fill="{FG}">total −dE/dt</text>')
parts[0].o.append(f'<line x1="210" x2="228" y1="50" y2="50" stroke="{FG}" stroke-width="2"/><text x="234" y="54" fill="{FG}">physical ν⟨|∇u|²⟩ (colour of the scheme); the gap is numerical</text>')
save("physical_vs_numerical.svg", 820, 360, parts)

# 3. energy vs time, all schemes and grids
c = Chart(820, 400, "Kinetic energy vs time", 0, 20, 0, 0.13, R=40, T=96)
c.axes([0, 5, 10, 15, 20], [0, 0.04, 0.08, 0.12], "time", "kinetic energy")
for s, lab, col in SCH:
    for N, dash in [("32", "2 3"), ("64", "6 4"), ("128", None)]:
        L = load(N, s)
        if not L: continue
        c.line([x[0] for x in L[0]], [x[1] for x in L[0]], col, w=1.6, dash=dash)
lx = c.X(0)
for s, lab, col in SCH:
    c.label(lx, 50, lab, col); lx += 30 + 7.4 * len(lab)
c.o.append(f'<text x="{lx+10:.1f}" y="54" fill="{MUTED}">dotted 32³ · dashed 64³ · solid 128³</text>')
save("energy.svg", 820, 400, [c])
print("ok")
