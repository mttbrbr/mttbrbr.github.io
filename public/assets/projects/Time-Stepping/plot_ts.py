"""Charts for the time-stepping article. Usage: python3 plot_ts.py datadir outdir"""
import sys, os, csv
_here = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [_here, os.path.join(_here, "..", "common")]
from svgchart import *

src, out = sys.argv[1], sys.argv[2]
set_out(out)
rd = lambda n: list(csv.DictReader(open(os.path.join(src, n))))
SCH = [("rk4", "RK4 (explicit)"), ("euler-implicit", "Euler implicit"), ("bdf2", "BDF2 (backward)"),
       ("crank-nicolson", "Crank-Nicolson"), ("theta-0.55", "theta = 0.55")]
COL = {k: SER[i] for i, (k, _) in enumerate(SCH)}

# 1. advection: peak and phase vs Courant
a = rd("advection.csv")
parts = []
for p, (key, title, ylab, y0, y1, yt, fy) in enumerate([
        ("peak", "Peak height after one period", "peak (exact = 1)", 0, 1.05, [0, 0.25, 0.5, 0.75, 1], lambda v: f"{v:g}"),
        ("phase_shift", "Phase lag after one period", "lag (fraction of domain)", 0, 0.085, [0, 0.02, 0.04, 0.06, 0.08], lambda v: f"{v:g}")]):
    c = Chart(410, 420, title, 0.08, 12, y0, y1, L=64, R=20, T=96, xlog=True, ox=p * 410)
    c.axes([0.1, 1, 10], yt, "Courant number (log)", ylab, fy=fy)
    for k, lab in SCH:
        d = [(float(r["courant"]), float(r[key])) for r in a if r["scheme"] == k and r[key] not in ("nan", "-nan")]
        if key == "phase_shift": d = [(x, -y) for x, y in d]
        if not d: continue
        c.line([x for x, _ in d], [y for _, y in d], COL[k]); c.dots([x for x, _ in d], [y for _, y in d], COL[k])
    parts.append(c)
lx = 74
for i, (k, lab) in enumerate(SCH):
    y = 50 if i < 3 else 70
    if i == 3: lx = 74
    parts[0].label(lx, y, lab, COL[k]); lx += 30 + 7.4 * len(lab)
save("advection.svg", 820, 420, parts)

# 2. advection profiles at Courant 2 (RK4 unstable there)
pr = rd("advection_profiles.csv")
c = Chart(820, 380, "Gaussian after one period, Courant number 2", 0.3, 0.7, -0.05, 1.05, R=40, T=70)
c.axes([0.3, 0.4, 0.5, 0.6, 0.7], [0, 0.25, 0.5, 0.75, 1], "x", "u")
import math
xs = [0.3 + i * 0.4 / 200 for i in range(201)]
c.line(xs, [math.exp(-0.5 * ((x - 0.5) / 0.04) ** 2) for x in xs], FG, dash="5 4", w=1.5)
lx = c.X(0.3)
c.o.append(f'<line x1="{lx}" x2="{lx+18}" y1="46" y2="46" stroke="{FG}" stroke-dasharray="5 4"/><text x="{lx+24}" y="50" fill="{FG}">exact</text>')
lx += 90
for k, lab in SCH[1:]:
    d = [(float(r["x"]), float(r["u"])) for r in pr if r["scheme"] == k and r["courant"] == "2" and 0.3 <= float(r["x"]) <= 0.7]
    c.line([x for x, _ in d], [y for _, y in d], COL[k])
    c.label(lx, 46, lab, COL[k]); lx += 30 + 7.4 * len(lab)
save("advection_profiles.svg", 820, 380, [c])

# 3. diffusion of a step: zoom on the left edge after 1 and 20 steps
dp = rd("diffusion_profiles.csv")
parts = []
for p, st in enumerate(["1", "20"]):
    c = Chart(410, 400, f"After {st} step{'s' if st != '1' else ''}", 0.2, 0.3, -0.1, 1.1, L=56, R=20, T=84, ox=p * 410)
    c.axes([0.2, 0.25, 0.3], [0, 0.25, 0.5, 0.75, 1], "x", "u")
    d = [(float(r["x"]), float(r["u"])) for r in dp if r["steps"] == st and r["scheme"] == "exact" and 0.2 <= float(r["x"]) <= 0.3]
    c.line([x for x, _ in d], [y for _, y in d], FG, dash="5 4", w=1.5)
    for k, lab in SCH[1:]:
        d = [(float(r["x"]), float(r["u"])) for r in dp if r["steps"] == st and r["scheme"] == k and 0.2 <= float(r["x"]) <= 0.3]
        c.line([x for x, _ in d], [y for _, y in d], COL[k], w=1.8); c.dots([x for x, _ in d], [y for _, y in d], COL[k])
    parts.append(c)
lx = 66
parts[0].o.append(f'<line x1="{lx}" x2="{lx+18}" y1="50" y2="50" stroke="{FG}" stroke-dasharray="5 4"/><text x="{lx+24}" y="54" fill="{FG}">exact</text>')
lx += 90
for k, lab in SCH[1:]:
    parts[0].label(lx, 50, lab, COL[k]); lx += 30 + 7.4 * len(lab)
save("diffusion.svg", 820, 400, parts)
print("ok")
