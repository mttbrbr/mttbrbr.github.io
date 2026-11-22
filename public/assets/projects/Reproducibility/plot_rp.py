"""Charts for the reproducibility article. Usage: python3 plot_rp.py datadir outdir"""
import sys, os, csv
_here = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [_here, os.path.join(_here, "..", "common")]
from svgchart import *

src, out = sys.argv[1], sys.argv[2]
set_out(out)
rd = lambda n: list(csv.DictReader(open(os.path.join(src, n))))
tur = rd("divergence_nu0.002_F1.csv")
lam = rd("divergence_nu0.05_F0.2.csv")

# 1. difference from the 1-process run vs time
c = Chart(820, 440, "Difference from the 1-process run (relative L2 norm of the vorticity)", 0, 80, 1e-17, 3, R=40, T=92, ylog=True)
c.axes([0, 10, 20, 30, 40, 50, 60, 70, 80], [1e-16, 1e-12, 1e-8, 1e-4, 1], "time", "relative difference (log)", fy=LOGF)
curves = [("diff_P2", "2 processes", SER[0]), ("diff_P4", "4 processes", SER[1]), ("diff_P8", "8 processes", SER[2]),
          ("diff_ulp", "1 process, one value changed by 1 ulp", SER[4])]
for key, lab, col in curves:
    d = [(float(r["t"]), float(r[key])) for r in tur if float(r["t"]) > 0 and float(r[key]) > 0]
    c.line([x for x, _ in d], [max(y, 1e-17) for _, y in d], col, w=1.8)
d = [(float(r["t"]), float(r["diff_P8"])) for r in lam if 0 < float(r["t"]) <= 60 and float(r["diff_P8"]) > 0]
c.line([x for x, _ in d], [max(y, 1e-17) for _, y in d], MUTED, dash="5 4", w=1.5)
lx = c.X(0)
for i, (key, lab, col) in enumerate(curves):
    if i == 3: lx = c.X(0)
    c.label(lx, 48 if i < 3 else 68, lab, col); lx += 30 + 7.4 * len(lab)
c.o.append(f'<line x1="{lx+10}" x2="{lx+28}" y1="68" y2="68" stroke="{MUTED}" stroke-dasharray="5 4"/><text x="{lx+34}" y="72" fill="{FG}">laminar case, 8 processes</text>')
save("divergence.svg", 820, 440, [c])

# 2. energy vs time, 1 vs 8 processes
c = Chart(820, 380, "Kinetic energy of the turbulent run", 0, 80, 0, 3.5, R=40, T=70)
c.axes([0, 10, 20, 30, 40, 50, 60, 70, 80], [0, 1, 2, 3], "time", "kinetic energy")
for key, lab, col in [("E_P1", "1 process", SER[0]), ("E_P8", "8 processes", SER[2])]:
    d = [(float(r["t"]), float(r[key])) for r in tur]
    c.line([x for x, _ in d], [y for _, y in d], col, w=1.6)
c.label(c.X(0), 46, "1 process", SER[0]); c.label(c.X(0) + 130, 46, "8 processes", SER[2])
save("energy.svg", 820, 380, [c])
# 3. final vorticity fields
for run, lab in [(0, "P1"), (3, "P8"), (4, "ulp")]:
    rows = [list(map(float, l.split(","))) for l in open(os.path.join(src, f"vorticity_nu0.002_F1_run{run}.csv"))]
    vals = [x for r in rows for x in r]
    field_png(vals, len(rows[0]), len(rows), os.path.join(out, f"vorticity_{lab}.png"), scale=4, vmin=-8, vmax=8, diverging=True)
print("ok")
