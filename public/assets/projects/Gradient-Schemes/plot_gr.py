"""Charts for the gradient-schemes article. Usage: python3 plot_gr.py datadir outdir"""
import sys, os, csv
_here = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [_here, os.path.join(_here, "..", "common")]
from svgchart import *
src, out = sys.argv[1], sys.argv[2]
set_out(out)
d = list(csv.DictReader(open(os.path.join(src, "gradients.csv"))))
M = [("uniform", SER[0]), ("skewed", SER[1]), ("random", SER[2]), ("stretched+curved", SER[3])]
parts = []
for q, (key, title) in enumerate([("gauss", "Gauss linear"), ("leastsquares", "leastSquares")]):
    c = Chart(410, 420, title, 12, 330, 3e-5, 1, L=64, R=20, T=96, xlog=True, ylog=True, ox=q * 410)
    c.axes([16, 32, 64, 128, 256], [1e-4, 1e-3, 1e-2, 0.1, 1], "cells per side (log)", "relative error (log)", fy=lambda v: f"{v*100:g}%")
    for name, col in M:
        p = [(int(r["N"]), float(r[key])) for r in d if r["mesh"] == name]
        c.line([x for x, _ in p], [y for _, y in p], col, w=1.6); c.dots([x for x, _ in p], [y for _, y in p], col)
    parts.append(c)
lx = 74
for name, col in M:
    parts[0].label(lx, 52, name, col); lx += 30 + 7.4 * len(name)
save("gradients.svg", 820, 420, parts)
# limiter effect, uniform and random meshes
c = Chart(820, 380, "leastSquares with and without cellLimited 1", 12, 330, 3e-5, 0.2, R=40, T=84, xlog=True, ylog=True)
c.axes([16, 32, 64, 128, 256], [1e-4, 1e-3, 1e-2, 0.1], "cells per side (log)", "relative error (log)", fy=lambda v: f"{v*100:g}%")
for name, col in [("uniform", SER[0]), ("random", SER[2])]:
    for key, dash in [("leastsquares", None), ("leastsquares_cellLimited", "6 4")]:
        p = [(int(r["N"]), float(r[key])) for r in d if r["mesh"] == name]
        c.line([x for x, _ in p], [y for _, y in p], col, w=1.6, dash=dash); c.dots([x for x, _ in p], [y for _, y in p], col)
c.label(c.X(12), 50, "uniform", SER[0]); c.label(c.X(12) + 110, 50, "random", SER[2])
c.o.append(f'<text x="{c.X(12)+220:.1f}" y="54" fill="{MUTED}">solid: no limiter · dashed: cellLimited 1</text>')
save("limiter.svg", 820, 380, [c])
print("ok")
