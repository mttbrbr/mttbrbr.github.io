"""Charts for the skewness article. Usage: python3 plot_sk.py datadir outdir"""
import sys, os, csv
_here = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [_here, os.path.join(_here, "..", "common")]
from svgchart import *
src, out = sys.argv[1], sys.argv[2]
set_out(out)
d = list(csv.DictReader(open(os.path.join(src, "skewness.csv"))))
SK = ["0", "0.05", "0.1", "0.2", "0.4"]
parts = []
for q, (key, title) in enumerate([("grad_linear", "Gauss gradient, linear interpolation"), ("grad_corrected", "Gauss gradient, skewCorrected")]):
    c = Chart(410, 420, title, 12, 700, 1e-5, 1, L=64, R=20, T=96, xlog=True, ylog=True, ox=q * 410)
    c.axes([16, 32, 64, 128, 256, 512], [1e-5, 1e-4, 1e-3, 1e-2, 0.1, 1], "cells per side (log)", "relative error (log)",
           fy=lambda v: f"{v*100:g}%")
    for k, s in enumerate(SK):
        p = [(int(r["N"]), float(r[key])) for r in d if r["skew"] == s]
        c.line([x for x, _ in p], [y for _, y in p], SER[k], w=1.6); c.dots([x for x, _ in p], [y for _, y in p], SER[k])
    parts.append(c)
lx = 74
for k, s in enumerate(SK):
    t = f"skewness {s}"
    parts[0].label(lx, 52, t, SER[k]); lx += 30 + 7.4 * len(t)
save("skewness.svg", 820, 420, parts)
print("ok")
