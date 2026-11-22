"""Charts for the non-orthogonality article. Usage: python3 plot_no.py datadir outdir"""
import sys, os, csv
_here = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [_here, os.path.join(_here, "..", "common")]
from svgchart import *

src, out = sys.argv[1], sys.argv[2]
set_out(out)
rd = lambda n: list(csv.DictReader(open(os.path.join(src, n))))

# 1. error vs angle
a = rd("angle.csv")
S5 = [("orthogonal", "orthogonal"), ("uncorrected", "uncorrected"), ("corrected_1", "corrected, 1 corrector"),
      ("corrected_2", "corrected, 2 correctors"), ("corrected_conv", "corrected, converged")]
c = Chart(820, 420, "Error vs non-orthogonality, 64 x 64 cells", 0, 85, 3e-4, 20, R=210, ylog=True)
c.axes([0, 10, 20, 30, 40, 50, 60, 70, 80], [1e-3, 1e-2, 0.1, 1, 10], "non-orthogonality angle (degrees)", "relative L2 error (log)",
       fy=lambda v: f"{v*100:g}%")
items = []
for k, (key, lab) in enumerate(S5):
    xs = [float(r["theta"]) for r in a]; ys = [float(r[key]) for r in a]
    c.line(xs, ys, SER[k]); c.dots(xs, ys, SER[k])
    items.append((c.Y(ys[-1]), lab, SER[k]))
legend_right(c, items)
save("error_angle.svg", 820, 420, [c])

# 2. error vs mesh size
g = rd("grid.csv")
c = Chart(820, 420, "Error vs mesh size", 12, 330, 3e-5, 5, R=210, xlog=True, ylog=True)
c.axes([16, 32, 64, 128, 256], [1e-4, 1e-3, 1e-2, 0.1, 1], "cells per side (log)", "relative L2 error (log)",
       fy=lambda v: f"{v*100:g}%")
items = []
for k, (th, key, lab, col) in enumerate([("60", "uncorrected", "60°, uncorrected", SER[1]), ("30", "uncorrected", "30°, uncorrected", SER[3]),
                                         ("60", "corrected_conv", "60°, corrected", SER[2]), ("0", "corrected_conv", "0° (orthogonal mesh)", SER[0])]):
    d = [(int(r["N"]), float(r[key])) for r in g if r["theta"] == th]
    c.line([x for x, _ in d], [y for _, y in d], col, dash="6 4" if th == "0" else None)
    c.dots([x for x, _ in d], [y for _, y in d], col)
    items.append((c.Y(d[-1][1]), lab, col))
legend_right(c, items)
save("error_grid.svg", 820, 420, [c])

# 3. error vs number of correctors
hst = rd("correctors.csv")
c = Chart(820, 400, "Error vs number of non-orthogonal correctors, 64 x 64 cells", 0, 70, 5e-4, 1, R=150, ylog=True)
c.axes([0, 10, 20, 30, 40, 50, 60, 70], [1e-3, 1e-2, 0.1, 1], "correctors", "relative L2 error (log)", fy=lambda v: f"{v*100:g}%")
items = []
for k, th in enumerate(["10", "30", "60", "80"]):
    d = [(int(r["corrector"]), float(r["error"])) for r in hst if r["theta"] == th and r["N"] == "64"]
    c.line([x for x, _ in d], [y for _, y in d], SER[k])
    items.append((c.Y(d[-1][1]), f"{th}°", SER[k]))
for y, tx, col in items: pass
lx = c.X(0) + 10
for k, th in enumerate(["10", "30", "60", "80"]):
    c.label(lx, 46, f"{th}°", SER[k]); lx += 70
save("error_correctors.svg", 820, 400, [c])
print("ok")
