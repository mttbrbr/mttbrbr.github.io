"""Charts for the multigrid article. Usage: python3 plot_mg.py datadir outdir"""
import sys, os, csv
_here = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [_here, os.path.join(_here, "..", "common")]   # svgchart.py next to this file or in ../common
from svgchart import *

src, out = sys.argv[1], sys.argv[2]
set_out(out)
rows = list(csv.DictReader(open(os.path.join(src, "scaling.csv"))))
M = [("gs", "Gauss-Seidel"), ("pcg", "PCG (Jacobi)"), ("mg", "multigrid"), ("mgpcg", "PCG + multigrid")]

for key_y, fname, title, ylab, y0, y1, yt in [
        ("iters", "iterations.svg", "Iterations to reduce the residual by 1e8", "iterations (log)", 1, 1e6, [1, 10, 100, 1000, 1e4, 1e5, 1e6]),
        ("seconds", "time.svg", "Time to solve, one core", "seconds (log)", 1e-4, 100, [1e-4, 1e-3, 1e-2, 0.1, 1, 10, 100])]:
    c = Chart(820, 420, title, 150, 1.6e6, y0, y1, R=170, xlog=True, ylog=True)
    c.axes([225, 3969, 65025, 1046529], yt, "unknowns (log)", ylab,
           fx=lambda v: {225: "225", 3969: "4k", 65025: "65k", 1046529: "1M"}[v], fy=LOGF)
    items = []
    for k, (m, lab) in enumerate(M):
        d = [(int(r["unknowns"]), max(float(r[key_y]), y0)) for r in rows if r["method"] == m]
        c.line([a for a, _ in d], [b for _, b in d], SER[k]); c.dots([a for a, _ in d], [b for _, b in d], SER[k])
        items.append((c.Y(d[-1][1]), lab, SER[k]))
    legend_right(c, items)
    save(fname, 820, 420, [c])

h = defaultdict(list)
for r in csv.DictReader(open(os.path.join(src, "history.csv"))): h[r["case"]].append((int(r["iter"]), float(r["residual"])))
c = Chart(820, 420, "Convergence on 127 x 127 unknowns", 1, 40000, 1e-10, 2, R=40, T=70, xlog=True, ylog=True)
c.axes([1, 10, 100, 1000, 10000], [1e-8, 1e-6, 1e-4, 1e-2, 1], "iteration (log)", "residual / initial (log)", fy=LOGF)
lx = c.X(1)
for k, (m, lab) in enumerate(M):
    d = h[m]
    c.line([a for a, _ in d], [max(b, 1e-10) for _, b in d], SER[k])
    t = f"{lab}: {d[-1][0]}"
    c.label(lx, 46, t, SER[k]); lx += 30 + 7.4 * len(t)
save("history.svg", 820, 420, [c])
print("ok")
