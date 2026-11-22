"""Charts for the relTol article. Usage: python3 plot_rt.py datadir outdir [N]"""
import sys, os, csv
_here = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [_here, os.path.join(_here, "..", "common")]
from svgchart import *

src, out = sys.argv[1], sys.argv[2]
N = sys.argv[3] if len(sys.argv) > 3 else "64"
set_out(out)
s = list(csv.DictReader(open(os.path.join(src, f"simple_reltol_N{N}.csv"))))
cp = list(csv.DictReader(open(os.path.join(src, f"coupled_reltol_N{N}.csv"))))

parts = []
for p, (key, title, ylab, y0, y1, yt) in enumerate([
        ("outer", "SIMPLE: outer iterations", "outer iterations (log)", 100, 20000, [100, 1000, 10000]),
        ("seconds", "SIMPLE: wall time", "seconds (log)", 1, 200, [1, 10, 100])]):
    c = Chart(410, 420, title, 5e-7, 1.5, y0, y1, L=64, R=20, T=84, xlog=True, ylog=True, ox=p * 410)
    c.axes([1e-6, 1e-4, 1e-2, 1], yt, "pressure relTol (log)", ylab, fy=LOGF, fx=LOGF)
    for k, ru in enumerate(["0.5", "0.1", "0.01"]):
        d = [(float(r["relP"]), float(r[key])) for r in s if r["relU"] == ru and r["converged"] == "1"]
        c.line([x for x, _ in d], [y for _, y in d], SER[k]); c.dots([x for x, _ in d], [y for _, y in d], SER[k])
    parts.append(c)
lx = 74
for k, ru in enumerate(["0.5", "0.1", "0.01"]):
    t = f"momentum relTol {ru}"
    parts[0].label(lx, 52, t, SER[k]); lx += 30 + 7.4 * len(t)

parts[0].o.append(f'<text x="{parts[0].X(0.9)-8:.1f}" y="{parts[0].Y(1500):.1f}" fill="{MUTED}" text-anchor="end">0.9 with momentum 0.01: did not converge</text>')
save("simple_reltol.svg", 820, 420, parts)

parts = []
for p, (key, title, ylab, y0, y1, yt) in enumerate([
        ("outer", "Coupled: outer iterations", "outer iterations (log)", 5, 300, [10, 100]),
        ("seconds", "Coupled: wall time", "seconds (log)", 0.5, 20, [1, 10])]):
    c = Chart(410, 380, title, 5e-7, 1.5, y0, y1, L=64, R=20, T=60, xlog=True, ylog=True, ox=p * 410)
    c.axes([1e-6, 1e-4, 1e-2, 1], yt, "FGMRES relTol (log)", ylab, fx=LOGF, fy=LOGF)
    d = [(float(r["relLin"]), float(r[key])) for r in cp if r["converged"] == "1"]
    c.line([x for x, _ in d], [y for _, y in d], SER[2]); c.dots([x for x, _ in d], [y for _, y in d], SER[2])
    parts.append(c)
save("coupled_reltol.svg", 820, 380, parts)
print("ok")
