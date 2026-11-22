"""Charts for the Amdahl article. Usage: python3 plot_am.py datadir outdir"""
import sys, os, csv, glob, re, statistics
_here = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [_here, os.path.join(_here, "..", "common")]
from svgchart import *

src, out = sys.argv[1], sys.argv[2]
set_out(out)
rows = list(csv.DictReader(open(os.path.join(src, "amdahl.csv"))))
def T(kind, s, P): return next(float(r["seconds"]) for r in rows if r["kind"] == kind and float(r["serial_fraction"]) == s and int(r["threads"]) == P)
SF = [0.0, 0.05, 0.2, 0.5]
T1c = min(T("compute", s, 1) for s in SF)      # same total work for every s: use the best single-thread time
T1m = min(T("memory", s, 1) for s in SF)

# 1. compute-bound: measured vs Amdahl
c = Chart(820, 440, "Compute-bound work: measured speedup (dots) vs Amdahl's law (dashed)", 1, 8, 0, 8.5, R=150, T=84)
c.axes([1, 2, 3, 4, 5, 6, 7, 8], [0, 2, 4, 6, 8], "threads (4 physical cores, 5-8 = hyperthreads)", "speedup")
c.o.append(f'<rect x="{c.X(4):.1f}" y="{c.Y(8.5):.1f}" width="{c.X(8)-c.X(4):.1f}" height="{c.Y(0)-c.Y(8.5):.1f}" fill="{FG}" opacity="0.04"/>')
items = []
for k, s in enumerate(SF):
    Ps = list(range(1, 9))
    c.line(Ps, [1 / (s + (1 - s) / P) for P in Ps], SER[k], dash="5 4", w=1.3)
    sp = [T1c / T("compute", s, P) for P in Ps]
    c.line(Ps, sp, SER[k], w=1.6); c.dots(Ps, sp, SER[k])
    items.append((c.Y(sp[-1]), f"serial {int(s*100)}%", SER[k]))
legend_right(c, items)
save("amdahl_compute.svg", 820, 440, [c])

# 2. compute vs memory, no serial part
c = Chart(820, 400, "No serial part at all: compute-bound vs memory-bound", 1, 8, 0, 8.5, R=150, T=84)
c.axes([1, 2, 3, 4, 5, 6, 7, 8], [0, 2, 4, 6, 8], "threads", "speedup")
Ps = list(range(1, 9))
c.line(Ps, Ps, MUTED, dash="2 3", w=1)
for kind, T1, col, lab in [("compute", T1c, SER[0], "compute-bound"), ("memory", T1m, SER[1], "memory-bound")]:
    sp = [T1 / T(kind, 0.0, P) for P in Ps]
    c.line(Ps, sp, col, w=1.8); c.dots(Ps, sp, col)
    c.label(c.X(8) + 10, c.Y(sp[-1]), lab, col)
save("amdahl_memory.svg", 820, 400, [c])

# 3. MPI strong and weak scaling (metis*, median of 3)
def parse(f):
    d = {}
    for l in open(f):
        if l.startswith(("metis", "simple")):
            p = l.strip().split(","); d[p[0]] = float(p[9]); d[p[0] + "_exch"] = float(p[6]); d[p[0] + "_sync"] = float(p[8])
        if l.startswith("Serial"): d["serial"] = float(re.findall(r"([\d.]+) ms", l)[0])
    return d
med = {}
for kind in ["strong", "weak"]:
    for n in [1, 2, 3, 4]:
        runs = [parse(f) for f in sorted(glob.glob(os.path.join(src, f"{kind}_np{n}_r*.txt")))]
        med[(kind, n)] = {k: statistics.median([r[k] for r in runs]) for k in runs[0]}
ser = statistics.median([med[("strong", n)]["serial"] for n in [1, 2, 3, 4]])
with open(os.path.join(out, "mpi_scaling.csv"), "w") as o:
    o.write("kind,processes,iteration_ms_metis,iteration_ms_simple,exchange_ms,sync_ms\n")
    for (kind, n), d in sorted(med.items()):
        o.write(f"{kind},{n},{d['metis*']:.4f},{d['simple']:.4f},{d['metis*_exch']:.4f},{d['metis*_sync']:.4f}\n")
parts = []
c = Chart(410, 400, "Strong scaling, 931k cells", 1, 4, 0, 4.2, L=56, R=20, T=84)
c.axes([1, 2, 3, 4], [0, 1, 2, 3, 4], "processes", "speedup vs serial code")
c.line([1, 4], [1, 4], MUTED, dash="2 3", w=1)
sp = [ser / med[("strong", n)]["metis*"] for n in [1, 2, 3, 4]]
c.line([1, 2, 3, 4], sp, SER[2], w=1.8); c.dots([1, 2, 3, 4], sp, SER[2])
parts.append(c)
c = Chart(410, 400, "Weak scaling, 233k cells per process", 1, 4, 0, 1.1, L=56, R=20, T=84, ox=410)
c.axes([1, 2, 3, 4], [0, 0.25, 0.5, 0.75, 1], "processes", "efficiency")
c.line([1, 4], [1, 1], MUTED, dash="2 3", w=1)
ef = [med[("weak", 1)]["metis*"] / med[("weak", n)]["metis*"] for n in [1, 2, 3, 4]]
c.line([1, 2, 3, 4], ef, SER[2], w=1.8); c.dots([1, 2, 3, 4], ef, SER[2])
parts.append(c)
parts[0].label(66, 52, "measured (metis*, median of 3)", SER[2])
parts[0].o.append(f'<line x1="330" x2="348" y1="52" y2="52" stroke="{MUTED}" stroke-dasharray="2 3"/><text x="354" y="56" fill="{FG}">ideal</text>')
save("mpi_scaling.svg", 820, 400, parts)
print("serial %.3f ms; strong speedups" % ser, [round(x, 2) for x in sp], "weak efficiency", [round(x, 2) for x in ef])
for n in [2, 3, 4]:
    S = sp[n - 1]; print("Karp-Flatt P=%d: %.2f" % (n, (1 / S - 1 / n) / (1 - 1 / n)))
print("compute T1 %.3f memory T1 %.3f" % (T1c, T1m))
