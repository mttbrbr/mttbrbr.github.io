"""Richardson extrapolation and Grid Convergence Index (Roache; Celik et al. 2008)
from grids.csv, plus a chart. Usage: python3 gci_analysis.py datadir outdir"""
import sys, os, csv, math
_here = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [_here, os.path.join(_here, "..", "common")]
from svgchart import *

src, out = sys.argv[1], sys.argv[2]
set_out(out)
g = list(csv.DictReader(open(os.path.join(src, "grids.csv"))))
Ns = [int(r["N"]) for r in g]
Q = [("u_centre", "u at the centre", -0.20581), ("v_centre", "v at the centre", 0.05454), ("kinetic_energy", "kinetic energy", None)]
REF = "Ghia et al. (1982), 129 x 129, second order"

rows = []
for key, lab, ref in Q:
    f = [float(r[key]) for r in g]
    for i in range(len(f) - 2):
        f3, f2, f1 = f[i], f[i + 1], f[i + 2]          # coarse, medium, fine
        rr = 2.0
        e32, e21 = f3 - f2, f2 - f1
        p = math.log(abs(e32 / e21)) / math.log(rr) if e21 != 0 and e32 / e21 > 0 else float("nan")
        fext = f1 + (f1 - f2) / (rr ** p - 1)
        gci = 1.25 * abs((f1 - f2) / f1) / (rr ** p - 1)
        rows.append((lab, Ns[i], Ns[i + 1], Ns[i + 2], f1, p, fext, gci, ref))
with open(os.path.join(out, "gci.csv"), "w") as o:
    o.write("quantity,N_coarse,N_medium,N_fine,value_fine,observed_order,extrapolated,GCI_fine,reference\n")
    for r in rows: o.write(",".join(str(x) for x in r) + "\n")
for r in rows:
    print(f"{r[0]:16s} grids {r[1]:4d}/{r[2]:4d}/{r[3]:4d}  fine {r[4]: .6f}  p = {r[5]:.3f}  extrapolated {r[6]: .6f}  GCI {100*r[7]:.2f}%"
          + (f"  ref {r[8]}" if r[8] is not None else ""))

# chart: u at the centre vs cell size, with the Richardson extrapolation and the reference
f = [float(r["u_centre"]) for r in g]; hs = [1.0 / n for n in Ns]
last = [r for r in rows if r[0] == "u at the centre"][-1]
p, fext = last[5], last[6]
c = Chart(820, 420, "u at the centre of the cavity vs cell size, Re = 100", 0, 0.068, -0.214, -0.16, R=40, T=84)
c.axes([0, 1/64, 1/32, 1/16], [-0.21, -0.20, -0.19, -0.18, -0.17, -0.16], "cell size h", "u(0.5, 0.5)",
       fx=lambda v: "0" if v == 0 else f"1/{round(1/v)}")
# fit through the finest grid: f = fext + C h^p
C = (f[-1] - fext) / hs[-1] ** p
hh = [i * 0.068 / 100 for i in range(101)]
c.line(hh, [fext + C * x ** p for x in hh], MUTED, dash="5 4", w=1.5)
c.dots(hs, f, SER[0]); c.line(hs, f, SER[0])
c.o.append(f'<circle cx="{c.X(0):.1f}" cy="{c.Y(fext):.1f}" r="6" fill="none" stroke="{SER[2]}" stroke-width="2"/>')
c.o.append(f'<rect x="{c.X(0)-5:.1f}" y="{c.Y(-0.20581)-5:.1f}" width="10" height="10" fill="{SER[1]}"/>')
c.label(c.X(0), 50, "computed (first-order upwind)", SER[0])
c.label(c.X(0) + 270, 50, f"Richardson extrapolation (p = {p:.2f})", SER[2])
c.label(c.X(0) + 270, 70, "Ghia et al. 1982", SER[1])
c.o.append(f'<line x1="{c.X(0)}" x2="{c.X(0)+18}" y1="70" y2="70" stroke="{MUTED}" stroke-dasharray="5 4"/><text x="{c.X(0)+24}" y="74" fill="{FG}">fit f = f_ext + C h^p</text>')
save("gci.svg", 820, 420, [c])
print("ok")
