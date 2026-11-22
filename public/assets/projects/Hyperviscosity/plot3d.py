"""Charts for the hyperviscosity (3D) article. Usage: python3 plot3d.py datadir outdir"""
import sys, os, csv, math
_here = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [_here, os.path.join(_here, "..", "common")]
from svgchart import *

src, out = sys.argv[1], sys.argv[2]
set_out(out)
EPS = 0.1
RUNS = [("64", "1", "p = 1 (Navier-Stokes), 64³", SER[0], None), ("64", "2", "p = 2, 64³", SER[1], None),
        ("64", "4", "p = 4, 64³", SER[2], None), ("64", "8", "p = 8, 64³", SER[3], None),
        ("96", "1", "p = 1 (Navier-Stokes), 96³", SER[4], "6 4")]
spec = {(N, p): [(int(r["k"]), float(r["E"]), float(r["D"])) for r in csv.DictReader(open(os.path.join(src, f"spectrum3d_N{N}_p{p}.csv")))]
        for N, p, *_ in RUNS}

def legend(c, x0, rows=((0, 3), (3, 5)), y0=50):
    for r, (a, b) in enumerate(rows):
        lx = x0
        for N, p, lab, col, dash in RUNS[a:b]:
            c.label(lx, y0 + 22 * r, lab, col); lx += 30 + 7.4 * len(lab)

# 1. energy spectra
c = Chart(820, 460, "Energy spectrum, averaged over the steady state", 1, 40, 1e-5, 0.3, R=40, T=110, xlog=True, ylog=True)
c.axes([1, 2, 5, 10, 20, 40], [1e-5, 1e-4, 1e-3, 1e-2, 0.1], "wavenumber k (log)", "E(k) (log)", fy=LOGF)
for N, p, lab, col, dash in RUNS:
    d = [(k, E) for k, E, _ in spec[(N, p)] if E > 0]
    c.line([k for k, _ in d], [E for _, E in d], col, w=1.7, dash=dash)
ks = [3 + i * 0.2 for i in range(60)]
c.line(ks, [1.6 * EPS ** (2 / 3) * k ** (-5 / 3) * 1.6 for k in ks], FG, dash="2 3", w=1)
c.o.append(f'<text x="{c.X(8):.1f}" y="{c.Y(1.6 * 1.6 * EPS ** (2/3) * 8 ** (-5/3)) - 10:.1f}" fill="{FG}">k^(-5/3)</text>')
c.o.append(f'<line x1="{c.X(64/3):.1f}" x2="{c.X(64/3):.1f}" y1="{c.Y(0.3)}" y2="{c.Y(1e-5)}" stroke="{MUTED}" stroke-dasharray="4 4"/>')
c.o.append(f'<text x="{c.X(64/3)-6:.1f}" y="{c.Y(2e-5):.1f}" fill="{MUTED}" text-anchor="end">k_max, 64³</text>')
legend(c, c.X(1))
save("spectra.svg", 820, 460, [c])

# 2. compensated spectra
c = Chart(820, 460, "Compensated spectrum: flat at the Kolmogorov constant if the cascade were ideal", 1, 34, 0, 5.2, R=40, T=110, xlog=True)
c.axes([1, 2, 5, 10, 20, 30], [0, 1, 2, 3, 4, 5], "wavenumber k (log)", "ε^(-2/3) k^(5/3) E(k)")
c.line([1, 34], [1.6, 1.6], FG, dash="2 3", w=1)
c.o.append(f'<text x="{c.X(1.1):.1f}" y="{c.Y(1.6) - 8:.1f}" fill="{FG}">Kolmogorov constant ≈ 1.6</text>')
for N, p, lab, col, dash in RUNS:
    d = [(k, EPS ** (-2 / 3) * k ** (5 / 3) * E) for k, E, _ in spec[(N, p)] if E > 0 and k <= int(int(N) / 3)]
    c.line([k for k, _ in d], [v for _, v in d], col, w=1.8, dash=dash)
legend(c, c.X(1))
save("compensated.svg", 820, 460, [c])

# 3. dissipation spectrum normalized
c = Chart(820, 440, "Where the energy is dissipated (dissipation spectrum, normalized)", 1, 34, 0, 1.05, R=40, T=110, xlog=True)
c.axes([1, 2, 5, 10, 20, 30], [0, 0.25, 0.5, 0.75, 1], "wavenumber k (log)", "D(k) / max D(k)")
for N, p, lab, col, dash in RUNS:
    d = [(k, D) for k, _, D in spec[(N, p)] if k <= int(int(N) / 3)]
    m = max(D for _, D in d)
    c.line([k for k, _ in d], [D / m for _, D in d], col, w=1.8, dash=dash)
legend(c, c.X(1))
save("dissipation.svg", 820, 440, [c])

# 4. energy vs time
c = Chart(820, 400, "Kinetic energy vs time (statistics averaged over 15 ≤ t ≤ 30)", 0, 30, 0, 0.8, R=40, T=110)
c.axes([0, 5, 10, 15, 20, 25, 30], [0, 0.2, 0.4, 0.6, 0.8], "time", "kinetic energy")
c.o.append(f'<rect x="{c.X(15):.1f}" y="{c.Y(0.8):.1f}" width="{c.X(30)-c.X(15):.1f}" height="{c.Y(0)-c.Y(0.8):.1f}" fill="{FG}" opacity="0.04"/>')
for N, p, lab, col, dash in RUNS:
    d = [(float(r["t"]), float(r["energy"])) for r in csv.DictReader(open(os.path.join(src, f"history3d_N{N}_p{p}.csv")))]
    c.line([t for t, _ in d], [e for _, e in d], col, w=1.5, dash=dash)
legend(c, c.X(0))
save("energy.svg", 820, 400, [c])
print("ok")
