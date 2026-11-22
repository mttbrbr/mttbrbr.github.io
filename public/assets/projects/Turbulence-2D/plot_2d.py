"""Charts for the 2D turbulence article. Usage: python3 plot_2d.py snapdir decaydir outdir
snapdir : 512^2 run with vorticity snapshots (nu = 1.25e-4); decaydir: 256^2 run (nu = 5e-4)"""
import sys, os, csv, math, subprocess
_here = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [_here, os.path.join(_here, "..", "common")]
from svgchart import *

snap, dec, out = sys.argv[1], sys.argv[2], sys.argv[3]
set_out(out)

# 1. vorticity snapshots
F = subprocess.run(["fc-match", "-f", "%{file}", "monospace"], capture_output=True, text=True).stdout
pngs = []
for t in ["00", "05", "10", "30"]:
    rows = [list(map(float, l.split(","))) for l in open(os.path.join(snap, f"vort_N512_p1_t{t}.csv"))]
    vals = [x for r in rows for x in r]; n = len(rows)
    m = sorted(abs(v) for v in vals)[int(0.995 * len(vals))]
    p = os.path.join(out, f"_w{t}.png")
    field_png(vals, n, n, p, scale=1, vmin=-m, vmax=m, diverging=True)
    lab = os.path.join(out, f"_l{t}.png")
    subprocess.run(["magick", p, "-resize", "300x300", "-background", BG, "-fill", FG, "-font", F, "-pointsize", "16",
                    "-gravity", "north", "-splice", "0x30", "-annotate", "+0+6", f"t = {int(t)}", "-bordercolor", BG, "-border", "6", "+repage", lab], check=True)
    pngs.append(lab); os.remove(p)
subprocess.run(["magick"] + pngs + ["+append", "+repage", os.path.join(out, "vorticity.png")], check=True)
for p in pngs: os.remove(p)

# 2. energy and enstrophy vs time, normalized
def hist(d, N):
    r = list(csv.DictReader(open(os.path.join(d, f"history_N{N}_p1.csv"))))
    return [(float(x["t"]), float(x["energy"]), float(x["enstrophy"])) for x in r]
H = [(hist(dec, 256), "ν = 5·10⁻⁴ (256²)", SER[0]), (hist(snap, 512), "ν = 1.25·10⁻⁴ (512²)", SER[2])]
parts = []
for q, (idx, title) in enumerate([(1, "Energy / initial energy"), (2, "Enstrophy / initial enstrophy")]):
    c = Chart(410, 400, title, 0, 30, 0, 1.05, L=56, R=20, T=84, ox=q * 410)
    c.axes([0, 10, 20, 30], [0, 0.25, 0.5, 0.75, 1], "time", "")
    for h, lab, col in H:
        c.line([x[0] for x in h], [x[idx] / h[0][idx] for x in h], col, w=1.8)
    parts.append(c)
lx = 66
for h, lab, col in H:
    parts[0].label(lx, 52, lab, col); lx += 30 + 7.4 * len(lab)
save("energy_enstrophy.svg", 820, 400, parts)

# 3. spectra at several times (512^2), with k^-3 and k^-5 references
r = list(csv.DictReader(open(os.path.join(snap, "spectra_N512_p1.csv"))))
c = Chart(820, 440, "Energy spectrum at different times (512², ν = 1.25·10⁻⁴)", 1, 170, 1e-12, 1, R=40, T=84, xlog=True, ylog=True)
c.axes([1, 10, 100], [1e-12, 1e-9, 1e-6, 1e-3, 1], "wavenumber k (log)", "E(k) (log)", fy=LOGF)
for k, (t, col) in enumerate([(5.0, SER[0]), (10.0, SER[1]), (20.0, SER[2]), (30.0, SER[3])]):
    d = [(int(x["k"]), float(x["E"])) for x in r if float(x["t"]) == t and float(x["E"]) > 0]
    c.line([a for a, _ in d], [b for _, b in d], col, w=1.6)
    c.label(c.X(1) + 8 + k * 90, 50, f"t = {int(t)}", col)
ks = [10 + i for i in range(60)]
c.line(ks, [2e-2 * (k / 10) ** -3 for k in ks], FG, dash="2 3", w=1)
c.o.append(f'<text x="{c.X(40):.1f}" y="{c.Y(2e-2 * 4 ** -3) - 8:.1f}" fill="{FG}">k⁻³ (theory)</text>')
c.line(ks, [3e-3 * (k / 10) ** -5 for k in ks], MUTED, dash="2 3", w=1)
c.o.append(f'<text x="{c.X(30):.1f}" y="{c.Y(3e-3 * 3 ** -5) + 18:.1f}" fill="{MUTED}">k⁻⁵</text>')
save("spectra.svg", 820, 440, [c])
print("ok")
