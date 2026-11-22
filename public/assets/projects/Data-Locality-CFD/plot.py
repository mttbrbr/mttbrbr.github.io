"""Turn results.csv into a dark-theme SVG (no dependencies). Median of repetitions.
Usage: python3 plot.py results.csv out.svg
"""
import csv, math, statistics, sys
from collections import defaultdict

src, out = sys.argv[1], sys.argv[2]
runs = defaultdict(list)
mb = {}
for r in csv.DictReader(open(src)):
    key = (r["layout"], r["order"])
    n = int(r["N"])
    runs[(key, n)].append(float(r["ns_per_cell_update"]))
    mb[n] = float(r["working_set_MB"])

series = [  # fixed categorical order (validated on #08090a)
    (("SoA", "i-j (contiguous)"), "SoA, contiguous", "#3987e5"),
    (("AoS", "i-j (contiguous)"), "AoS, contiguous", "#d95926"),
    (("AoS", "j-i (strided)"), "AoS, strided", "#199e70"),
    (("SoA", "j-i (strided)"), "SoA, strided", "#c98500"),
]
Ns = sorted(mb)
med = {(k, n): statistics.median(runs[(k, n)]) for (k, n) in runs}

# print a markdown table for the article
print("| Grid | Working set | " + " | ".join(s[1] for s in series) + " |")
print("|---|---|" + "---|" * len(series))
for n in Ns:
    print(f"| {n}² | {mb[n]:.1f} MB | " + " | ".join(f"{med[(s[0], n)]:.2f}" for s in series) + " |")

W, H = 820, 480
L, R, T, B = 64, 150, 48, 56
xmin, xmax = math.log10(min(mb.values()) / 1.3), math.log10(max(mb.values()) * 1.3)
ymax = max(med.values()) * 1.08
ymax = math.ceil(ymax / 2) * 2
X = lambda v: L + (math.log10(v) - xmin) / (xmax - xmin) * (W - L - R)
Y = lambda v: H - B - v / ymax * (H - T - B)

FG, MUTED, GRID, BG = "#e6e6e6", "#8a8f94", "#23272b", "#08090a"
o = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" font-family="ui-monospace, SFMono-Regular, Menlo, monospace" font-size="12">',
     f'<rect width="{W}" height="{H}" fill="{BG}"/>',
     f'<text x="{L}" y="24" fill="{FG}" font-size="14">Time per cell update vs. working set (i7-1165G7, gcc -O3 -march=native)</text>']
# y grid
step = 2 if ymax <= 20 else 5
for v in range(0, int(ymax) + 1, step):
    o.append(f'<line x1="{L}" x2="{W-R}" y1="{Y(v):.1f}" y2="{Y(v):.1f}" stroke="{GRID}" stroke-width="1"/>')
    o.append(f'<text x="{L-8}" y="{Y(v)+4:.1f}" fill="{MUTED}" text-anchor="end">{v}</text>')
o.append(f'<text x="16" y="{(T+H-B)/2:.0f}" fill="{MUTED}" transform="rotate(-90 16 {(T+H-B)/2:.0f})" text-anchor="middle">ns per cell update</text>')
# x ticks
for v in [0.1, 1, 10, 100, 1000]:
    if math.log10(v) < xmin or math.log10(v) > xmax: continue
    o.append(f'<text x="{X(v):.1f}" y="{H-B+18}" fill="{MUTED}" text-anchor="middle">{v:g} MB</text>')
o.append(f'<text x="{(L+W-R)/2:.0f}" y="{H-12}" fill="{MUTED}" text-anchor="middle">working set (6 fields × N² doubles), log scale</text>')
# cache boundaries
for v, lab in [(1.25, "L2 1.25 MB"), (12, "L3 12 MB")]:
    o.append(f'<line x1="{X(v):.1f}" x2="{X(v):.1f}" y1="{T}" y2="{H-B}" stroke="{MUTED}" stroke-width="1" stroke-dasharray="4 4"/>')
    o.append(f'<text x="{X(v)+6:.1f}" y="{T+12}" fill="{MUTED}">{lab}</text>')
# series
ends = []
for key, name, col in series:
    pts = [(X(mb[n]), Y(med[(key, n)])) for n in Ns]
    d = " ".join(f"{x:.1f},{y:.1f}" for x, y in pts)
    o.append(f'<polyline points="{d}" fill="none" stroke="{col}" stroke-width="2" stroke-linejoin="round"/>')
    for x, y in pts:
        o.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="4" fill="{col}" stroke="{BG}" stroke-width="2"/>')
    ends.append([pts[-1][1], name, col, pts[-1][0]])
# annotate the power-of-two spikes on the worst series
worst = series[3][0]
for n in Ns:
    if n >= 1024 and n & (n - 1) == 0:
        o.append(f'<text x="{X(mb[n]):.1f}" y="{Y(med[(worst, n)])-12:.1f}" fill="{FG}" text-anchor="middle">N={n}</text>')
# direct labels, de-collided
ends.sort()
for i in range(1, len(ends)):
    if ends[i][0] - ends[i-1][0] < 16: ends[i][0] = ends[i-1][0] + 16
for y, name, col, x in ends:
    o.append(f'<rect x="{x+12}" y="{y-5}" width="10" height="10" rx="2" fill="{col}"/>')
    o.append(f'<text x="{x+27}" y="{y+4}" fill="{FG}">{name}</text>')
o.append("</svg>")
open(out, "w").write("\n".join(o))
