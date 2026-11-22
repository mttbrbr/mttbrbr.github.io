"""Dark-theme SVG charts for the core-communication article (no dependencies).
Usage: python3 plot_comm.py pp_cores.csv halo_runs.csv outdir
"""
import csv, math, statistics, sys
from collections import defaultdict

pp_csv, halo_csv, out = sys.argv[1], sys.argv[2], sys.argv[3]
FG, MUTED, GRID, BG = "#e6e6e6", "#8a8f94", "#23272b", "#08090a"
S1, S2, S3 = "#3987e5", "#d95926", "#199e70"          # validated on #08090a
FONT = 'font-family="ui-monospace, SFMono-Regular, Menlo, monospace" font-size="12"'

def human(b):
    for u, k in (("MB", 1 << 20), ("KB", 1 << 10)):
        if b >= k: return f"{b // k} {u}"
    return f"{b} B"

pp = [(int(r["bytes"]), float(r["one_way_us"]), float(r["bandwidth_GBps"])) for r in csv.DictReader(open(pp_csv))]

def line_chart(fname, title, ylab, ys, ylog, marks):
    W, H, L, R, T, B = 820, 400, 70, 30, 48, 56
    xs = [math.log2(b) for b, *_ in pp]
    x0, x1 = xs[0] - 0.5, xs[-1] + 0.5
    if ylog:
        y0, y1 = math.floor(math.log10(min(ys))), math.ceil(math.log10(max(ys)))
        Y = lambda v: H - B - (math.log10(v) - y0) / (y1 - y0) * (H - T - B)
        ticks = [10 ** e for e in range(y0, y1 + 1)]
    else:
        y0, y1 = 0, math.ceil(max(ys) / 5) * 5
        Y = lambda v: H - B - (v - y0) / (y1 - y0) * (H - T - B)
        ticks = list(range(0, int(y1) + 1, 5))
    X = lambda v: L + (v - x0) / (x1 - x0) * (W - L - R)
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" {FONT}>',
         f'<rect width="{W}" height="{H}" fill="{BG}"/>',
         f'<text x="{L}" y="24" fill="{FG}" font-size="14">{title}</text>']
    for t in ticks:
        o.append(f'<line x1="{L}" x2="{W-R}" y1="{Y(t):.1f}" y2="{Y(t):.1f}" stroke="{GRID}"/>')
        o.append(f'<text x="{L-8}" y="{Y(t)+4:.1f}" fill="{MUTED}" text-anchor="end">{t:g}</text>')
    o.append(f'<text x="18" y="{(T+H-B)/2:.0f}" fill="{MUTED}" transform="rotate(-90 18 {(T+H-B)/2:.0f})" text-anchor="middle">{ylab}</text>')
    for e in range(3, 27, 3):
        if x0 <= e <= x1:
            o.append(f'<text x="{X(e):.1f}" y="{H-B+18}" fill="{MUTED}" text-anchor="middle">{human(1 << e)}</text>')
    o.append(f'<text x="{(L+W-R)/2:.0f}" y="{H-12}" fill="{MUTED}" text-anchor="middle">message size (log scale)</text>')
    for b, lab in marks:
        x = X(math.log2(b))
        o.append(f'<line x1="{x:.1f}" x2="{x:.1f}" y1="{T}" y2="{H-B}" stroke="{MUTED}" stroke-dasharray="4 4"/>')
        o.append(f'<text x="{x+6:.1f}" y="{T+12}" fill="{MUTED}">{lab}</text>')
    pts = " ".join(f"{X(x):.1f},{Y(y):.1f}" for x, y in zip(xs, ys))
    o.append(f'<polyline points="{pts}" fill="none" stroke="{S1}" stroke-width="2" stroke-linejoin="round"/>')
    for x, y in zip(xs, ys):
        o.append(f'<circle cx="{X(x):.1f}" cy="{Y(y):.1f}" r="4" fill="{S1}" stroke="{BG}" stroke-width="2"/>')
    o.append("</svg>")
    open(f"{out}/{fname}", "w").write("\n".join(o))

line_chart("pingpong_latency.svg", "One-way time to send a message between two cores (MPI, shared memory)",
           "time (µs), log scale", [p[1] for p in pp], True, [(4096, "protocol switch 4 KB"), (12 << 20, "L3 12 MB")])
line_chart("pingpong_bandwidth.svg", "Bandwidth between two cores vs message size",
           "bandwidth (GB/s)", [p[2] for p in pp], False, [(4096, "protocol switch 4 KB"), (12 << 20, "L3 12 MB")])

# ---- halo breakdown: stacked horizontal bars, one panel per mesh ----
cols = "case,method,cells_max,halo_cells_total,msgs_per_iter,MB_per_iter,compute_ms_max,exchange_ms_avg,exchange_ms_max,sync_ms_avg,iter_ms,speedup,check".split(",")
d = defaultdict(list)
for r in csv.reader(open(halo_csv)):
    row = dict(zip(cols, r)); d[(row["case"], row["method"])].append(row)
med = lambda k, c: statistics.median(float(x[c]) for x in d[k])
methods = ["random", "simple", "scotch*", "metis*"]

print("| Mesh | Method | Compute (ms) | Exchange (ms) | Sync/wait (ms) | Iteration (ms) | Speedup |")
print("|---|---|---|---|---|---|---|")
for case in ("small", "big"):
    for m in methods:
        k = (case, m)
        print(f"| {case} | {m} | {med(k,'compute_ms_max'):.3f} | {med(k,'exchange_ms_avg'):.3f} | {med(k,'sync_ms_avg'):.3f} | {med(k,'iter_ms'):.3f} | {med(k,'speedup'):.2f}× |")

def panel(fname, case, title, serial):
    W, L, R, T, rowh = 820, 90, 40, 70, 46
    H = T + rowh * len(methods) + 50
    xmax = max(med((case, m), "iter_ms") for m in methods + []) * 1.08
    xmax = max(xmax, serial * 1.15)
    X = lambda v: L + v / xmax * (W - L - R)
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" {FONT}>',
         f'<rect width="{W}" height="{H}" fill="{BG}"/>',
         f'<text x="{L}" y="24" fill="{FG}" font-size="14">{title}</text>']
    # legend
    lx = L
    for col, lab in ((S1, "compute"), (S2, "halo exchange"), (S3, "sync / waiting")):
        o.append(f'<rect x="{lx}" y="38" width="10" height="10" rx="2" fill="{col}"/>')
        o.append(f'<text x="{lx+15}" y="47" fill="{FG}">{lab}</text>')
        lx += 30 + 8 * len(lab)
    step = 0.5 if xmax < 3 else (1 if xmax < 6 else 2)
    v = 0.0
    while v <= xmax:
        o.append(f'<line x1="{X(v):.1f}" x2="{X(v):.1f}" y1="{T-6}" y2="{H-44}" stroke="{GRID}"/>')
        o.append(f'<text x="{X(v):.1f}" y="{H-28}" fill="{MUTED}" text-anchor="middle">{v:g}</text>')
        v += step
    o.append(f'<text x="{(L+W-R)/2:.0f}" y="{H-8}" fill="{MUTED}" text-anchor="middle">ms per iteration (slowest rank), lower is better</text>')
    for i, m in enumerate(methods):
        k = (case, m); y = T + i * rowh
        comp, exch, sync = med(k, "compute_ms_max"), med(k, "exchange_ms_avg"), med(k, "sync_ms_avg")
        o.append(f'<text x="{L-10}" y="{y+21}" fill="{FG}" text-anchor="end">{m}</text>')
        x = L
        for val, col in ((comp, S1), (exch, S2), (sync, S3)):
            w = max(val / xmax * (W - L - R) - 2, 0.5)
            o.append(f'<rect x="{x:.1f}" y="{y+6}" width="{w:.1f}" height="22" rx="2" fill="{col}"/>')
            x += w + 2
        o.append(f'<text x="{x+6:.1f}" y="{y+21}" fill="{FG}">{med(k,"speedup"):.2f}×</text>')
    sx = X(serial)
    o.append(f'<line x1="{sx:.1f}" x2="{sx:.1f}" y1="{T-6}" y2="{H-44}" stroke="{FG}" stroke-dasharray="5 4"/>')
    o.append(f'<text x="{sx+6:.1f}" y="{T-10}" fill="{FG}">1 core: {serial:.2f} ms</text>')
    o.append("</svg>")
    open(f"{out}/{fname}", "w").write("\n".join(o))

serial = lambda case: statistics.median(float(x["iter_ms"]) * float(x["speedup"]) for m in methods for x in d[(case, m)])
panel("halo_small.svg", "small", "4 cores, 233k cells: every subdomain fits in cache", serial("small"))
panel("halo_big.svg", "big", "4 cores, 931k cells: limited by shared memory bandwidth", serial("big"))
