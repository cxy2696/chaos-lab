#!/usr/bin/env python3
"""Render all figures in docs/figures from the committed data files.

Run with the bundled Python (Pillow is required):
    /path/to/python scripts/plot_figures.py
"""
import csv
import math
import os

from PIL import Image

from miniplot import Canvas, hstack, vstack

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
FIGS = os.path.join(ROOT, "docs", "figures")
os.makedirs(FIGS, exist_ok=True)

NAVY = (31, 58, 147)
RED = (192, 57, 43)
ORANGE = (211, 84, 0)
PURPLE = (142, 68, 173)
TEAL = (22, 160, 133)
GRAY = (110, 118, 130)


def read_csv(name):
    with open(os.path.join(DATA, name), newline="") as fh:
        rows = list(csv.reader(fh))
    header = rows[0]
    body = []
    for row in rows[1:]:
        if not row:
            continue
        body.append([float(v) if v.strip() else float("nan") for v in row])
    return header, body


def col(rows, i):
    return [r[i] for r in rows]


def logistic_bifurcation():
    _, bif = read_csv("logistic_bifurcation.csv")
    _, lyap = read_csv("logistic_lyapunov.csv")

    top = Canvas(1500, 780, margins=(100, 40, 70, 80))
    top.set_limits(2.5, 4.0, 0.0, 1.0, pad=0.005)
    top.scatter(col(bif, 0), col(bif, 1), color=NAVY, radius=0.35)
    top.axes("r", "x", "Logistic map: bifurcation diagram",
             "x(n+1) = r x(n) (1 - x(n)), 120000 points")
    p1 = top.save(os.path.join(FIGS, "logistic_bifurcation.png"))

    bottom = Canvas(1500, 620, margins=(100, 40, 70, 80))
    xs, ys = col(lyap, 0), col(lyap, 1)
    bottom.set_limits(2.5, 4.0, -2.2, 0.8, pad=0.02)
    bottom.horizontal(0.0, color=(80, 88, 100), width=1.0, dash=8)
    bottom.vertical(3.569945672, color=(150, 150, 150), width=0.9, dash=8)
    bottom.line(xs, ys, color=RED, width=1.2)
    bottom.text(3.575, 0.55, "r_inf = 3.5699457", color=GRAY, size=12)
    bottom.axes("r", "lambda(r)", "Lyapunov exponent of the logistic map",
                "lambda = ln 2 at r = 4; positive means chaos")
    p2 = bottom.save(os.path.join(FIGS, "logistic_lyapunov.png"))
    vstack([p1, p2], os.path.join(FIGS, "logistic.png"))
    print("logistic.png")


def lorenz_attractor():
    _, rows = read_csv("lorenz_trajectory.csv")
    x, y, z = col(rows, 1), col(rows, 2), col(rows, 3)

    a = Canvas(980, 760, margins=(90, 40, 65, 75))
    a.set_limits(min(x), max(x), min(z), max(z), pad=0.03)
    a.gradient_line(x, z, c0=(31, 58, 147), c1=(192, 57, 43), width=0.75)
    a.scatter([x[0]], [z[0]], color=(39, 174, 96), radius=2.8)
    a.axes("x", "z", "Lorenz attractor (x-z projection)", "RK4, h = 0.01")
    pa = a.save(os.path.join(FIGS, "lorenz_xz.png"))

    b = Canvas(980, 760, margins=(90, 40, 65, 75))
    b.set_limits(min(y), max(y), min(z), max(z), pad=0.03)
    b.gradient_line(y, z, c0=(31, 58, 147), c1=(211, 84, 0), width=0.75)
    b.scatter([y[0]], [z[0]], color=(39, 174, 96), radius=2.8)
    b.axes("y", "z", "Lorenz attractor (y-z projection)", "green dot: initial condition")
    pb = b.save(os.path.join(FIGS, "lorenz_yz.png"))
    hstack([pa, pb], os.path.join(FIGS, "lorenz_attractor.png"))
    print("lorenz_attractor.png")


def lorenz_series():
    _, rows = read_csv("lorenz_series.csv")
    t, x = col(rows, 0), col(rows, 1)
    c = Canvas(1500, 600, margins=(100, 40, 65, 80))
    c.set_limits(min(t), max(t), min(x), max(x), pad=0.03)
    c.line(t, x, color=NAVY, width=0.9)
    c.axes("t", "x(t)", "Lorenz x(t): the scalar series used for embedding",
           "t in [50, 130], dt = 0.01")
    c.save(os.path.join(FIGS, "lorenz_series.png"))
    print("lorenz_series.png")


def henon_attractor():
    _, rows = read_csv("henon_orbit.csv")
    c = Canvas(900, 780, margins=(90, 40, 65, 80))
    c.set_limits(-1.4, 1.35, -0.42, 0.42, pad=0.02)
    c.scatter(col(rows, 1), col(rows, 2), color=PURPLE, radius=0.45)
    c.axes("x", "y", "Henon map attractor", "a = 1.4, b = 0.3, 5000 iterations")
    c.save(os.path.join(FIGS, "henon_attractor.png"))
    print("henon_attractor.png")


def rossler_attractor():
    _, rows = read_csv("rossler_trajectory.csv")
    x, y = col(rows, 1), col(rows, 2)
    c = Canvas(900, 780, margins=(90, 40, 65, 80))
    c.set_limits(min(x), max(x), min(y), max(y), pad=0.03)
    c.gradient_line(x, y, c0=(22, 160, 133), c1=(192, 57, 43), width=0.8)
    c.axes("x", "y", "Rossler attractor (x-y projection)", "a = b = 0.2, c = 5.7")
    c.save(os.path.join(FIGS, "rossler_attractor.png"))
    print("rossler_attractor.png")


def chua_attractor():
    _, rows = read_csv("chua_trajectory.csv")
    x, z = col(rows, 1), col(rows, 3)
    c = Canvas(950, 720, margins=(90, 40, 65, 80))
    c.set_limits(min(x), max(x), min(z), max(z), pad=0.03)
    c.gradient_line(x, z, c0=(31, 58, 147), c1=(192, 57, 43), width=0.55)
    c.axes("x", "z", "Chua circuit: double scroll", "alpha=15.6, beta=28, m0=-1.143, m1=-0.714")
    c.save(os.path.join(FIGS, "chua_attractor.png"))
    print("chua_attractor.png")


def lorenz_poincare():
    _, rows = read_csv("lorenz_poincare.csv")
    c = Canvas(900, 780, margins=(90, 40, 65, 80))
    c.set_limits(min(col(rows, 0)), max(col(rows, 0)),
                 min(col(rows, 1)), max(col(rows, 1)), pad=0.05)
    c.scatter(col(rows, 0), col(rows, 1), color=ORANGE, radius=2.0)
    c.axes("x", "y", "Lorenz Poincare section", "z = 27, downward crossings, 240 points")
    c.save(os.path.join(FIGS, "lorenz_poincare.png"))
    print("lorenz_poincare.png")


def lorenz_gp():
    _, rows = read_csv("lorenz_gp.csv")
    lr, lc = col(rows, 2), col(rows, 3)
    slope = [r[4] if len(r) > 4 else float("nan") for r in rows]
    valid = [(a, b) for a, b in zip(lr, slope) if not math.isnan(b) and b > 0.5]
    lr_s, slope_s = [a for a, _ in valid], [b for _, b in valid]

    a = Canvas(1100, 700, margins=(100, 40, 65, 80))
    a.set_limits(min(lr), max(lr), min(lc), max(lc), pad=0.03)
    a.line(lr, lc, color=NAVY, width=1.8)
    a.scatter(lr[::6], lc[::6], color=RED, radius=2.2)
    a.axes("log10 r", "log10 C(r)", "Correlation integral, Lorenz x(t)",
           "m = 3, tau = 17, Theiler = 50")
    pa = a.save(os.path.join(FIGS, "lorenz_gp_curve.png"))

    b = Canvas(1100, 560, margins=(100, 40, 65, 80))
    b.set_limits(min(lr_s), max(lr_s), min(slope_s) - 0.1, max(slope_s) + 0.1, pad=0.02)
    b.horizontal(2.0466, color=RED, width=1.2, dash=9)
    b.line(lr_s, slope_s, color=TEAL, width=1.6)
    b.text(lr_s[int(len(lr_s) * 0.62)], 2.10, "fitted D2 = 2.047", color=RED, size=13)
    b.axes("log10 r", "local slope", "Local scaling exponent",
           "plateau near 2.0; shaded guide is the fitted D2")
    pb = b.save(os.path.join(FIGS, "lorenz_gp_slope.png"))
    vstack([pa, pb], os.path.join(FIGS, "lorenz_gp.png"))
    print("lorenz_gp.png")


def lorenz_rosenstein():
    _, rows = read_csv("lorenz_rosenstein.csv")
    k, d = col(rows, 0), col(rows, 1)
    c = Canvas(1200, 650, margins=(100, 40, 65, 80))
    c.set_limits(min(k), max(k), min(d) - 0.1, max(d) + 0.1, pad=0.02)
    c.line(k, d, color=GRAY, width=1.3)
    fit = [(kk, dd) for kk, dd in zip(k, d) if 148 <= kk <= 297]
    c.line([p[0] for p in fit], [p[1] for p in fit], color=RED, width=3.0)
    c.text(210, -1.55, "fit: k in [148, 297]", color=RED, size=13)
    c.axes("k (samples)", "mean ln divergence", "Rosenstein divergence curve",
           "Lorenz x(t), m = 3, tau = 17, theiler = 100; lambda = 0.855")
    c.save(os.path.join(FIGS, "lorenz_rosenstein.png"))
    print("lorenz_rosenstein.png")


def lorenz_rqa():
    with open(os.path.join(DATA, "lorenz_rqa_matrix.csv"), newline="") as fh:
        rows = [[int(float(v)) for v in row] for row in csv.reader(fh) if row]
    n = len(rows)
    c = Canvas(1000, 980, margins=(100, 55, 70, 90))
    c.set_limits(0, n, 0, n, pad=0.0)
    c.imshow(rows, cmap=[(255, 255, 255), (31, 58, 147)])
    ticks = [0, n // 4, n // 2, 3 * n // 4, n]
    c.axes("i", "j", "Lorenz recurrence plot",
           f"n = {n}, RR = 2.96%, DET = 0.69, LAM = 0.29",
           xticks=ticks, yticks=ticks)
    c.save(os.path.join(FIGS, "lorenz_rqa.png"))
    print("lorenz_rqa.png")


def three_body():
    _, rows = read_csv("three_body_figure_eight.csv")
    c = Canvas(1000, 760, margins=(90, 40, 65, 80))
    bx1, by1 = col(rows, 1), col(rows, 2)
    bx2, by2 = col(rows, 3), col(rows, 4)
    bx3, by3 = col(rows, 5), col(rows, 6)
    allx = bx1 + bx2 + bx3
    ally = by1 + by2 + by3
    c.set_limits(min(allx), max(allx), min(ally), max(ally), pad=0.06)
    # All three bodies trace the same curve with a phase shift; draw it once
    # in grey, then mark the three positions at one time slice.
    c.line(bx1, by1, color=(150, 155, 165), width=1.3)
    k = min(len(bx1) - 1, int(2.0 / 0.001))
    c.scatter([bx1[k], bx2[k], bx3[k]], [by1[k], by2[k], by3[k]],
              color=(30, 30, 30), radius=2.2)
    c.scatter([bx1[k]], [by1[k]], color=NAVY, radius=4.2)
    c.scatter([bx2[k]], [by2[k]], color=RED, radius=4.2)
    c.scatter([bx3[k]], [by3[k]], color=TEAL, radius=4.2)
    c.axes("x", "y", "Figure-eight three-body choreography",
           "equal masses, G = 1; dots: bodies at t = 2.0")
    c.legend([("body 1", NAVY, "dot"), ("body 2", RED, "dot"), ("body 3", TEAL, "dot")])
    c.save(os.path.join(FIGS, "three_body_figure_eight.png"))
    print("three_body_figure_eight.png")


def double_pendulum():
    _, rows = read_csv("double_pendulum.csv")
    t = col(rows, 0)
    c = Canvas(1300, 620, margins=(100, 40, 65, 80))
    th1, th2 = col(rows, 1), col(rows, 2)
    c.set_limits(min(t), max(t), min(min(th1), min(th2)), max(max(th1), max(th2)), pad=0.04)
    c.line(t, th1, color=NAVY, width=1.2)
    c.line(t, th2, color=ORANGE, width=1.2)
    c.axes("t", "angle (rad)", "Double pendulum: sensitive dependence",
           "m1=m2=l1=l2=1, g=9.81; initial theta1=1, theta2=2")
    c.legend([("theta1", NAVY, "line"), ("theta2", ORANGE, "line")])
    c.save(os.path.join(FIGS, "double_pendulum.png"))
    print("double_pendulum.png")


def lorenz_spectrum():
    _, rows = read_csv("lorenz_spectrum.csv")
    f, p = col(rows, 0), col(rows, 1)
    keep = [(ff, pp) for ff, pp in zip(f, p) if ff <= 5.0 and pp > 0]
    c = Canvas(1200, 600, margins=(100, 40, 65, 80))
    c.set_limits(0, 5, math.log10(min(q[1] for q in keep)) - 0.2,
                 math.log10(max(q[1] for q in keep)) + 0.2, pad=0.02)
    c.line([q[0] for q in keep], [math.log10(q[1]) for q in keep], color=NAVY, width=1.0)
    c.axes("frequency", "log10 power", "Lorenz power spectrum",
           "broadband continuum, no isolated spectral line")
    c.save(os.path.join(FIGS, "lorenz_spectrum.png"))
    print("lorenz_spectrum.png")


def sync_error():
    _, rows = read_csv("lorenz_sync_error.csv")
    t, e = col(rows, 0), col(rows, 1)
    keep = [(tt, max(ee, 1e-16)) for tt, ee in zip(t, e)]
    c = Canvas(1200, 600, margins=(100, 40, 65, 80))
    c.set_limits(0, max(t), -16, 1, pad=0.02)
    c.line([q[0] for q in keep], [math.log10(q[1]) for q in keep], color=NAVY, width=1.3)
    c.horizontal(-6, color=RED, width=1.0, dash=8)
    c.text(21, -5.4, "1e-6", color=RED, size=12)
    c.axes("t", "log10 |error|", "Pecora-Carroll synchronization error",
           "drive x, response (y,z); error reaches machine precision")
    c.save(os.path.join(FIGS, "lorenz_sync_error.png"))
    print("lorenz_sync_error.png")


def feigenbaum():
    _, rows = read_csv("feigenbaum.csv")
    pts = [(r[0], r[2]) for r in rows if len(r) > 2 and r[2] != ""]
    c = Canvas(1000, 620, margins=(100, 40, 65, 80))
    c.set_limits(0, max(p[0] for p in pts) * 1.1, 4.64, 4.72, pad=0.04)
    c.horizontal(4.669201609, color=RED, width=1.2, dash=9)
    c.line([p[0] for p in pts], [p[1] for p in pts], color=NAVY, width=1.8)
    c.scatter([p[0] for p in pts], [p[1] for p in pts], color=TEAL, radius=3.0)
    c.text(20, 4.6715, "delta = 4.669201609", color=RED, size=13)
    c.axes("period", "delta estimate", "Feigenbaum delta from superstable cycles",
           "logistic map; period 2^k, superstable parameters s_k")
    c.save(os.path.join(FIGS, "feigenbaum.png"))
    print("feigenbaum.png")


if __name__ == "__main__":
    logistic_bifurcation()
    lorenz_attractor()
    lorenz_series()
    henon_attractor()
    rossler_attractor()
    chua_attractor()
    lorenz_poincare()
    lorenz_gp()
    lorenz_rosenstein()
    lorenz_rqa()
    three_body()
    double_pendulum()
    lorenz_spectrum()
    sync_error()
    feigenbaum()
