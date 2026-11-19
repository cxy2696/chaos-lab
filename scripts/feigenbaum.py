#!/usr/bin/env python3
"""Estimate Feigenbaum's delta from superstable cycles of the logistic map.

A period-p cycle is superstable when x=1/2 lies on it, i.e. f^p(1/2)=1/2.
The sequence of superstable parameters converges to r_infinity with the same
universal ratio delta = 4.669201609... as the period-doubling points.

The intervals below isolate successive period-2^k windows; they are wide
enough to be located by a coarse scan but narrow enough to exclude the
previous window.  The root is then refined by bisection.

Output: data/feigenbaum.csv and a summary on stdout.
"""
import csv
import os

def f(x, r):
    return r * x * (1.0 - x)

def iterate(x, r, n):
    for _ in range(n):
        x = f(x, r)
    return x

def superstable_residual(r, period):
    return iterate(0.5, r, period) - 0.5

def find_superstable(lo, hi, period, iterations=120):
    flo = superstable_residual(lo, period)
    fhi = superstable_residual(hi, period)
    if flo == 0.0:
        return lo
    if fhi == 0.0:
        return hi
    if flo * fhi > 0.0:
        raise RuntimeError(f"bracket does not straddle a root for period {period}: {lo},{hi}")
    for _ in range(iterations):
        mid = 0.5 * (lo + hi)
        fm = superstable_residual(mid, period)
        if flo * fm <= 0.0:
            hi, fhi = mid, fm
        else:
            lo, flo = mid, fm
    return 0.5 * (lo + hi)

def main():
    # Period 1, 2, 4, ..., 64.  The intervals shrink geometrically.
    brackets = [
        (1.9, 2.1),
        (3.20, 3.28),
        (3.49, 3.51),
        (3.55, 3.56),
        (3.566, 3.5672),
        (3.5691, 3.5694),
        (3.56975, 3.56988),
        (3.569905, 3.569925),
    ]
    rows = []
    for k, (lo, hi) in enumerate(brackets):
        period = 2 ** k
        r = find_superstable(lo, hi, period)
        rows.append((period, r))
    print("period  r_superstable        delta_estimate")
    deltas = []
    for i, (p, r) in enumerate(rows):
        if i >= 2:
            d = (rows[i - 1][1] - rows[i - 2][1]) / (rows[i][1] - rows[i - 1][1])
            deltas.append((p, d))
            print(f"{p:6d}  {r:.15f}  {d:.12f}")
        else:
            print(f"{p:6d}  {r:.15f}")
    os.makedirs("data", exist_ok=True)
    with open("data/feigenbaum.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["period", "r_superstable", "delta_estimate"])
        for i, (p, r) in enumerate(rows):
            d = ""
            if i >= 2:
                d = (rows[i - 1][1] - rows[i - 2][1]) / (rows[i][1] - rows[i - 1][1])
            w.writerow([p, f"{r:.15f}", f"{d:.12f}" if d != "" else ""])
    if deltas:
        print(f"\nlast delta estimate: {deltas[-1][1]:.9f}  (reference 4.669201609)")

if __name__ == "__main__":
    main()
