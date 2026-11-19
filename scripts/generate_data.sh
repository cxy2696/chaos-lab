#!/usr/bin/env bash
# Regenerate every committed data file from the current library build.
# Usage: ./scripts/generate_data.sh [path-to-chaos_cli]
set -euo pipefail
cd "$(dirname "$0")/.."
CLI="${1:-./build/chaos_cli}"
if [[ ! -x "$CLI" ]]; then
  echo "chaos_cli not found at $CLI; build first: cmake -S . -B build && cmake --build build -j" >&2
  exit 1
fi

mkdir -p data

echo "[1/15] Lorenz trajectory"
"$CLI" simulate --system lorenz --n 6000 --h 0.01 --transient 5 --out data/lorenz_trajectory.csv

echo "[2/15] Lorenz scalar series"
"$CLI" simulate --system lorenz --n 8000 --h 0.01 --transient 50 --out data/lorenz_series.csv

echo "[3/15] Rossler trajectory"
"$CLI" simulate --system rossler --n 8000 --h 0.01 --transient 20 --out data/rossler_trajectory.csv

echo "[4/15] Chua trajectory"
"$CLI" simulate --system chua --n 10000 --h 0.01 --transient 20 --out data/chua_trajectory.csv

echo "[5/15] Figure-eight three-body orbit"
"$CLI" simulate --system three_body --n 7000 --h 0.001 --out data/three_body_figure_eight.csv

echo "[6/15] Double pendulum"
"$CLI" simulate --system double_pendulum --n 5000 --h 0.005 --out data/double_pendulum.csv

echo "[7/15] Henon orbit"
"$CLI" simulate --map henon --n 5000 --transient 1000 --out data/henon_orbit.csv

echo "[8/15] Logistic orbit"
"$CLI" simulate --map logistic --n 2000 --transient 1000 --out data/logistic_orbit.csv

echo "[9/15] Logistic bifurcation diagram"
"$CLI" bifurcation --map logistic --pmin 2.5 --pmax 4.0 --nparams 1500 --transient 1000 \
  --keep 80 --stride 1 --out data/logistic_bifurcation.csv

echo "[10/15] Logistic Lyapunov scan"
"$CLI" lyapunov-scan --map logistic --pmin 2.5 --pmax 4.0 --nparams 1500 \
  --transient 1000 --n 20000 --out data/logistic_lyapunov.csv

echo "[11/15] Lorenz Poincare section (z = 27, downward)"
"$CLI" poincare --system lorenz --coord 2 --value 27 --direction -1 --h 0.005 \
  --transient 20 --t 200 --out data/lorenz_poincare.csv

echo "[12/15] Lorenz correlation dimension"
"$CLI" dimension --input data/lorenz_series.csv --column 1 --m 3 --tau 17 \
  --theiler 50 --nr 60 --out data/lorenz_gp.csv

echo "[13/15] Lorenz Rosenstein curve and spectrum"
"$CLI" rosenstein --input data/lorenz_series.csv --column 1 --m 3 --tau 17 \
  --theiler 100 --max-k 300 --dt 0.01 --out data/lorenz_rosenstein.csv
"$CLI" spectrum --input data/lorenz_series.csv --column 1 --dt 0.01 \
  --out data/lorenz_spectrum.csv

echo "[14/15] Lorenz recurrence matrix and 0-1 test"
"$CLI" rqa --input data/lorenz_series.csv --column 1 --stride 10 --m 3 --tau 3 \
  --theiler 5 --rate 0.03 --matrix-out data/lorenz_rqa_matrix.csv
"$CLI" zeroone --input data/lorenz_series.csv --column 1 --stride 10 \
  --out data/lorenz_zeroone.csv

echo "[15/15] Lorenz SINDy coefficients and synchronization error"
"$CLI" sindy --system lorenz --n 5000 --transient 5000 --h 0.002 --degree 2 \
  --threshold 0.1 --out data/lorenz_sindy.csv
"$CLI" sync --t 40 --h 0.005 --out data/lorenz_sync_error.csv

echo "Data written under data/"
