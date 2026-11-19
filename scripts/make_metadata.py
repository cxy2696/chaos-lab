#!/usr/bin/env python3
"""Write data/checksums.sha256 and data/metrics.json from committed data.

The numerical values are the outputs of the commands recorded in
docs/data-provenance.md.  Re-run after regenerating data.
"""
import csv
import hashlib
import json
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")

FILES = [
    "lorenz_trajectory.csv",
    "lorenz_series.csv",
    "rossler_trajectory.csv",
    "chua_trajectory.csv",
    "three_body_figure_eight.csv",
    "double_pendulum.csv",
    "henon_orbit.csv",
    "logistic_orbit.csv",
    "logistic_bifurcation.csv",
    "logistic_lyapunov.csv",
    "lorenz_poincare.csv",
    "lorenz_gp.csv",
    "lorenz_rosenstein.csv",
    "lorenz_spectrum.csv",
    "lorenz_rqa_matrix.csv",
    "lorenz_zeroone.csv",
    "lorenz_sindy.csv",
    "lorenz_sync_error.csv",
    "feigenbaum.csv",
]

METRICS = {
    "schema": 1,
    "library_version": "0.1.0",
    "generated": "2026-11-19",
    "notes": "Values measured with the commands listed in docs/data-provenance.md; they are not literature values.",
    "measurements": {
        "logistic_lambda_r4": 0.693149,
        "feigenbaum_delta_last": 4.669157182,
        "feigenbaum_r_infinity_reference": 3.569945672,
        "lorenz_lyapunov_spectrum": [0.910064, -0.004271, -14.572357],
        "lorenz_lyapunov_sum": -13.666564,
        "lorenz_divergence": -13.666667,
        "lorenz_correlation_dimension": 2.04662,
        "lorenz_gp_fit_r2": 0.999931,
        "lorenz_gp_radius_range": [0.105461, 0.820631],
        "lorenz_rosenstein_lambda": 0.855017,
        "lorenz_zero_one_k_stride10": 0.990549,
        "lorenz_rqa": {
            "epsilon": 3.20195,
            "rr": 0.0296232,
            "det": 0.691293,
            "lam": 0.293427,
            "entropy_nats": 2.03944,
            "trapping_time": 2.04445,
            "longest_diagonal": 73,
        },
        "henon_lyapunov_spectrum": [0.420636, -1.624609],
        "henon_lyapunov_sum": -1.203973,
        "henon_log_b": -1.203973,
        "rossler_lyapunov_spectrum": [0.079038, 0.000406, -5.417134],
        "chua_lyapunov_spectrum": [0.424696, -0.007094, -4.217896],
        "three_body_period_return_distance": 2.026e-4,
        "synchronization": {
            "final_error": 0.0,
            "max_error_after_transient": 1.25347e-4,
            "settling_time": 7.58,
        },
        "sindy_terms": 10,
        "sindy_active_terms": 7,
    },
    "data_files": [],
}


def main():
    with open(os.path.join(DATA, "checksums.sha256"), "w") as out:
        for name in FILES:
            path = os.path.join(DATA, name)
            if not os.path.exists(path):
                raise SystemExit(f"missing data file: {name}")
            h = hashlib.sha256()
            with open(path, "rb") as fh:
                for chunk in iter(lambda: fh.read(1 << 20), b""):
                    h.update(chunk)
            digest = h.hexdigest()
            with open(path, newline="") as fh:
                # Only the recurrence matrix is headerless.
                skip = 0 if name.endswith("rqa_matrix.csv") else 1
                rows = max(0, sum(1 for _ in fh) - skip)
            out.write(f"{digest}  data/{name}\n")
            METRICS["data_files"].append({"name": name, "rows": rows, "sha256": digest})
    with open(os.path.join(DATA, "metrics.json"), "w") as out:
        json.dump(METRICS, out, indent=2)
        out.write("\n")
    print(f"wrote checksums for {len(FILES)} files and data/metrics.json")


if __name__ == "__main__":
    main()
