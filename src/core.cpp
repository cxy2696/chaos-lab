// SPDX-License-Identifier: MIT
#include "chaos/core.hpp"

#include <cmath>
#include <stdexcept>

namespace chaos {

Matrix numerical_jacobian(const OdeSystem& sys, const State& x, double t, double eps) {
    if (!sys.rhs) throw std::invalid_argument("numerical_jacobian: OdeSystem has no rhs");
    const int n = sys.dim;
    require_dim(x, n, "numerical_jacobian");
    Matrix jac(n, n);
    State xp = x;
    State xm = x;
    for (int j = 0; j < n; ++j) {
        const double h = eps * std::max(1.0, std::fabs(x[j]));
        xp[j] = x[j] + h;
        xm[j] = x[j] - h;
        const State fp = sys.rhs(xp, t);
        const State fm = sys.rhs(xm, t);
        if (static_cast<int>(fp.size()) != n || static_cast<int>(fm.size()) != n) {
            throw std::runtime_error("numerical_jacobian: rhs returned wrong dimension");
        }
        for (int i = 0; i < n; ++i) jac(i, j) = (fp[i] - fm[i]) / (2.0 * h);
        xp[j] = x[j];
        xm[j] = x[j];
    }
    return jac;
}

Matrix numerical_jacobian(const MapSystem& map, const State& x, double eps) {
    if (!map.step) throw std::invalid_argument("numerical_jacobian: MapSystem has no step");
    const int n = map.dim;
    require_dim(x, n, "numerical_jacobian");
    Matrix jac(n, n);
    State xp = x;
    State xm = x;
    for (int j = 0; j < n; ++j) {
        const double h = eps * std::max(1.0, std::fabs(x[j]));
        xp[j] = x[j] + h;
        xm[j] = x[j] - h;
        const State fp = map.step(xp);
        const State fm = map.step(xm);
        if (static_cast<int>(fp.size()) != n || static_cast<int>(fm.size()) != n) {
            throw std::runtime_error("numerical_jacobian: map returned wrong dimension");
        }
        for (int i = 0; i < n; ++i) jac(i, j) = (fp[i] - fm[i]) / (2.0 * h);
        xp[j] = x[j];
        xm[j] = x[j];
    }
    return jac;
}

std::vector<double> solve_linear(Matrix A, std::vector<double> b) {
    const std::size_t n = A.rows;
    if (A.cols != n || b.size() != n) {
        throw std::invalid_argument("solve_linear: dimension mismatch");
    }
    for (std::size_t k = 0; k < n; ++k) {
        std::size_t pivot = k;
        double best = std::fabs(A(k, k));
        for (std::size_t i = k + 1; i < n; ++i) {
            const double v = std::fabs(A(i, k));
            if (v > best) {
                best = v;
                pivot = i;
            }
        }
        if (best < 1e-15) throw std::runtime_error("solve_linear: singular matrix");
        if (pivot != k) {
            for (std::size_t j = 0; j < n; ++j) std::swap(A(k, j), A(pivot, j));
            std::swap(b[k], b[pivot]);
        }
        const double d = A(k, k);
        for (std::size_t i = k + 1; i < n; ++i) {
            const double f = A(i, k) / d;
            if (f == 0.0) continue;
            for (std::size_t j = k; j < n; ++j) A(i, j) -= f * A(k, j);
            b[i] -= f * b[k];
        }
    }
    std::vector<double> x(n);
    for (std::size_t ii = n; ii-- > 0;) {
        double s = b[ii];
        for (std::size_t j = ii + 1; j < n; ++j) s -= A(ii, j) * x[j];
        x[ii] = s / A(ii, ii);
    }
    return x;
}

double relative_error(double a, double b) {
    const double denom = std::max({std::fabs(a), std::fabs(b), 1e-300});
    return std::fabs(a - b) / denom;
}

}  // namespace chaos
