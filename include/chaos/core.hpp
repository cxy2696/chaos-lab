// chaos-lab: shared types and small linear-algebra helpers.
// SPDX-License-Identifier: MIT
#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace chaos {

inline constexpr double kPi = 3.141592653589793238462643383279502884;
inline constexpr double kTwoPi = 2.0 * kPi;

using State = std::vector<double>;

// Row-major dense matrix.  Small dimensions only; this is a teaching library,
// not a replacement for BLAS/LAPACK.
struct Matrix {
    std::size_t rows = 0;
    std::size_t cols = 0;
    std::vector<double> data;

    Matrix() = default;
    Matrix(std::size_t r, std::size_t c, double value = 0.0)
        : rows(r), cols(c), data(r * c, value) {}

    static Matrix identity(std::size_t n) {
        Matrix m(n, n);
        for (std::size_t i = 0; i < n; ++i) m(i, i) = 1.0;
        return m;
    }

    double& operator()(std::size_t i, std::size_t j) { return data[i * cols + j]; }
    double operator()(std::size_t i, std::size_t j) const { return data[i * cols + j]; }
    double* row(std::size_t i) { return data.data() + i * cols; }
    const double* row(std::size_t i) const { return data.data() + i * cols; }
};

// A continuous-time autonomous or non-autonomous ODE dx/dt = f(x, t).
// `jacobian` is optional; the Lyapunov routines fall back to finite differences.
struct OdeSystem {
    std::string name;
    std::string description;
    int dim = 0;
    std::vector<std::string> param_names;
    std::vector<double> params;
    std::function<State(const State&, double)> rhs;
    std::function<Matrix(const State&, double)> jacobian;
};

// A discrete map x_{n+1} = f(x_n).
struct MapSystem {
    std::string name;
    std::string description;
    int dim = 0;
    std::vector<std::string> param_names;
    std::vector<double> params;
    std::function<State(const State&)> step;
    std::function<Matrix(const State&)> jacobian;
};

inline double dot(const State& a, const State& b) {
    double s = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
    return s;
}

inline double norm2(const State& a) { return std::sqrt(dot(a, a)); }

inline double distance(const State& a, const State& b) {
    double s = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        const double d = a[i] - b[i];
        s += d * d;
    }
    return std::sqrt(s);
}

inline State axpy(const State& x, const State& y, double a) {
    State r(x.size());
    for (std::size_t i = 0; i < x.size(); ++i) r[i] = x[i] + a * y[i];
    return r;
}

inline State scale(const State& x, double a) {
    State r(x.size());
    for (std::size_t i = 0; i < x.size(); ++i) r[i] = a * x[i];
    return r;
}

inline double mean(const std::vector<double>& x) {
    if (x.empty()) return 0.0;
    double s = 0.0;
    for (double v : x) s += v;
    return s / static_cast<double>(x.size());
}

inline std::vector<double> linspace(double a, double b, std::size_t n) {
    std::vector<double> out;
    if (n == 0) return out;
    out.reserve(n);
    if (n == 1) {
        out.push_back(a);
        return out;
    }
    const double h = (b - a) / static_cast<double>(n - 1);
    for (std::size_t i = 0; i < n; ++i) out.push_back(a + h * static_cast<double>(i));
    return out;
}

inline void require_dim(const State& x, int dim, const char* where) {
    if (static_cast<int>(x.size()) != dim) {
        throw std::invalid_argument(std::string(where) + ": expected state dimension " +
                                    std::to_string(dim) + ", got " +
                                    std::to_string(x.size()));
    }
}

// Central finite-difference Jacobian of an ODE right-hand side.
Matrix numerical_jacobian(const OdeSystem& sys, const State& x, double t,
                          double eps = 1e-7);
Matrix numerical_jacobian(const MapSystem& map, const State& x, double eps = 1e-7);

// Small dense solvers used by SINDy / least-squares code.
std::vector<double> solve_linear(Matrix A, std::vector<double> b);
double relative_error(double a, double b);

inline bool is_finite_state(const State& x) {
    for (double v : x) {
        if (!std::isfinite(v)) return false;
    }
    return true;
}

}  // namespace chaos
