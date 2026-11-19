// SPDX-License-Identifier: MIT
#include "chaos/lyapunov.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "chaos/integrators.hpp"

namespace chaos {

namespace {

void set_identity_basis(State& y, int d, int k) {
    for (int i = 0; i < d; ++i) {
        for (int j = 0; j < k; ++j) y[d + i * k + j] = (i == j) ? 1.0 : 0.0;
    }
}

// Modified Gram-Schmidt with one reorthogonalization pass.  `y` holds the
// d x k tangent basis column-major inside the augmented state.  Log stretch
// factors are appended to `logs`.
void reorthonormalize(State& y, int d, int k, std::vector<double>& logs) {
    std::vector<double> v(d);
    for (int j = 0; j < k; ++j) {
        for (int i = 0; i < d; ++i) v[i] = y[d + i * k + j];
        for (int pass = 0; pass < 2; ++pass) {
            for (int q = 0; q < j; ++q) {
                double r = 0.0;
                for (int i = 0; i < d; ++i) r += y[d + i * k + q] * v[i];
                for (int i = 0; i < d; ++i) v[i] -= r * y[d + i * k + q];
            }
        }
        double n = 0.0;
        for (int i = 0; i < d; ++i) n += v[i] * v[i];
        n = std::sqrt(n);
        if (!(n > 1e-300) || !std::isfinite(n)) {
            // Degenerate basis: replace with a coordinate axis orthogonalized
            // against the previous columns.  This is a rare numerical fallback.
            for (int i = 0; i < d; ++i) v[i] = (i == (j % d)) ? 1.0 : 0.0;
            for (int q = 0; q < j; ++q) {
                double r = 0.0;
                for (int i = 0; i < d; ++i) r += y[d + i * k + q] * v[i];
                for (int i = 0; i < d; ++i) v[i] -= r * y[d + i * k + q];
            }
            n = 0.0;
            for (int i = 0; i < d; ++i) n += v[i] * v[i];
            n = std::sqrt(n);
            if (n < 1e-300) n = 1e-300;
        }
        logs.push_back(std::log(n));
        const double inv = 1.0 / n;
        for (int i = 0; i < d; ++i) y[d + i * k + j] = v[i] * inv;
    }
}

State tangent_rhs(const OdeSystem& sys, const State& y, double t, int d, int k) {
    State x(y.begin(), y.begin() + d);
    const State dx = sys.rhs(x, t);
    const Matrix J = sys.jacobian ? sys.jacobian(x, t) : numerical_jacobian(sys, x, t);
    State dy(d + d * k, 0.0);
    for (int i = 0; i < d; ++i) dy[i] = dx[i];
    for (int i = 0; i < d; ++i) {
        for (int j = 0; j < k; ++j) {
            double s = 0.0;
            for (int l = 0; l < d; ++l) s += J(i, l) * y[d + l * k + j];
            dy[d + i * k + j] = s;
        }
    }
    return dy;
}

}  // namespace

LyapunovSpectrum lyapunov_spectrum(const OdeSystem& sys, State x0,
                                   const LyapunovOptions& options) {
    if (!sys.rhs) throw std::invalid_argument("lyapunov_spectrum: OdeSystem has no rhs");
    if (sys.dim <= 0) throw std::invalid_argument("lyapunov_spectrum: bad dimension");
    require_dim(x0, sys.dim, "lyapunov_spectrum");
    if (!(options.h > 0.0)) throw std::invalid_argument("lyapunov_spectrum: h must be > 0");
    if (!(options.duration > 0.0)) {
        throw std::invalid_argument("lyapunov_spectrum: duration must be > 0");
    }

    const int d = sys.dim;
    const int k = (options.n_exponents > 0) ? std::min(options.n_exponents, d) : d;
    const std::size_t reorth = std::max<std::size_t>(1, options.reorthonormalize_every);

    double t = 0.0;
    State x = x0;
    std::size_t ntrans = static_cast<std::size_t>(std::llround(std::max(0.0, options.transient) / options.h));
    for (std::size_t i = 0; i < ntrans; ++i) {
        x = rk4_step(sys, x, t, options.h);
        t += options.h;
        if (!is_finite_state(x)) throw std::runtime_error("lyapunov_spectrum: transient diverged");
    }

    State y(d + d * k, 0.0);
    for (int i = 0; i < d; ++i) y[i] = x[i];
    set_identity_basis(y, d, k);

    const std::size_t nsteps =
        static_cast<std::size_t>(std::llround(options.duration / options.h));
    if (nsteps == 0) {
        throw std::invalid_argument("lyapunov_spectrum: duration is shorter than h");
    }
    std::vector<double> logs;
    logs.reserve(nsteps / reorth * k + k);
    std::size_t since = 0;
    std::size_t reorth_count = 0;
    for (std::size_t s = 0; s < nsteps; ++s) {
        y = rk4_step_t([&](const State& yy, double tt) { return tangent_rhs(sys, yy, tt, d, k); },
                       y, t, options.h);
        t += options.h;
        ++since;
        if (since == reorth && s + 1 < nsteps) {
            reorthonormalize(y, d, k, logs);
            since = 0;
            ++reorth_count;
        }
    }
    reorthonormalize(y, d, k, logs);
    ++reorth_count;

    const double total_time = static_cast<double>(nsteps) * options.h;
    LyapunovSpectrum result;
    result.exponents.assign(k, 0.0);
    const std::size_t nblocks = logs.size() / static_cast<std::size_t>(k);
    for (std::size_t b = 0; b < nblocks; ++b) {
        for (int j = 0; j < k; ++j) {
            result.exponents[j] += logs[b * k + j];
        }
    }
    for (double& v : result.exponents) v /= total_time;
    std::sort(result.exponents.rbegin(), result.exponents.rend());
    result.integration_time = total_time;
    result.steps = nsteps;
    result.reorthonormalizations = reorth_count;
    return result;
}

LyapunovSpectrum lyapunov_spectrum_map(const MapSystem& map, State x0,
                                       std::size_t transient, std::size_t n,
                                       std::size_t reorthonormalize_every,
                                       int n_exponents) {
    if (!map.step) throw std::invalid_argument("lyapunov_spectrum_map: no step function");
    require_dim(x0, map.dim, "lyapunov_spectrum_map");
    const int d = map.dim;
    const int k = (n_exponents > 0) ? std::min(n_exponents, d) : d;
    const std::size_t reorth = std::max<std::size_t>(1, reorthonormalize_every);
    if (n == 0) throw std::invalid_argument("lyapunov_spectrum_map: n must be positive");

    State x = x0;
    for (std::size_t i = 0; i < transient; ++i) x = map.step(x);

    State y(d + d * k, 0.0);
    for (int i = 0; i < d; ++i) y[i] = x[i];
    set_identity_basis(y, d, k);

    std::vector<double> logs;
    logs.reserve((n / reorth + 1) * k);
    std::size_t since = 0;
    std::vector<double> new_q(static_cast<std::size_t>(d) * k);
    for (std::size_t s = 0; s < n; ++s) {
        const Matrix J = map.jacobian ? map.jacobian(x) : numerical_jacobian(map, x);
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < k; ++j) {
                double sum = 0.0;
                for (int l = 0; l < d; ++l) sum += J(i, l) * y[d + l * k + j];
                new_q[static_cast<std::size_t>(i) * k + j] = sum;
            }
        }
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < k; ++j) {
                y[d + i * k + j] = new_q[static_cast<std::size_t>(i) * k + j];
            }
        }
        x = map.step(x);
        for (int i = 0; i < d; ++i) y[i] = x[i];
        ++since;
        if (since == reorth && s + 1 < n) {
            reorthonormalize(y, d, k, logs);
            since = 0;
        }
    }
    reorthonormalize(y, d, k, logs);

    LyapunovSpectrum result;
    result.exponents.assign(k, 0.0);
    const std::size_t nblocks = logs.size() / static_cast<std::size_t>(k);
    for (std::size_t b = 0; b < nblocks; ++b) {
        for (int j = 0; j < k; ++j) result.exponents[j] += logs[b * k + j];
    }
    for (double& v : result.exponents) v /= static_cast<double>(n);
    std::sort(result.exponents.rbegin(), result.exponents.rend());
    result.integration_time = static_cast<double>(n);
    result.steps = n;
    result.reorthonormalizations = nblocks;
    return result;
}

double lyapunov_max_two_trajectories(const OdeSystem& sys, State x0, double d0,
                                     double h, double transient, double duration) {
    if (!(d0 > 0.0) || !(h > 0.0)) {
        throw std::invalid_argument("two-trajectory Lyapunov: d0 and h must be positive");
    }
    if (!sys.rhs || sys.dim <= 0) {
        throw std::invalid_argument("two-trajectory Lyapunov: invalid system");
    }
    require_dim(x0, sys.dim, "two-trajectory Lyapunov");
    double t = 0.0;
    State a = x0;
    const std::size_t ntrans = static_cast<std::size_t>(std::llround(transient / h));
    for (std::size_t i = 0; i < ntrans; ++i) {
        a = rk4_step(sys, a, t, h);
        t += h;
    }
    State b = a;
    b[0] += d0;
    const std::size_t n = static_cast<std::size_t>(std::llround(duration / h));
    if (n == 0) throw std::invalid_argument("two-trajectory Lyapunov: duration < h");
    double sum = 0.0;
    std::size_t completed = 0;
    for (std::size_t i = 0; i < n; ++i) {
        a = rk4_step(sys, a, t, h);
        b = rk4_step(sys, b, t, h);
        t += h;
        const double d = distance(a, b);
        if (!(d > 0.0)) break;
        sum += std::log(d / d0);
        const double scale = d0 / d;
        for (std::size_t j = 0; j < b.size(); ++j) b[j] = a[j] + (b[j] - a[j]) * scale;
        ++completed;
    }
    if (completed == 0) return 0.0;
    return sum / (static_cast<double>(completed) * h);
}

namespace {

void linear_fit(const std::vector<double>& x, const std::vector<double>& y, double& slope,
                double& intercept, double& r2) {
    const std::size_t n = x.size();
    double sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0, syy = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        sx += x[i];
        sy += y[i];
        sxx += x[i] * x[i];
        sxy += x[i] * y[i];
        syy += y[i] * y[i];
    }
    const double denom = static_cast<double>(n) * sxx - sx * sx;
    if (std::fabs(denom) < 1e-300) {
        slope = 0.0;
        intercept = 0.0;
        r2 = 0.0;
        return;
    }
    slope = (static_cast<double>(n) * sxy - sx * sy) / denom;
    intercept = (sy - slope * sx) / static_cast<double>(n);
    const double ss_tot = syy - sy * sy / static_cast<double>(n);
    double ss_res = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double e = y[i] - (intercept + slope * x[i]);
        ss_res += e * e;
    }
    r2 = ss_tot > 0.0 ? 1.0 - ss_res / ss_tot : 1.0;
}

}  // namespace

RosensteinResult rosenstein_lyapunov(const std::vector<double>& series, int m, int tau,
                                     int theiler, int max_k, double dt, int fit_begin,
                                     int fit_end) {
    if (m < 1 || tau < 1 || theiler < 0 || max_k < 1 || !(dt > 0.0)) {
        throw std::invalid_argument("rosenstein_lyapunov: invalid m/tau/theiler/max_k/dt");
    }
    const int n = static_cast<int>(series.size());
    const int M = n - (m - 1) * tau;
    if (M <= max_k + 2 || M < 20) {
        throw std::invalid_argument("rosenstein_lyapunov: series too short");
    }
    const int dim = m;
    std::vector<double> emb(static_cast<std::size_t>(M) * dim);
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < dim; ++j) {
            emb[static_cast<std::size_t>(i) * dim + j] = series[i + j * tau];
        }
    }

    std::vector<double> mean_log(static_cast<std::size_t>(max_k) + 1, 0.0);
    std::vector<std::size_t> count(static_cast<std::size_t>(max_k) + 1, 0);
    for (int i = 0; i < M - max_k; ++i) {
        int best = -1;
        double best_d2 = std::numeric_limits<double>::max();
        for (int j = 0; j < M - max_k; ++j) {
            if (std::abs(i - j) <= theiler) continue;
            double d2 = 0.0;
            for (int q = 0; q < dim; ++q) {
                const double diff = emb[static_cast<std::size_t>(i) * dim + q] -
                                    emb[static_cast<std::size_t>(j) * dim + q];
                d2 += diff * diff;
            }
            if (d2 < best_d2) {
                best_d2 = d2;
                best = j;
            }
        }
        if (best < 0 || best_d2 <= 0.0) continue;
        for (int kk = 0; kk <= max_k && i + kk < M && best + kk < M; ++kk) {
            double d2 = 0.0;
            for (int q = 0; q < dim; ++q) {
                const double diff = emb[static_cast<std::size_t>(i + kk) * dim + q] -
                                    emb[static_cast<std::size_t>(best + kk) * dim + q];
                d2 += diff * diff;
            }
            mean_log[static_cast<std::size_t>(kk)] +=
                0.5 * std::log(std::max(d2, 1e-300));
            ++count[static_cast<std::size_t>(kk)];
        }
    }

    RosensteinResult result;
    result.dt = dt;
    for (int kk = 0; kk <= max_k; ++kk) {
        if (count[static_cast<std::size_t>(kk)] == 0) continue;
        result.k.push_back(static_cast<double>(kk));
        result.divergence.push_back(
            mean_log[static_cast<std::size_t>(kk)] /
            static_cast<double>(count[static_cast<std::size_t>(kk)]));
    }
    if (result.k.size() < 3) throw std::runtime_error("rosenstein_lyapunov: no neighbours");

    const std::size_t K = result.k.size();
    std::size_t a = 0;
    std::size_t b = K - 1;
    if (fit_begin >= 0 && fit_end >= 0 &&
        static_cast<std::size_t>(fit_end) < K &&
        fit_end > fit_begin) {
        a = static_cast<std::size_t>(fit_begin);
        b = static_cast<std::size_t>(fit_end);
    } else if (K >= 10) {
        // Choose a long, nearly linear, positive-slope window.  The score
        // rewards both fit quality and window length, which separates the
        // asymptotic growth regime from the initial transient and saturation.
        std::vector<double> adj;
        adj.reserve(K - 1);
        for (std::size_t i = 0; i + 1 < K; ++i) {
            const double dx = result.k[i + 1] - result.k[i];
            if (dx > 0.0) adj.push_back((result.divergence[i + 1] - result.divergence[i]) / dx);
        }
        std::sort(adj.begin(), adj.end());
        const double smax = adj.empty() ? 0.0 : adj[adj.size() * 9 / 10];
        const double cutoff = 0.2 * smax;
        double best_score = -1.0;
        const std::size_t min_len = std::max<std::size_t>(8, K / 20);
        std::vector<std::size_t> lengths;
        for (std::size_t L : {min_len, K / 8, K / 5, K / 3, K / 2}) {
            if (L >= min_len && L < K) lengths.push_back(L);
        }
        std::sort(lengths.begin(), lengths.end());
        lengths.erase(std::unique(lengths.begin(), lengths.end()), lengths.end());
        for (std::size_t L : lengths) {
            const std::size_t step = std::max<std::size_t>(1, L / 4);
            for (std::size_t start = 0; start + L <= K; start += step) {
                double sl = 0.0, ic = 0.0, r2 = 0.0;
                linear_fit(std::vector<double>(result.k.begin() + static_cast<std::ptrdiff_t>(start),
                                               result.k.begin() + static_cast<std::ptrdiff_t>(start + L)),
                           std::vector<double>(result.divergence.begin() + static_cast<std::ptrdiff_t>(start),
                                               result.divergence.begin() + static_cast<std::ptrdiff_t>(start + L)),
                           sl, ic, r2);
                if (sl <= cutoff || r2 <= 0.0) continue;
                const double score = r2 * std::sqrt(static_cast<double>(L) / static_cast<double>(K));
                if (score > best_score) {
                    best_score = score;
                    a = start;
                    b = start + L - 1;
                }
            }
        }
        if (best_score < 0.0) {
            a = 0;
            b = K - 1;
        }
    }

    result.fit_begin = a;
    result.fit_end = b;
    double slope = 0.0;
    linear_fit(std::vector<double>(result.k.begin() + static_cast<std::ptrdiff_t>(a),
                                   result.k.begin() + static_cast<std::ptrdiff_t>(b + 1)),
               std::vector<double>(result.divergence.begin() + static_cast<std::ptrdiff_t>(a),
                                   result.divergence.begin() + static_cast<std::ptrdiff_t>(b + 1)),
               slope, result.intercept, result.r2);
    // Consecutive embedding vectors are separated by one sample, dt, so the
    // divergence slope is a per-sample exponent.  Note: do NOT divide by tau.
    result.lambda = slope / dt;
    return result;
}

}  // namespace chaos
