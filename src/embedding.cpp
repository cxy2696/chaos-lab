// SPDX-License-Identifier: MIT
#include "chaos/embedding.hpp"

#include "chaos/core.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace chaos {

Embedding delay_embed(const std::vector<double>& series, int m, int tau) {
    if (m < 1 || tau < 1) throw std::invalid_argument("delay_embed: m,tau must be >= 1");
    const int n = static_cast<int>(series.size());
    const int rows = n - (m - 1) * tau;
    if (rows < 2) throw std::invalid_argument("delay_embed: series too short");
    Embedding emb;
    emb.m = m;
    emb.tau = tau;
    emb.n = static_cast<std::size_t>(rows);
    emb.data.resize(static_cast<std::size_t>(rows) * m);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < m; ++j) {
            emb.data[static_cast<std::size_t>(i) * m + j] = series[i + j * tau];
        }
    }
    return emb;
}

std::vector<double> average_mutual_information(const std::vector<double>& series,
                                               int max_tau, int bins) {
    if (series.size() < 4) throw std::invalid_argument("AMI: series too short");
    if (max_tau < 1) return {};
    bins = std::max(4, bins);
    double lo = *std::min_element(series.begin(), series.end());
    double hi = *std::max_element(series.begin(), series.end());
    if (!(hi > lo)) throw std::invalid_argument("AMI: constant series");
    const double width = (hi - lo) / static_cast<double>(bins);

    auto bin_of = [&](double v) {
        int b = static_cast<int>((v - lo) / width);
        return std::clamp(b, 0, bins - 1);
    };

    std::vector<double> out(static_cast<std::size_t>(max_tau) + 1, 0.0);
    std::vector<double> px(static_cast<std::size_t>(bins));
    std::vector<double> py(static_cast<std::size_t>(bins));
    std::vector<double> pxy(static_cast<std::size_t>(bins) * bins);

    for (int tau = 1; tau <= max_tau; ++tau) {
        const std::size_t n = series.size() - static_cast<std::size_t>(tau);
        std::fill(px.begin(), px.end(), 0.0);
        std::fill(py.begin(), py.end(), 0.0);
        std::fill(pxy.begin(), pxy.end(), 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            const int a = bin_of(series[i]);
            const int b = bin_of(series[i + static_cast<std::size_t>(tau)]);
            px[static_cast<std::size_t>(a)] += 1.0;
            py[static_cast<std::size_t>(b)] += 1.0;
            pxy[static_cast<std::size_t>(a) * bins + static_cast<std::size_t>(b)] += 1.0;
        }
        const double inv = 1.0 / static_cast<double>(n);
        for (double& v : px) v *= inv;
        for (double& v : py) v *= inv;
        for (double& v : pxy) v *= inv;
        double mi = 0.0;
        for (int a = 0; a < bins; ++a) {
            if (px[static_cast<std::size_t>(a)] <= 0.0) continue;
            for (int b = 0; b < bins; ++b) {
                const double p = pxy[static_cast<std::size_t>(a) * bins + b];
                if (p <= 0.0 || py[static_cast<std::size_t>(b)] <= 0.0) continue;
                mi += p * std::log(p / (px[static_cast<std::size_t>(a)] *
                                        py[static_cast<std::size_t>(b)]));
            }
        }
        out[static_cast<std::size_t>(tau)] = mi;
    }
    // tau = 0: entropy of the marginal distribution.
    std::vector<double> p(static_cast<std::size_t>(bins), 0.0);
    for (double v : series) p[static_cast<std::size_t>(bin_of(v))] += 1.0;
    const double inv = 1.0 / static_cast<double>(series.size());
    double h = 0.0;
    for (double q : p) {
        q *= inv;
        if (q > 0.0) h -= q * std::log(q);
    }
    out[0] = h;
    return out;
}

int first_ami_minimum(const std::vector<double>& ami) {
    for (std::size_t tau = 1; tau + 1 < ami.size(); ++tau) {
        if (ami[tau] < ami[tau - 1] && ami[tau] <= ami[tau + 1]) {
            return static_cast<int>(tau);
        }
    }
    return -1;
}

std::vector<double> autocorrelation(const std::vector<double>& series, int max_lag) {
    std::vector<double> out;
    if (series.empty() || max_lag < 0) return out;
    const double mu = mean(series);
    double var = 0.0;
    for (double v : series) var += (v - mu) * (v - mu);
    if (!(var > 0.0)) throw std::invalid_argument("autocorrelation: constant series");
    const int n = static_cast<int>(series.size());
    max_lag = std::min(max_lag, n - 1);
    out.resize(static_cast<std::size_t>(max_lag) + 1);
    for (int lag = 0; lag <= max_lag; ++lag) {
        double s = 0.0;
        for (int i = 0; i + lag < n; ++i) s += (series[i] - mu) * (series[i + lag] - mu);
        out[static_cast<std::size_t>(lag)] = s / var;
    }
    return out;
}

int autocorrelation_1e(const std::vector<double>& ac) {
    constexpr double target = 1.0 / 2.71828182845904523536;
    for (std::size_t lag = 1; lag < ac.size(); ++lag) {
        if (std::fabs(ac[lag]) < target) return static_cast<int>(lag);
    }
    return -1;
}

FnnResult false_nearest_neighbors(const std::vector<double>& series, int max_m, int tau,
                                  int theiler, double rtol, double atol) {
    if (max_m < 1 || tau < 1) throw std::invalid_argument("FNN: invalid m/tau");
    const int n = static_cast<int>(series.size());
    FnnResult result;
    result.false_ratio.assign(static_cast<std::size_t>(max_m) + 1, 0.0);
    result.valid_points.assign(static_cast<std::size_t>(max_m) + 1, 0);

    for (int m = 1; m < max_m; ++m) {
        const int M = n - (m - 1) * tau;
        // Need one extra coordinate for the (m+1)-th component.
        const int usable = n - m * tau;
        if (usable <= 1) break;
        std::size_t false_count = 0;
        std::size_t total = 0;
        for (int i = 0; i < usable; ++i) {
            int best = -1;
            double best_d2 = std::numeric_limits<double>::max();
            for (int j = 0; j < usable; ++j) {
                if (std::abs(i - j) <= theiler) continue;
                double d2 = 0.0;
                for (int q = 0; q < m; ++q) {
                    const double diff = series[i + q * tau] - series[j + q * tau];
                    d2 += diff * diff;
                }
                if (d2 < best_d2) {
                    best_d2 = d2;
                    best = j;
                }
            }
            if (best < 0) continue;
            ++total;
            const double d = std::sqrt(best_d2);
            const double extra =
                std::fabs(series[i + m * tau] - series[best + m * tau]);
            bool is_false = false;
            if (d <= atol) {
                is_false = extra > 0.0;
            } else {
                is_false = extra / d > rtol;
            }
            if (is_false) ++false_count;
            (void)M;
        }
        result.false_ratio[static_cast<std::size_t>(m)] =
            total ? static_cast<double>(false_count) / static_cast<double>(total) : 1.0;
        result.valid_points[static_cast<std::size_t>(m)] = static_cast<int>(total);
    }

    result.suggested_m = max_m;
    for (int m = 1; m <= max_m; ++m) {
        if (result.valid_points[static_cast<std::size_t>(m)] > 0 &&
            result.false_ratio[static_cast<std::size_t>(m)] < 0.01) {
            result.suggested_m = m;
            break;
        }
    }
    return result;
}

namespace {

double quantile(std::vector<double>& v, double q) {
    if (v.empty()) return 0.0;
    const std::size_t idx = std::min(v.size() - 1,
                                     static_cast<std::size_t>(q * static_cast<double>(v.size())));
    std::nth_element(v.begin(), v.begin() + static_cast<std::ptrdiff_t>(idx), v.end());
    return v[idx];
}

}  // namespace

CorrelationResult correlation_dimension(const Embedding& emb, double r_min, double r_max,
                                        int n_r, int theiler, double fit_lo, double fit_hi) {
    if (emb.n < 10 || emb.m < 1) throw std::invalid_argument("correlation_dimension: bad embedding");
    if (n_r < 3) throw std::invalid_argument("correlation_dimension: n_r must be >= 3");
    if (fit_lo < 0.0 || fit_hi > 1.0 || fit_lo >= fit_hi) {
        throw std::invalid_argument("correlation_dimension: require 0 <= fit_lo < fit_hi <= 1");
    }
    const std::size_t N = emb.n;
    const int m = emb.m;

    // Collect distances for the automatic radius range.  O(N^2) memory is
    // acceptable for the intended teaching-scale data sets; explicit radii
    // avoid storing the whole distance list.
    std::vector<double> sample;
    if (!(r_min > 0.0) || !(r_max > 0.0)) {
        sample.reserve(N * (N - 1) / 2);
        for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t j = i + 1; j < N; ++j) {
                if (std::abs(static_cast<long long>(i) - static_cast<long long>(j)) <= theiler) {
                    continue;
                }
                double d2 = 0.0;
                for (int q = 0; q < m; ++q) {
                    const double diff = emb.at(i, q) - emb.at(j, q);
                    d2 += diff * diff;
                }
                sample.push_back(std::sqrt(d2));
            }
        }
        if (sample.empty()) throw std::runtime_error("correlation_dimension: no pairs");
        // Correlation sums between ~0.005% and 0.3% of all pairs usually sit
        // inside the scaling region for low-dimensional attractors; the upper
        // end avoids the saturation that dominates at large radii.
        const double q_lo = std::max(5e-5, 100.0 / static_cast<double>(sample.size()));
        if (!(r_max > 0.0)) r_max = quantile(sample, 0.003);
        if (!(r_min > 0.0)) r_min = quantile(sample, q_lo);
        if (!(r_min > 0.0) || !(r_max > r_min)) {
            throw std::runtime_error("correlation_dimension: degenerate radius range");
        }
    }

    CorrelationResult out;
    out.r_min = r_min;
    out.r_max = r_max;
    out.radii = linspace(r_min, r_max, static_cast<std::size_t>(n_r));
    out.log_r.resize(static_cast<std::size_t>(n_r));
    out.log_c.resize(static_cast<std::size_t>(n_r), -std::numeric_limits<double>::infinity());

    std::vector<double> counts(static_cast<std::size_t>(n_r), 0.0);
    double pairs = 0.0;
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = i + 1; j < N; ++j) {
            if (std::abs(static_cast<long long>(i) - static_cast<long long>(j)) <= theiler) {
                continue;
            }
            double d2 = 0.0;
            for (int q = 0; q < m; ++q) {
                const double diff = emb.at(i, q) - emb.at(j, q);
                d2 += diff * diff;
            }
            const double d = std::sqrt(d2);
            ++pairs;
            for (int k = 0; k < n_r; ++k) {
                if (d <= out.radii[static_cast<std::size_t>(k)]) {
                    counts[static_cast<std::size_t>(k)] += 1.0;
                }
            }
        }
    }
    if (pairs <= 0.0) throw std::runtime_error("correlation_dimension: no pairs");

    for (int k = 0; k < n_r; ++k) {
        const double c = counts[static_cast<std::size_t>(k)] / pairs;
        out.log_r[static_cast<std::size_t>(k)] = std::log(out.radii[static_cast<std::size_t>(k)]);
        if (c > 0.0) out.log_c[static_cast<std::size_t>(k)] = std::log(c);
    }
    out.local_slope.assign(static_cast<std::size_t>(n_r),
                           std::numeric_limits<double>::quiet_NaN());
    for (int k = 0; k + 1 < n_r; ++k) {
        const double dy = out.log_c[static_cast<std::size_t>(k + 1)] -
                          out.log_c[static_cast<std::size_t>(k)];
        const double dx = out.log_r[static_cast<std::size_t>(k + 1)] -
                          out.log_r[static_cast<std::size_t>(k)];
        if (std::isfinite(dy) && dx > 0.0) {
            out.local_slope[static_cast<std::size_t>(k)] = dy / dx;
        }
    }

    const std::size_t begin = static_cast<std::size_t>(
        std::clamp(fit_lo, 0.0, 0.9) * static_cast<double>(n_r - 1));
    const std::size_t end = static_cast<std::size_t>(
        std::clamp(fit_hi, fit_lo + 0.05, 1.0) * static_cast<double>(n_r - 1));
    out.fit_begin = begin;
    out.fit_end = end;

    // Linear regression on log C versus log r.
    double sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0, syy = 0.0;
    std::size_t nfit = 0;
    for (std::size_t k = begin; k <= end && k < out.log_r.size(); ++k) {
        if (!std::isfinite(out.log_c[k])) continue;
        sx += out.log_r[k];
        sy += out.log_c[k];
        sxx += out.log_r[k] * out.log_r[k];
        sxy += out.log_r[k] * out.log_c[k];
        syy += out.log_c[k] * out.log_c[k];
        ++nfit;
    }
    if (nfit >= 3) {
        const double denom = static_cast<double>(nfit) * sxx - sx * sx;
        if (std::fabs(denom) > 1e-300) {
            out.d2 = (static_cast<double>(nfit) * sxy - sx * sy) / denom;
            const double intercept = (sy - out.d2 * sx) / static_cast<double>(nfit);
            const double ss_tot = syy - sy * sy / static_cast<double>(nfit);
            double ss_res = 0.0;
            for (std::size_t k = begin; k <= end && k < out.log_r.size(); ++k) {
                if (!std::isfinite(out.log_c[k])) continue;
                const double e = out.log_c[k] - (intercept + out.d2 * out.log_r[k]);
                ss_res += e * e;
            }
            out.r2 = ss_tot > 0.0 ? 1.0 - ss_res / ss_tot : 1.0;
        }
    }
    return out;
}

CorrelationResult correlation_dimension(const std::vector<double>& series, int m, int tau,
                                        double r_min, double r_max, int n_r, int theiler,
                                        double fit_lo, double fit_hi) {
    return correlation_dimension(delay_embed(series, m, tau), r_min, r_max, n_r, theiler,
                                 fit_lo, fit_hi);
}

}  // namespace chaos
