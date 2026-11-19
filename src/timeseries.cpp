// SPDX-License-Identifier: MIT
#include "chaos/timeseries.hpp"

#include "chaos/core.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <numeric>
#include <stdexcept>

namespace chaos {

namespace {

using Complex = std::complex<double>;

void fft(std::vector<Complex>& a, bool inverse) {
    const std::size_t n = a.size();
    for (std::size_t i = 1, j = 0; i < n; ++i) {
        std::size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(a[i], a[j]);
    }
    for (std::size_t len = 2; len <= n; len <<= 1) {
        const double angle = (inverse ? 2.0 : -2.0) * kPi / static_cast<double>(len);
        const Complex wlen(std::cos(angle), std::sin(angle));
        for (std::size_t i = 0; i < n; i += len) {
            Complex w(1.0, 0.0);
            for (std::size_t j = 0; j < len / 2; ++j) {
                const Complex u = a[i + j];
                const Complex v = a[i + j + len / 2] * w;
                a[i + j] = u + v;
                a[i + j + len / 2] = u - v;
                w *= wlen;
            }
        }
    }
    if (inverse) {
        for (Complex& v : a) v /= static_cast<double>(n);
    }
}

std::size_t next_power_of_two(std::size_t n) {
    std::size_t p = 1;
    while (p < n) p <<= 1;
    return p;
}

// O(n^2) DFT for lengths that are not powers of two.  Used only by the
// surrogate generator, where n is typically a few thousand.  A table of roots
// of unity avoids trigonometric calls in the inner loop.
void dft_any(std::vector<Complex>& a, bool inverse) {
    const std::size_t n = a.size();
    std::vector<Complex> roots(n);
    const double sign = inverse ? 1.0 : -1.0;
    for (std::size_t k = 0; k < n; ++k) {
        const double angle = sign * kTwoPi * static_cast<double>(k) / static_cast<double>(n);
        roots[k] = Complex(std::cos(angle), std::sin(angle));
    }
    std::vector<Complex> out(n, Complex(0.0, 0.0));
    for (std::size_t k = 0; k < n; ++k) {
        Complex sum(0.0, 0.0);
        std::size_t idx = 0;
        for (std::size_t j = 0; j < n; ++j) {
            sum += a[j] * roots[idx];
            idx += k;
            if (idx >= n) idx -= n;
        }
        out[k] = inverse ? sum / static_cast<double>(n) : sum;
    }
    a.swap(out);
}

}  // namespace

Spectrum power_spectrum(const std::vector<double>& series, double dt, bool detrend,
                        bool hann, bool zero_pad) {
    if (series.size() < 4) throw std::invalid_argument("power_spectrum: series too short");
    if (!(dt > 0.0)) throw std::invalid_argument("power_spectrum: dt must be positive");
    const std::size_t n = series.size();
    const std::size_t nfft = zero_pad ? next_power_of_two(n) : n;
    std::vector<Complex> a(nfft, Complex(0.0, 0.0));
    const double mu = detrend ? mean(series) : 0.0;
    double window_sum = 0.0;
    double window_sq_sum = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        double w = 1.0;
        if (hann && n > 1) {
            w = 0.5 * (1.0 - std::cos(2.0 * kPi * static_cast<double>(i) /
                                      static_cast<double>(n - 1)));
        }
        a[i] = Complex((series[i] - mu) * w, 0.0);
        window_sum += w;
        window_sq_sum += w * w;
    }
    if (nfft == n && (n & (n - 1)) != 0) {
        dft_any(a, false);
    } else {
        fft(a, false);
    }
    const std::size_t half = nfft / 2;
    Spectrum out;
    out.nfft = nfft;
    out.df = 1.0 / (static_cast<double>(nfft) * dt);
    out.frequency.resize(half + 1);
    out.power.resize(half + 1, 0.0);
    const double norm = 1.0 / (static_cast<double>(nfft) * std::max(window_sq_sum, 1e-300));
    for (std::size_t k = 0; k <= half; ++k) {
        out.frequency[k] = static_cast<double>(k) * out.df;
        double p = std::norm(a[k]) * norm;
        if (k != 0 && k != half) p *= 2.0;
        out.power[k] = p;
    }
    (void)window_sum;
    return out;
}

std::vector<double> phase_randomized_surrogate(const std::vector<double>& series, Rng& rng) {
    if (series.size() < 4) throw std::invalid_argument("surrogate: series too short");
    const std::size_t n = series.size();
    std::vector<Complex> a(n);
    for (std::size_t i = 0; i < n; ++i) a[i] = Complex(series[i], 0.0);
    const bool pow2 = (n & (n - 1)) == 0;
    if (pow2) {
        fft(a, false);
    } else {
        dft_any(a, false);
    }
    // Keep DC and Nyquist real; randomize the phase of conjugate pairs.
    for (std::size_t k = 1; k <= (n - 1) / 2; ++k) {
        const double phase = rng.uniform(0.0, kTwoPi);
        const Complex unit(std::cos(phase), std::sin(phase));
        a[k] *= unit;
        if (k != n - k) a[n - k] = std::conj(a[k]);
    }
    if (n % 2 == 0) a[n / 2] = Complex(a[n / 2].real(), 0.0);
    if (pow2) {
        fft(a, true);
    } else {
        dft_any(a, true);
    }
    std::vector<double> out(n);
    for (std::size_t i = 0; i < n; ++i) out[i] = a[i].real();
    return out;
}

ZeroOneResult zero_one_test(const std::vector<double>& series, int n_c, std::uint64_t seed) {
    constexpr int kMaxC = 100000;
    if (series.size() < 100) throw std::invalid_argument("0-1 test: series too short");
    n_c = std::clamp(n_c, 1, kMaxC);
    const std::size_t n = series.size();
    const std::size_t n_max = std::max<std::size_t>(2, n / 10);
    const double mu = mean(series);
    std::vector<double> phi(n);
    for (std::size_t i = 0; i < n; ++i) phi[i] = series[i] - mu;

    Rng rng(seed);
    ZeroOneResult out;
    out.K_values.reserve(static_cast<std::size_t>(n_c));
    std::vector<double> p(n), q(n), m(n_max + 1), idx(n_max);
    for (std::size_t i = 0; i < n_max; ++i) idx[i] = static_cast<double>(i + 1);

    for (int trial = 0; trial < n_c; ++trial) {
        const double c = rng.uniform(kPi / 5.0, 4.0 * kPi / 5.0);
        double pc = 0.0;
        double qc = 0.0;
        for (std::size_t j = 0; j < n; ++j) {
            pc += phi[j] * std::cos(static_cast<double>(j) * c);
            qc += phi[j] * std::sin(static_cast<double>(j) * c);
            p[j] = pc;
            q[j] = qc;
        }
        for (std::size_t nn = 1; nn <= n_max; ++nn) {
            double s = 0.0;
            const std::size_t count = n - nn;
            for (std::size_t j = 0; j < count; ++j) {
                const double dp = p[j + nn] - p[j];
                const double dq = q[j + nn] - q[j];
                s += dp * dp + dq * dq;
            }
            // Standard 0-1 definition uses 1/N, not 1/(N-n); see the
            // asymptotic formula in the documentation.
            m[nn] = s / static_cast<double>(n);
        }
        // Correlation of M(n) with n; use indices 1..n_max.
        double sm = 0.0;
        double sn = 0.0;
        double smm = 0.0;
        double snn = 0.0;
        double smn = 0.0;
        const std::size_t cnt = n_max;
        for (std::size_t i = 1; i <= n_max; ++i) {
            const double x = static_cast<double>(i);
            const double y = m[i];
            sm += y;
            sn += x;
            smm += y * y;
            snn += x * x;
            smn += x * y;
        }
        const double denom = std::sqrt((static_cast<double>(cnt) * smm - sm * sm) *
                                       (static_cast<double>(cnt) * snn - sn * sn));
        const double k = denom > 0.0
                             ? (static_cast<double>(cnt) * smn - sm * sn) / denom
                             : 0.0;
        out.K_values.push_back(k);
    }
    std::vector<double> sorted = out.K_values;
    std::sort(sorted.begin(), sorted.end());
    out.K = sorted[sorted.size() / 2];
    return out;
}

std::vector<std::uint8_t> recurrence_matrix(const Embedding& emb, double epsilon,
                                            int theiler) {
    if (!(epsilon > 0.0)) throw std::invalid_argument("recurrence_matrix: epsilon > 0");
    const std::size_t n = emb.n;
    std::vector<std::uint8_t> R(n * n, 0);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            if (std::abs(static_cast<long long>(i) - static_cast<long long>(j)) <= theiler) {
                continue;
            }
            double d2 = 0.0;
            for (int q = 0; q < emb.m; ++q) {
                const double diff = emb.at(i, q) - emb.at(j, q);
                d2 += diff * diff;
            }
            if (std::sqrt(d2) <= epsilon) R[i * n + j] = 1;
        }
    }
    return R;
}

double epsilon_for_recurrence_rate(const Embedding& emb, double rate, int theiler) {
    if (!(rate > 0.0 && rate < 1.0)) {
        throw std::invalid_argument("epsilon_for_recurrence_rate: rate in (0,1)");
    }
    std::vector<double> d;
    const std::size_t n = emb.n;
    d.reserve(n * (n - 1) / 2);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i + 1; j < n; ++j) {
            if (std::abs(static_cast<long long>(i) - static_cast<long long>(j)) <= theiler) {
                continue;
            }
            double d2 = 0.0;
            for (int q = 0; q < emb.m; ++q) {
                const double diff = emb.at(i, q) - emb.at(j, q);
                d2 += diff * diff;
            }
            d.push_back(std::sqrt(d2));
        }
    }
    if (d.empty()) throw std::runtime_error("epsilon_for_recurrence_rate: no pairs");
    const std::size_t idx = std::min(d.size() - 1, static_cast<std::size_t>(rate * d.size()));
    std::nth_element(d.begin(), d.begin() + static_cast<std::ptrdiff_t>(idx), d.end());
    return d[idx];
}

RqaResult rqa_from_matrix(const std::vector<std::uint8_t>& matrix, std::size_t n,
                          double epsilon, int min_line) {
    if (matrix.size() != n * n) throw std::invalid_argument("rqa: matrix size mismatch");
    min_line = std::max(2, min_line);
    RqaResult out;
    out.epsilon = epsilon;

    auto at = [&](std::size_t i, std::size_t j) -> bool {
        return matrix[i * n + j] != 0;
    };

    std::size_t recurrent = 0;
    std::size_t total = 0;
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            if (i == j) continue;
            ++total;
            if (at(i, j)) ++recurrent;
        }
    }
    out.points = recurrent;
    out.recurrence_rate = total ? static_cast<double>(recurrent) / static_cast<double>(total) : 0.0;

    // Diagonal lines (excluding the main diagonal).
    std::size_t diag_points = 0;
    std::size_t diag_lines = 0;
    std::size_t longest = 0;
    std::size_t sum_len = 0;
    std::vector<std::size_t> diag_lengths;
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            if (i == j || !at(i, j)) continue;
            if (i > 0 && j > 0 && at(i - 1, j - 1)) continue;  // not a line start
            std::size_t len = 0;
            std::size_t ii = i;
            std::size_t jj = j;
            while (ii < n && jj < n && ii != jj && at(ii, jj)) {
                ++len;
                ++ii;
                ++jj;
            }
            if (len >= static_cast<std::size_t>(min_line)) {
                diag_points += len;
                ++diag_lines;
                longest = std::max(longest, len);
                sum_len += len;
                diag_lengths.push_back(len);
            }
        }
    }
    out.determinism = recurrent ? static_cast<double>(diag_points) / static_cast<double>(recurrent) : 0.0;
    out.longest_diagonal = static_cast<double>(longest);
    out.mean_diagonal = diag_lines ? static_cast<double>(sum_len) / static_cast<double>(diag_lines) : 0.0;
    if (!diag_lengths.empty()) {
        std::vector<std::size_t> uniq = diag_lengths;
        std::sort(uniq.begin(), uniq.end());
        uniq.erase(std::unique(uniq.begin(), uniq.end()), uniq.end());
        double h = 0.0;
        for (std::size_t len : uniq) {
            const std::size_t count = static_cast<std::size_t>(
                std::count(diag_lengths.begin(), diag_lengths.end(), len));
            const double p = static_cast<double>(count) / static_cast<double>(diag_lengths.size());
            h -= p * std::log(p);
        }
        out.entropy = h;
    }

    // Vertical lines (excluding the main diagonal).
    std::size_t vert_points = 0;
    std::size_t vert_lines = 0;
    std::size_t vert_sum = 0;
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            if (i == j || !at(i, j)) continue;
            if (i > 0 && at(i - 1, j)) continue;
            std::size_t len = 0;
            std::size_t ii = i;
            while (ii < n && ii != j && at(ii, j)) {
                ++len;
                ++ii;
            }
            if (len >= static_cast<std::size_t>(min_line)) {
                vert_points += len;
                ++vert_lines;
                vert_sum += len;
            }
        }
    }
    out.laminarity = recurrent ? static_cast<double>(vert_points) / static_cast<double>(recurrent) : 0.0;
    out.trapping_time = vert_lines ? static_cast<double>(vert_sum) / static_cast<double>(vert_lines) : 0.0;
    return out;
}

RqaResult rqa(const Embedding& emb, double epsilon, int theiler, int min_line) {
    return rqa_from_matrix(recurrence_matrix(emb, epsilon, theiler), emb.n, epsilon, min_line);
}

}  // namespace chaos
