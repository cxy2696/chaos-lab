// Time-series, embedding, spectral and recurrence tests.
// SPDX-License-Identifier: MIT
#include <algorithm>
#include <cmath>
#include <vector>

#include "chaos/embedding.hpp"
#include "chaos/integrators.hpp"
#include "chaos/lyapunov.hpp"
#include "chaos/systems.hpp"
#include "chaos/timeseries.hpp"
#include "chaos/tools.hpp"
#include "test_framework.hpp"

using namespace chaos;

namespace {

const std::vector<double>& lorenz_x_series() {
    static std::vector<double> series;
    if (series.empty()) {
        Lorenz l;
        OdeSystem sys = l.ode();
        State x{1.0, 1.0, 1.0};
        double t = 0.0;
        for (int i = 0; i < 5000; ++i) {
            x = rk4_step(sys, x, t, 0.01);
            t += 0.01;
        }
        series.reserve(8000);
        for (int i = 0; i < 8000; ++i) {
            series.push_back(x[0]);
            x = rk4_step(sys, x, t, 0.01);
            t += 0.01;
        }
    }
    return series;
}

const std::vector<double>& henon_x_series() {
    static std::vector<double> series;
    if (series.empty()) {
        Henon h;
        MapSystem map = h.map();
        State x{0.0, 0.0};
        for (int i = 0; i < 1000; ++i) x = map.step(x);
        series.reserve(6000);
        for (int i = 0; i < 6000; ++i) {
            x = map.step(x);
            series.push_back(x[0]);
        }
    }
    return series;
}

double power_correlation(const Spectrum& a, const Spectrum& b) {
    const std::size_t n = std::min(a.power.size(), b.power.size());
    double ma = 0.0, mb = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        ma += a.power[i];
        mb += b.power[i];
    }
    ma /= static_cast<double>(n);
    mb /= static_cast<double>(n);
    double sa = 0.0, sb = 0.0, sab = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double da = a.power[i] - ma;
        const double db = b.power[i] - mb;
        sa += da * da;
        sb += db * db;
        sab += da * db;
    }
    return sab / std::sqrt(sa * sb);
}

}  // namespace

TEST_CASE(ami_and_fnn_recover_lorenz_embedding) {
    const auto& s = lorenz_x_series();
    const auto ami = average_mutual_information(s, 60, 32);
    const int tau = first_ami_minimum(ami);
    CHECK_IN(tau, 10, 24);
    const auto ac = autocorrelation(s, 120);
    const int ac_tau = autocorrelation_1e(ac);
    CHECK_IN(ac_tau, 20, 40);
    const FnnResult fnn = false_nearest_neighbors(s, 8, tau, 20, 15.0);
    CHECK(fnn.suggested_m == 3);
}

TEST_CASE(grassberger_procaccia_lorenz) {
    const auto& s = lorenz_x_series();
    const CorrelationResult r = correlation_dimension(s, 3, 17, 0.0, 0.0, 60, 50, 0.2, 0.8);
    CHECK_IN(r.d2, 1.98, 2.12);
    CHECK(r.r2 > 0.99);
}

TEST_CASE(grassberger_procaccia_henon) {
    const auto& s = henon_x_series();
    const CorrelationResult r = correlation_dimension(s, 2, 1, 0.0, 0.0, 50, 10, 0.2, 0.8);
    CHECK_IN(r.d2, 1.10, 1.25);
    CHECK(r.r2 > 0.99);
}

TEST_CASE(grassberger_procaccia_henon_manual_scaling) {
    const auto& s = henon_x_series();
    const CorrelationResult r =
        correlation_dimension(s, 2, 1, 0.0005, 0.5, 50, 10, 0.2, 0.8);
    CHECK_IN(r.d2, 1.18, 1.25);
    CHECK(r.r2 > 0.99);
}

TEST_CASE(rosenstein_henon_explicit_window) {
    const auto& s = henon_x_series();
    const RosensteinResult r = rosenstein_lyapunov(s, 2, 1, 10, 40, 1.0, 0, 10);
    CHECK_IN(r.lambda, 0.39, 0.44);
    CHECK(r.r2 > 0.99);
}

TEST_CASE(rosenstein_lorenz_automatic_window) {
    const auto& s = lorenz_x_series();
    const RosensteinResult r = rosenstein_lyapunov(s, 3, 17, 100, 300, 0.01);
    CHECK_IN(r.lambda, 0.78, 0.95);
    CHECK(r.fit_end > r.fit_begin);
}

TEST_CASE(zero_one_test_separates_regular_and_chaotic) {
    std::vector<double> sine(4000);
    for (int i = 0; i < 4000; ++i) {
        sine[static_cast<std::size_t>(i)] = std::sin(0.05 * kTwoPi * i);
    }
    CHECK(zero_one_test(sine, 32).K < 0.1);

    std::vector<double> logistic(4000);
    double x = 0.123456789;
    for (int i = 0; i < 4000; ++i) {
        x = 4.0 * x * (1.0 - x);
        logistic[static_cast<std::size_t>(i)] = x;
    }
    CHECK(zero_one_test(logistic, 32).K > 0.9);

    // Oversampled flows can mask the diffusion; stride to a sampling interval
    // comparable to the correlation time before applying the 0-1 test.
    const auto& s = lorenz_x_series();
    std::vector<double> sub;
    for (std::size_t i = 0; i < s.size(); i += 10) sub.push_back(s[i]);
    CHECK(zero_one_test(sub, 32).K > 0.9);
}

TEST_CASE(fft_finds_sinusoid) {
    const std::size_t n = 4096;
    const double dt = 0.001;
    const double f0 = 50.0;
    std::vector<double> x(n);
    for (std::size_t i = 0; i < n; ++i) {
        x[i] = std::sin(kTwoPi * f0 * static_cast<double>(i) * dt);
    }
    const Spectrum sp = power_spectrum(x, dt);
    std::size_t peak = 0;
    for (std::size_t i = 1; i < sp.power.size(); ++i) {
        if (sp.power[i] > sp.power[peak]) peak = i;
    }
    CHECK_NEAR(sp.frequency[peak], f0, 0.5);
}

TEST_CASE(phase_randomized_surrogate_preserves_spectrum) {
    const auto& s = lorenz_x_series();
    Rng rng(42);
    const std::vector<double> surrogate = phase_randomized_surrogate(s, rng);
    CHECK(surrogate.size() == s.size());
    CHECK_NEAR(mean(surrogate), mean(s), 0.05);
    // Compare the exact (unwindowed) periodograms: phase randomization
    // preserves |DFT| by construction.  A Hann window would mix phases.
    const Spectrum a = power_spectrum(s, 0.01, false, false, false);
    const Spectrum b = power_spectrum(surrogate, 0.01, false, false, false);
    CHECK(power_correlation(a, b) > 0.9999);
}

TEST_CASE(rqa_regular_vs_noise) {
    std::vector<double> sine(3000), noise(3000);
    Rng rng(7);
    for (std::size_t i = 0; i < sine.size(); ++i) {
        sine[i] = std::sin(0.05 * kTwoPi * static_cast<double>(i));
        noise[i] = rng.normal();
    }
    const Embedding es = delay_embed(sine, 2, 10);
    const Embedding en = delay_embed(noise, 2, 10);
    const double eps_s = epsilon_for_recurrence_rate(es, 0.05, 10);
    const double eps_n = epsilon_for_recurrence_rate(en, 0.05, 10);
    const RqaResult rs = rqa(es, eps_s, 10, 2);
    const RqaResult rn = rqa(en, eps_n, 10, 2);
    CHECK_IN(rs.recurrence_rate, 0.04, 0.06);
    CHECK(rs.determinism > 0.45);
    CHECK(rn.determinism < 0.25);
}

TEST_CASE(poincare_section_of_lorenz) {
    Lorenz l;
    PoincareOptions po;
    po.coordinate = 2;
    po.value = 27.0;
    po.direction = -1;
    po.h = 0.005;
    po.t_transient = 20.0;
    po.t_max = 100.0;
    const auto pts = poincare_section(l.ode(), {1.0, 1.0, 1.0}, po);
    CHECK(pts.size() > 80);
    CHECK(pts.size() < 200);
    for (const auto& p : pts) CHECK_NEAR(p[2], 27.0, 1e-12);
}

TEST_CASE(bifurcation_scan_stable_fixed_point) {
    auto f = [](double x, double r) { return r * x * (1.0 - x); };
    const auto pts = bifurcation_scan_1d(f, 0.123456789, 2.5, 4.0, 50, 500, 10);
    CHECK(pts.size() == 500);
    CHECK_NEAR(pts.front().parameter, 2.5, 1e-12);
    CHECK_NEAR(pts.front().value, 0.6, 1e-9);
    CHECK_NEAR(pts.back().parameter, 4.0, 1e-12);
}
