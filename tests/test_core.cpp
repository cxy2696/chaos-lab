// Integrator, model and Jacobian tests.
// SPDX-License-Identifier: MIT
#include <cmath>
#include <vector>

#include "chaos/integrators.hpp"
#include "chaos/lyapunov.hpp"
#include "chaos/systems.hpp"
#include "test_framework.hpp"

using namespace chaos;

namespace {

double max_jacobian_error(const Matrix& a, const Matrix& b) {
    double e = 0.0;
    for (std::size_t i = 0; i < a.rows; ++i) {
        for (std::size_t j = 0; j < a.cols; ++j) {
            e = std::max(e, std::fabs(a(i, j) - b(i, j)) /
                                (1.0 + std::fabs(b(i, j))));
        }
    }
    return e;
}

}  // namespace

TEST_CASE(rk4_exponential_decay) {
    OdeSystem sys;
    sys.name = "decay";
    sys.dim = 1;
    sys.rhs = [](const State& x, double) { return State{-x[0]}; };
    State x{1.0};
    const double h = 0.01;
    for (int i = 0; i < 100; ++i) x = rk4_step(sys, x, i * h, h);
    CHECK_NEAR(x[0], std::exp(-1.0), 1e-9);
}

TEST_CASE(rk45_adaptive_accuracy) {
    Lorenz l;
    OdeSystem sys = l.ode();
    State a{1.0, 1.0, 1.0};
    const double h = 1e-3;
    for (int i = 0; i < 1000; ++i) a = rk4_step(sys, a, i * h, h);
    const Rk45Result r = rk45_integrate(sys, {1.0, 1.0, 1.0}, 0.0, 1.0, 1e-3,
                                        Rk45Options{1e-10, 1e-12});
    CHECK(r.ok);
    CHECK(r.accepted > 10);
    CHECK(distance(a, r.x) < 1e-6);
}

TEST_CASE(analytic_jacobians_match_finite_differences) {
    {
        Lorenz l;
        const State x{1.2, -0.7, 2.4};
        CHECK(max_jacobian_error(l.jacobian(x), numerical_jacobian(l.ode(), x, 0.0)) < 1e-6);
    }
    {
        Rossler r;
        const State x{1.2, -0.7, 2.4};
        CHECK(max_jacobian_error(r.jacobian(x), numerical_jacobian(r.ode(), x, 0.0)) < 1e-6);
    }
    {
        Chua c;
        const State x{0.4, -0.2, 0.1};
        CHECK(max_jacobian_error(c.jacobian(x), numerical_jacobian(c.ode(), x, 0.0)) < 1e-6);
    }
    {
        Henon h;
        const State x{0.3, -0.2};
        CHECK(max_jacobian_error(h.jacobian(x), numerical_jacobian(h.map(), x)) < 1e-6);
    }
}

TEST_CASE(lorenz_lyapunov_spectrum) {
    Lorenz l;
    LyapunovOptions o;
    o.h = 0.01;
    o.transient = 30.0;
    o.duration = 300.0;
    const LyapunovSpectrum sp = lyapunov_spectrum(l.ode(), {1.0, 1.0, 1.0}, o);
    CHECK(sp.exponents.size() == 3);
    CHECK_IN(sp.exponents[0], 0.88, 0.94);
    CHECK_NEAR(sp.exponents[1], 0.0, 0.05);
    CHECK_IN(sp.exponents[2], -14.65, -14.50);
    CHECK_NEAR(sp.exponents[0] + sp.exponents[1] + sp.exponents[2], l.divergence(), 1e-2);
}

TEST_CASE(lorenz_two_trajectory_max_exponent) {
    Lorenz l;
    const double lambda = lyapunov_max_two_trajectories(l.ode(), {1.0, 1.0, 1.0}, 1e-9,
                                                        0.01, 30.0, 300.0);
    CHECK_IN(lambda, 0.88, 0.94);
}

TEST_CASE(henon_lyapunov_spectrum) {
    Henon h;
    const LyapunovSpectrum sp =
        lyapunov_spectrum_map(h.map(), {0.0, 0.0}, 1000, 100000, 1);
    CHECK_IN(sp.exponents[0], 0.40, 0.44);
    CHECK_IN(sp.exponents[1], -1.65, -1.60);
    CHECK_NEAR(sp.exponents[0] + sp.exponents[1], std::log(h.b), 5e-3);
}

TEST_CASE(cat_map_lyapunov_exponent) {
    CatMap cat;
    const LyapunovSpectrum sp = lyapunov_spectrum_map(cat.map(), {0.1, 0.2}, 10, 1000, 1);
    const double expected = std::log((3.0 + std::sqrt(5.0)) / 2.0);
    CHECK_NEAR(sp.exponents[0], expected, 1e-3);
}

TEST_CASE(logistic_lyapunov_known_values) {
    CHECK_NEAR(Logistic{4.0}.lyapunov(0.123456789, 1000, 1000000), std::log(2.0), 1e-3);
    CHECK(Logistic{3.5}.lyapunov(0.123456789, 1000, 100000) < 0.0);
    // Period-3 window ("period three implies chaos" is about a different
    // parameter regime; here the stable 3-cycle has a negative exponent).
    CHECK(Logistic{3.83}.lyapunov(0.5, 1000, 100000) < 0.0);
}

TEST_CASE(symplectic_integrators_conserve_energy) {
    auto accel = [](const State& q) { return State{-q[0], -q[1]}; };
    auto potential = [](const State& q) { return 0.5 * (q[0] * q[0] + q[1] * q[1]); };
    State q1{1.0, 0.0}, v1{0.0, 1.0};
    const double e0 = hamiltonian_energy(q1, v1, potential);
    for (int i = 0; i < 10000; ++i) verlet_step(accel, q1, v1, 0.01);
    CHECK(std::fabs(hamiltonian_energy(q1, v1, potential) - e0) < 1e-6);

    State q2{1.0, 0.0}, v2{0.0, 1.0};
    for (int i = 0; i < 10000; ++i) yoshida4_step(accel, q2, v2, 0.01);
    CHECK(std::fabs(hamiltonian_energy(q2, v2, potential) - e0) < 1e-9);
}

TEST_CASE(three_body_figure_eight_periodic) {
    ThreeBody tb;
    OdeSystem sys = tb.ode();
    State x = ThreeBody::figure_eight_initial_state();
    const State x0 = x;
    const double e0 = tb.energy(x);
    const double h = 0.002;
    const double period = 6.32591398;
    const int n = static_cast<int>(std::llround(period / h));
    double t = 0.0;
    for (int i = 0; i < n; ++i) {
        x = rk4_step(sys, x, t, h);
        t += h;
    }
    CHECK(distance(x, x0) < 1e-2);
    CHECK(std::fabs(tb.energy(x) - e0) < 1e-9);
}

TEST_CASE(double_pendulum_conserves_energy) {
    DoublePendulum dp;
    OdeSystem sys = dp.ode();
    State x{1.0, 2.0, 0.0, 0.0};
    const double e0 = dp.energy(x);
    const double h = 0.001;
    double t = 0.0;
    for (int i = 0; i < 5000; ++i) {
        x = rk4_step(sys, x, t, h);
        t += h;
    }
    CHECK(std::fabs(dp.energy(x) - e0) < 1e-6);
}

TEST_CASE(chua_stays_bounded) {
    Chua c;
    OdeSystem sys = c.ode();
    State x{0.1, 0.0, 0.0};
    double t = 0.0;
    for (int i = 0; i < 20000; ++i) {
        x = rk4_step(sys, x, t, 0.01);
        t += 0.01;
    }
    CHECK(is_finite_state(x));
    CHECK(std::fabs(x[0]) < 10.0);
    CHECK(std::fabs(x[1]) < 10.0);
    CHECK(std::fabs(x[2]) < 20.0);
}

TEST_CASE(rossler_lyapunov_spectrum) {
    Rossler r;
    LyapunovOptions o;
    o.h = 0.01;
    o.transient = 50.0;
    o.duration = 500.0;
    const LyapunovSpectrum sp = lyapunov_spectrum(r.ode(), {1.0, 1.0, 1.0}, o);
    CHECK_IN(sp.exponents[0], 0.07, 0.09);
    CHECK_IN(sp.exponents[2], -5.5, -5.3);
    CHECK(sp.exponents[0] + sp.exponents[1] + sp.exponents[2] < 0.0);
}

TEST_CASE(chua_lyapunov_spectrum) {
    Chua c;
    LyapunovOptions o;
    o.h = 0.01;
    o.transient = 50.0;
    o.duration = 500.0;
    const LyapunovSpectrum sp = lyapunov_spectrum(c.ode(), {0.1, 0.0, 0.0}, o);
    CHECK_IN(sp.exponents[0], 0.35, 0.50);
    CHECK_IN(sp.exponents[2], -4.3, -4.1);
    CHECK(sp.exponents[0] + sp.exponents[1] + sp.exponents[2] < 0.0);
}
