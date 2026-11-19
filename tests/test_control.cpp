// Control, synchronization and SINDy tests.
// SPDX-License-Identifier: MIT
#include <cmath>
#include <string>
#include <vector>

#include "chaos/control.hpp"
#include "chaos/integrators.hpp"
#include "chaos/sindy.hpp"
#include "chaos/systems.hpp"
#include "test_framework.hpp"

using namespace chaos;

TEST_CASE(ogy_stabilizes_logistic_fixed_point) {
    const double r = 3.8;
    const double x_star = 1.0 - 1.0 / r;
    OgyOptions o;
    o.x_star = x_star;
    o.p0 = r;
    o.lambda_u = 2.0 - r;
    o.g = x_star * (1.0 - x_star);
    o.tolerance = 0.2;
    o.dp_max = 0.02;
    o.steps = 200000;
    const OgyResult res =
        ogy_control_1d([](double x, double p) { return p * x * (1.0 - x); }, 0.2, o);
    CHECK(res.locked);
    CHECK(res.activations > 1000);
    CHECK(res.max_error_tail < 1e-8);
}

TEST_CASE(pyragas_delayed_feedback_stabilizes) {
    const double r = 3.8;
    const double x_star = 1.0 - 1.0 / r;
    const auto x = pyragas_map_control(
        [r](double v) { return r * v * (1.0 - v); }, x_star + 1e-6, 1, -0.5, 5000);
    CHECK(std::fabs(x.back() - x_star) < 1e-10);
}

TEST_CASE(pecora_carroll_synchronization) {
    const SyncResult r =
        pecora_carroll_lorenz(Lorenz{}, {1.0, 1.0, 1.0}, {0.0, 0.0}, 0.005, 5.0, 40.0);
    CHECK(r.synchronized);
    CHECK(r.final_error < 1e-8);
    CHECK(r.settling_time > 0.0);
    CHECK(r.max_error_after_transient < 5e-4);
}

TEST_CASE(sindy_recovers_lorenz_from_exact_derivatives) {
    Lorenz l;
    OdeSystem sys = l.ode();
    State x{1.0, 1.0, 1.0};
    double t = 0.0;
    for (int i = 0; i < 5000; ++i) {
        x = rk4_step(sys, x, t, 0.002);
        t += 0.002;
    }
    std::vector<State> X, dX;
    for (int i = 0; i < 4000; ++i) {
        X.push_back(x);
        dX.push_back(sys.rhs(x, t));
        x = rk4_step(sys, x, t, 0.002);
        t += 0.002;
    }
    SindyOptions o;
    o.degree = 2;
    o.threshold = 0.1;
    const SindyModel model = sindy_fit(X, dX, o);
    auto coeff = [&](const std::string& term, int dim) {
        for (std::size_t i = 0; i < model.terms.size(); ++i) {
            if (model.terms[i].name == term) return model.coefficients(i, static_cast<std::size_t>(dim));
        }
        throw std::runtime_error("missing SINDy term " + term);
    };
    CHECK_NEAR(coeff("x0", 0), -10.0, 1e-6);
    CHECK_NEAR(coeff("x1", 0), 10.0, 1e-6);
    CHECK_NEAR(coeff("x0", 1), 28.0, 1e-6);
    CHECK_NEAR(coeff("x1", 1), -1.0, 1e-6);
    CHECK_NEAR(coeff("x0*x2", 1), -1.0, 1e-6);
    CHECK_NEAR(coeff("x2", 2), -8.0 / 3.0, 1e-6);
    CHECK_NEAR(coeff("x0*x1", 2), 1.0, 1e-6);
}

TEST_CASE(sindy_recovers_lorenz_from_finite_differences) {
    Lorenz l;
    OdeSystem sys = l.ode();
    const double h = 0.001;
    State x{1.0, 1.0, 1.0};
    double t = 0.0;
    for (int i = 0; i < 2000; ++i) {
        x = rk4_step(sys, x, t, h);
        t += h;
    }
    std::vector<State> X;
    for (int i = 0; i < 4000; ++i) {
        X.push_back(x);
        x = rk4_step(sys, x, t, h);
        t += h;
    }
    const std::vector<State> dX = finite_difference_derivatives(X, h);
    SindyOptions o;
    o.degree = 2;
    o.threshold = 0.2;
    const SindyModel model = sindy_fit(X, dX, o);
    auto coeff = [&](const std::string& term, int dim) {
        for (std::size_t i = 0; i < model.terms.size(); ++i) {
            if (model.terms[i].name == term) return model.coefficients(i, static_cast<std::size_t>(dim));
        }
        throw std::runtime_error("missing SINDy term " + term);
    };
    CHECK_NEAR(coeff("x0", 0), -10.0, 2e-3);
    CHECK_NEAR(coeff("x1", 0), 10.0, 2e-3);
    CHECK_NEAR(coeff("x0", 1), 28.0, 5e-3);
    CHECK_NEAR(coeff("x0*x2", 1), -1.0, 2e-3);
    CHECK_NEAR(coeff("x2", 2), -8.0 / 3.0, 2e-3);
    CHECK_NEAR(coeff("x0*x1", 2), 1.0, 2e-3);
}

TEST_CASE(sindy_cubic_library_includes_mixed_terms) {
    // dx/dt = -x + 0.5 x^2 y, dy/dt = -y.  The cubic library must contain
    // x0^2*x1, not just pure powers (a regression test for the generator).
    OdeSystem sys;
    sys.name = "cubic_test";
    sys.dim = 2;
    sys.rhs = [](const State& s, double) {
        return State{-s[0] + 0.5 * s[0] * s[0] * s[1], -s[1]};
    };
    // One decaying trajectory is not enough to identify the coefficients:
    // x, y and x^2 y are strongly collinear along it.  Use several initial
    // conditions, which is the usual SINDy requirement of sufficiently rich
    // excitation.
    const double h = 0.001;
    std::vector<State> X, dX;
    for (double x0 : {0.4, 0.8, 1.2, 1.6, 2.0}) {
        for (double y0 : {0.2, 0.6, 1.0, 1.4}) {
            State x{x0, y0};
            double t = 0.0;
            for (int i = 0; i < 500; ++i) {
                X.push_back(x);
                dX.push_back(sys.rhs(x, t));
                x = rk4_step(sys, x, t, h);
                t += h;
            }
        }
    }
    SindyOptions o;
    o.degree = 3;
    o.threshold = 0.05;
    const SindyModel model = sindy_fit(X, dX, o);
    auto coeff = [&](const std::string& term, int dim) {
        for (std::size_t i = 0; i < model.terms.size(); ++i) {
            if (model.terms[i].name == term) {
                return model.coefficients(i, static_cast<std::size_t>(dim));
            }
        }
        throw std::runtime_error("missing SINDy term " + term);
    };
    CHECK_NEAR(coeff("x0", 0), -1.0, 1e-6);
    CHECK_NEAR(coeff("x0*x0*x1", 0), 0.5, 1e-6);
    CHECK_NEAR(coeff("x1", 1), -1.0, 1e-6);
}

TEST_CASE(sindy_quartic_library_contains_intermediate_degrees) {
    // Degree 4 must include cubic and quartic terms, not only the top degree.
    std::vector<State> X(40, State{0.5, 0.2});
    std::vector<State> dX(40, State{0.0, 0.0});
    SindyOptions o;
    o.degree = 4;
    const SindyModel model = sindy_fit(X, dX, o);
    bool has_cubic = false;
    bool has_quartic = false;
    for (const auto& term : model.terms) {
        if (term.name == "x0*x0*x1") has_cubic = true;
        if (term.name == "x0*x0*x0*x1") has_quartic = true;
    }
    CHECK(has_cubic);
    CHECK(has_quartic);
}
