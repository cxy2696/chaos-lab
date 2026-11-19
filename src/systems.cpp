// SPDX-License-Identifier: MIT
#include "chaos/systems.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace chaos {

namespace {

void expect_dim(const State& s, std::size_t n, const char* who) {
    if (s.size() != n) {
        throw std::invalid_argument(std::string(who) + ": bad state dimension");
    }
}

double wrap_unit(double x) { return x - std::floor(x); }

}  // namespace

// ---------------------------------------------------------------------------
// Logistic
// ---------------------------------------------------------------------------

double Logistic::lyapunov(double x0, std::size_t transient, std::size_t n) const {
    double x = x0;
    for (std::size_t i = 0; i < transient; ++i) x = step(x);
    double sum = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double d = std::fabs(derivative(x));
        sum += std::log(std::max(d, 1e-300));
        x = step(x);
    }
    return n ? sum / static_cast<double>(n) : 0.0;
}

MapSystem Logistic::map() const {
    MapSystem m;
    m.name = "logistic";
    m.description = "Logistic map x' = r x (1-x)";
    m.dim = 1;
    m.param_names = {"r"};
    m.params = {r};
    const double rr = r;
    m.step = [rr](const State& s) {
        expect_dim(s, 1, "logistic");
        return State{rr * s[0] * (1.0 - s[0])};
    };
    m.jacobian = [rr](const State& s) {
        expect_dim(s, 1, "logistic");
        return Matrix(1, 1, rr * (1.0 - 2.0 * s[0]));
    };
    return m;
}

// ---------------------------------------------------------------------------
// Henon
// ---------------------------------------------------------------------------

State Henon::step(const State& s) const {
    expect_dim(s, 2, "henon");
    return State{1.0 - a * s[0] * s[0] + s[1], b * s[0]};
}

Matrix Henon::jacobian(const State& s) const {
    expect_dim(s, 2, "henon");
    Matrix j(2, 2);
    j(0, 0) = -2.0 * a * s[0];
    j(0, 1) = 1.0;
    j(1, 0) = b;
    j(1, 1) = 0.0;
    return j;
}

MapSystem Henon::map() const {
    MapSystem m;
    m.name = "henon";
    m.description = "Henon map";
    m.dim = 2;
    m.param_names = {"a", "b"};
    m.params = {a, b};
    const Henon self = *this;
    m.step = [self](const State& s) { return self.step(s); };
    m.jacobian = [self](const State& s) { return self.jacobian(s); };
    return m;
}

// ---------------------------------------------------------------------------
// Cat map
// ---------------------------------------------------------------------------

State CatMap::step(const State& s) const {
    expect_dim(s, 2, "cat");
    return State{wrap_unit(s[0] + s[1]), wrap_unit(s[0] + 2.0 * s[1])};
}

Matrix CatMap::jacobian(const State&) const {
    Matrix j(2, 2);
    j(0, 0) = 1.0;
    j(0, 1) = 1.0;
    j(1, 0) = 1.0;
    j(1, 1) = 2.0;
    return j;
}

MapSystem CatMap::map() const {
    MapSystem m;
    m.name = "cat";
    m.description = "Arnold cat map on the unit torus";
    m.dim = 2;
    const CatMap self = *this;
    m.step = [self](const State& s) { return self.step(s); };
    m.jacobian = [self](const State& s) { return self.jacobian(s); };
    return m;
}

// ---------------------------------------------------------------------------
// Lorenz
// ---------------------------------------------------------------------------

State Lorenz::f(const State& s, double) const {
    expect_dim(s, 3, "lorenz");
    return State{sigma * (s[1] - s[0]), s[0] * (rho - s[2]) - s[1],
                 s[0] * s[1] - beta * s[2]};
}

Matrix Lorenz::jacobian(const State& s, double) const {
    expect_dim(s, 3, "lorenz");
    Matrix j(3, 3);
    j(0, 0) = -sigma;
    j(0, 1) = sigma;
    j(1, 0) = rho - s[2];
    j(1, 1) = -1.0;
    j(1, 2) = -s[0];
    j(2, 0) = s[1];
    j(2, 1) = s[0];
    j(2, 2) = -beta;
    return j;
}

OdeSystem Lorenz::ode() const {
    OdeSystem sys;
    sys.name = "lorenz";
    sys.description = "Lorenz (1963) convection model";
    sys.dim = 3;
    sys.param_names = {"sigma", "rho", "beta"};
    sys.params = {sigma, rho, beta};
    const Lorenz self = *this;
    sys.rhs = [self](const State& s, double t) { return self.f(s, t); };
    sys.jacobian = [self](const State& s, double t) { return self.jacobian(s, t); };
    return sys;
}

// ---------------------------------------------------------------------------
// Rossler
// ---------------------------------------------------------------------------

State Rossler::f(const State& s, double) const {
    expect_dim(s, 3, "rossler");
    return State{-s[1] - s[2], s[0] + a * s[1], b + s[2] * (s[0] - c)};
}

Matrix Rossler::jacobian(const State& s, double) const {
    expect_dim(s, 3, "rossler");
    Matrix j(3, 3);
    j(0, 0) = 0.0;
    j(0, 1) = -1.0;
    j(0, 2) = -1.0;
    j(1, 0) = 1.0;
    j(1, 1) = a;
    j(1, 2) = 0.0;
    j(2, 0) = s[2];
    j(2, 1) = 0.0;
    j(2, 2) = s[0] - c;
    return j;
}

OdeSystem Rossler::ode() const {
    OdeSystem sys;
    sys.name = "rossler";
    sys.description = "Rossler (1976) funnel attractor";
    sys.dim = 3;
    sys.param_names = {"a", "b", "c"};
    sys.params = {a, b, c};
    const Rossler self = *this;
    sys.rhs = [self](const State& s, double t) { return self.f(s, t); };
    sys.jacobian = [self](const State& s, double t) { return self.jacobian(s, t); };
    return sys;
}

// ---------------------------------------------------------------------------
// Chua
// ---------------------------------------------------------------------------

double Chua::nonlinearity(double x) const {
    return m1 * x + 0.5 * (m0 - m1) * (std::fabs(x + 1.0) - std::fabs(x - 1.0));
}

double Chua::nonlinearity_derivative(double x) const {
    if (x > 1.0 || x < -1.0) return m1;
    return m0;
}

State Chua::f(const State& s, double) const {
    expect_dim(s, 3, "chua");
    return State{alpha * (s[1] - s[0] - nonlinearity(s[0])), s[0] - s[1] + s[2],
                 -beta * s[1]};
}

Matrix Chua::jacobian(const State& s, double) const {
    expect_dim(s, 3, "chua");
    Matrix j(3, 3);
    j(0, 0) = alpha * (-1.0 - nonlinearity_derivative(s[0]));
    j(0, 1) = alpha;
    j(1, 0) = 1.0;
    j(1, 1) = -1.0;
    j(1, 2) = 1.0;
    j(2, 1) = -beta;
    return j;
}

OdeSystem Chua::ode() const {
    OdeSystem sys;
    sys.name = "chua";
    sys.description = "Chua circuit (double scroll)";
    sys.dim = 3;
    sys.param_names = {"alpha", "beta", "m0", "m1"};
    sys.params = {alpha, beta, m0, m1};
    const Chua self = *this;
    sys.rhs = [self](const State& s, double t) { return self.f(s, t); };
    sys.jacobian = [self](const State& s, double t) { return self.jacobian(s, t); };
    return sys;
}

// ---------------------------------------------------------------------------
// Double pendulum
// ---------------------------------------------------------------------------

State DoublePendulum::f(const State& s, double) const {
    expect_dim(s, 4, "double_pendulum");
    const double th1 = s[0];
    const double th2 = s[1];
    const double w1 = s[2];
    const double w2 = s[3];
    const double delta = th1 - th2;
    const double denom = 2.0 * m1 + m2 - m2 * std::cos(2.0 * delta);
    if (std::fabs(denom) < 1e-14) return State{0, 0, 0, 0};

    const double a1 =
        (-g * (2.0 * m1 + m2) * std::sin(th1) - m2 * g * std::sin(th1 - 2.0 * th2) -
         2.0 * std::sin(delta) * m2 * (w2 * w2 * l2 + w1 * w1 * l1 * std::cos(delta))) /
        (l1 * denom);
    const double a2 =
        (2.0 * std::sin(delta) *
         (w1 * w1 * l1 * (m1 + m2) + g * (m1 + m2) * std::cos(th1) +
          w2 * w2 * l2 * m2 * std::cos(delta))) /
        (l2 * denom);
    return State{w1, w2, a1, a2};
}

OdeSystem DoublePendulum::ode() const {
    OdeSystem sys;
    sys.name = "double_pendulum";
    sys.description = "Planar double pendulum";
    sys.dim = 4;
    sys.param_names = {"m1", "m2", "l1", "l2", "g"};
    sys.params = {m1, m2, l1, l2, g};
    const DoublePendulum self = *this;
    sys.rhs = [self](const State& s, double t) { return self.f(s, t); };
    return sys;
}

double DoublePendulum::energy(const State& s) const {
    expect_dim(s, 4, "double_pendulum");
    const double th1 = s[0];
    const double th2 = s[1];
    const double w1 = s[2];
    const double w2 = s[3];
    const double kinetic =
        0.5 * (m1 + m2) * l1 * l1 * w1 * w1 + 0.5 * m2 * l2 * l2 * w2 * w2 +
        m2 * l1 * l2 * w1 * w2 * std::cos(th1 - th2);
    const double potential = -(m1 + m2) * g * l1 * std::cos(th1) -
                             m2 * g * l2 * std::cos(th2);
    return kinetic + potential;
}

// ---------------------------------------------------------------------------
// Three-body
// ---------------------------------------------------------------------------

State ThreeBody::f(const State& s, double) const {
    expect_dim(s, 12, "three_body");
    const std::array<double, 3> m{m1, m2, m3};
    State out(12, 0.0);
    for (int i = 0; i < 3; ++i) {
        out[2 * i] = s[6 + 2 * i];
        out[2 * i + 1] = s[7 + 2 * i];
    }
    for (int i = 0; i < 3; ++i) {
        double ax = 0.0;
        double ay = 0.0;
        for (int j = 0; j < 3; ++j) {
            if (i == j) continue;
            const double dx = s[2 * j] - s[2 * i];
            const double dy = s[2 * j + 1] - s[2 * i + 1];
            const double r2 = dx * dx + dy * dy;
            const double inv_r3 = 1.0 / (r2 * std::sqrt(r2) + 1e-300);
            ax += G * m[j] * dx * inv_r3;
            ay += G * m[j] * dy * inv_r3;
        }
        out[6 + 2 * i] = ax;
        out[7 + 2 * i] = ay;
    }
    return out;
}

OdeSystem ThreeBody::ode() const {
    OdeSystem sys;
    sys.name = "three_body";
    sys.description = "Planar three-body problem (figure-eight choreography)";
    sys.dim = 12;
    sys.param_names = {"m1", "m2", "m3", "G"};
    sys.params = {m1, m2, m3, G};
    const ThreeBody self = *this;
    sys.rhs = [self](const State& s, double t) { return self.f(s, t); };
    return sys;
}

double ThreeBody::energy(const State& s) const {
    expect_dim(s, 12, "three_body");
    const std::array<double, 3> m{m1, m2, m3};
    double kinetic = 0.0;
    for (int i = 0; i < 3; ++i) {
        const double vx = s[6 + 2 * i];
        const double vy = s[7 + 2 * i];
        kinetic += 0.5 * m[i] * (vx * vx + vy * vy);
    }
    double potential = 0.0;
    for (int i = 0; i < 3; ++i) {
        for (int j = i + 1; j < 3; ++j) {
            const double dx = s[2 * j] - s[2 * i];
            const double dy = s[2 * j + 1] - s[2 * i + 1];
            potential -= G * m[i] * m[j] / std::sqrt(dx * dx + dy * dy);
        }
    }
    return kinetic + potential;
}

State ThreeBody::figure_eight_initial_state() {
    return State{-0.97000436, 0.24308753, 0.97000436, -0.24308753, 0.0, 0.0,
                 0.4662036850, 0.4323657300, 0.4662036850, 0.4323657300,
                 -0.93240737, -0.86473146};
}

// ---------------------------------------------------------------------------
// Catalogue
// ---------------------------------------------------------------------------

const std::vector<OdeSystemInfo>& ode_catalog() {
    static const std::vector<OdeSystemInfo> catalog = {
        {"lorenz",
         "Lorenz (1963) convection model; sigma=10, rho=28, beta=8/3",
         {"sigma", "rho", "beta"},
         {10.0, 28.0, 8.0 / 3.0},
         [](const std::vector<double>& p) {
             return Lorenz{p[0], p[1], p[2]}.ode();
         }},
        {"rossler",
         "Rossler (1976) funnel attractor; a=b=0.2, c=5.7",
         {"a", "b", "c"},
         {0.2, 0.2, 5.7},
         [](const std::vector<double>& p) { return Rossler{p[0], p[1], p[2]}.ode(); }},
        {"chua",
         "Chua circuit (double-scroll); alpha=15.6, beta=28",
         {"alpha", "beta", "m0", "m1"},
         {15.6, 28.0, -1.143, -0.714},
         [](const std::vector<double>& p) {
             return Chua{p[0], p[1], p[2], p[3]}.ode();
         }},
        {"double_pendulum",
         "Planar double pendulum; state=[theta1,theta2,omega1,omega2]",
         {"m1", "m2", "l1", "l2", "g"},
         {1.0, 1.0, 1.0, 1.0, 9.81},
         [](const std::vector<double>& p) {
             return DoublePendulum{p[0], p[1], p[2], p[3], p[4]}.ode();
         }},
        {"three_body",
         "Planar three-body problem; state=[r1,r2,r3,v1,v2,v3]",
         {"m1", "m2", "m3", "G"},
         {1.0, 1.0, 1.0, 1.0},
         [](const std::vector<double>& p) {
             return ThreeBody{p[0], p[1], p[2], p[3]}.ode();
         }},
    };
    return catalog;
}

const std::vector<MapSystemInfo>& map_catalog() {
    static const std::vector<MapSystemInfo> catalog = {
        {"logistic",
         "Logistic map x' = r x (1-x), r in [0,4]",
         {"r"},
         {4.0},
         [](const std::vector<double>& p) { return Logistic{p[0]}.map(); }},
        {"henon",
         "Henon map x' = 1-a x^2+y, y'=b x; a=1.4, b=0.3",
         {"a", "b"},
         {1.4, 0.3},
         [](const std::vector<double>& p) { return Henon{p[0], p[1]}.map(); }},
        {"cat",
         "Arnold cat map on the unit torus",
         {},
         {},
         [](const std::vector<double>&) { return CatMap{}.map(); }},
    };
    return catalog;
}

namespace {

std::vector<double> checked_params(const std::string& name,
                                   const std::vector<std::string>& names,
                                   const std::vector<double>& defaults,
                                   const std::vector<double>& given) {
    if (given.empty()) return defaults;
    if (given.size() != defaults.size()) {
        std::ostringstream oss;
        oss << name << ": expected " << defaults.size() << " parameters (";
        for (std::size_t i = 0; i < names.size(); ++i) {
            oss << (i ? "," : "") << names[i];
        }
        oss << "), got " << given.size();
        throw std::invalid_argument(oss.str());
    }
    return given;
}

}  // namespace

OdeSystem make_ode_system(const std::string& name, const std::vector<double>& params) {
    for (const auto& info : ode_catalog()) {
        if (info.name == name) {
            return info.build(checked_params(name, info.param_names, info.defaults, params));
        }
    }
    throw std::invalid_argument("unknown ODE system: " + name);
}

MapSystem make_map_system(const std::string& name, const std::vector<double>& params) {
    for (const auto& info : map_catalog()) {
        if (info.name == name) {
            return info.build(checked_params(name, info.param_names, info.defaults, params));
        }
    }
    throw std::invalid_argument("unknown map: " + name);
}

State default_initial_state(const std::string& name) {
    if (name == "logistic") return State{0.123456789};
    if (name == "henon") return State{0.0, 0.0};
    if (name == "cat") return State{0.1, 0.2};
    if (name == "lorenz") return State{1.0, 1.0, 1.0};
    if (name == "rossler") return State{1.0, 1.0, 1.0};
    if (name == "chua") return State{0.1, 0.0, 0.0};
    if (name == "double_pendulum") return State{1.0, 2.0, 0.0, 0.0};
    if (name == "three_body") return ThreeBody::figure_eight_initial_state();
    throw std::invalid_argument("no default initial state for " + name);
}

State parse_state(const std::string& text) {
    State out;
    std::stringstream ss(text);
    std::string token;
    while (std::getline(ss, token, ',')) {
        if (token.empty()) throw std::invalid_argument("empty coordinate in state: " + text);
        std::size_t used = 0;
        const double v = std::stod(token, &used);
        if (used != token.size()) throw std::invalid_argument("bad number in state: " + token);
        out.push_back(v);
    }
    if (out.empty()) throw std::invalid_argument("empty state");
    return out;
}

}  // namespace chaos
