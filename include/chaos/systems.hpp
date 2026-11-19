// Classical chaotic systems and maps used in the library, examples and tests.
// SPDX-License-Identifier: MIT
#pragma once

#include <array>
#include <functional>
#include <string>
#include <vector>

#include "chaos/core.hpp"

namespace chaos {

// ---------------------------------------------------------------------------
// One-dimensional maps
// ---------------------------------------------------------------------------

// x_{n+1} = r x_n (1 - x_n), r in [0, 4].
struct Logistic {
    double r = 4.0;
    static constexpr int dimension = 1;

    double step(double x) const { return r * x * (1.0 - x); }
    double derivative(double x) const { return r * (1.0 - 2.0 * x); }
    // Exact Lyapunov exponent of an orbit, averaged in log space.
    double lyapunov(double x0 = 0.123456789, std::size_t transient = 1000,
                    std::size_t n = 100000) const;
    MapSystem map() const;
};

// x_{n+1} = 1 - a x_n^2 + y_n, y_{n+1} = b x_n.
struct Henon {
    double a = 1.4;
    double b = 0.3;
    static constexpr int dimension = 2;

    State step(const State& s) const;
    Matrix jacobian(const State& s) const;
    MapSystem map() const;
};

// Arnold's cat map on the unit torus: (x, y) -> (x+y, x+2y) mod 1.
struct CatMap {
    static constexpr int dimension = 2;
    State step(const State& s) const;
    Matrix jacobian(const State&) const;
    MapSystem map() const;
};

// ---------------------------------------------------------------------------
// Continuous-time flows
// ---------------------------------------------------------------------------

// Lorenz (1963): x' = sigma(y-x), y' = x(rho-z)-y, z' = xy-beta z.
struct Lorenz {
    double sigma = 10.0;
    double rho = 28.0;
    double beta = 8.0 / 3.0;
    static constexpr int dimension = 3;

    State f(const State& s, double t = 0.0) const;
    Matrix jacobian(const State& s, double t = 0.0) const;
    double divergence() const { return -(sigma + 1.0 + beta); }
    OdeSystem ode() const;
};

// Rossler (1976): x' = -y-z, y' = x+a y, z' = b + z(x-c).
struct Rossler {
    double a = 0.2;
    double b = 0.2;
    double c = 5.7;
    static constexpr int dimension = 3;

    State f(const State& s, double t = 0.0) const;
    Matrix jacobian(const State& s, double t = 0.0) const;
    OdeSystem ode() const;
};

// Chua's circuit in dimensionless form (Chua, Komuro & Matsumoto 1986):
// x' = alpha(y - x - f(x)), y' = x - y + z, z' = -beta y,
// f(x) = m1 x + (m0-m1)(|x+1|-|x-1|)/2.
struct Chua {
    double alpha = 15.6;
    double beta = 28.0;
    double m0 = -1.143;
    double m1 = -0.714;
    static constexpr int dimension = 3;

    double nonlinearity(double x) const;
    double nonlinearity_derivative(double x) const;
    State f(const State& s, double t = 0.0) const;
    Matrix jacobian(const State& s, double t = 0.0) const;
    OdeSystem ode() const;
};

// Planar double pendulum.  State = [theta1, theta2, omega1, omega2].
struct DoublePendulum {
    double m1 = 1.0;
    double m2 = 1.0;
    double l1 = 1.0;
    double l2 = 1.0;
    double g = 9.81;
    static constexpr int dimension = 4;

    State f(const State& s, double t = 0.0) const;
    OdeSystem ode() const;
    double energy(const State& s) const;
};

// Planar equal-mass three-body problem (12-dimensional state).
// State = [x1,y1,x2,y2,x3,y3, vx1,vy1,vx2,vy2,vx3,vy3].
struct ThreeBody {
    double m1 = 1.0;
    double m2 = 1.0;
    double m3 = 1.0;
    double G = 1.0;
    static constexpr int dimension = 12;

    State f(const State& s, double t = 0.0) const;
    OdeSystem ode() const;
    double energy(const State& s) const;
    // Chenciner-Montgomery figure-eight choreography (period 6.32591398...).
    static State figure_eight_initial_state();
};

// ---------------------------------------------------------------------------
// Catalogue used by the command-line tool
// ---------------------------------------------------------------------------

struct OdeSystemInfo {
    std::string name;
    std::string summary;
    std::vector<std::string> param_names;
    std::vector<double> defaults;
    std::function<OdeSystem(const std::vector<double>&)> build;
};

struct MapSystemInfo {
    std::string name;
    std::string summary;
    std::vector<std::string> param_names;
    std::vector<double> defaults;
    std::function<MapSystem(const std::vector<double>&)> build;
};

const std::vector<OdeSystemInfo>& ode_catalog();
const std::vector<MapSystemInfo>& map_catalog();

OdeSystem make_ode_system(const std::string& name, const std::vector<double>& params = {});
MapSystem make_map_system(const std::string& name, const std::vector<double>& params = {});
State default_initial_state(const std::string& name);
State parse_state(const std::string& text);

}  // namespace chaos
