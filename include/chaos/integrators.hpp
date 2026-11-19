// Fixed and adaptive integrators.
// SPDX-License-Identifier: MIT
#pragma once

#include <cstddef>
#include <functional>
#include <string>
#include <vector>

#include "chaos/core.hpp"

namespace chaos {

// Classic fourth-order Runge-Kutta.  Template version avoids std::function
// overhead for tight loops (e.g. Poincare sections).
template <class RHS>
State rk4_step_t(const RHS& rhs, const State& x, double t, double h) {
    const std::size_t n = x.size();
    const State k1 = rhs(x, t);
    State tmp(n);
    for (std::size_t i = 0; i < n; ++i) tmp[i] = x[i] + 0.5 * h * k1[i];
    const State k2 = rhs(tmp, t + 0.5 * h);
    for (std::size_t i = 0; i < n; ++i) tmp[i] = x[i] + 0.5 * h * k2[i];
    const State k3 = rhs(tmp, t + 0.5 * h);
    for (std::size_t i = 0; i < n; ++i) tmp[i] = x[i] + h * k3[i];
    const State k4 = rhs(tmp, t + h);

    State out(n);
    for (std::size_t i = 0; i < n; ++i) {
        out[i] = x[i] + h / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    }
    return out;
}

State rk4_step(const OdeSystem& sys, const State& x, double t, double h);

// Dormand-Prince 5(4) adaptive integrator.
struct Rk45Options {
    double rtol = 1e-8;
    double atol = 1e-10;
    double h_min = 1e-14;
    double h_max = 0.0;  // 0 => no explicit cap beyond the remaining interval
    std::size_t max_steps = 20000000;
};

struct Rk45Result {
    bool ok = false;
    std::string message;
    State x;
    double t = 0.0;
    double h_last = 0.0;
    std::size_t accepted = 0;
    std::size_t rejected = 0;
};

Rk45Result rk45_integrate(const OdeSystem& sys, State x, double t0, double t1,
                          double h0, const Rk45Options& options = {});

// Velocity Verlet for q'' = a(q).  Second order, symplectic for separable
// Hamiltonians; intended for double pendulums, N-body problems, etc.
using Accel = std::function<State(const State&)>;
void verlet_step(const Accel& accel, State& q, State& v, double h);

// Fourth-order Yoshida composition of velocity Verlet steps.
void yoshida4_step(const Accel& accel, State& q, State& v, double h);

// Energy of a unit-mass separable Hamiltonian H = |v|^2/2 + V(q), where
// `potential` evaluates V(q).  Handy in tests and examples.
double hamiltonian_energy(const State& q, const State& v,
                          const std::function<double(const State&)>& potential);

}  // namespace chaos
