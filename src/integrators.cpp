// SPDX-License-Identifier: MIT
#include "chaos/integrators.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace chaos {

State rk4_step(const OdeSystem& sys, const State& x, double t, double h) {
    if (!sys.rhs) throw std::invalid_argument("rk4_step: OdeSystem has no rhs");
    require_dim(x, sys.dim, "rk4_step");
    return rk4_step_t(sys.rhs, x, t, h);
}

namespace {

struct DpStep {
    State x5;
    State error;
    bool ok = true;
};

// Dormand-Prince 5(4) coefficients (Dormand & Prince 1980).
DpStep dormand_prince_step(const OdeSystem& sys, const State& x, double t, double h) {
    const std::size_t n = x.size();
    const State k1 = sys.rhs(x, t);
    State tmp(n);
    for (std::size_t i = 0; i < n; ++i) tmp[i] = x[i] + h * (1.0 / 5.0) * k1[i];
    const State k2 = sys.rhs(tmp, t + h * (1.0 / 5.0));

    for (std::size_t i = 0; i < n; ++i) {
        tmp[i] = x[i] + h * ((3.0 / 40.0) * k1[i] + (9.0 / 40.0) * k2[i]);
    }
    const State k3 = sys.rhs(tmp, t + h * (3.0 / 10.0));

    for (std::size_t i = 0; i < n; ++i) {
        tmp[i] = x[i] + h * ((44.0 / 45.0) * k1[i] - (56.0 / 15.0) * k2[i] +
                             (32.0 / 9.0) * k3[i]);
    }
    const State k4 = sys.rhs(tmp, t + h * (4.0 / 5.0));

    for (std::size_t i = 0; i < n; ++i) {
        tmp[i] = x[i] + h * ((19372.0 / 6561.0) * k1[i] - (25360.0 / 2187.0) * k2[i] +
                             (64448.0 / 6561.0) * k3[i] - (212.0 / 729.0) * k4[i]);
    }
    const State k5 = sys.rhs(tmp, t + h * (8.0 / 9.0));

    for (std::size_t i = 0; i < n; ++i) {
        tmp[i] = x[i] + h * ((9017.0 / 3168.0) * k1[i] - (355.0 / 33.0) * k2[i] +
                             (46732.0 / 5247.0) * k3[i] + (49.0 / 176.0) * k4[i] -
                             (5103.0 / 18656.0) * k5[i]);
    }
    const State k6 = sys.rhs(tmp, t + h);

    for (std::size_t i = 0; i < n; ++i) {
        tmp[i] = x[i] + h * ((35.0 / 384.0) * k1[i] + (500.0 / 1113.0) * k3[i] +
                             (125.0 / 192.0) * k4[i] - (2187.0 / 6784.0) * k5[i] +
                             (11.0 / 84.0) * k6[i]);
    }
    const State k7 = sys.rhs(tmp, t + h);

    DpStep out;
    out.x5.resize(n);
    out.error.resize(n);
    const double b5[7] = {35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0,
                          -2187.0 / 6784.0, 11.0 / 84.0, 0.0};
    const double b4[7] = {5179.0 / 57600.0, 0.0, 7571.0 / 16695.0, 393.0 / 640.0,
                          -92097.0 / 339200.0, 187.0 / 2100.0, 1.0 / 40.0};
    const State* ks[7] = {&k1, &k2, &k3, &k4, &k5, &k6, &k7};
    for (std::size_t i = 0; i < n; ++i) {
        double hi = 0.0;
        double lo = 0.0;
        for (int j = 0; j < 7; ++j) {
            hi += b5[j] * (*ks[j])[i];
            lo += b4[j] * (*ks[j])[i];
        }
        out.x5[i] = x[i] + h * hi;
        out.error[i] = h * (hi - lo);
    }
    return out;
}

}  // namespace

Rk45Result rk45_integrate(const OdeSystem& sys, State x, double t0, double t1,
                          double h0, const Rk45Options& options) {
    Rk45Result result;
    if (!sys.rhs) {
        result.message = "rk45_integrate: OdeSystem has no rhs";
        return result;
    }
    if (sys.dim <= 0 || static_cast<int>(x.size()) != sys.dim) {
        result.message = "rk45_integrate: state dimension mismatch";
        return result;
    }
    if (!(t1 > t0)) {
        result.message = "rk45_integrate: require t1 > t0";
        return result;
    }
    if (!(h0 > 0.0)) {
        result.message = "rk45_integrate: require h0 > 0";
        return result;
    }

    double t = t0;
    double h = std::min(h0, t1 - t0);
    const double span = t1 - t0;
    h = std::min(h, span);
    const double h_max = options.h_max > 0.0 ? std::min(options.h_max, span) : span;

    while (t < t1 - 1e-15 * std::max(1.0, std::fabs(t1))) {
        if (result.accepted + result.rejected > options.max_steps) {
            result.message = "rk45_integrate: step limit exceeded";
            result.x = x;
            result.t = t;
            return result;
        }
        if (h > t1 - t) h = t1 - t;
        if (h < options.h_min) {
            result.message = "rk45_integrate: step size underflow";
            result.x = x;
            result.t = t;
            return result;
        }

        const DpStep step = dormand_prince_step(sys, x, t, h);
        double err = 0.0;
        for (std::size_t i = 0; i < x.size(); ++i) {
            const double sc = options.atol + options.rtol * std::max(std::fabs(x[i]),
                                                                     std::fabs(step.x5[i]));
            const double e = step.error[i] / (sc > 0.0 ? sc : 1.0);
            err += e * e;
        }
        err = std::sqrt(err / static_cast<double>(x.size()));

        if (err <= 1.0) {
            t += h;
            x = step.x5;
            ++result.accepted;
            result.h_last = h;
            if (!is_finite_state(x)) {
                result.message = "rk45_integrate: non-finite state";
                result.x = x;
                result.t = t;
                return result;
            }
        } else {
            ++result.rejected;
        }

        double factor = 0.9 * std::pow(err > 1e-16 ? 1.0 / err : 1e16, 0.2);
        factor = std::clamp(factor, 0.2, 5.0);
        if (err > 1.0) factor = std::min(factor, 0.5);
        h *= factor;
        if (h > h_max) h = h_max;
        if (h < options.h_min && err <= 1.0) {
            // The remaining interval is smaller than h_min; finish exactly.
            if (t1 - t <= options.h_min) break;
        }
    }

    if (t < t1) {
        // Snap to the endpoint when the residual is at round-off scale.
        if (t1 - t > 1e-12 * std::max(1.0, span)) {
            result.message = "rk45_integrate: did not reach t1";
            result.x = x;
            result.t = t;
            return result;
        }
        t = t1;
    }
    result.ok = true;
    result.x = std::move(x);
    result.t = t;
    result.message = "ok";
    return result;
}

void verlet_step(const Accel& accel, State& q, State& v, double h) {
    if (!accel) throw std::invalid_argument("verlet_step: accel is empty");
    if (q.size() != v.size()) throw std::invalid_argument("verlet_step: q/v size mismatch");
    const State a0 = accel(q);
    for (std::size_t i = 0; i < q.size(); ++i) {
        v[i] += 0.5 * h * a0[i];
        q[i] += h * v[i];
    }
    const State a1 = accel(q);
    for (std::size_t i = 0; i < q.size(); ++i) v[i] += 0.5 * h * a1[i];
}

void yoshida4_step(const Accel& accel, State& q, State& v, double h) {
    // Yoshida (1990): S4(h) = S2(w1 h) S2(w0 h) S2(w1 h)
    const double cbrt2 = std::cbrt(2.0);
    const double w1 = 1.0 / (2.0 - cbrt2);
    const double w0 = -cbrt2 * w1;
    verlet_step(accel, q, v, w1 * h);
    verlet_step(accel, q, v, w0 * h);
    verlet_step(accel, q, v, w1 * h);
}

double hamiltonian_energy(const State& q, const State& v,
                          const std::function<double(const State&)>& potential) {
    double kinetic = 0.0;
    for (double x : v) kinetic += 0.5 * x * x;
    return kinetic + (potential ? potential(q) : 0.0);
}

}  // namespace chaos
