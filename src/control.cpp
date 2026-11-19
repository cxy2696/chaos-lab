// SPDX-License-Identifier: MIT
#include "chaos/control.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "chaos/integrators.hpp"

namespace chaos {

OgyResult ogy_control_1d(const std::function<double(double, double)>& map, double x0,
                         const OgyOptions& options) {
    if (!map) throw std::invalid_argument("ogy_control_1d: empty map");
    if (!(options.steps > 0)) throw std::invalid_argument("ogy: steps must be positive");
    OgyResult result;
    result.x.reserve(options.steps);
    result.p.reserve(options.steps);
    double x = x0;
    for (std::size_t n = 0; n < options.steps; ++n) {
        const double err = x - options.x_star;
        double p = options.p0;
        if (std::fabs(err) < options.tolerance && std::fabs(options.g) > 1e-300) {
            double dp = -options.lambda_u * err / options.g;
            dp = std::clamp(dp, -options.dp_max, options.dp_max);
            p += dp;
            ++result.activations;
        }
        x = map(x, p);
        result.x.push_back(x);
        result.p.push_back(p);
    }
    result.final_error = std::fabs(x - options.x_star);
    const std::size_t tail = std::max<std::size_t>(1, options.steps / 10);
    double mx = 0.0;
    for (std::size_t i = options.steps - tail; i < options.steps; ++i) {
        mx = std::max(mx, std::fabs(result.x[i] - options.x_star));
    }
    result.max_error_tail = mx;
    result.locked = result.max_error_tail < 1e-5;
    return result;
}

std::vector<double> pyragas_map_control(const std::function<double(double)>& map,
                                        double x0, int delay, double gain,
                                        std::size_t steps) {
    if (!map) throw std::invalid_argument("pyragas_map_control: empty map");
    if (delay < 1) throw std::invalid_argument("pyragas: delay must be >= 1");
    std::vector<double> x;
    x.reserve(steps + 1);
    x.push_back(x0);
    for (std::size_t n = 0; n < steps; ++n) {
        double next = map(x[n]);
        if (n >= static_cast<std::size_t>(delay)) {
            next += gain * (x[n - static_cast<std::size_t>(delay)] - x[n]);
        }
        x.push_back(next);
    }
    return x;
}

SyncResult pecora_carroll_lorenz(const Lorenz& lorenz, const State& drive0,
                                 const State& response0, double h, double t_transient,
                                 double t_total) {
    if (drive0.size() != 3 || response0.size() != 2) {
        throw std::invalid_argument("pecora_carroll_lorenz: need drive0(3), response0(2)");
    }
    if (!(h > 0.0) || !(t_total > t_transient)) {
        throw std::invalid_argument("pecora_carroll_lorenz: bad times");
    }
    State s = {drive0[0], drive0[1], drive0[2], response0[0], response0[1]};
    const Lorenz p = lorenz;
    auto rhs = [p](const State& y, double) {
        State d(5);
        d[0] = p.sigma * (y[1] - y[0]);
        d[1] = y[0] * (p.rho - y[2]) - y[1];
        d[2] = y[0] * y[1] - p.beta * y[2];
        d[3] = y[0] * (p.rho - y[4]) - y[3];
        d[4] = y[0] * y[3] - p.beta * y[4];
        return d;
    };

    SyncResult out;
    const std::size_t nsteps =
        static_cast<std::size_t>(std::llround(t_total / h));
    out.t.reserve(nsteps + 1);
    out.error.reserve(nsteps + 1);
    double t = 0.0;
    double settle_candidate = -1.0;
    bool settled = false;
    for (std::size_t i = 0; i <= nsteps; ++i) {
        const double err = std::hypot(s[3] - s[1], s[4] - s[2]);
        out.t.push_back(t);
        out.error.push_back(err);
        if (t >= t_transient) {
            out.max_error_after_transient =
                std::max(out.max_error_after_transient, err);
            if (err < 1e-6) {
                if (!settled) {
                    settled = true;
                    settle_candidate = t;
                }
            } else {
                settled = false;
                settle_candidate = -1.0;
            }
        }
        if (i == nsteps) break;
        s = rk4_step_t(rhs, s, t, h);
        t += h;
        if (!is_finite_state(s)) break;
    }
    out.final_error = out.error.empty() ? 0.0 : out.error.back();
    out.synchronized = out.final_error < 1e-6;
    out.settling_time = settled ? settle_candidate : -1.0;
    return out;
}

}  // namespace chaos
