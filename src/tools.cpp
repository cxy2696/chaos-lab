// SPDX-License-Identifier: MIT
#include "chaos/tools.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "chaos/integrators.hpp"

namespace chaos {

std::vector<BifurcationPoint> bifurcation_scan_1d(
    const std::function<double(double, double)>& map, double x0, double p_min,
    double p_max, int n_params, std::size_t transient, std::size_t keep,
    std::size_t stride) {
    if (!map) throw std::invalid_argument("bifurcation_scan_1d: empty map");
    if (n_params < 2 || !(p_max > p_min)) {
        throw std::invalid_argument("bifurcation_scan_1d: bad parameter range");
    }
    stride = std::max<std::size_t>(1, stride);
    std::vector<BifurcationPoint> out;
    out.reserve(static_cast<std::size_t>(n_params) * keep / stride + 1);
    for (int k = 0; k < n_params; ++k) {
        const double p = p_min + (p_max - p_min) * static_cast<double>(k) /
                                     static_cast<double>(n_params - 1);
        double x = x0;
        for (std::size_t i = 0; i < transient; ++i) {
            x = map(x, p);
            if (!std::isfinite(x)) break;
        }
        if (!std::isfinite(x)) continue;
        for (std::size_t i = 0; i < keep; ++i) {
            x = map(x, p);
            if (!std::isfinite(x)) break;
            if (i % stride == 0) out.push_back({p, x});
        }
    }
    return out;
}

std::vector<State> poincare_section(const OdeSystem& sys, State x0,
                                    const PoincareOptions& options) {
    if (!(options.h > 0.0)) throw std::invalid_argument("poincare_section: h must be positive");
    if (!sys.rhs || sys.dim <= 0) throw std::invalid_argument("poincare_section: invalid system");
    if (options.coordinate < 0 || options.coordinate >= sys.dim) {
        throw std::invalid_argument("poincare_section: coordinate out of range");
    }
    if (options.direction < -1 || options.direction > 1) {
        throw std::invalid_argument("poincare_section: direction must be -1, 0 or 1");
    }
    require_dim(x0, sys.dim, "poincare_section");
    std::vector<State> out;
    State prev = std::move(x0);
    double t = 0.0;
    double g_prev = prev[static_cast<std::size_t>(options.coordinate)] - options.value;
    while (t < options.t_max && out.size() < options.max_points) {
        State next = rk4_step(sys, prev, t, options.h);
        const double t_next = t + options.h;
        const double g_next = next[static_cast<std::size_t>(options.coordinate)] - options.value;
        bool crossing = false;
        if (options.direction > 0) {
            crossing = (g_prev < 0.0 && g_next >= 0.0);
        } else if (options.direction < 0) {
            crossing = (g_prev > 0.0 && g_next <= 0.0);
        } else {
            crossing = ((g_prev < 0.0 && g_next >= 0.0) || (g_prev > 0.0 && g_next <= 0.0));
        }
        if (crossing && t_next >= options.t_transient) {
            const double denom = g_prev - g_next;
            const double frac = std::fabs(denom) > 0.0 ? g_prev / denom : 0.0;
            State cross(prev.size());
            for (std::size_t i = 0; i < prev.size(); ++i) {
                cross[i] = prev[i] + frac * (next[i] - prev[i]);
            }
            cross[static_cast<std::size_t>(options.coordinate)] = options.value;
            out.push_back(std::move(cross));
        }
        prev = std::move(next);
        t = t_next;
        g_prev = g_next;
        if (!is_finite_state(prev)) break;
    }
    return out;
}

}  // namespace chaos
