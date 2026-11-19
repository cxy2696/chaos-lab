// Bifurcation scans and Poincare sections.
// SPDX-License-Identifier: MIT
#pragma once

#include <cstddef>
#include <functional>
#include <vector>

#include "chaos/core.hpp"

namespace chaos {

struct BifurcationPoint {
    double parameter = 0.0;
    double value = 0.0;
};

// Bifurcation diagram for a one-dimensional parametric map f(x, p).
// `keep` late-iterate values are recorded per parameter value (stride applies).
std::vector<BifurcationPoint> bifurcation_scan_1d(
    const std::function<double(double, double)>& map, double x0, double p_min,
    double p_max, int n_params, std::size_t transient, std::size_t keep,
    std::size_t stride = 1);

// Poincare section of a flow on the hyperplane x[coordinate] = value.
// direction: +1 upward crossings, -1 downward crossings, 0 both.
// Crossings are located by linear interpolation within one RK4 step.
struct PoincareOptions {
    int coordinate = 2;
    double value = 0.0;
    int direction = 0;
    double h = 0.005;
    double t_transient = 0.0;
    double t_max = 1000.0;
    std::size_t max_points = 20000;
};

std::vector<State> poincare_section(const OdeSystem& sys, State x0,
                                    const PoincareOptions& options);

}  // namespace chaos
