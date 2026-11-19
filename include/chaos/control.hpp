// Chaos control and synchronization.
// SPDX-License-Identifier: MIT
#pragma once

#include <cstddef>
#include <functional>
#include <vector>

#include "chaos/core.hpp"
#include "chaos/systems.hpp"

namespace chaos {

struct OgyOptions {
    double x_star = 0.0;       // target fixed point at p0
    double p0 = 0.0;           // nominal parameter
    double lambda_u = 0.0;     // unstable multiplier of the map at (x_star, p0)
    double g = 0.0;            // d f / d p at (x_star, p0)
    double tolerance = 0.1;    // only control inside this neighbourhood
    double dp_max = 0.05;      // actuator limit
    std::size_t steps = 100000;
};

struct OgyResult {
    bool locked = false;
    std::size_t activations = 0;
    double final_error = 0.0;
    double max_error_tail = 0.0;
    std::vector<double> x;
    std::vector<double> p;
};

// One-dimensional OGY targeting: delta p = -lambda_u (x-x*) / g.
// The perturbation is recomputed every step and reset to p0 otherwise.
OgyResult ogy_control_1d(const std::function<double(double, double)>& map, double x0,
                         const OgyOptions& options);

// Pyragas delayed feedback for a map: x_{n+1} = f(x_n) + K (x_{n-d} - x_n).
std::vector<double> pyragas_map_control(const std::function<double(double)>& map,
                                        double x0, int delay, double gain,
                                        std::size_t steps);

struct SyncResult {
    std::vector<double> t;
    std::vector<double> error;
    double final_error = 0.0;
    double max_error_after_transient = 0.0;
    double settling_time = -1.0;
    bool synchronized = false;
};

// Pecora-Carroll drive-response synchronization of two identical Lorenz
// systems.  x drives the (y,z) slave; the slave starts elsewhere.
SyncResult pecora_carroll_lorenz(const Lorenz& lorenz, const State& drive0,
                                 const State& response0, double h = 0.005,
                                 double t_transient = 5.0, double t_total = 40.0);

}  // namespace chaos
