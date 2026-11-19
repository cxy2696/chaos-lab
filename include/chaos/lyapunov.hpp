// Lyapunov exponents: Benettin/Gram-Schmidt for flows and maps, two-trajectory
// max exponent, and Rosenstein's method for scalar time series.
// SPDX-License-Identifier: MIT
#pragma once

#include <cstddef>
#include <vector>

#include "chaos/core.hpp"

namespace chaos {

struct LyapunovOptions {
    double h = 0.01;
    double transient = 100.0;   // discarded before tangent-space averaging
    double duration = 1000.0;   // averaging window
    std::size_t reorthonormalize_every = 1;
    int n_exponents = 0;        // 0 => full spectrum
};

struct LyapunovSpectrum {
    std::vector<double> exponents;
    double integration_time = 0.0;
    std::size_t steps = 0;
    std::size_t reorthonormalizations = 0;
};

// Benettin et al. (1980): integrate the state and a tangent basis, periodically
// re-orthonormalize with modified Gram-Schmidt and average log stretching
// factors.  If the system has no analytic Jacobian, central differences are used.
LyapunovSpectrum lyapunov_spectrum(const OdeSystem& sys, State x0,
                                   const LyapunovOptions& options = {});

// Discrete analogue, exponents per map iteration.
LyapunovSpectrum lyapunov_spectrum_map(const MapSystem& map, State x0,
                                       std::size_t transient, std::size_t n,
                                       std::size_t reorthonormalize_every = 1,
                                       int n_exponents = 0);

// Largest exponent from two nearby trajectories with periodic rescaling.
// Returns a per-unit-time exponent.
double lyapunov_max_two_trajectories(const OdeSystem& sys, State x0, double d0,
                                     double h, double transient, double duration);

struct RosensteinResult {
    double lambda = 0.0;          // per unit time (slope of ln divergence / dt)
    double intercept = 0.0;
    double r2 = 0.0;
    double dt = 1.0;
    std::size_t fit_begin = 0;    // index into k/divergence
    std::size_t fit_end = 0;
    std::vector<double> k;        // embedding-step indices
    std::vector<double> divergence;  // <ln d_i(k)>
};

// Rosenstein, Collins & De Luca (1994), for a scalar, evenly sampled series.
// Theiler window and tau are measured in samples; dt is the sample interval.
// fit_begin/fit_end select the linear part of the divergence curve; leave
// both negative to use the automatic plateau search.
RosensteinResult rosenstein_lyapunov(const std::vector<double>& series, int m,
                                     int tau, int theiler, int max_k,
                                     double dt = 1.0, int fit_begin = -1,
                                     int fit_end = -1);

}  // namespace chaos
