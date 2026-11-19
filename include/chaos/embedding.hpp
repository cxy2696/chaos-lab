// Delay-coordinate embedding and nonlinear time-series diagnostics.
// SPDX-License-Identifier: MIT
#pragma once

#include <cstddef>
#include <vector>

namespace chaos {

struct Embedding {
    int m = 1;
    int tau = 1;
    std::size_t n = 0;
    std::vector<double> data;  // n x m, row-major

    double at(std::size_t i, int j) const {
        return data[i * static_cast<std::size_t>(m) + static_cast<std::size_t>(j)];
    }
    const double* row(std::size_t i) const {
        return data.data() + i * static_cast<std::size_t>(m);
    }
};

Embedding delay_embed(const std::vector<double>& series, int m, int tau);

// Fraser & Swinney (1986) average mutual information, estimated with
// equiprobable-in-width bins.  Returns values for tau = 0..max_tau.
std::vector<double> average_mutual_information(const std::vector<double>& series,
                                               int max_tau, int bins = 32);

// First interior local minimum of the AMI curve, or -1 if none.
int first_ami_minimum(const std::vector<double>& ami);

// Normalized autocorrelation for lags 0..max_lag.
std::vector<double> autocorrelation(const std::vector<double>& series, int max_lag);

// First lag at which |autocorrelation| drops below 1/e; -1 if none.
int autocorrelation_1e(const std::vector<double>& ac);

struct FnnResult {
    int suggested_m = 1;
    std::vector<double> false_ratio;  // index m: fraction of false neighbours
    std::vector<int> valid_points;
};

// Kennel, Brown & Abarbanel (1992) false-nearest-neighbours method.
// The threshold is on |x_{i+m tau} - x_{j+m tau}| / ||X_i - X_j||_m.
FnnResult false_nearest_neighbors(const std::vector<double>& series, int max_m, int tau,
                                  int theiler = 0, double rtol = 15.0,
                                  double atol = 0.0);

struct CorrelationResult {
    double d2 = 0.0;
    double r2 = 0.0;
    double r_min = 0.0;
    double r_max = 0.0;
    std::vector<double> radii;
    std::vector<double> log_r;
    std::vector<double> log_c;
    std::vector<double> local_slope;
    std::size_t fit_begin = 0;
    std::size_t fit_end = 0;
};

// Grassberger-Procaccia (1983) correlation dimension.  If r_min/r_max are
// non-positive they are selected from the distribution of pairwise distances.
// The slope is fitted over the middle part of the scaling range.
CorrelationResult correlation_dimension(const Embedding& emb, double r_min = 0.0,
                                        double r_max = 0.0, int n_r = 40,
                                        int theiler = 0, double fit_lo = 0.25,
                                        double fit_hi = 0.75);

CorrelationResult correlation_dimension(const std::vector<double>& series, int m, int tau,
                                        double r_min = 0.0, double r_max = 0.0,
                                        int n_r = 40, int theiler = 0,
                                        double fit_lo = 0.25, double fit_hi = 0.75);

}  // namespace chaos
