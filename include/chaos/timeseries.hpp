// Spectral, recurrence and nonlinearity diagnostics for scalar time series.
// SPDX-License-Identifier: MIT
#pragma once

#include <cstdint>
#include <vector>

#include "chaos/embedding.hpp"
#include "chaos/rng.hpp"

namespace chaos {

struct Spectrum {
    std::vector<double> frequency;
    std::vector<double> power;
    double df = 0.0;
    std::size_t nfft = 0;
};

// One-sided periodogram.  The mean is removed and a Hann window is applied by
// default; the signal is zero-padded to the next power of two.
Spectrum power_spectrum(const std::vector<double>& series, double dt,
                        bool detrend = true, bool hann = true, bool zero_pad = true);

// Fourier phase-randomized surrogate (Theiler et al. 1992).  This preserves the
// power spectrum and the mean, but not the amplitude distribution; it is a
// first-line test against a linear Gaussian null hypothesis, not a full IAAFT.
std::vector<double> phase_randomized_surrogate(const std::vector<double>& series, Rng& rng);

struct ZeroOneResult {
    double K = 0.0;  // median over c
    std::vector<double> K_values;
};

// Gottwald-Melbourne 0-1 test.  K ~ 0 for regular dynamics and K ~ 1 for
// chaotic dynamics.  The scalar series is demeaned internally.
ZeroOneResult zero_one_test(const std::vector<double>& series, int n_c = 64,
                            std::uint64_t seed = 0x123456789abcdef0ULL);

// Recurrence quantification analysis (Eckmann et al. 1987; Marwan et al. 2007).
struct RqaResult {
    double epsilon = 0.0;
    double recurrence_rate = 0.0;
    double determinism = 0.0;
    double laminarity = 0.0;
    double entropy = 0.0;      // Shannon entropy of diagonal line lengths (nats)
    double trapping_time = 0.0;  // mean vertical line length
    double longest_diagonal = 0.0;
    double mean_diagonal = 0.0;
    std::size_t points = 0;
};

std::vector<std::uint8_t> recurrence_matrix(const Embedding& emb, double epsilon,
                                           int theiler = 0);

// Pick epsilon so that the off-diagonal recurrence rate is approximately `rate`.
double epsilon_for_recurrence_rate(const Embedding& emb, double rate, int theiler = 0);

RqaResult rqa_from_matrix(const std::vector<std::uint8_t>& matrix, std::size_t n,
                          double epsilon, int min_line = 2);

RqaResult rqa(const Embedding& emb, double epsilon, int theiler = 0, int min_line = 2);

}  // namespace chaos
