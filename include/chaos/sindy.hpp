// Sparse identification of nonlinear dynamics (SINDy), Brunton et al. (2016).
// A compact STLSQ implementation with a polynomial/trigonometric library.
// SPDX-License-Identifier: MIT
#pragma once

#include <functional>
#include <string>
#include <vector>

#include "chaos/core.hpp"

namespace chaos {

struct SindyOptions {
    int degree = 2;
    bool include_constant = true;
    bool include_sine = false;
    bool include_cosine = false;
    double threshold = 0.1;
    double ridge = 1e-10;
    int max_iter = 20;
};

struct SindyTerm {
    std::string name;
    std::function<double(const State&)> evaluate;
};

struct SindyModel {
    std::vector<SindyTerm> terms;
    Matrix coefficients;  // terms x state dimension
    State predict_derivative(const State& x) const;
    std::string pretty() const;
};

// X: samples x state dimension, dX: time derivatives (same shape).
SindyModel sindy_fit(const std::vector<State>& X, const std::vector<State>& dX,
                     const SindyOptions& options = {});

// Central differences from evenly sampled trajectory data (sample interval dt).
std::vector<State> finite_difference_derivatives(const std::vector<State>& X, double dt);

}  // namespace chaos
