// SPDX-License-Identifier: MIT
#include "chaos/sindy.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <sstream>
#include <stdexcept>

namespace chaos {

namespace {

struct Term {
    std::string name;
    std::vector<int> powers;
    bool sine = false;
    bool cosine = false;
    int trig_index = 0;
};

double eval_term(const Term& t, const State& x) {
    if (t.sine) return std::sin(x[static_cast<std::size_t>(t.trig_index)]);
    if (t.cosine) return std::cos(x[static_cast<std::size_t>(t.trig_index)]);
    double v = 1.0;
    for (std::size_t i = 0; i < t.powers.size(); ++i) {
        if (t.powers[i] == 0) continue;
        v *= std::pow(x[i], t.powers[i]);
    }
    return v;
}

std::string power_name(int var) { return "x" + std::to_string(var); }

void add_monomials(int dim, int degree, int var, int remaining,
                   std::vector<int>& powers, std::vector<Term>& terms) {
    if (var == dim) {
        if (remaining != 0) return;
        std::string name;
        for (int v = 0; v < dim; ++v) {
            for (int k = 0; k < powers[static_cast<std::size_t>(v)]; ++k) {
                if (!name.empty()) name += "*";
                name += power_name(v);
            }
        }
        terms.push_back({name, powers, false, false, 0});
        return;
    }
    for (int p = 0; p <= remaining; ++p) {
        powers[static_cast<std::size_t>(var)] = p;
        add_monomials(dim, degree, var + 1, remaining - p, powers, terms);
    }
    powers[static_cast<std::size_t>(var)] = 0;
}

std::vector<Term> make_library(int dim, const SindyOptions& opt) {
    std::vector<Term> terms;
    if (opt.include_constant) terms.push_back({"1", std::vector<int>(dim, 0), false, false, 0});
    for (int v = 0; v < dim; ++v) {
        std::vector<int> p(dim, 0);
        p[static_cast<std::size_t>(v)] = 1;
        terms.push_back({power_name(v), p, false, false, 0});
    }
    if (opt.degree >= 2) {
        for (int i = 0; i < dim; ++i) {
            for (int j = i; j < dim; ++j) {
                std::vector<int> p(dim, 0);
                p[static_cast<std::size_t>(i)] += 1;
                p[static_cast<std::size_t>(j)] += 1;
                std::string name = power_name(i);
                if (i == j) {
                    name += "^2";
                } else {
                    name += "*" + power_name(j);
                }
                terms.push_back({name, p, false, false, 0});
            }
        }
    }
    if (opt.degree >= 3) {
        std::vector<int> powers(dim, 0);
        for (int degree = 3; degree <= opt.degree; ++degree) {
            add_monomials(dim, degree, 0, degree, powers, terms);
        }
    }
    if (opt.include_sine) {
        for (int v = 0; v < dim; ++v) {
            terms.push_back({"sin(" + power_name(v) + ")", std::vector<int>(dim, 0), true,
                             false, v});
        }
    }
    if (opt.include_cosine) {
        for (int v = 0; v < dim; ++v) {
            terms.push_back({"cos(" + power_name(v) + ")", std::vector<int>(dim, 0), false,
                             true, v});
        }
    }
    return terms;
}

}  // namespace

State SindyModel::predict_derivative(const State& x) const {
    State out(coefficients.cols, 0.0);
    for (std::size_t i = 0; i < terms.size(); ++i) {
        const double v = terms[i].evaluate(x);
        for (std::size_t j = 0; j < out.size(); ++j) out[j] += coefficients(i, j) * v;
    }
    return out;
}

std::string SindyModel::pretty() const {
    std::ostringstream oss;
    for (std::size_t j = 0; j < coefficients.cols; ++j) {
        oss << "d x" << j << "/dt = ";
        bool first = true;
        for (std::size_t i = 0; i < terms.size(); ++i) {
            const double c = coefficients(i, j);
            if (std::fabs(c) < 1e-12) continue;
            if (!first) oss << (c >= 0 ? " + " : " - ");
            if (!first) {
                oss << std::fabs(c) << "*" << terms[i].name;
            } else {
                oss << c << "*" << terms[i].name;
                first = false;
            }
            first = false;
        }
        if (first) oss << "0";
        oss << "\n";
    }
    return oss.str();
}

SindyModel sindy_fit(const std::vector<State>& X, const std::vector<State>& dX,
                     const SindyOptions& options) {
    if (X.empty() || X.size() != dX.size()) {
        throw std::invalid_argument("sindy_fit: X/dX size mismatch");
    }
    const std::size_t n = X.size();
    const std::size_t dim = X[0].size();
    if (dim == 0 || dX[0].size() != dim) {
        throw std::invalid_argument("sindy_fit: bad state dimension");
    }
    if (n < 10) throw std::invalid_argument("sindy_fit: need at least 10 samples");

    std::vector<Term> terms = make_library(static_cast<int>(dim), options);
    const std::size_t p = terms.size();
    Matrix theta(n, p);
    for (std::size_t i = 0; i < n; ++i) {
        if (X[i].size() != dim || dX[i].size() != dim) {
            throw std::invalid_argument("sindy_fit: inconsistent state dimension");
        }
        for (std::size_t j = 0; j < p; ++j) theta(i, j) = eval_term(terms[j], X[i]);
    }

    Matrix coef(p, dim, 0.0);
    std::vector<double> scale(p, 1.0);
    for (std::size_t j = 0; j < p; ++j) {
        double s = 0.0;
        for (std::size_t i = 0; i < n; ++i) s += theta(i, j) * theta(i, j);
        if (s > 1e-24) scale[j] = std::sqrt(s);
    }

    for (std::size_t d = 0; d < dim; ++d) {
        std::vector<std::size_t> active;
        for (std::size_t j = 0; j < p; ++j) {
            if (scale[j] > 1e-12) active.push_back(j);
        }
        std::vector<double> c(p, 0.0);
        for (int iter = 0; iter < options.max_iter; ++iter) {
            const std::size_t q = active.size();
            if (q == 0) break;
            Matrix A(q, q, 0.0);
            std::vector<double> b(q, 0.0);
            for (std::size_t a = 0; a < q; ++a) {
                const std::size_t ja = active[a];
                for (std::size_t bb = 0; bb < q; ++bb) {
                    const std::size_t jb = active[bb];
                    double s = 0.0;
                    for (std::size_t i = 0; i < n; ++i) {
                        s += (theta(i, ja) / scale[ja]) * (theta(i, jb) / scale[jb]);
                    }
                    A(a, bb) = s;
                }
                A(a, a) += options.ridge;
                double s = 0.0;
                for (std::size_t i = 0; i < n; ++i) {
                    s += (theta(i, ja) / scale[ja]) * dX[i][d];
                }
                b[a] = s;
            }
            const std::vector<double> w = solve_linear(A, b);
            std::fill(c.begin(), c.end(), 0.0);
            for (std::size_t a = 0; a < q; ++a) {
                c[active[a]] = w[a] / scale[active[a]];
            }
            std::vector<std::size_t> next;
            for (std::size_t j : active) {
                if (std::fabs(c[j]) >= options.threshold) next.push_back(j);
            }
            if (next.size() == active.size()) break;
            active.swap(next);
        }
        for (std::size_t j = 0; j < p; ++j) coef(j, d) = c[j];
    }

    SindyModel model;
    model.coefficients = std::move(coef);
    for (const Term& t : terms) {
        model.terms.push_back(
            {t.name, [t](const State& x) { return eval_term(t, x); }});
    }
    return model;
}

std::vector<State> finite_difference_derivatives(const std::vector<State>& X, double dt) {
    if (X.size() < 3 || !(dt > 0.0)) {
        throw std::invalid_argument("finite_difference_derivatives: need >=3 samples, dt>0");
    }
    std::vector<State> out(X.size(), State(X[0].size(), 0.0));
    for (std::size_t i = 1; i + 1 < X.size(); ++i) {
        for (std::size_t j = 0; j < X[i].size(); ++j) {
            out[i][j] = (X[i + 1][j] - X[i - 1][j]) / (2.0 * dt);
        }
    }
    out.front() = out[1];
    out.back() = out[X.size() - 2];
    return out;
}

}  // namespace chaos
