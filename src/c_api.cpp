// SPDX-License-Identifier: MIT
#include "chaos/chaos_c.h"

#include <cmath>
#include <exception>
#include <string>
#include <vector>

#include "chaos/embedding.hpp"
#include "chaos/integrators.hpp"
#include "chaos/lyapunov.hpp"
#include "chaos/systems.hpp"
#include "chaos/timeseries.hpp"
#include "chaos/version.hpp"

namespace {
thread_local std::string g_error;
void set_error(const char* msg) { g_error = msg ? msg : "unknown error"; }
}  // namespace

extern "C" {

const char* chaos_version(void) { return CHAOS_VERSION_STRING; }

const char* chaos_last_error(void) { return g_error.c_str(); }

int chaos_lorenz_trajectory(double sigma, double rho, double beta, const double x0[3],
                            double h, size_t n, double* out_t, double* out_xyz) {
    try {
        if (!x0 || !out_xyz || !(h > 0.0)) {
            set_error("invalid argument");
            return -1;
        }
        chaos::Lorenz model{sigma, rho, beta};
        chaos::OdeSystem sys = model.ode();
        chaos::State x{x0[0], x0[1], x0[2]};
        double t = 0.0;
        for (size_t i = 0; i < n; ++i) {
            out_xyz[i * 3 + 0] = x[0];
            out_xyz[i * 3 + 1] = x[1];
            out_xyz[i * 3 + 2] = x[2];
            if (out_t) out_t[i] = t;
            x = chaos::rk4_step(sys, x, t, h);
            t += h;
        }
        return 0;
    } catch (const std::exception& e) {
        set_error(e.what());
        return -1;
    }
}

int chaos_lorenz_lyapunov(double sigma, double rho, double beta, const double x0[3],
                          double h, double transient, double duration,
                          double out_exponents[3]) {
    try {
        if (!x0 || !out_exponents) {
            set_error("invalid argument");
            return -1;
        }
        chaos::Lorenz model{sigma, rho, beta};
        chaos::LyapunovOptions opt;
        opt.h = h;
        opt.transient = transient;
        opt.duration = duration;
        const chaos::LyapunovSpectrum sp =
            chaos::lyapunov_spectrum(model.ode(), chaos::State{x0[0], x0[1], x0[2]}, opt);
        if (sp.exponents.size() != 3) {
            set_error("unexpected spectrum size");
            return -1;
        }
        for (int i = 0; i < 3; ++i) out_exponents[i] = sp.exponents[static_cast<size_t>(i)];
        return 0;
    } catch (const std::exception& e) {
        set_error(e.what());
        return -1;
    }
}

double chaos_logistic_lyapunov(double r, size_t transient, size_t n) {
    try {
        return chaos::Logistic{r}.lyapunov(0.123456789, transient, n);
    } catch (const std::exception& e) {
        set_error(e.what());
        return std::nan("");
    }
}

int chaos_correlation_dimension(const double* series, size_t n, int m, int tau,
                                int theiler, double* out_d2) {
    try {
        if (!series || !out_d2 || n < 10) {
            set_error("invalid argument");
            return -1;
        }
        const std::vector<double> x(series, series + n);
        const chaos::CorrelationResult r = chaos::correlation_dimension(
            chaos::delay_embed(x, m, tau), 0.0, 0.0, 40, theiler);
        *out_d2 = r.d2;
        return 0;
    } catch (const std::exception& e) {
        set_error(e.what());
        return -1;
    }
}

int chaos_zero_one_test(const double* series, size_t n, int n_c, double* out_k) {
    try {
        if (!series || !out_k || n < 100) {
            set_error("invalid argument");
            return -1;
        }
        const std::vector<double> x(series, series + n);
        const chaos::ZeroOneResult r = chaos::zero_one_test(x, n_c, 12345);
        *out_k = r.K;
        return 0;
    } catch (const std::exception& e) {
        set_error(e.what());
        return -1;
    }
}

int chaos_rosenstein_lambda(const double* series, size_t n, int m, int tau, int theiler,
                            int max_k, double dt, double* out_lambda) {
    try {
        if (!series || !out_lambda || n < 20) {
            set_error("invalid argument");
            return -1;
        }
        const std::vector<double> x(series, series + n);
        const chaos::RosensteinResult r =
            chaos::rosenstein_lyapunov(x, m, tau, theiler, max_k, dt);
        *out_lambda = r.lambda;
        return 0;
    } catch (const std::exception& e) {
        set_error(e.what());
        return -1;
    }
}

}  // extern "C"
