/* C interface for chaos-lab.  Link against the C++ library (libchaos).
 * SPDX-License-Identifier: MIT */
#ifndef CHAOS_C_H
#define CHAOS_C_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

const char* chaos_version(void);
const char* chaos_last_error(void);

/* Integrate Lorenz with RK4.  out must hold n*3 doubles. */
int chaos_lorenz_trajectory(double sigma, double rho, double beta, const double x0[3],
                            double h, size_t n, double* out_t, double* out_xyz);

/* Full Lorenz Lyapunov spectrum via Benettin/Gram-Schmidt. out_exponents[3]. */
int chaos_lorenz_lyapunov(double sigma, double rho, double beta, const double x0[3],
                          double h, double transient, double duration,
                          double out_exponents[3]);

/* Logistic map Lyapunov exponent (natural log per iteration). */
double chaos_logistic_lyapunov(double r, size_t transient, size_t n);

/* Grassberger-Procaccia correlation dimension of a scalar series. */
int chaos_correlation_dimension(const double* series, size_t n, int m, int tau,
                                int theiler, double* out_d2);

/* Gottwald-Melbourne 0-1 test. */
int chaos_zero_one_test(const double* series, size_t n, int n_c, double* out_k);

/* Rosenstein maximum Lyapunov exponent from a scalar series. */
int chaos_rosenstein_lambda(const double* series, size_t n, int m, int tau, int theiler,
                            int max_k, double dt, double* out_lambda);

#ifdef __cplusplus
}
#endif

#endif /* CHAOS_C_H */
