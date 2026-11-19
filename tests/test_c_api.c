/* C API smoke test: verifies that the public C header works from a C
 * translation unit and that the numerical entry points agree with theory. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "chaos/chaos_c.h"

static int failures = 0;
#define CHECK(cond)                                                        \
    do {                                                                   \
        if (!(cond)) {                                                     \
            fprintf(stderr, "FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond); \
            failures++;                                                    \
        }                                                                  \
    } while (0)

int main(void) {
    CHECK(strcmp(chaos_version(), "0.1.0") == 0);

    double xyz[30 * 3];
    double t[30];
    const double x0[3] = {1.0, 1.0, 1.0};
    CHECK(chaos_lorenz_trajectory(10.0, 28.0, 8.0 / 3.0, x0, 0.001, 30, t, xyz) == 0);
    for (int i = 0; i < 30; ++i) {
        CHECK(isfinite(xyz[i * 3]));
        CHECK(isfinite(t[i]));
    }

    double spectrum[3];
    CHECK(chaos_lorenz_lyapunov(10.0, 28.0, 8.0 / 3.0, x0, 0.01, 30.0, 120.0,
                                spectrum) == 0);
    CHECK(spectrum[0] > 0.8);
    CHECK(spectrum[2] < -14.0);
    CHECK(spectrum[0] + spectrum[1] + spectrum[2] < -13.0);

    const double lambda = chaos_logistic_lyapunov(4.0, 1000, 200000);
    CHECK(fabs(lambda - log(2.0)) < 1e-2);
    CHECK(fabs(chaos_logistic_lyapunov(3.5, 1000, 20000)) > 0.1);

    const int n = 4000;
    double* logistic = (double*)malloc(sizeof(double) * n);
    CHECK(logistic != NULL);
    double v = 0.123456789;
    for (int i = 0; i < n; ++i) {
        v = 4.0 * v * (1.0 - v);
        logistic[i] = v;
    }
    double d2 = 0.0;
    CHECK(chaos_correlation_dimension(logistic, n, 2, 1, 10, &d2) == 0);
    CHECK(d2 > 0.85 && d2 < 1.25);

    double k = 0.0;
    CHECK(chaos_zero_one_test(logistic, n, 32, &k) == 0);
    CHECK(k > 0.85);

    double rosen = 0.0;
    CHECK(chaos_rosenstein_lambda(logistic, n, 2, 1, 10, 30, 1.0, &rosen) == 0);
    CHECK(rosen > 0.5 && rosen < 0.9);

    free(logistic);

    double dummy = 0.0;
    CHECK(chaos_zero_one_test(NULL, 10, 32, &dummy) != 0);
    CHECK(chaos_last_error()[0] != '\0');

    if (failures == 0) {
        printf("C API smoke test passed\n");
        return 0;
    }
    fprintf(stderr, "%d C API checks failed\n", failures);
    return 1;
}
