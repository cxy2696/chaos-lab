// Minimal test harness: no external dependencies.
// SPDX-License-Identifier: MIT
#pragma once

#include <cmath>
#include <exception>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace chaos {
namespace test {

struct TestCase {
    std::string name;
    void (*fn)();
};

inline std::vector<TestCase>& registry() {
    static std::vector<TestCase> r;
    return r;
}

struct Registrar {
    Registrar(const char* name, void (*fn)()) { registry().push_back({name, fn}); }
};

struct Failure : std::runtime_error {
    explicit Failure(const std::string& what) : std::runtime_error(what) {}
};

inline void check(bool ok, const char* expr, const char* file, int line) {
    if (!ok) {
        std::ostringstream oss;
        oss << file << ":" << line << ": check failed: " << expr;
        throw Failure(oss.str());
    }
}

inline void check_near(double a, double b, double tol, const char* expr, const char* file,
                       int line) {
    if (!(std::fabs(a - b) <= tol)) {
        std::ostringstream oss;
        oss << file << ":" << line << ": " << expr << " expected " << b << " +/- " << tol
            << ", got " << a;
        throw Failure(oss.str());
    }
}

inline void check_in(double a, double lo, double hi, const char* expr, const char* file,
                     int line) {
    if (!(a >= lo && a <= hi)) {
        std::ostringstream oss;
        oss << file << ":" << line << ": " << expr << " expected in [" << lo << ", " << hi
            << "], got " << a;
        throw Failure(oss.str());
    }
}

}  // namespace test
}  // namespace chaos

#define TEST_CASE(name)                                                     \
    static void name();                                                     \
    static ::chaos::test::Registrar test_registrar_##name(#name, name);     \
    static void name()

#define CHECK(expr) ::chaos::test::check((expr), #expr, __FILE__, __LINE__)
#define CHECK_NEAR(a, b, tol) \
    ::chaos::test::check_near((a), (b), (tol), #a " ~= " #b, __FILE__, __LINE__)
#define CHECK_IN(a, lo, hi) \
    ::chaos::test::check_in((a), (lo), (hi), #a, __FILE__, __LINE__)
