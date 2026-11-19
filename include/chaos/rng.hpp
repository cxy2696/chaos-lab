// Small deterministic RNG (xoshiro256**), so that stochastic diagnostics are
// reproducible without pulling in a dependency.
// SPDX-License-Identifier: MIT
#pragma once

#include <cmath>
#include <cstdint>

namespace chaos {

class Rng {
public:
    explicit Rng(std::uint64_t s = 0x9e3779b97f4a7c15ULL) { seed(s); }

    void seed(std::uint64_t s) {
        for (auto& v : s_) v = splitmix64(s);
        have_normal_ = false;
    }

    std::uint64_t next_u64() {
        const std::uint64_t result = rotl(s_[1] * 5ULL, 7) * 9ULL;
        const std::uint64_t t = s_[1] << 17;
        s_[2] ^= s_[0];
        s_[3] ^= s_[1];
        s_[1] ^= s_[2];
        s_[0] ^= s_[3];
        s_[2] ^= t;
        s_[3] = rotl(s_[3], 45);
        return result;
    }

    // Uniform in [0, 1).
    double uniform() {
        return static_cast<double>(next_u64() >> 11) * 0x1.0p-53;
    }

    double uniform(double a, double b) { return a + (b - a) * uniform(); }

    int uniform_int(int n) {
        if (n <= 0) return 0;
        return static_cast<int>(uniform() * static_cast<double>(n));
    }

    // Standard normal via Box-Muller (cached second value).
    double normal() {
        if (have_normal_) {
            have_normal_ = false;
            return normal_cache_;
        }
        double u1 = 0.0;
        do {
            u1 = uniform();
        } while (u1 <= 1e-300);
        const double u2 = uniform();
        const double r = std::sqrt(-2.0 * std::log(u1));
        constexpr double kTwoPi = 6.283185307179586476925286766559;
        normal_cache_ = r * std::sin(kTwoPi * u2);
        have_normal_ = true;
        return r * std::cos(kTwoPi * u2);
    }

private:
    static std::uint64_t rotl(std::uint64_t x, int k) {
        return (x << k) | (x >> (64 - k));
    }

    static std::uint64_t splitmix64(std::uint64_t& x) {
        x += 0x9e3779b97f4a7c15ULL;
        std::uint64_t z = x;
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        return z ^ (z >> 31);
    }

    std::uint64_t s_[4]{};
    bool have_normal_ = false;
    double normal_cache_ = 0.0;
};

}  // namespace chaos
