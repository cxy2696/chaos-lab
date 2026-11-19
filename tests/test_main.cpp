// SPDX-License-Identifier: MIT
#include <algorithm>
#include <chrono>
#include <cstring>
#include <iostream>
#include <string>

#include "test_framework.hpp"

int main(int argc, char** argv) {
    const std::string filter = argc > 1 ? argv[1] : "";
    int passed = 0;
    int failed = 0;
    const auto start = std::chrono::steady_clock::now();
    for (const auto& test : chaos::test::registry()) {
        if (!filter.empty() && test.name.find(filter) == std::string::npos) continue;
        try {
            test.fn();
            std::cout << "[ PASS ] " << test.name << "\n";
            ++passed;
        } catch (const std::exception& e) {
            std::cout << "[ FAIL ] " << test.name << "\n         " << e.what() << "\n";
            ++failed;
        } catch (...) {
            std::cout << "[ FAIL ] " << test.name << "\n         unknown exception\n";
            ++failed;
        }
    }
    const double seconds = std::chrono::duration<double>(
                               std::chrono::steady_clock::now() - start)
                               .count();
    std::cout << "\n" << passed << " passed, " << failed << " failed, " << seconds << " s\n";
    return failed == 0 ? 0 : 1;
}
