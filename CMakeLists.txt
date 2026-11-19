cmake_minimum_required(VERSION 3.20)
project(chaos_lab VERSION 0.1.0 LANGUAGES C CXX)

include(GNUInstallDirs)

option(CHAOS_BUILD_TESTS "Build the test suite" ON)
option(CHAOS_BUILD_CLI "Build the chaos_cli command-line tool" ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CHAOS_SOURCES
    src/core.cpp
    src/integrators.cpp
    src/systems.cpp
    src/lyapunov.cpp
    src/embedding.cpp
    src/timeseries.cpp
    src/tools.cpp
    src/control.cpp
    src/sindy.cpp
    src/c_api.cpp
)

add_library(chaos ${CHAOS_SOURCES})
add_library(chaos::chaos ALIAS chaos)
target_compile_features(chaos PUBLIC cxx_std_17)
target_include_directories(chaos
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
set_target_properties(chaos PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
    POSITION_INDEPENDENT_CODE ON
)
if(MSVC)
    target_compile_options(chaos PRIVATE /W4 /permissive-)
else()
    target_compile_options(chaos PRIVATE -Wall -Wextra -Wpedantic)
endif()

if(CHAOS_BUILD_CLI)
    add_executable(chaos_cli apps/chaos_cli.cpp)
    target_link_libraries(chaos_cli PRIVATE chaos::chaos)
    install(TARGETS chaos_cli RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
    target_compile_features(chaos_cli PRIVATE cxx_std_17)
    if(MSVC)
        target_compile_options(chaos_cli PRIVATE /W4 /permissive-)
    else()
        target_compile_options(chaos_cli PRIVATE -Wall -Wextra -Wpedantic)
    endif()
endif()

if(CHAOS_BUILD_TESTS)
    enable_testing()
    add_executable(chaos_tests
        tests/test_main.cpp
        tests/test_core.cpp
        tests/test_analysis.cpp
        tests/test_control.cpp
    )
    target_link_libraries(chaos_tests PRIVATE chaos::chaos)
    target_include_directories(chaos_tests PRIVATE tests)
    target_compile_features(chaos_tests PRIVATE cxx_std_17)
    if(NOT MSVC)
        target_compile_options(chaos_tests PRIVATE -Wall -Wextra -Wpedantic)
    endif()
    add_test(NAME chaos_tests COMMAND chaos_tests)
    set_tests_properties(chaos_tests PROPERTIES TIMEOUT 600)

    add_executable(chaos_c_smoke tests/test_c_api.c)
    target_link_libraries(chaos_c_smoke PRIVATE chaos::chaos)
    set_target_properties(chaos_c_smoke PROPERTIES LINKER_LANGUAGE CXX)
    add_test(NAME chaos_c_api COMMAND chaos_c_smoke)
endif()

install(TARGETS chaos
    EXPORT chaosTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(EXPORT chaosTargets
    FILE chaosTargets.cmake
    NAMESPACE chaos::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/chaos
)
include(CMakePackageConfigHelpers)
configure_package_config_file(
    cmake/chaosConfig.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/chaosConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/chaos
)
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/chaosConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)
install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/chaosConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/chaosConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/chaos
)
