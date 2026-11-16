cmake_minimum_required(VERSION 3.28)
set(PO_VERSION 0.1.0 CACHE STRING "PortOpt version (x.y.z)")
project(PortOpt VERSION ${PO_VERSION} LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
# No C++20 modules: without this, Ninja adds @modmap flags to compile_commands.json that clang-tidy cannot read.
set(CMAKE_CXX_SCAN_FOR_MODULES OFF)
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)
set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_C_VISIBILITY_PRESET hidden)
set(CMAKE_VISIBILITY_INLINES_HIDDEN ON)
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS OFF)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    set(PO_IS_TOP_LEVEL ON)
else()
    set(PO_IS_TOP_LEVEL OFF)
endif()
option(PO_BUILD_TESTS "Build the PortOpt tests" ${PO_IS_TOP_LEVEL})
option(PO_BUILD_BENCHMARKS "Build the PortOpt benchmarks" OFF)

if(MSVC)
    set(PO_WARNING_FLAGS /W4 /WX)
else()
    set(PO_WARNING_FLAGS -Wall -Wextra -Wpedantic -Werror)
endif()
# Every return of a named local must be elided (STYLEGUIDE "Return values"). Only GCC 14+ can report a missed
# elision; GCC and Clang already reject pessimizing and redundant moves through -Wall and -Wextra.
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 14)
    list(APPEND PO_WARNING_FLAGS $<$<COMPILE_LANGUAGE:CXX>:-Wnrvo>)
endif()

include(FetchContent)

add_library(portopt_core STATIC
    src/linalg/cholesky.cpp
    src/linalg/qr.cpp
    src/qp/goldfarb_idnani.cpp
    src/core/data.cpp
    src/problems/constraints.cpp
    src/problems/homogeneous.cpp
    src/problems/max_return.cpp
    src/problems/parametric_path.cpp
    src/problems/mean_variance.cpp
    src/problems/tracking_error.cpp
    src/riskparity/risk_parity.cpp
    src/riskparity/spinu.cpp
    src/solve/requirements.cpp
    src/solve/solve.cpp
    src/metrics/portfolio.cpp)
set_target_properties(portopt_core PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(portopt_core PUBLIC src include)
target_compile_options(portopt_core PRIVATE ${PO_WARNING_FLAGS})

# The C entry points are compiled once and shared by the shared library and the tests, so that po_last_error() and
# po::capi::guard() in the tests see the same thread-local error as the entry points they call.
add_library(portopt_capi OBJECT
    src/capi/error.cpp
    src/capi/version_api.cpp
    src/capi/spec_api.cpp
    src/capi/data_api.cpp
    src/capi/solve_api.cpp
    src/capi/result_api.cpp
    src/capi/metrics_api.cpp)
set_target_properties(portopt_capi PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_compile_definitions(portopt_capi PRIVATE PO_BUILDING PO_VERSION_STRING="${PROJECT_VERSION}")
target_compile_options(portopt_capi PRIVATE ${PO_WARNING_FLAGS})
target_link_libraries(portopt_capi PUBLIC portopt_core)

add_library(portopt SHARED $<TARGET_OBJECTS:portopt_capi>)
target_include_directories(portopt PUBLIC include)
target_link_libraries(portopt PRIVATE portopt_core)
if(MINGW)
    target_link_options(portopt PRIVATE
        -static-libgcc -static-libstdc++
        -Wl,-Bstatic,--whole-archive -lwinpthread -Wl,--no-whole-archive,-Bdynamic
        -Wl,--output-def,${CMAKE_CURRENT_BINARY_DIR}/portopt.def)
elseif(NOT WIN32 AND NOT APPLE)
    # libstdc++.so.6 and libgcc_s.so.1 are not on every distribution (musl, minimal images).
    target_link_options(portopt PRIVATE -static-libstdc++ -static-libgcc)
endif()

if(PO_BUILD_TESTS)
    enable_testing()
    set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_Declare(googletest
        URL https://github.com/google/googletest/archive/refs/tags/v1.17.0.tar.gz)
    FetchContent_MakeAvailable(googletest)

    add_library(portopt_c_header_check OBJECT tests/unit/capi/c_header_check.c)
    target_link_libraries(portopt_c_header_check PRIVATE portopt)
    target_compile_options(portopt_c_header_check PRIVATE ${PO_WARNING_FLAGS})

    add_executable(portopt_tests
        tests/unit/capi/abi_test.cpp
        tests/unit/capi/spec_test.cpp
        tests/unit/capi/data_test.cpp
        tests/unit/capi/solve_test.cpp
        tests/unit/capi/metrics_test.cpp
        tests/unit/linalg/cholesky_test.cpp
        tests/unit/linalg/qr_test.cpp
        tests/unit/qp/goldfarb_idnani_test.cpp
        tests/unit/core/data_test.cpp
        tests/unit/fuzz/leverage_fuzz_test.cpp
        tests/unit/fuzz/max_return_fuzz_test.cpp
        tests/unit/fuzz/ratio_fuzz_test.cpp
        tests/unit/fuzz/risk_parity/risk_parity_reference_test.cpp
        tests/unit/problems/leverage_test.cpp
        tests/unit/problems/max_return_test.cpp
        tests/unit/problems/max_sharpe_test.cpp
        tests/unit/problems/min_variance_test.cpp
        tests/unit/problems/utility_test.cpp
        tests/unit/problems/tracking_error_test.cpp
        tests/unit/riskparity/risk_parity_test.cpp
        tests/unit/riskparity/spinu_test.cpp
        tests/unit/solve/requirements_test.cpp
        tests/unit/solve/solve_test.cpp
        tests/unit/metrics/portfolio_test.cpp)
    target_include_directories(portopt_tests PRIVATE tests/unit)
    target_compile_definitions(portopt_tests PRIVATE PO_BUILDING)
    target_link_libraries(portopt_tests PRIVATE portopt_capi GTest::gtest_main)
    target_compile_options(portopt_tests PRIVATE ${PO_WARNING_FLAGS})
    if(MINGW)
        # An older libstdc++-6.dll earlier on PATH crashes std::format in the test executable.
        target_link_options(portopt_tests PRIVATE -static-libgcc -static-libstdc++)
    endif()
    include(GoogleTest)
    gtest_discover_tests(portopt_tests)

    find_package(Python3 REQUIRED COMPONENTS Interpreter)
    add_test(NAME style
        COMMAND Python3::Interpreter ${CMAKE_SOURCE_DIR}/tests/style/check_style.py ${CMAKE_SOURCE_DIR})
endif()
