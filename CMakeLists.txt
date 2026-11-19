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
    src/core/portable_math.cpp
    src/problems/constraints.cpp
    src/problems/homogeneous.cpp
    src/problems/max_return.cpp
    src/problems/parametric_path.cpp
    src/problems/mean_variance.cpp
    src/problems/tracking_error.cpp
    src/riskparity/augmented_lagrangian.cpp
    src/riskparity/kkt.cpp
    src/riskparity/newton_polish.cpp
    src/riskparity/objective.cpp
    src/riskparity/polyhedral_descent.cpp
    src/riskparity/risk_parity.cpp
    src/riskparity/spinu.cpp
    src/solve/requirements.cpp
    src/solve/solve.cpp
    src/metrics/portfolio.cpp)
set_target_properties(portopt_core PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(portopt_core PUBLIC src include)
target_compile_options(portopt_core PRIVATE ${PO_WARNING_FLAGS})
# The compensated products need every a * b + c rounded as written. GCC and Clang contract it into an fma by
# default when the target has one (-march=native), even under -std=c++23; MSVC's default /fp:precise does not.
if(NOT MSVC)
    target_compile_options(portopt_core PRIVATE -ffp-contract=off)
endif()

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
# The entry points do no arithmetic of their own, but they include core headers: inline core code (and, under LTO,
# any core function) compiled here is contracted into fma unless this unit says otherwise.
if(NOT MSVC)
    target_compile_options(portopt_capi PRIVATE -ffp-contract=off)
endif()
target_link_libraries(portopt_capi PUBLIC portopt_core)

add_library(portopt SHARED $<TARGET_OBJECTS:portopt_capi>)
target_include_directories(portopt PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(portopt PRIVATE portopt_core)
if(MINGW)
    # No --output-def: ld 2.47 fails at random with it and writes an empty file. The export list of the release .lib
    # comes from the ABI contract (scripts/abi_contract.py def), checked against the DLL's own export table.
    target_link_options(portopt PRIVATE
        -static-libgcc -static-libstdc++
        -Wl,-Bstatic,--whole-archive -lwinpthread -Wl,--no-whole-archive,-Bdynamic)
elseif(NOT WIN32 AND NOT APPLE)
    # libstdc++.so.6 and libgcc_s.so.1 are not on every distribution (musl, minimal images).
    # --exclude-libs keeps the symbols of the static libstdc++ out of the export table: the library exports the C
    # interface (po_*) and nothing else, so it cannot interpose on the libstdc++ of the program that loads it, and the
    # ABI dump lists only the interface.
    target_link_options(portopt PRIVATE -static-libstdc++ -static-libgcc -Wl,--exclude-libs,ALL)
elseif(APPLE)
    # libc++ marks the typeinfo of its templates visible whatever -fvisibility says, and they end up in the export
    # table; only the C interface is exported.
    target_link_options(portopt PRIVATE "LINKER:-exported_symbol,_po_*")
endif()

# The header-only C++ wrapper (cpp/include/portopt/portopt.hpp) over the shared library.
add_library(portopt_cpp INTERFACE)
target_include_directories(portopt_cpp INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/cpp/include>)
target_compile_features(portopt_cpp INTERFACE cxx_std_23)
target_link_libraries(portopt_cpp INTERFACE portopt)
set_target_properties(portopt_cpp PROPERTIES EXPORT_NAME cpp)
add_library(PortOpt::portopt ALIAS portopt)
add_library(PortOpt::cpp ALIAS portopt_cpp)

option(PO_INSTALL "Generate the install rules and the find_package(PortOpt) package" ${PO_IS_TOP_LEVEL})
if(PO_INSTALL)
    include(GNUInstallDirs)
    include(CMakePackageConfigHelpers)
    set(PO_CONFIG_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/PortOpt")
    install(TARGETS portopt portopt_cpp EXPORT PortOptTargets
        RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
        LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
        ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
        INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")
    install(FILES include/portopt/portopt.h cpp/include/portopt/portopt.hpp
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/portopt")
    install(EXPORT PortOptTargets NAMESPACE PortOpt:: DESTINATION "${PO_CONFIG_DIR}")
    configure_package_config_file(cmake/PortOptConfig.cmake.in "${CMAKE_CURRENT_BINARY_DIR}/PortOptConfig.cmake"
        INSTALL_DESTINATION "${PO_CONFIG_DIR}")
    # Until 1.0 the ABI may break between minor versions (PO_ABI_VERSION 0): only the same x.y is compatible.
    write_basic_package_version_file("${CMAKE_CURRENT_BINARY_DIR}/PortOptConfigVersion.cmake"
        VERSION "${PROJECT_VERSION}" COMPATIBILITY SameMinorVersion)
    install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}/PortOptConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/PortOptConfigVersion.cmake"
        DESTINATION "${PO_CONFIG_DIR}")
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
        tests/unit/qp/active_set_test.cpp
        tests/unit/qp/goldfarb_idnani_test.cpp
        tests/unit/core/data_test.cpp
        tests/unit/core/portable_math_test.cpp
        tests/unit/fuzz/leverage_fuzz_test.cpp
        tests/unit/fuzz/max_return_fuzz_test.cpp
        tests/unit/fuzz/ratio_fuzz_test.cpp
        tests/unit/fuzz/risk_parity/risk_parity_descent_fuzz_test.cpp
        tests/unit/fuzz/risk_parity/risk_parity_fuzz_test.cpp
        tests/unit/fuzz/risk_parity/risk_parity_multiplier_measure_test.cpp
        tests/unit/fuzz/risk_parity/risk_parity_reference_test.cpp
        tests/unit/fuzz/risk_parity/risk_parity_solve_sweep_test.cpp
        tests/unit/problems/leverage_test.cpp
        tests/unit/problems/max_return_test.cpp
        tests/unit/problems/max_sharpe_test.cpp
        tests/unit/problems/min_variance_test.cpp
        tests/unit/problems/parametric_path_test.cpp
        tests/unit/problems/utility_test.cpp
        tests/unit/problems/tracking_error_test.cpp
        tests/unit/riskparity/augmented_lagrangian_test.cpp
        tests/unit/riskparity/newton_polish_test.cpp
        tests/unit/riskparity/objective_test.cpp
        tests/unit/riskparity/polyhedral_descent_test.cpp
        tests/unit/riskparity/risk_parity_test.cpp
        tests/unit/riskparity/spinu_test.cpp
        tests/unit/solve/requirements_test.cpp
        tests/unit/solve/solve_test.cpp
        tests/unit/metrics/portfolio_test.cpp)
    target_include_directories(portopt_tests PRIVATE tests/unit)
    target_compile_definitions(portopt_tests PRIVATE PO_BUILDING)
    target_link_libraries(portopt_tests PRIVATE portopt_capi GTest::gtest_main)
    target_compile_options(portopt_tests PRIVATE ${PO_WARNING_FLAGS})
    # The test oracle (compensated_dot_of in tests/unit/fuzz/risk_parity) needs the same rounding as the core.
    if(NOT MSVC)
        target_compile_options(portopt_tests PRIVATE -ffp-contract=off)
    endif()
    if(MSVC)
        # The measurement tests read PORTOPT_* variables with getenv and write files with fopen.
        target_compile_definitions(portopt_tests PRIVATE _CRT_SECURE_NO_WARNINGS)
    endif()
    if(MINGW)
        # An older libstdc++-6.dll earlier on PATH crashes std::format in the test executable.
        target_link_options(portopt_tests PRIVATE -static-libgcc -static-libstdc++)
    endif()
    include(GoogleTest)
    gtest_discover_tests(portopt_tests)

    find_package(Python3 REQUIRED COMPONENTS Interpreter)
    add_test(NAME style
        COMMAND Python3::Interpreter ${CMAKE_SOURCE_DIR}/tests/style/check_style.py ${CMAKE_SOURCE_DIR})
    # include/portopt/portopt.abi is the contract of portopt.h; a change to the header is a change to that file.
    add_test(NAME abi_contract
        COMMAND Python3::Interpreter ${CMAKE_SOURCE_DIR}/scripts/abi_contract.py check)
    add_test(NAME abi_contract_unit
        COMMAND Python3::Interpreter -m unittest discover -s ${CMAKE_SOURCE_DIR}/tests/abi)

    add_subdirectory(cpp/tests)
endif()

if(PO_BUILD_BENCHMARKS)
    set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
    set(BENCHMARK_ENABLE_GTEST_TESTS OFF CACHE BOOL "" FORCE)
    set(BENCHMARK_ENABLE_INSTALL OFF CACHE BOOL "" FORCE)
    set(BENCHMARK_ENABLE_WERROR OFF CACHE BOOL "" FORCE)
    FetchContent_Declare(benchmark
        URL https://github.com/google/benchmark/archive/refs/tags/v1.9.1.tar.gz
        URL_HASH SHA256=32131c08ee31eeff2c8968d7e874f3cb648034377dfc32a4c377fa8796d84981
        SYSTEM)
    FetchContent_MakeAvailable(benchmark)

    # Both executables replace the global operator new to count allocations (alloc_counter.cpp); nothing else does.
    # They call the core directly (and the C entry points, objects of portopt_capi, for BM_CApi and `--via capi`)
    # and compile its headers, so they need the same rounding as the core.
    add_executable(portopt_benchmarks
        tests/benchmarks/alloc_counter.cpp
        tests/benchmarks/bench_main.cpp
        tests/benchmarks/bm_kernel.cpp
        tests/benchmarks/capi/bm_capi.cpp
        tests/benchmarks/bm_solve.cpp)
    add_executable(portopt_profile
        tests/benchmarks/alloc_counter.cpp
        tests/benchmarks/profile.cpp)
    target_link_libraries(portopt_benchmarks PRIVATE portopt_capi benchmark::benchmark)
    target_link_libraries(portopt_profile PRIVATE portopt_capi)
    foreach(target portopt_benchmarks portopt_profile)
        target_include_directories(${target} PRIVATE tests/unit tests/benchmarks)
        target_compile_definitions(${target} PRIVATE PO_BUILDING)
        target_compile_options(${target} PRIVATE ${PO_WARNING_FLAGS})
        if(NOT MSVC)
            target_compile_options(${target} PRIVATE -ffp-contract=off)
        endif()
        if(MINGW)
            target_link_options(${target} PRIVATE -static-libgcc -static-libstdc++)
        endif()
    endforeach()

    enable_testing()
    find_package(Python3 REQUIRED COMPONENTS Interpreter)
    add_test(NAME benchmark_allocation_selftest COMMAND portopt_benchmarks --selftest)
    add_test(NAME benchmark_quick_selftest
        COMMAND Python3::Interpreter ${CMAKE_SOURCE_DIR}/tests/benchmarks/quick.py selftest)
endif()
