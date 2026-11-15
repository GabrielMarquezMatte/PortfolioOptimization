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
    src/capi/error.cpp)
set_target_properties(portopt_core PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(portopt_core PUBLIC src include)
target_compile_options(portopt_core PRIVATE ${PO_WARNING_FLAGS})

add_library(portopt SHARED
    src/capi/version_api.cpp)
target_include_directories(portopt PUBLIC include)
target_compile_definitions(portopt PRIVATE PO_BUILDING PO_VERSION_STRING="${PROJECT_VERSION}")
target_compile_options(portopt PRIVATE ${PO_WARNING_FLAGS})
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

    # The C entry points are compiled into the test binary next to portopt_core so that po_last_error() and
    # po::capi::guard() share one thread-local error; a separately loaded portopt would hold its own copy.
    add_executable(portopt_tests
        src/capi/version_api.cpp
        tests/unit/capi/abi_test.cpp)
    target_include_directories(portopt_tests PRIVATE tests/unit)
    target_compile_definitions(portopt_tests PRIVATE PO_BUILDING PO_VERSION_STRING="${PROJECT_VERSION}")
    target_link_libraries(portopt_tests PRIVATE portopt_core GTest::gtest_main)
    target_compile_options(portopt_tests PRIVATE ${PO_WARNING_FLAGS})
    include(GoogleTest)
    gtest_discover_tests(portopt_tests)

    find_package(Python3 REQUIRED COMPONENTS Interpreter)
    add_test(NAME style
        COMMAND Python3::Interpreter ${CMAKE_SOURCE_DIR}/tests/style/check_style.py ${CMAKE_SOURCE_DIR})
endif()
