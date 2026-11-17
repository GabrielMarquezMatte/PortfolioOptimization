cmake_minimum_required(VERSION 3.28)

# Consumes the C++ wrapper the way a downstream project does, in one of two modes:
#   add_subdirectory  this repository as a subdirectory of the consumer
#   find_package      an install prefix made by `cmake --install` (pass -DCMAKE_PREFIX_PATH=<prefix>)
# It is configured as its own project, so a package that only works inside the build tree of PortOpt fails here.
project(portopt_package_consumer LANGUAGES CXX)

set(PORTOPT_CONSUMER_MODE "find_package" CACHE STRING "add_subdirectory or find_package")
set(PORTOPT_CONSUMER_VERSION "" CACHE STRING "Version find_package asks for (empty: any)")
set(PORTOPT_REPOSITORY "${CMAKE_CURRENT_SOURCE_DIR}/../../.." CACHE PATH "Root of the PortOpt repository")

if(PORTOPT_CONSUMER_MODE STREQUAL "add_subdirectory")
    add_subdirectory("${PORTOPT_REPOSITORY}" portopt)
elseif(PORTOPT_CONSUMER_MODE STREQUAL "find_package")
    find_package(PortOpt ${PORTOPT_CONSUMER_VERSION} REQUIRED)
else()
    message(FATAL_ERROR "PORTOPT_CONSUMER_MODE must be add_subdirectory or find_package")
endif()

add_executable(package_consumer main.cpp)
target_link_libraries(package_consumer PRIVATE PortOpt::cpp)
if(WIN32)
    # A DLL that is neither beside the executable nor on PATH is not found.
    add_custom_command(TARGET package_consumer POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_RUNTIME_DLLS:package_consumer> $<TARGET_FILE_DIR:package_consumer>
        COMMAND_EXPAND_LISTS)
endif()

enable_testing()
add_test(NAME package_consumer COMMAND package_consumer)
