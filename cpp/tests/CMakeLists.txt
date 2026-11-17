# The wrapper tests link the shared library, as a user of the wrapper does (the core tests link the C entry points as
# objects instead).
add_executable(portopt_cpp_tests wrapper_test.cpp)
target_link_libraries(portopt_cpp_tests PRIVATE PortOpt::cpp GTest::gtest_main)
target_compile_options(portopt_cpp_tests PRIVATE ${PO_WARNING_FLAGS})
if(MINGW)
    # An older libstdc++-6.dll earlier on PATH crashes the test executable.
    target_link_options(portopt_cpp_tests PRIVATE -static-libgcc -static-libstdc++)
endif()
if(WIN32)
    # Before test discovery, which runs the executable: the DLL has to be beside it.
    add_custom_command(TARGET portopt_cpp_tests POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_RUNTIME_DLLS:portopt_cpp_tests> $<TARGET_FILE_DIR:portopt_cpp_tests>
        COMMAND_EXPAND_LISTS)
endif()
gtest_discover_tests(portopt_cpp_tests)

# Every PO_API function of the C header has to be called by the wrapper.
add_test(NAME abi_coverage
    COMMAND Python3::Interpreter ${CMAKE_CURRENT_SOURCE_DIR}/../scripts/check_abi_coverage.py)
