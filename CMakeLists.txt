# ============================================================================
# Calculator CMake Build Configuration
# ============================================================================

cmake_minimum_required(VERSION 3.16)
project(Calculator VERSION 1.0.0 LANGUAGES CXX)

# C++ standard
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_BUILD_TYPE Release)

# Compiler flags
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-Wall -Wextra -pedantic)

    # Debug build flags
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
        add_compile_options(-g -O0)
    elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
        add_compile_options(-O2)
    endif()

    # Architecture-specific optimizations for high-precision math
    include(CheckCXXCompilerFlag)

    # AVX2 support for SIMD optimizations
    check_cxx_compiler_flag("-mavx2" COMPILER_SUPPORTS_AVX2)
    if(COMPILER_SUPPORTS_AVX2)
        add_compile_options(-mavx2)
    endif()

    # SSE4.2 support
    check_cxx_compiler_flag("-msse4.2" COMPILER_SUPPORTS_SSE42)
    if(COMPILER_SUPPORTS_SSE42)
        add_compile_options(-msse4.2)
    endif()

    # Native architecture optimization (optional, can be disabled)
    option(USE_NATIVE_ARCH "Use -march=native for maximum performance" OFF)
    if(USE_NATIVE_ARCH)
        add_compile_options(-march=native)
    endif()
endif()

# Source directories
set(SRC_DIR ${CMAKE_SOURCE_DIR}/src)

# Collect all subdirectories under src/ for includes
file(GLOB SRC_SUBDIRS LIST_DIRECTORIES true ${SRC_DIR}/*)
foreach(subdir ${SRC_SUBDIRS})
    if(IS_DIRECTORY ${subdir})
        list(APPEND INCLUDE_DIRS ${subdir})
        # Also add subdirectories (like math/helpers)
        file(GLOB SUBSUBDIRS LIST_DIRECTORIES true ${subdir}/*)
        foreach(subsubdir ${SUBSUBDIRS})
            if(IS_DIRECTORY ${subsubdir})
                list(APPEND INCLUDE_DIRS ${subsubdir})
            endif()
        endforeach()
    endif()
endforeach()
list(APPEND INCLUDE_DIRS ${SRC_DIR} ${CMAKE_SOURCE_DIR}/test)
list(REMOVE_DUPLICATES INCLUDE_DIRS)

# Collect all source files using wildcards
file(GLOB_RECURSE COMMON_SOURCES ${SRC_DIR}/*.cpp)

# Main executable (exclude main.cpp from common sources)
list(FILTER COMMON_SOURCES EXCLUDE REGEX ".*/app/main\\.cpp$")

add_executable(calculator
    ${SRC_DIR}/app/main.cpp
    ${COMMON_SOURCES}
)

target_include_directories(calculator PRIVATE ${INCLUDE_DIRS})

# Test executable (optional)
option(BUILD_TESTS "Build test executable" ON)

if(BUILD_TESTS)
    # Collect test source files
    file(GLOB TEST_SUITE_SOURCES ${CMAKE_SOURCE_DIR}/test/suites/*.cpp)

    add_executable(calculator_tests
        ${CMAKE_SOURCE_DIR}/test/main.cpp
        ${TEST_SUITE_SOURCES}
        ${COMMON_SOURCES}
    )

    target_include_directories(calculator_tests PRIVATE ${INCLUDE_DIRS})

    # Enable testing
    enable_testing()
    add_test(NAME CalculatorTests COMMAND calculator_tests)
endif()

# Install targets (optional)
option(INSTALL_CALCULATOR "Install calculator binary" OFF)

if(INSTALL_CALCULATOR)
    install(TARGETS calculator DESTINATION bin)
endif()

# Print build info
message(STATUS "")
message(STATUS "Calculator Build Configuration:")
message(STATUS "  Version:      ${PROJECT_VERSION}")
message(STATUS "  C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "  Build Type:   ${CMAKE_BUILD_TYPE}")
message(STATUS "  Compiler:     ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "  Build Tests:  ${BUILD_TESTS}")
message(STATUS "  Sources:      ${CMAKE_CURRENT_SOURCE_DIR}/src/**/*.cpp")
message(STATUS "")