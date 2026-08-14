# ============================================================================
# Calculator CMake Build Configuration (Modular Architecture)
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
        file(GLOB SUBSUBDIRS LIST_DIRECTORIES true ${subdir}/*)
        foreach(subsubdir ${SUBSUBDIRS})
            if(IS_DIRECTORY ${subsubdir})
                list(APPEND INCLUDE_DIRS ${subsubdir})
            endif()
        endforeach()
    endif()
endforeach()
list(APPEND INCLUDE_DIRS ${SRC_DIR})
list(REMOVE_DUPLICATES INCLUDE_DIRS)

# ============================================================================
# Modular Component Sources
# ============================================================================

file(GLOB_RECURSE TYPES_SOURCES ${SRC_DIR}/types/*.cpp)
file(GLOB_RECURSE MATH_SOURCES ${SRC_DIR}/math/*.cpp)
file(GLOB_RECURSE MATRIX_SOURCES ${SRC_DIR}/matrix/*.cpp)
file(GLOB_RECURSE POLYNOMIAL_SOURCES ${SRC_DIR}/polynomial/*.cpp)
file(GLOB_RECURSE STATISTICS_SOURCES ${SRC_DIR}/statistics/*.cpp)
file(GLOB_RECURSE DSP_SOURCES ${SRC_DIR}/dsp/*.cpp)
file(GLOB_RECURSE SYMBOLIC_SOURCES ${SRC_DIR}/symbolic/*.cpp)
file(GLOB_RECURSE ANALYSIS_SOURCES ${SRC_DIR}/analysis/*.cpp)
file(GLOB_RECURSE PARSER_SOURCES ${SRC_DIR}/parser/*.cpp)
file(GLOB_RECURSE EXECUTION_SOURCES ${SRC_DIR}/execution/*.cpp)
file(GLOB_RECURSE IO_SOURCES ${SRC_DIR}/io/*.cpp)
file(GLOB_RECURSE PLOT_SOURCES ${SRC_DIR}/plot/*.cpp)
file(GLOB_RECURSE TIME_SOURCES ${SRC_DIR}/time/*.cpp)
file(GLOB_RECURSE MODULE_SOURCES ${SRC_DIR}/module/*.cpp)
file(GLOB_RECURSE CORE_SOURCES ${SRC_DIR}/core/*.cpp)

# Common Calculator Engine Library (Modular Static Target)
add_library(calc_core STATIC
    ${TYPES_SOURCES}
    ${MATH_SOURCES}
    ${MATRIX_SOURCES}
    ${POLYNOMIAL_SOURCES}
    ${STATISTICS_SOURCES}
    ${DSP_SOURCES}
    ${SYMBOLIC_SOURCES}
    ${ANALYSIS_SOURCES}
    ${PARSER_SOURCES}
    ${EXECUTION_SOURCES}
    ${IO_SOURCES}
    ${PLOT_SOURCES}
    ${TIME_SOURCES}
    ${MODULE_SOURCES}
    ${CORE_SOURCES}
)

target_include_directories(calc_core PUBLIC ${INCLUDE_DIRS})

# ============================================================================
# Applications
# ============================================================================

file(GLOB APP_SOURCES ${SRC_DIR}/app/*.cpp)
list(FILTER APP_SOURCES EXCLUDE REGEX ".*/app/main\\.cpp$")

add_executable(calculator
    ${SRC_DIR}/app/main.cpp
    ${APP_SOURCES}
)

target_include_directories(calculator PRIVATE ${INCLUDE_DIRS})
target_link_libraries(calculator PRIVATE calc_core)

# ============================================================================
# Tests
# ============================================================================

option(BUILD_TESTS "Build test executable" ON)

if(BUILD_TESTS)
    file(GLOB TEST_SUITE_SOURCES ${CMAKE_SOURCE_DIR}/test/suites/*.cpp)

    add_executable(calculator_tests
        ${CMAKE_SOURCE_DIR}/test/main.cpp
        ${TEST_SUITE_SOURCES}
        ${APP_SOURCES}
    )

    target_include_directories(calculator_tests PRIVATE ${INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/test)
    target_link_libraries(calculator_tests PRIVATE calc_core)

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
message(STATUS "")