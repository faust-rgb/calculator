#include "suites/test_symbolic.h"
#include "calculator.h"
#include "test_helpers.h"
#include "symbolic_expression.h"
#include <iostream>
#include <vector>
#include <string>
#include <exception>
#include <fstream>
#include <filesystem>

namespace test_suites {

int run_symbolic_tests(int& passed, int& failed) {
    Calculator calculator;
    using namespace test_helpers;
    const std::vector<DisplayCase> matrix_display_cases = {
        {"vec(1, 2, 3)", false, "[1, 2, 3]"},
        {"[1, 2, 3]", false, "[1, 2, 3]"},
        {"[1, 2; 3, 4]", false, "[[1, 2], [3, 4]]"},
        {"[1, 2; 3]", false, "[[1, 2], [3, 0]]"},
        {"[1, , 3; 4]", false, "[[1, 0, 3], [4, 0, 0]]"},
        {"mat(2, 3, 1, 2, 3, 4, 5, 6)", false, "[[1, 2, 3], [4, 5, 6]]"},
        {"zeros(2, 2)", false, "[[0, 0], [0, 0]]"},
        {"eye(3)", false, "[[1, 0, 0], [0, 1, 0], [0, 0, 1]]"},
        {"resize(mat(2, 2, 1, 2, 3, 4), 3, 3)", false, "[[1, 2, 0], [3, 4, 0], [0, 0, 0]]"},
        {"resize(mat(2, 3, 1, 2, 3, 4, 5, 6), 1, 2)", false, "[1, 2]"},
        {"append_row(mat(1, 2, 1, 2), 3, 4)", false, "[[1, 2], [3, 4]]"},
        {"append_row([1, 2], 3)", false, "[[1, 2], [3, 0]]"},
        {"append_row([1, 2], 3, 4, 5)", false, "[[1, 2, 0], [3, 4, 5]]"},
        {"append_col(mat(2, 1, 1, 2), 3, 4)", false, "[[1, 3], [2, 4]]"},
        {"append_col([1; 2], 3)", false, "[[1, 3], [2, 0]]"},
        {"append_col([1; 2], 3, 4, 5)", false, "[[1, 3], [2, 4], [0, 5]]"},
        {"mat(2, 2, 1, 2, 3, 4) + eye(2)", false, "[[2, 2], [3, 5]]"},
        {"eye(2) + 2", false, "[[3, 2], [2, 3]]"},
        {"2 + eye(2)", false, "[[3, 2], [2, 3]]"},
        {"mat(2, 2, 5, 6, 7, 8) - eye(2)", false, "[[4, 6], [7, 7]]"},
        {"10 - mat(2, 2, 1, 2, 3, 4)", false, "[[9, 8], [7, 6]]"},
        {"mat(2, 2, 1, 2, 3, 4) * 2", false, "[[2, 4], [6, 8]]"},
        {"2 * mat(2, 2, 1, 2, 3, 4)", false, "[[2, 4], [6, 8]]"},
        {"mat(2, 3, 1, 2, 3, 4, 5, 6) * mat(3, 1, 7, 8, 9)", false, "[[50], [122]]"},
        {"mat(1, 3, 1e16, 1, -1e16) * mat(3, 1, 1, 1, 1)", false, "[1]"},
        {"mat(2, 2, 2, 4, 6, 8) / 2", false, "[[1, 2], [3, 4]]"},
        {"mat(2, 2, 1, 1, 0, 1) ^ 3", false, "[[1, 3], [0, 1]]"},
        {"mat(2, 2, 1, 2, 3, 4) ^ -1", false, "[[-2, 1], [1.5, -0.5]]"},
        {"mat(2, 2, 1, 1, 0, 1) ^ -2", false, "[[1, -2], [0, 1]]"},
        {"transpose(mat(2, 3, 1, 2, 3, 4, 5, 6))", false, "[[1, 4], [2, 5], [3, 6]]"},
        {"inverse(mat(2, 2, 1, 2, 3, 4))", false, "[[-2, 1], [1.5, -0.5]]"},
        {"inverse(mat(3, 3, 1, 2, 3, 0, 1, 4, 5, 6, 0))", false, "[[-24, 18, 5], [20, -15, -4], [-5, 4, 1]]"},
        {"pinv(mat(2, 2, 1, 2, 3, 4))", false, "[[-2, 1], [1.5, -0.5]]"},
        {"dot(vec(1, 2, 3), vec(4, 5, 6))", false, "32"},
        {"dot(vec(1e16, 1, -1e16), vec(1, 1, 1))", false, "1"},
        {"outer(vec(1, 2), vec(3, 4, 5))", false, "[[3, 4, 5], [6, 8, 10]]"},
        {"kron(mat(2, 2, 1, 2, 3, 4), mat(2, 1, 5, 6))", false, "[[5, 10], [6, 12], [15, 20], [18, 24]]"},
        {"hadamard(mat(2, 2, 1, 2, 3, 4), mat(2, 2, 5, 6, 7, 8))", false, "[[5, 12], [21, 32]]"},
        {"null(mat(2, 3, 1, 2, 3, 2, 4, 6))", false, "[[-2, -3], [1, 0], [0, 1]]"},
        {"least_squares(mat(2, 1, 1, 1), vec(2, 4))", false, "[3]"},
        {"least_squares(mat(2, 3, 1, 0, 0, 0, 1, 0), vec(1, 2))", false, "[[1], [2], [0]]"},
        {"qr_q(mat(2, 2, 2, 0, 0, 3))", false, "[[1, 0], [0, 1]]"},
        {"qr_r(mat(2, 2, 2, 0, 0, 3))", false, "[[2, 0], [0, 3]]"},
        {"qr_q(mat(2, 3, 1, 0, 0, 0, 1, 0))", false, "[[1, 0], [0, 1]]"},
        {"qr_r(mat(2, 3, 1, 0, 0, 0, 1, 0))", false, "[[1, 0, 0], [0, 1, 0]]"},
        {"lu_l(mat(2, 2, 4, 3, 6, 3))", false, "[[1, 0], [0.666666666667, 1]]"},
        {"lu_u(mat(2, 2, 4, 3, 6, 3))", false, "[[6, 3], [0, 1]]"},
        {"lu_l(mat(3, 3, 2, 1, 1, 4, -6, 0, -2, 7, 2))", false, "[[1, 0, 0], [0.5, 1, 0], [-0.5, 1, 1]]"},
        {"lu_u(mat(3, 3, 2, 1, 1, 4, -6, 0, -2, 7, 2))", false, "[[4, -6, 0], [0, 4, 1], [0, 0, 1]]"},
        {"svd_u(mat(3, 2, 3, 0, 0, 2, 0, 0))", false, "[[1, 0], [0, 1], [0, 0]]"},
        {"svd_s(mat(3, 2, 3, 0, 0, 2, 0, 0))", false, "[[3, 0], [0, 2]]"},
        {"svd_vt(mat(3, 2, 3, 0, 0, 2, 0, 0))", false, "[[1, 0], [0, 1]]"},
        {"solve(mat(2, 2, 2, 1, 5, 3), vec(1, 2))", false, "[[1], [-1]]"},
        {"solve(mat(2, 2, 2, 1, 5, 3), mat(2, 1, 1, 2))", false, "[[1], [-1]]"},
        {"solve(mat(3, 3, 3, 2, -1, 2, -2, 4, -1, 0.5, -1), vec(1, -2, 0))", false, "[[1], [-2], [-2]]"},
        {"get([1, 2; 3, 4], 1, 0)", false, "3"},
        {"get(mat(2, 2, 1, 2, 3, 4), 1, 0)", false, "3"},
        {"get(vec(5, 6, 7), 2)", false, "7"},
        {"set([1, 2; 3], 1, 2, 9)", false, "[[1, 2, 0], [3, 0, 9]]"},
        {"set(mat(2, 2, 1, 2, 3, 4), 0, 1, 9)", false, "[[1, 9], [3, 4]]"},
        {"set([1, 2], 4, 9)", false, "[1, 2, 0, 0, 9]"},
        {"set(vec(5, 6, 7), 1, 42)", false, "[5, 42, 7]"},
        {"norm(vec(3, 4))", false, "5"},
        {"norm(mat(2, 2, 1, 2, 3, 4))", false, "5.47722557505"},
        {"cond(mat(2, 2, 1, 0, 0, 2))", false, "2"},
        {"pinv(mat(2, 2, 1, 2, 2, 4))", false, "[[0.04, 0.08], [0.08, 0.16]]"},
        {"trace(mat(2, 2, 1, 2, 3, 4))", false, "5"},
        {"det(mat(2, 2, 1, 2, 3, 4))", false, "-2"},
        {"det(mat(3, 3, 1, 2, 3, 0, 1, 4, 5, 6, 0))", false, "1"},
        {"rank(mat(2, 2, 1, 2, 2, 4))", false, "1"},
        {"rref(mat(2, 3, 1, 2, 3, 2, 4, 6))", false, "[[1, 2, 3], [0, 0, 0]]"},
        {"rref(mat(3, 4, 1, 2, -1, -4, 2, 3, -1, -11, -2, 0, -3, 22))", false, "[[1, 0, 0, -8], [0, 1, 0, 1], [0, 0, 1, -2]]"},
        {"eigvals(mat(2, 2, 2, 0, 0, 3))", false, "[3, 2]"},
        {"eigvecs(mat(2, 2, 2, 0, 0, 3))", false, "[[0, 1], [1, 0]]"},
        {"diag(vec(1, 2, 3))", false, "[[1, 0, 0], [0, 2, 0], [0, 0, 3]]"},
        {"diag(mat(2, 2, 1, 2, 3, 4))", false, "[[1], [4]]"},
        {"reshape(mat(2, 2, 1, 2, 3, 4), 4, 1)", false, "[[1], [2], [3], [4]]"},
        {"vec(mat(2, 2, 1, 2, 3, 4))", false, "[[1], [3], [2], [4]]"},
        {"cholesky(mat(2, 2, 4, 2, 2, 3))", false, "[[2, 0], [1, 1.41421356237]]"},
        {"poly_eval(vec(1, 2, 3), 2)", false, "17"},
        {"poly_deriv(vec(1, 2, 3))", false, "[2, 6]"},
        {"poly_integ(vec(2, 6))", false, "[0, 2, 3]"},
        {"poly_compose(vec(1, 1), vec(0, 1, 1))", false, "[1, 1, 1]"},
        {"poly_gcd(vec(-1, 0, 1), vec(-1, 1))", false, "[-1, 1]"},
        {"poly_fit(vec(0, 1, 2), vec(1, 2, 5), 2)", false, "[1, 0, 1]"},
        {"polynomial_fit(vec(0, 1, 2), vec(1, 2, 5), 2)", false, "[1, 0, 1]"},
        {"poly_fit(vec(1000, 1001, 1002), vec(1, 4, 9), 2)",
         false,
         "[998001, -1998, 1]"},
        {"lagrange(vec(0, 1, 2), vec(1, 2, 5), 1.5)", false, "3.25"},
        {"linear_regression(vec(0, 1, 2), vec(1, 3, 5))", false, "[2, 1]"},
        {"dft(mat(1, 4, 1, 0, 0, 0))", false, "[[1, 0], [1, 0], [1, 0], [1, 0]]"},
        {"idft(mat(4, 2, 1, 0, 1, 0, 1, 0, 1, 0))", false, "[1, 0, 0, 0]"},
        {"convolve(mat(1, 2, 1, 2), mat(1, 3, 3, 4, 5))", false, "[3, 10, 13, 10]"},
        {"complex(3, 4)", false, "complex(3, 4)"},
        {"real(complex(3, 4))", false, "3"},
        {"imag(complex(3, 4))", false, "4"},
        {"abs(complex(3, 4))", false, "5"},
        {"arg(complex(0, 1))", false, "1.57079632679"},
        {"conj(complex(3, 4))", false, "complex(3, -4)"},
        {"polar(2, pi / 2)", false, "complex(0, 2)"},
        {"complex(1, 2) * complex(3, 4)", false, "complex(-5, 10)"},
        {"1 / complex(0, 1)", false, "complex(0, -1)"},
        {"eigvals(mat(2, 2, 0, -1, 1, 0))", false, "[[0, 1], [0, -1]]"},
        {"mean(vec(1, 2, 3))", false, "2"},
        {"mode(vec(1, 2, 2, 3))", false, "2"},
        {"var(vec(1, 2, 3))", false, "0.666666666667"},
        {"std(vec(1, 2, 3))", false, "0.816496580928"},
        {"skewness(vec(1, 2, 3))", false, "0"},
        {"kurtosis(vec(1, 2, 3))", false, "-1.5"},
        {"cov(vec(1, 2, 3), vec(2, 4, 6))", false, "1.33333333333"},
        {"corr(vec(1, 2, 3), vec(2, 4, 6))", false, "1"},
        {"hann(5)", false, "[0, 0.5, 1, 0.5, 0]"},
        {"hamming(3)", false, "[0.08, 1, 0.08]"},
        {"blackman(3)", false, "[0, 1, 0]"},
        {"divisors(12)", false, "[1, 2, 3, 4, 6, 12]"},
        {"extended_gcd(240, 46)", false, "[2, -9, 47]"},
    };

    for (const auto& test : matrix_display_cases) {
        try {
            const std::string actual =
                calculator.evaluate_for_display(test.expression, test.exact_mode);
            if (actual == test.expected) {
                ++passed;
            } else {
                ++failed;
                std::cout << "FAIL: matrix display " << test.expression
                          << " expected " << test.expected << " got "
                          << actual << '\n';
            }
        } catch (const std::exception& ex) {
            ++failed;
            std::cout << "FAIL: matrix display " << test.expression
                      << " threw unexpected error: " << ex.what() << '\n';
        }
    }

    try {
        const double actual =
            calculator.evaluate("det(mat(2, 2, 100001, 100000, 100000, 99999))");
        if (nearly_equal(actual, -1.0, 1e-5)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: cancellation-sensitive determinant expected about -1 got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: cancellation-sensitive determinant threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string actual =
            calculator.evaluate_for_display("null(mat(2, 2, 1, 2, 3, 4))", false);
        if (actual == "[]") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: full-rank null space expected [] got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: full-rank null space threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const double actual =
            calculator.evaluate("cond(mat(2, 2, 1, 1, 1, 1))");
        if (!mymath::isfinite(actual) && actual > 0.0) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: singular cond expected inf got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: singular cond threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string assigned = calculator.process_line("m = mat(2, 2, 1, 2, 3, 4)", false);
        if (assigned == "m = [[1, 2], [3, 4]]") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: matrix assignment expected matrix display got "
                      << assigned << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: matrix assignment threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string assigned = calculator.process_line("b = [1, 2; 3]", false);
        if (assigned == "b = [[1, 2], [3, 0]]") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: bracket matrix assignment expected [[1, 2], [3, 0]] got "
                      << assigned << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: bracket matrix assignment threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string assigned = calculator.process_line("v = vec(5, 6)", true);
        if (assigned == "v = [5, 6]") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: vector assignment expected vector display got "
                      << assigned << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: vector assignment threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string assigned =
            calculator.process_line("n = m + eye(2)", false);
        if (assigned == "n = [[2, 2], [3, 5]]") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: matrix operation assignment expected [[2, 2], [3, 5]] got "
                      << assigned << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: matrix operation assignment threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string assigned =
            calculator.process_line("m = set(m, 1, 0, 8)", false);
        if (assigned == "m = [[1, 2], [8, 4]]") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: matrix set assignment expected [[1, 2], [8, 4]] got "
                      << assigned << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: matrix set assignment threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string assigned =
            calculator.process_line("b = set(b, 2, 3, 7)", false);
        if (assigned == "b = [[1, 2, 0, 0], [3, 0, 0, 0], [0, 0, 0, 7]]") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: bracket matrix set expansion expected [[1, 2, 0, 0], [3, 0, 0, 0], [0, 0, 0, 7]] got "
                      << assigned << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: bracket matrix set expansion threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string element =
            calculator.evaluate_for_display("get(m, 1, 0)", false);
        if (element == "8") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: get(matrix variable) expected 8 got "
                      << element << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: get(matrix variable) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string assigned =
            calculator.process_line("v = set(v, 0, -3)", false);
        if (assigned == "v = [-3, 6]") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: vector set assignment expected [-3, 6] got "
                      << assigned << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: vector set assignment threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string element =
            calculator.evaluate_for_display("get(v, 1)", false);
        if (element == "6") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: get(vector variable) expected 6 got "
                      << element << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: get(vector variable) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string resized =
            calculator.evaluate_for_display("resize(m, 3, 1)", false);
        if (resized == "[[1], [8], [0]]") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: resize(matrix variable) expected [[1], [8], [0]] got "
                      << resized << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: resize(matrix variable) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string vars_output = calculator.list_variables();
        if (vars_output.find("b = [[1, 2, 0, 0], [3, 0, 0, 0], [0, 0, 0, 7]]") != std::string::npos &&
            vars_output.find("m = [[1, 2], [8, 4]]") != std::string::npos &&
            vars_output.find("n = [[2, 2], [3, 5]]") != std::string::npos &&
            vars_output.find("v = [-3, 6]") != std::string::npos) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: matrix vars expected m/v listing got "
                      << vars_output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: matrix list_variables threw unexpected error: "
                  << ex.what() << '\n';
    }

    const std::vector<ErrorCase> matrix_error_cases = {
        {"mat(2, 2, 1, 2, 3)"},
        {"resize(3, 2, 2)"},
        {"transpose(3)"},
        {"inverse(mat(2, 3, 1, 2, 3, 4, 5, 6))"},
        {"inverse(mat(2, 2, 1, 2, 2, 4))"},
        {"dot(mat(2, 2, 1, 2, 3, 4), vec(1, 2))"},
        {"dot(vec(1, 2), vec(1, 2, 3))"},
        {"outer(mat(2, 2, 1, 2, 3, 4), vec(1, 2))"},
        {"null(3)"},
        {"least_squares(mat(2, 1, 1, 1), mat(2, 2, 1, 2, 3, 4))"},
        {"least_squares(mat(2, 2, 1, 2, 3, 4), vec(1, 2, 3))"},
        {"lu_l(3)"},
        {"lu_u(3)"},
        {"lu_l(mat(2, 3, 1, 2, 3, 4, 5, 6))"},
        {"svd_u(3)"},
        {"svd_s(3)"},
        {"svd_vt(3)"},
        {"solve(3, vec(1, 2))"},
        {"zeros(2.5, 2)"},
        {"mat(2, 2, 1, 2, 3, 4) + mat(1, 2, 1, 2)"},
        {"mat(2, 2, 1, 2, 3, 4) * mat(1, 2, 1, 2)"},
        {"mat(2, 2, 1, 2, 3, 4) / eye(2)"},
        {"mat(2, 3, 1, 2, 3, 4, 5, 6) ^ 2"},
        {"mat(2, 2, 1, 2, 3, 4) ^ 1.5"},
        {"mat(2, 2, 1, 2, 2, 4) ^ -1"},
        {"solve(mat(2, 2, 1, 2, 2, 4), vec(1, 2))"},
        {"solve(mat(2, 3, 1, 2, 3, 4, 5, 6), vec(1, 2))"},
        {"solve(mat(2, 2, 1, 2, 3, 4), vec(1, 2, 3))"},
        {"solve(mat(2, 2, 1, 2, 3, 4), mat(2, 2, 1, 2, 3, 4))"},
        {"get(mat(2, 2, 1, 2, 3, 4), 2, 0)"},
        {"get(mat(2, 2, 1, 2, 3, 4), 1)"},
        {"set(vec(1, 2, 3), -1, 9)"},
        {"norm(3)"},
        {"cond(3)"},
        {"pinv(3)"},
        {"kron(3, eye(2))"},
        {"hadamard(mat(2, 2, 1, 2, 3, 4), mat(1, 2, 1, 2))"},
        {"trace(mat(2, 3, 1, 2, 3, 4, 5, 6))"},
        {"det(mat(2, 3, 1, 2, 3, 4, 5, 6))"},
        {"rank(3)"},
        {"rref(3)"},
        {"eigvals(mat(2, 3, 1, 2, 3, 4, 5, 6))"},
        {"eigvecs(mat(2, 2, 0, -1, 1, 0))"},
        {"eigvecs(mat(2, 3, 1, 2, 3, 4, 5, 6))"},
        {"ode_system(vec(y2), 0, vec(1, 2), 1)"},
        {"lp_max(vec(1, 2), mat(1, 3, 1, 2, 3), vec(4), vec(0, 0), vec(1, 1))"},
        {"ilp_max(vec(1), mat(1, 1, 1), vec(2), vec(0.5), vec(3))"},
        {"milp_max(vec(1, 2), mat(1, 2, 1, 1), mat(1, 1, 3), vec(0, 0), vec(1, 1), mat(1, 1, 1))"},
        {"reshape(mat(2, 2, 1, 2, 3, 4), 3, 2)"},
        {"diag(3)"},
        {"cholesky(mat(2, 2, 1, 2, 2, 1))"},
        {"real(vec(1, 2, 3))"},
        {"arg(vec(1, 2, 3))"},
        {"poly_eval(3, 2)"},
        {"poly_fit(vec(0, 1), vec(1), 1)"},
        {"lagrange(vec(0, 1), vec(1), 0.5)"},
        {"linear_regression(vec(1, 1), vec(2, 3))"},
        {"dft(mat(2, 3, 1, 2, 3, 4, 5, 6))"},
        {"convolve(mat(2, 3, 1, 2, 3, 4, 5, 6), mat(1, 2, 1, 2))"},
    };

    for (const auto& test : matrix_error_cases) {
        try {
            (void)calculator.evaluate_for_display(test.expression, false);
            ++failed;
            std::cout << "FAIL: matrix error " << test.expression
                      << " expected an error but succeeded\n";
        } catch (const std::exception&) {
            ++passed;
        }
    }

    try {
        const double value = calculator.evaluate("rand()");
        if (value >= 0.0 && value < 1.0) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: rand() expected [0, 1) got " << value << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: rand() threw unexpected error: " << ex.what() << '\n';
    }

    try {
        const double value = calculator.evaluate("randint(2, 4)");
        if (value == 2.0 || value == 3.0 || value == 4.0) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: randint(2, 4) expected integer in range got "
                      << value << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: randint(2, 4) threw unexpected error: " << ex.what() << '\n';
    }

    const std::vector<DisplayCase> command_display_cases = {
        {"solve(x^2 - 2, 1)", false, "1.41421356237"},
        {"bisect(x^2 - 2, 1, 2)", false, "1.41421356238"},
        {"secant(x^2 - 2, 1, 2)", false, "1.41421356237"},
        {"fixed_point(cos(x), 0.5)", false, "0.73908513325"},
        {"pade(exp(x), 0, 2, 2)", false, "(1/12 * x ^ 2 + 1/2 * x + 1) / (1/12 * x ^ 2 - 1/2 * x + 1)"},
        {"puiseux((1 + x) ^ (1 / 2), 0, 4, 2)", false, "1 + 1/2 * x - 1/8 * x ^ 2"},
        {"puiseux(sqrt(x), 0, 4, 2)", false, "x ^ (1 / 2)"},
        {"puiseux(sqrt(sin(x)), 0, 4, 2)", false, "x ^ (1 / 2)"},
        {"series_sum(n^2, n, 1, N)", false, "1/3 * (N ^ 3 + 3/2 * N ^ 2 + 1/2 * N)"},
        {"summation(0.5^n, n, 0, inf)", false, "2"},
        {"laplace(step(t))", false, "1 / s"},
        {"ilaplace(1 / s)", false, "step(t)"},
        {"ilaplace(1 / (s + 2), s, t)", false, "exp(-2 * t) * step(t)"},
        {"ilaplace(3 / (2 * s + 4), s, t)", false, "3/2 * exp(-2 * t) * step(t)"},
        {"fourier(exp(-2 * t) * step(t), t, w)", false, "1 / (i * w + 2)"},
        {"fourier(exp(-2 * t) * step(t) + 3 * exp(-4 * t) * step(t), t, w)", false, "1 / (i * w + 2) + 3 * 1 / (i * w + 4)"},
        {"fourier(delta(t - 2))", false, "exp(-2 * i * w)"},
        {"ifourier(delta(w - 3))", false, "0.159154943092 * exp(3 * i * t)"},
        {"ztrans(step(n - 2))", false, "z ^ -1 / (z - 1)"},
        {"iztrans(z ^ -2)", false, "delta(n - 2)"},
        {"iztrans(z / (z - 1), z, n)", false, "step(n)"},
        {"iztrans(3 * z / (z - 1), z, n)", false, "3 * step(n)"},
        {"iztrans(2 * z / (z - 1) ^ 2, z, n)", false, "2 * step(n) * n"},
        {"iztrans(5 * z / (z - 3) - 2 * z ^ -2 + z / (z - 1), z, n)", false, "5 * step(n) * 3 ^ n - 2 * delta(n - 2) + step(n)"},
    };

    for (const auto& test : command_display_cases) {
        try {
            std::string output;
            const bool handled =
                calculator.try_process_function_command(test.expression, &output);
            const bool z_transform_equivalent =
                test.expression == "ztrans(step(n - 2))" &&
                (output == "z ^ -1 / (z - 1)" ||
                 output == "1 / (z * (z - 1))");
            const bool series_sum_equivalent =
                test.expression == "series_sum(n^2, n, 1, N)" &&
                (output == "N * (N + 1) * (2 * N + 1) / 6" ||
                 output == "N * (2 * N + 1) * (N + 1) / 6");
            if (handled && (output == test.expected ||
                            z_transform_equivalent ||
                            series_sum_equivalent)) {
                ++passed;
            } else {
                ++failed;
                std::cout << "FAIL: function command " << test.expression
                          << " expected " << test.expected << " got "
                          << output << '\n';
            }
        } catch (const std::exception& ex) {
            ++failed;
            std::cout << "FAIL: function command " << test.expression
                      << " threw unexpected error: " << ex.what() << '\n';
        }
    }

    try {
        std::string output;
        const bool handled = calculator.try_process_function_command(
            "svd(mat(2, 2, 1, 0, 0, 2))",
            &output);
        if (handled &&
            output == "U: [[0, 1], [1, 0]]\nS: [[2, 0], [0, 1]]\nVt: [[0, 1], [1, 0]]") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: svd(...) expected formatted factors got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: svd(...) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        std::string output;
        const bool handled = calculator.try_process_function_command(
            "eig(mat(2, 2, 0, -1, 1, 0))",
            &output);
        if (handled &&
            output == "values: [complex(0, 1), complex(0, -1)]\nvectors: unavailable for complex eigenvalues") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: eig(...) expected complex eigenvalue display got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: eig(...) threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator script_calculator;
        const std::string output = script_calculator.execute_script(
            "greeting = \"hello\"\n"
            "print(greeting, \"world\", 2 + 3)\n",
            false);
        if (output == "hello world 5") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: script print expected hello world 5 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: script print threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator script_calculator;
        const std::string output = script_calculator.execute_script(
            "p(x) = x^2 - 1;\n"
            "q(x) = x - 1;\n"
            "r(x) = x + 1;\n"
            "poly_mul(poly_add(p, q), r);\n",
            false);
        if (output == "x ^ 3 + 2 * x ^ 2 - x - 2") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: script nested poly command expected x ^ 3 + 2 * x ^ 2 - x - 2 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: script nested poly command threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator script_calculator;
        const std::string output = script_calculator.execute_script(
            "p(x) = x^2 - 1;\n"
            "q(x) = x - 1;\n"
            "poly_mul(poly_div(poly_mul(p, q), q), q);\n",
            false);
        if (output == "x ^ 3 - x ^ 2 - x + 1") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: script nested poly_div/poly_mul expected x ^ 3 - x ^ 2 - x + 1 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: script nested poly_div/poly_mul threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator script_calculator;
        const std::string output = script_calculator.execute_script(
            "message = \"done\"\n"
            "message\n",
            false);
        if (output == "\"done\"") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: script string expression expected \"done\" got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: script string expression threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator script_calculator;
        const std::string output = script_calculator.execute_script(
            "x = 0\n"
            "sum = 0\n"
            "while x < 5:\n"
            "  sum = sum + x\n"
            "  x = x + 1\n"
            "sum\n",
            false);
        if (output == "10") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: script while expected 10 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: script while threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator script_calculator;
        const std::string output = script_calculator.execute_script(
            "flag = 0\n"
            "if 2 > 3:\n"
            "  flag = 1\n"
            "elif 3 > 4:\n"
            "  flag = 2\n"
            "else:\n"
            "  flag = 7\n"
            "flag\n",
            false);
        if (output == "7") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: script if/else expected 7 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: script if/else threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator script_calculator;
        const std::string output = script_calculator.execute_script(
            "sum = 0\n"
            "for i in range(0, 6):\n"
            "  if i == 2:\n"
            "    continue\n"
            "  if i == 5:\n"
            "    break\n"
            "  sum = sum + i\n"
            "sum\n",
            false);
        if (output == "8") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: script for/break/continue expected 8 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: script for/break/continue threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator script_calculator;
        const std::string output = script_calculator.execute_script(
            "def fact(n):\n"
            "  if n <= 1:\n"
            "    return 1\n"
            "  else:\n"
            "    return n * fact(n - 1)\n"
            "fact(5)\n",
            false);
        if (output == "120") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: script recursive function expected 120 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: script recursive function threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator script_calculator;
        const std::string output = script_calculator.execute_script(
            "# leading comment\n"
            "def add(a, b):\n"
            "  # comment inside function body\n"
            "  return a + b\n"
            "# comment before print\n"
            "print(add(2, 5))\n",
            false);
        if (output == "7") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: script hash comment expected 7 got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: script hash comment threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator script_calculator;
        const std::string output = script_calculator.execute_script(
            "x = 10;\n"
            "fn f(n) {\n"
            "  temp = n + 1;\n"
            "  return temp;\n"
            "}\n"
            "y = f(4);\n",
            false);
        const std::string vars_output = script_calculator.list_variables();
        if (output == "y = 5" && vars_output == "x = 10\ny = 5") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: script function scope expected output y = 5 and vars x/y only got output "
                      << output << " vars " << vars_output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: script function scope threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator script_calculator;
        (void)script_calculator.execute_script(
            "fn add(a, b) {\n"
            "  return a + b;\n"
            "}\n",
            false);
        std::string output;
        const bool handled = script_calculator.try_process_function_command(":funcs", &output);
        const bool ok =
            handled &&
            output.find("add(a, b) = { ... }") != std::string::npos;
        if (ok) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: :funcs should list script functions, got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: :funcs script listing threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator script_calculator;
        (void)script_calculator.execute_script(
            "fn add(a, b) {\n"
            "  return a + b;\n"
            "}\n",
            false);
        std::string output;
        const bool handled =
            script_calculator.try_process_function_command(":clearfunc add", &output);
        if (handled && output == "Cleared custom function: add") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: :clearfunc add for script function returned "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: :clearfunc script function threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator script_calculator;
        const std::string output = script_calculator.execute_script(
            "m = [1, 2; 3];\n"
            "for (i = 0; i < 2; i = i + 1) {\n"
            "  m = set(m, i, 2, i + 7);\n"
            "}\n"
            "m;\n",
            false);
        if (output == "[[1, 2, 7], [3, 0, 8]]") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: script matrix loop expected [[1, 2, 7], [3, 0, 8]] got "
                      << output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: script matrix loop threw unexpected error: "
                  << ex.what() << '\n';
    }

    const std::string matrix_save_path =
        make_test_path("calculator_matrix_state_test.txt").string();

    try {
        Calculator matrix_save;
        (void)matrix_save.process_line("a = vec(1, 2)", false);
        (void)matrix_save.process_line("m = mat(2, 2, 1, 2, 3, 4)", false);
        (void)matrix_save.save_state(matrix_save_path);
        Calculator matrix_loaded;
        (void)matrix_loaded.load_state(matrix_save_path);
        const std::string vars = matrix_loaded.list_variables();
        if (vars == "a = [1, 2]\nm = [[1, 2], [3, 4]]") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: matrix save/load expected persisted variables got "
                      << vars << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: matrix save/load threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string cleared = calculator.clear_all_variables();
        if (cleared == "Cleared all variables.") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: clear_all_variables after matrix tests expected confirmation got "
                      << cleared << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: clear_all_variables after matrix tests threw unexpected error: "
                  << ex.what() << '\n';
    }

    const std::string save_path =
        make_test_path("calculator_state_test.txt").string();
    try {
        (void)calculator.process_line("a = 5/6", true);
        (void)calculator.process_line("b = max(4, 9)", false);
        (void)calculator.process_line("prec = 0.12345678901234567890123456789", false);
        (void)calculator.execute_script(
            "label = \"persisted\";\n"
            "fn add(a, b) {\n"
            "  return a + b;\n"
            "}\n",
            false);
        const std::string saved = calculator.save_state(save_path);
        if (saved == "Saved variables to: " + save_path) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: save_state expected confirmation got "
                      << saved << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: save_state threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator loaded;
        const std::string loaded_message = loaded.load_state(save_path);
        if (loaded_message == "Loaded variables from: " + save_path) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: load_state expected confirmation got "
                      << loaded_message << '\n';
        }

        const std::string vars_output = loaded.list_variables();
        if (vars_output == "a = 5/6\nb = 9\nlabel = \"persisted\"\nprec = 0.12345678901234567890123456789") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: loaded vars expected a = 5/6\\nb = 9\\nlabel = \"persisted\"\\nprec = 0.12345678901234567890123456789 got "
                      << vars_output << '\n';
        }

        const std::string loaded_script_output =
            loaded.execute_script("print(add(2, 3), label);\n", false);
        if (loaded_script_output == "5 persisted") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: loaded script function/string expected 5 persisted got "
                      << loaded_script_output << '\n';
        }

        const std::string loaded_precise_sum = loaded.evaluate_for_display(
            "prec + 0.00000000000000000000000000001", false);
        if (loaded_precise_sum == "0.1234567890123456789012345679") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: loaded precise decimal expected 0.1234567890123456789012345679 got "
                      << loaded_precise_sum << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: load_state threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        Calculator missing;
        (void)missing.load_state(make_test_path("calculator_no_such_state_file.txt").string());
        ++failed;
        std::cout << "FAIL: load missing file expected an error but succeeded\n";
    } catch (const std::exception&) {
        ++passed;
    }

    const std::string bad_save_path =
        make_test_path("calculator_bad_state_test.txt").string();
    {
        std::ofstream bad_file(bad_save_path);
        bad_file << "not\ta\tvalid\tstate\n";
    }

    try {
        Calculator bad_loaded;
        (void)bad_loaded.load_state(bad_save_path);
        ++failed;
        std::cout << "FAIL: load invalid file expected an error but succeeded\n";
    } catch (const std::exception&) {
        ++passed;
    }

    try {
        Calculator symbolic_persisted;
        (void)symbolic_persisted.set_symbolic_constants_mode(true);
        (void)symbolic_persisted.process_line("sym = pi / 2", false);
        const std::string symbolic_path =
            make_test_path("calculator_symbolic_state_test.txt").string();
        (void)symbolic_persisted.save_state(symbolic_path);

        Calculator symbolic_loaded;
        (void)symbolic_loaded.set_symbolic_constants_mode(true);
        (void)symbolic_loaded.load_state(symbolic_path);
        const std::string vars_output = symbolic_loaded.list_variables();
        if (vars_output == "sym = pi / 2") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: symbolic persisted vars expected sym = pi / 2 got "
                      << vars_output << '\n';
        }

        std::error_code cleanup_error;
        std::filesystem::remove(symbolic_path, cleanup_error);
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: symbolic persistence threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const double actual = calculator.evaluate(
            "ode(-1000*(y-sin(x))+cos(x), 0, 0, 0.1, 20)");
        if (mymath::isfinite(actual) && nearly_equal(actual, mymath::sin(0.1), 1e-6)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: adaptive stiff ODE expected approximately sin(0.1) got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: adaptive stiff ODE threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string actual =
            calculator.evaluate_for_display(
                "solve(mat(2,2,1,1,1,1.000000000001), vec(2,2.000000000001))",
                false);
        if (actual == "[[1], [1]]") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: near-singular solve expected [[1], [1]] got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: near-singular solve threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const double actual =
            calculator.evaluate("rank(mat(2,2,1,1,1,1.000000000001))");
        if (nearly_equal(actual, 2.0)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: near-singular rank expected 2 got "
                      << actual << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: near-singular rank threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        (void)calculator.evaluate("eigvals(mat(3,3,0,-1,0,1,0,0,0,0,2))");
        ++failed;
        std::cout << "FAIL: complex eigenvalue case expected an error but succeeded\n";
    } catch (const std::exception&) {
        ++passed;
    }

    std::error_code cleanup_error;
    std::filesystem::remove(matrix_save_path, cleanup_error);
    cleanup_error.clear();
    std::filesystem::remove(save_path, cleanup_error);
    cleanup_error.clear();
    std::filesystem::remove(bad_save_path, cleanup_error);

    // 统计扩展测试
    try {
        Calculator calc;
        const double c = calc.evaluate("cov(vec(1, 2, 3), vec(4, 5, 6))");
        const double r = calc.evaluate("corr(vec(1, 2, 3), vec(4, 5, 6))");
        if (nearly_equal(c, 2.0/3.0) && nearly_equal(r, 1.0)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: cov/corr expected 0.666..., 1.0 got " << c << ", " << r << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: cov/corr threw: " << ex.what() << '\n';
    }

    // 整数扩展测试
    try {
        Calculator calc;
        const std::string x1 = calc.process_line("xgcd(12, 8)", false);
        const std::string x2 = calc.process_line("extended_gcd(12, 8)", false);
        const std::string d = calc.process_line("divisors(12)", false);
        if (x1 == "[4, 1, -1]" && x2 == "[4, 1, -1]" && d == "[1, 2, 3, 4, 6, 12]") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: xgcd/divisors got " << x1 << ", " << d << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: xgcd/divisors threw: " << ex.what() << '\n';
    }

    // 信号窗口测试
    try {
        Calculator calc;
        const std::string h = calc.process_line("hanning(3)", false);
        if (h.find("0") != std::string::npos && h.find("1") != std::string::npos) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: hanning(3) got " << h << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: hanning threw: " << ex.what() << '\n';
    }

    // Schur 分解测试 (矩阵)
    try {
        Calculator calc;
        const std::string s = calc.process_line("schur([1, 2; 3, 4])", false);
        if (s.find("[[") != std::string::npos && s.find("]]") != std::string::npos) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: schur got " << s << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: schur threw: " << ex.what() << '\n';
    }

    std::cout << "Passed: " << passed << '\n';
    std::cout << "Failed: " << failed << '\n';
    return failed == 0 ? 0 : 1;

    return 0;
}

} // namespace test_suites
