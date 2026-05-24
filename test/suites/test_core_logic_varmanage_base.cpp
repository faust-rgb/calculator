/**
 * @file test_core_logic_varmanage_base.cpp
 * @brief 变量管理与进制转换测试
 *
 * 测试内容：
 * - 变量名列表
 * - 变量清除
 * - 进制转换
 */

#include "suites/test_core_logic.h"
#include "calculator.h"
#include "test_helpers.h"
#include <iostream>
#include <vector>
#include <string>
#include <exception>

namespace test_suites {

int run_logic_varmanage_base_tests(int& passed, int& failed) {
    Calculator calculator;
    using namespace test_helpers;

    // 设置一些变量用于测试
    calculator.process_line("x = 1/3 + 1/4", true);
    calculator.process_line("agg = sum(1, 2, 3)", true);

    // ========== 变量管理测试 ==========
    // 测试变量名列表
    try {
        const std::vector<std::string> variable_names = calculator.variable_names();
        bool has_x = false;
        bool has_agg = false;
        for (const std::string& name : variable_names) {
            if (name == "x") {
                has_x = true;
            }
            if (name == "agg") {
                has_agg = true;
            }
        }
        if (has_x && has_agg) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: variable_names missing expected variables\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: variable_names threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string cleared = calculator.clear_variable("x");
        if (cleared == "Cleared variable: x") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: clear_variable expected confirmation got "
                      << cleared << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: clear_variable threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        (void)calculator.clear_variable("missing");
        ++failed;
        std::cout << "FAIL: clear_variable missing expected an error but succeeded\n";
    } catch (const std::exception&) {
        ++passed;
    }

    try {
        (void)calculator.process_line("1bad = 3", false);
        ++failed;
        std::cout << "FAIL: invalid variable name expected an error but succeeded\n";
    } catch (const std::exception&) {
        ++passed;
    }

    try {
        (void)calculator.process_line("good =", false);
        ++failed;
        std::cout << "FAIL: empty assignment expected an error but succeeded\n";
    } catch (const std::exception&) {
        ++passed;
    }

    try {
        const std::string cleared_all = calculator.clear_all_variables();
        if (cleared_all == "Cleared all variables.") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: clear_all_variables expected confirmation got "
                      << cleared_all << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: clear_all_variables threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string vars_output = calculator.list_variables();
        if (vars_output == "No variables defined.") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: empty :vars expected no variables got "
                      << vars_output << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: empty list_variables threw unexpected error: "
                  << ex.what() << '\n';
    }

    // ========== 进制转换测试 ==========
    // 测试因数分解
    try {
        const std::string fact = calculator.factor_expression("factor(360)");
        if (fact == "2^3 * 3^2 * 5") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: factor expected 2^3 * 3^2 * 5 got "
                      << fact << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: factor_expression threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string converted = calculator.base_conversion_expression("hex(255)");
        if (converted == "FF") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: hex(255) expected FF got " << converted << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: base_conversion_expression threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string converted = calculator.base_conversion_expression("base(-31, 16)");
        if (converted == "-1F") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: base(-31, 16) expected -1F got " << converted << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: base_conversion_expression negative threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string status = calculator.set_hex_prefix_mode(true);
        const std::string converted = calculator.base_conversion_expression("hex(255)");
        if (status == "Hex prefix mode: ON" && converted == "0xFF") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: hex prefix mode expected 0xFF got "
                      << converted << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: hex prefix mode threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const std::string status = calculator.set_hex_uppercase_mode(false);
        const std::string converted = calculator.base_conversion_expression("base(-31, 16)");
        if (status == "Hex letter case: LOWER" && converted == "-0x1f") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: lowercase hex mode expected -0x1f got "
                      << converted << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: lowercase hex mode threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        const bool ok = calculator.hex_prefix_mode() && !calculator.hex_uppercase_mode();
        calculator.set_hex_prefix_mode(false);
        calculator.set_hex_uppercase_mode(true);
        if (ok) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: hex format getter state mismatch\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: hex format getters threw unexpected error: "
                  << ex.what() << '\n';
    }

    try {
        (void)calculator.base_conversion_expression("base(10, 1)");
        ++failed;
        std::cout << "FAIL: base invalid radix expected error but succeeded\n";
    } catch (const std::exception&) {
        ++passed;
    }

    try {
        (void)calculator.factor_expression("factor(3.5)");
        ++failed;
        std::cout << "FAIL: factor_expression non-integer expected error but succeeded\n";
    } catch (const std::exception&) {
        ++passed;
    }

    try {
        const std::string fact = calculator.factor_expression("factor(0)");
        if (fact == "0") {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: factor(0) expected 0 got " << fact << '\n';
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: factor(0) threw unexpected error: "
                  << ex.what() << '\n';
    }

    return 0;
}

}  // namespace test_suites
