/**
 * @file test_core_suite.cpp
 * @brief 核心功能统一测试套件（基础运算、结果显示、变量管理、解析器回归与上下文隔离）
 */

#include "test_suites.h"
#include "test_framework.h"
#include "calculator.h"
#include "core/value/stored_value.h"
#include "parser/ast/expression_ast.h"
#include "parser/grammars/command_parser.h"
#include "parser/infra/syntax_validator.h"
#include "math/mymath.h"
#include <iostream>
#include <stdexcept>
#include <vector>
#include <string>

namespace test_suites {

// 1. 上下文与多实例隔离测试
void run_context_isolation_tests(int& passed, int& failed) {
    try {
        Calculator first;
        Calculator second;
        first.set_display_precision(5);
        second.set_display_precision(20);
        if (first.display_precision() != 5 || second.display_precision() != 20) {
            throw std::runtime_error("display precision leaked between calculators");
        }

        first.process_line(":assume isolated positive", false);
        const std::string second_assumptions = second.process_line(":assume", false);
        if (second_assumptions.find("isolated") != std::string::npos) {
            throw std::runtime_error("symbolic assumptions leaked between calculators");
        }

        StoredValue rational(Rational(1, 3));
        if (!rational.is_scalar() || rational.as_scalar() == Scalar(0)) {
            throw std::runtime_error("rational scalar access is invalid");
        }
        bool rejected = false;
        try {
            StoredValue nil;
            (void)nil.as_scalar();
        } catch (const std::runtime_error&) {
            rejected = true;
        }
        if (!rejected) throw std::runtime_error("nil scalar access was not rejected");
        ++passed;
    } catch (const std::exception& error) {
        ++failed;
        std::cerr << "Context isolation test failed: " << error.what() << "\n";
    }
}

// 2. 解析器回归测试
void run_parser_regression_tests(int& passed, int& failed) {
    try {
        Calculator calculator;
        calculator.process_line("x = 3", false);
        if (calculator.evaluate("2x^3") != 54) throw std::runtime_error("implicit multiplication precedence");
        if (calculator.evaluate("-2^2") != -4) throw std::runtime_error("power/unary precedence");
        if (calculator.evaluate("50%") != Scalar(0.5L)) throw std::runtime_error("percent postfix");
        if (calculator.evaluate("50% + 1") != Scalar(1.5L)) throw std::runtime_error("percent postfix precedence");
        if (calculator.evaluate("50 % 2") != Scalar(0.0L)) throw std::runtime_error("modulo precedence");

        const auto imaginary_product = compile_expression_ast("2i");
        if (!imaginary_product || imaginary_product->kind != ExprKind::kBinaryOp ||
            imaginary_product->op_char != '*' || imaginary_product->children.size() != 2 ||
            imaginary_product->children[1]->kind != ExprKind::kVariable ||
            imaginary_product->children[1]->identifier != "i") {
            throw std::runtime_error("imaginary suffix tokenization");
        }
        if (compile_expression_ast("1e+")) throw std::runtime_error("malformed exponent accepted");

        calculator.process_line("A = [1, 2]", false);
        calculator.process_line("A[0] = 9", false);
        if (calculator.evaluate("A[0]") != 9) throw std::runtime_error("interactive index assignment");

        const auto named = parse_command("print(1, color=2)");
        const auto* named_call = named.as_function_call();
        if (!named_call || named_call->arguments.size() != 1 ||
            named_call->named_args.size() != 1 || named_call->named_args[0].name != "color" ||
            named_call->named_args[0].value.text != "2") {
            throw std::runtime_error("named argument parsing");
        }

        calculator.process_line("M = [1, 2; 3, 4]", false);
        if (calculator.evaluate_for_display("M'", false).find("[[1, 3], [2, 4]]") == std::string::npos) {
            throw std::runtime_error("transpose postfix tokenization");
        }

        for (const char* invalid : {"f(1,)", "f(,1)", "f(1,,2)", "A[] = 1", "A[0] ="}) {
            bool rejected = false;
            try {
                (void)parse_command(invalid);
            } catch (const std::exception&) {
                rejected = true;
            }
            if (!rejected) throw std::runtime_error("empty command argument accepted");
        }

        std::string output;
        calculator.try_process_function_command("f(x) = x^2 + 1", &output);
        calculator.try_process_function_command("g(x) = f(x) + f(x + 1)", &output);
        if (calculator.evaluate("g(2)") != 15) throw std::runtime_error("nested AST function binding");

        const auto matrix = calculator.evaluate_for_display("[1, 2; 3, 4]'", false);
        if (matrix.find("[[1, 3], [2, 4]]") == std::string::npos) {
            throw std::runtime_error("matrix transpose postfix");
        }

        const auto imaginary = compile_expression_ast("i");
        if (!imaginary || imaginary->kind != ExprKind::kVariable || imaginary->identifier != "i") {
            throw std::runtime_error("imaginary AST variable");
        }
        const auto dictionary = compile_expression_ast("{\"a\": 1, \"b\": 2}");
        if (!dictionary || dictionary->kind != ExprKind::kDictLiteral || dictionary->dict_entries.size() != 2) {
            throw std::runtime_error("dictionary AST node");
        }
        const auto sliced = compile_expression_ast("items[1:4:2]");
        if (!sliced || sliced->kind != ExprKind::kIndexAccess || sliced->children.size() != 2 ||
            sliced->children[1]->kind != ExprKind::kSlice) {
            throw std::runtime_error("slice AST node");
        }
        const auto list = compile_list_expression_ast("[1, 2, 3]");
        if (!list || list->kind != ExprKind::kListLiteral || list->children.size() != 3) {
            throw std::runtime_error("list AST node");
        }

        calculator.process_line("A = [1, 2]", false);
        bool rejected_fractional_index = false;
        try {
            (void)calculator.evaluate_for_display("A[1.5]", false);
        } catch (const std::exception&) {
            rejected_fractional_index = true;
        }
        if (!rejected_fractional_index) throw std::runtime_error("fractional matrix index accepted");

        SyntaxValidator validator;
        if (!validator.has_errors("1 +")) throw std::runtime_error("canonical syntax validation");
        if (validator.has_errors("{\"a\": 1}")) throw std::runtime_error("dictionary syntax validation");

        bool diagnostic = false;
        try {
            (void)compile_expression_ast_diagnostic("1 + (2\n");
        } catch (const std::exception& error) {
            const std::string message = error.what();
            diagnostic = message.find("line") != std::string::npos &&
                         message.find("column") != std::string::npos;
        }
        if (!diagnostic) throw std::runtime_error("syntax diagnostic location");
        ++passed;
    } catch (const std::exception& error) {
        ++failed;
        std::cerr << "Parser regression test failed: " << error.what() << "\n";
    }
}

// 外部子函数前向声明（来自整合后的实现）
int run_core_basic_tests(int& passed, int& failed);
int run_core_display_tests(int& passed, int& failed);
int run_core_logic_tests(int& passed, int& failed);

void run_core_suite(int& passed, int& failed) {
    run_context_isolation_tests(passed, failed);
    run_parser_regression_tests(passed, failed);
    run_core_basic_tests(passed, failed);
    run_core_display_tests(passed, failed);
    run_core_logic_tests(passed, failed);
}

} // namespace test_suites
