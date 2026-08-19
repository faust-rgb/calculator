#include "calculator.h"
#include "parser/ast/expression_ast.h"
#include "parser/infra/syntax_validator.h"
#include <iostream>
#include <stdexcept>

namespace test_suites {

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
        bool malformed_exponent_rejected = false;
        try {
            (void)compile_expression_ast("1e+");
        } catch (const std::exception&) {
            malformed_exponent_rejected = true;
        }
        if (!malformed_exponent_rejected) throw std::runtime_error("malformed exponent accepted");

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

} // namespace test_suites
