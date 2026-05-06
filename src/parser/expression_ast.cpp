/**
 * @file expression_ast.cpp
 * @brief 编译表达式 AST 实现
 *
 * 提供表达式的编译、求值和分析功能，解决循环中重复解析问题。
 */

#include "expression_ast.h"
#include "parser/base_parser.h"
#include "parser/unified_parser_factory.h"
#include "calculator_exceptions.h"
#include "calculator_internal_types.h"
#include "execution/builtin_constants.h"
#include "mymath.h"
#include "variable_resolver.h"
#include "lazy_token_stream.h"

#include <algorithm>
#include <functional>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

// ============================================================================
// AST 编译器
// ============================================================================

/**
 * @class ASTCompiler
 * @brief 将表达式字符串编译为 AST
 */
class ASTCompiler {
	public:
    ASTCompiler(const std::string& source) : tokens_(source) {
    }

    std::unique_ptr<ExpressionAST> compile() {
        if (tokens_.is_at_end()) {
            return nullptr;
        }
        auto ast = parse_logical();
        if (tokens_.peek().kind != TokenKind::kEnd) {
            throw_syntax_error("unexpected token near: " + std::string(tokens_.peek().text));
        }
        return ast;
    }

private:
    LazyTokenStream tokens_;

    bool is_at_end() {
        return tokens_.peek().kind == TokenKind::kEnd;
    }

    const Token& peek_token() {
        return tokens_.peek();
    }

    Token advance_token() {
        return tokens_.advance();
    }

    bool match_operator(std::string_view op) {
        if (!tokens_.is_at_end() && tokens_.peek().kind == TokenKind::kOperator && tokens_.peek().text == op) {
            tokens_.advance();
            return true;
        }
        return false;
    }

    bool match_kind(TokenKind kind) {
        if (!tokens_.is_at_end() && tokens_.peek().kind == kind) {
            tokens_.advance();
            return true;
        }
        return false;
    }

    void expect_kind(TokenKind kind, const std::string& msg) {
        if (!match_kind(kind)) {
            throw_syntax_error(msg);
        }
    }

    [[noreturn]] void throw_syntax_error(const std::string& message) {
        std::size_t error_pos = tokens_.is_at_end() ? 0 : tokens_.peek().position;
        std::ostringstream oss;
        oss << message << " at position " << error_pos;
        throw SyntaxError(oss.str());
    }

    std::unique_ptr<ExpressionAST> parse_logical() {
        auto left = parse_comparison();

        while (true) {
            std::string op;
            if (match_operator("&&")) op = "&&";
            else if (match_operator("||")) op = "||";
            else break;

            auto right = parse_comparison();
            auto node = std::make_unique<ExpressionAST>(ExprKind::kLogicalOp);
            node->comparison_op = op;  // 复用 comparison_op 存储逻辑运算符
            node->position = left ? left->position : 0;
            node->children.push_back(std::move(left));
            node->children.push_back(std::move(right));
            left = std::move(node);
        }
        return left;
    }

    std::unique_ptr<ExpressionAST> parse_comparison() {
        auto left = parse_expression();

        while (true) {
            std::string op;
            if (match_operator("==")) op = "==";
            else if (match_operator("!=")) op = "!=";
            else if (match_operator("<=")) op = "<=";
            else if (match_operator(">=")) op = ">=";
            else if (match_operator("<")) op = "<";
            else if (match_operator(">")) op = ">";
            else break;

            auto right = parse_expression();
            auto node = std::make_unique<ExpressionAST>(ExprKind::kComparison);
            node->comparison_op = op;
            node->position = left ? left->position : 0;
            node->children.push_back(std::move(left));
            node->children.push_back(std::move(right));
            left = std::move(node);
        }
        return left;
    }

    std::unique_ptr<ExpressionAST> parse_expression() {
        auto left = parse_term();

        while (true) {
            char op = '\0';
            if (match_operator("+")) op = '+';
            else if (match_operator("-")) op = '-';
            else break;

            auto right = parse_term();
            auto node = std::make_unique<ExpressionAST>(ExprKind::kBinaryOp);
            node->op_char = op;
            node->position = left ? left->position : 0;
            node->children.push_back(std::move(left));
            node->children.push_back(std::move(right));
            left = std::move(node);
        }
        return left;
    }

    std::unique_ptr<ExpressionAST> parse_term() {
        auto left = parse_unary();

        while (true) {
            char op = '\0';
            if (match_operator("*")) op = '*';
            else if (match_operator("/")) op = '/';
            else if (match_operator("%")) op = '%';
            else break;

            auto right = parse_unary();
            auto node = std::make_unique<ExpressionAST>(ExprKind::kBinaryOp);
            node->op_char = op;
            node->position = left ? left->position : 0;
            node->children.push_back(std::move(left));
            node->children.push_back(std::move(right));
            left = std::move(node);
        }
        return left;
    }

    std::unique_ptr<ExpressionAST> parse_power() {
        auto base = parse_primary();

        if (match_operator("^")) {
            auto exponent = parse_unary();
            auto node = std::make_unique<ExpressionAST>(ExprKind::kBinaryOp);
            node->op_char = '^';
            node->position = base ? base->position : 0;
            node->children.push_back(std::move(base));
            node->children.push_back(std::move(exponent));
            return node;
        }
        return base;
    }

    std::unique_ptr<ExpressionAST> parse_unary() {
        if (match_operator("+")) {
            return parse_unary();
        }
        if (!is_at_end() && peek_token().kind == TokenKind::kOperator && peek_token().text == "-") {
            std::size_t pos = peek_token().position;
            advance_token();
            auto operand = parse_unary();
            auto node = std::make_unique<ExpressionAST>(ExprKind::kUnaryOp);
            node->op_char = '-';
            node->position = pos;
            node->children.push_back(std::move(operand));
            return node;
        }
        return parse_power();
    }

    std::unique_ptr<ExpressionAST> parse_primary() {
        if (is_at_end()) throw_syntax_error("expected expression");
        const auto& tok = peek_token();

        // 括号表达式
        if (match_kind(TokenKind::kLParen)) {
            auto expr = parse_logical();
            expect_kind(TokenKind::kRParen, "expected ')' after expression");
            return expr;
        }

        // 标识符或函数调用
        if (tok.kind == TokenKind::kIdentifier) {
            std::string name(tok.text);
            std::size_t pos = tok.position;
            advance_token();

            // 函数调用
            if (match_kind(TokenKind::kLParen)) {
                std::vector<std::unique_ptr<ExpressionAST>> args;
                if (!is_at_end() && peek_token().kind != TokenKind::kRParen) {
                    while (true) {
                        args.push_back(parse_logical());
                        if (!match_kind(TokenKind::kComma)) break;
                    }
                }
                expect_kind(TokenKind::kRParen, "expected ')' after arguments");

                auto node = std::make_unique<ExpressionAST>(ExprKind::kFunctionCall);
                node->identifier = name;
                node->position = pos;
                node->children = std::move(args);
                return node;
            }

            // 变量引用
            auto node = std::make_unique<ExpressionAST>(ExprKind::kVariable);
            node->identifier = name;
            node->position = pos;
            return node;
        }

        // 数字字面量
        if (tok.kind == TokenKind::kNumber) {
            auto node = std::make_unique<ExpressionAST>(ExprKind::kNumber);
            node->number_value = tok.number_value;
            node->string_value = std::string(tok.text);
            node->position = tok.position;
            advance_token();
            return node;
        }

        // 字符串字面量
        if (tok.kind == TokenKind::kString) {
            auto node = std::make_unique<ExpressionAST>(ExprKind::kString);
            node->string_value = tok.string_value;
            node->position = tok.position;
            advance_token();
            return node;
        }

        throw_syntax_error("unexpected token: " + std::string(tok.text));
    }
};

// ============================================================================
// AST 求值器
// ============================================================================

/**
 * @brief 求值编译后的 AST
 */
    template <typename ExceptionType = SyntaxError>
    [[noreturn]] static void throw_ast_error(const std::string& message, std::size_t pos) {
        std::ostringstream oss;
        oss << message << " at position " << pos;
        throw ExceptionType(oss.str());
    }

double evaluate_ast(const ExpressionAST* ast,
                    const VariableResolver& variables,
                    const std::map<std::string, CustomFunction>* functions,
                    const std::map<std::string, std::function<double(const std::vector<double>&)>>* scalar_functions,
                    const HasScriptFunctionCallback& has_script_function,
                    const InvokeScriptFunctionDecimalCallback& invoke_script_function) {
    if (!ast) {
        throw MathError("null AST node");
    }

    switch (ast->kind) {
        case ExprKind::kNumber:
            return ast->number_value;

        case ExprKind::kString:
            // 字符串不能作为标量求值，但在比较中可以使用
            // 返回非零值表示字符串存在
            return ast->string_value.empty() ? 0.0 : 1.0;

        case ExprKind::kVariable: {
            // 超快路径：已绑定到槽位索引
            if (ast->is_slot_bound && ast->variable_slot_index >= 0) {
                const StoredValue* found = variables.lookup_by_slot(ast->variable_slot_index);
                if (found) {
                    if (found->is_matrix || found->is_complex ||
                        found->is_string || found->has_symbolic_text) {
                        throw_ast_error<MathError>("unsupported variable type: " + ast->identifier, ast->position);
                    }
                    return found->exact ? rational_to_double(found->rational) : found->decimal;
                }
            }

            // 快速路径：内置常量
            if (ast->is_builtin_constant) {
                return ast->number_value;
            }

            const StoredValue* found = variables.lookup(ast->identifier);
            if (found) {
                if (found->is_matrix || found->is_complex ||
                    found->is_string || found->has_symbolic_text) {
                    throw_ast_error<MathError>("unsupported variable type: " + ast->identifier, ast->position);
                }
                return found->exact ? rational_to_double(found->rational) : found->decimal;
            }
            throw_ast_error<UndefinedError>("unknown variable: " + ast->identifier, ast->position);
        }

        case ExprKind::kBinaryOp: {
            if (ast->children.size() != 2) {
                throw_ast_error<MathError>("invalid binary operation", ast->position);
            }
            double left = evaluate_ast(ast->children[0].get(), variables,
                                       functions, scalar_functions,
                                       has_script_function, invoke_script_function);
            double right = evaluate_ast(ast->children[1].get(), variables,
                                        functions, scalar_functions,
                                        has_script_function, invoke_script_function);

            switch (ast->op_char) {
                case '+': return left + right;
                case '-': return left - right;
                case '*': return left * right;
                case '/':
                    if (right == 0.0) throw_ast_error<MathError>("division by zero", ast->position);
                    return left / right;
                case '%':
                    if (right == 0.0) throw_ast_error<MathError>("modulo by zero", ast->position);
                    return mymath::fmod(left, right);
                case '^': return mymath::pow(left, right);
                default:
                    throw_ast_error<MathError>("unknown operator", ast->position);
            }
        }

        case ExprKind::kUnaryOp: {
            if (ast->children.size() != 1) {
                throw_ast_error<MathError>("invalid unary operation", ast->position);
            }
            double operand = evaluate_ast(ast->children[0].get(), variables,
                                          functions, scalar_functions,
                                          has_script_function, invoke_script_function);
            switch (ast->op_char) {
                case '-': return -operand;
                case '+': return operand;
                default:
                    throw_ast_error<MathError>("unknown unary operator", ast->position);
            }
        }

        case ExprKind::kComparison: {
            if (ast->children.size() != 2) {
                throw_ast_error<MathError>("invalid comparison", ast->position);
            }

            // 检查是否是字符串比较
            auto* left_child = ast->children[0].get();
            auto* right_child = ast->children[1].get();

            if (left_child->kind == ExprKind::kString || right_child->kind == ExprKind::kString) {
                // 至少有一个是字符串，进行字符串比较
                std::string left_str, right_str;

                if (left_child->kind == ExprKind::kString) {
                    left_str = left_child->string_value;
                } else if (left_child->kind == ExprKind::kVariable) {
                    const StoredValue* val = variables.lookup(left_child->identifier);
                    if (val && val->is_string) {
                        left_str = val->string_value;
                    } else {
                        throw_ast_error<MathError>("cannot compare string with non-string", ast->position);
                    }
                } else {
                    throw_ast_error<MathError>("cannot compare string with non-string", ast->position);
                }

                if (right_child->kind == ExprKind::kString) {
                    right_str = right_child->string_value;
                } else if (right_child->kind == ExprKind::kVariable) {
                    const StoredValue* val = variables.lookup(right_child->identifier);
                    if (val && val->is_string) {
                        right_str = val->string_value;
                    } else {
                        throw_ast_error<MathError>("cannot compare string with non-string", ast->position);
                    }
                } else {
                    throw_ast_error<MathError>("cannot compare string with non-string", ast->position);
                }

                if (ast->comparison_op == "==") return left_str == right_str ? 1.0 : 0.0;
                if (ast->comparison_op == "!=") return left_str != right_str ? 1.0 : 0.0;
                if (ast->comparison_op == "<") return left_str < right_str ? 1.0 : 0.0;
                if (ast->comparison_op == ">") return left_str > right_str ? 1.0 : 0.0;
                if (ast->comparison_op == "<=") return left_str <= right_str ? 1.0 : 0.0;
                if (ast->comparison_op == ">=") return left_str >= right_str ? 1.0 : 0.0;

                throw_ast_error<MathError>("unknown comparison operator: " + ast->comparison_op, ast->position);
            }

            // 标量比较
            double left = evaluate_ast(left_child, variables,
                                       functions, scalar_functions,
                                       has_script_function, invoke_script_function);
            double right = evaluate_ast(right_child, variables,
                                        functions, scalar_functions,
                                        has_script_function, invoke_script_function);

            if (ast->comparison_op == "==") return mymath::is_near_zero(left - right, 1e-10) ? 1.0 : 0.0;
            if (ast->comparison_op == "!=") return mymath::is_near_zero(left - right, 1e-10) ? 0.0 : 1.0;
            if (ast->comparison_op == "<") return left < right ? 1.0 : 0.0;
            if (ast->comparison_op == ">") return left > right ? 1.0 : 0.0;
            if (ast->comparison_op == "<=") return left <= right ? 1.0 : 0.0;
            if (ast->comparison_op == ">=") return left >= right ? 1.0 : 0.0;

            throw_ast_error<MathError>("unknown comparison operator: " + ast->comparison_op, ast->position);
        }

        case ExprKind::kLogicalOp: {
            if (ast->children.size() != 2) {
                throw_ast_error<MathError>("invalid logical operation", ast->position);
            }

            // 逻辑运算符支持短路求值
            double left = evaluate_ast(ast->children[0].get(), variables,
                                       functions, scalar_functions,
                                       has_script_function, invoke_script_function);

            if (ast->comparison_op == "&&") {
                // &&: 如果左边为假，直接返回 0，不计算右边
                if (left == 0.0) return 0.0;
                double right = evaluate_ast(ast->children[1].get(), variables,
                                            functions, scalar_functions,
                                            has_script_function, invoke_script_function);
                return (right != 0.0) ? 1.0 : 0.0;
            }

            if (ast->comparison_op == "||") {
                // ||: 如果左边为真，直接返回 1，不计算右边
                if (left != 0.0) return 1.0;
                double right = evaluate_ast(ast->children[1].get(), variables,
                                            functions, scalar_functions,
                                            has_script_function, invoke_script_function);
                return (right != 0.0) ? 1.0 : 0.0;
            }

            throw_ast_error<MathError>("unknown logical operator: " + ast->comparison_op, ast->position);
        }

        case ExprKind::kFunctionCall: {
            std::vector<double> args;
            args.reserve(ast->children.size());
            for (const auto& child : ast->children) {
                args.push_back(evaluate_ast(child.get(), variables,
                                            functions, scalar_functions,
                                            has_script_function, invoke_script_function));
            }

            // 自定义函数
            if (functions) {
                auto it = functions->find(ast->identifier);
                if (it != functions->end()) {
                    if (args.size() != it->second.parameter_names.size()) {
                        throw_ast_error<MathError>("custom function " + ast->identifier + " expects " +
                                        std::to_string(it->second.parameter_names.size()) + " arguments", ast->position);
                    }
                    std::map<std::string, StoredValue> snapshot = variables.snapshot();
                    for (std::size_t i = 0; i < args.size(); ++i) {
                        StoredValue arg_value;
                        arg_value.decimal = args[i];
                        snapshot[it->second.parameter_names[i]] = arg_value;
                    }
                    // 使用缓存的 AST，避免每次调用重新编译
                    VariableResolver custom_vars(&snapshot, nullptr);
                    auto compiled = it->second.get_or_compile_ast();
                    if (!compiled) {
                        throw_ast_error<MathError>("failed to compile custom function " + ast->identifier, ast->position);
                    }
                    return evaluate_ast(compiled.get(), custom_vars, functions, scalar_functions, has_script_function, invoke_script_function);
                }
            }

            // 脚本函数
            if (has_script_function && has_script_function(ast->identifier)) {
                return invoke_script_function(ast->identifier, args);
            }

            // 标量函数
            if (scalar_functions) {
                auto it = scalar_functions->find(ast->identifier);
                if (it != scalar_functions->end()) {
                    return it->second(args);
                }
            }

            throw_ast_error<UndefinedError>("unknown function: " + ast->identifier, ast->position);
        }

        case ExprKind::kConditional: {
            if (ast->children.size() != 3) {
                throw_ast_error<MathError>("invalid conditional", ast->position);
            }
            double cond = evaluate_ast(ast->children[0].get(), variables,
                                       functions, scalar_functions,
                                       has_script_function, invoke_script_function);
            if (cond != 0.0) {
                return evaluate_ast(ast->children[1].get(), variables,
                                   functions, scalar_functions,
                                   has_script_function, invoke_script_function);
            }
            return evaluate_ast(ast->children[2].get(), variables,
                               functions, scalar_functions,
                               has_script_function, invoke_script_function);
        }

        default:
            throw_ast_error<MathError>("unknown AST node kind", ast->position);
    }
}

// ============================================================================
// 变量槽位绑定
// ============================================================================

/**
 * @brief 绑定 AST 中的变量到槽位索引
 * @param ast AST 节点
 * @param variables 变量解析器
 * @return 是否所有变量都成功绑定
 */
bool bind_variable_slots(ExpressionAST* ast, const VariableResolver& variables) {
    if (!ast) return false;

    switch (ast->kind) {
        case ExprKind::kNumber:
        case ExprKind::kString:
            return true;

        case ExprKind::kVariable: {
            // Flat-scope slot indices are unstable across recursive calls and
            // nested script frames, so variables must remain dynamically
            // resolved. Built-in constants are safe to fold.
            // 检查是否为内置常量
            double value = 0.0;
            if (lookup_builtin_constant(ast->identifier, &value)) {
                ast->is_builtin_constant = true;
                ast->number_value = value;
                return true;
            }
            // 无法绑定，保持动态查找
            return true;
        }

        case ExprKind::kBinaryOp:
        case ExprKind::kUnaryOp:
        case ExprKind::kComparison:
        case ExprKind::kLogicalOp:
            for (auto& child : ast->children) {
                bind_variable_slots(child.get(), variables);
            }
            return true;

        case ExprKind::kFunctionCall:
            for (auto& child : ast->children) {
                bind_variable_slots(child.get(), variables);
            }
            return true;

        case ExprKind::kConditional:
            for (auto& child : ast->children) {
                bind_variable_slots(child.get(), variables);
            }
            return true;

        default:
            return true;
    }
}

// ============================================================================
// 公共 API
// ============================================================================

bool can_compile_to_ast(const std::string& expression) {
    if (expression.empty()) return false;
    return true;
}

std::unique_ptr<ExpressionAST> compile_expression_ast(const std::string& expression) {
    if (!can_compile_to_ast(expression)) {
        return nullptr;
    }

    ASTCompiler compiler(expression);
    try {
        return compiler.compile();
    } catch (...) {
        // Will throw SyntaxError if not matching
        throw;
    }
}

double evaluate_compiled_ast(const ExpressionAST* ast,
                             const VariableResolver& variables,
                             const std::map<std::string, CustomFunction>* functions,
                             const std::map<std::string, std::function<double(const std::vector<double>&)>>* scalar_functions,
                             const HasScriptFunctionCallback& has_script_function,
                             const InvokeScriptFunctionDecimalCallback& invoke_script_function) {
    return evaluate_ast(ast, variables, functions, scalar_functions,
                        has_script_function, invoke_script_function);
}

// ============================================================================
// 表达式分析函数（委托给 UnifiedParserFactory）
// ============================================================================

namespace {
    UnifiedParserFactory& get_global_factory() {
        static UnifiedParserFactory factory;
        return factory;
    }
}

ExpressionFeature analyze_expression_features(const std::string& expression) {
    return get_global_factory().analyze_features(expression);
}

ExpressionHint analyze_expression_hint(const std::string& expression) {
    return get_global_factory().analyze(expression).hint;
}
