/**
 * @file expression_ast.cpp
 * @brief 编译表达式 AST 实现
 *
 * 提供表达式的编译和求值功能，解决循环中重复解析问题。
 */

#include "expression_ast.h"
#include "base_parser.h"
#include "calculator_exceptions.h"
#include "calculator_internal_types.h"
#include "decimal_parser.h"
#include "mymath.h"
#include "variable_resolver.h"

#include <algorithm>
#include <cmath>
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
class ASTCompiler : public BaseParser {
public:
    ASTCompiler(std::string_view source)
        : BaseParser(source) {}

    std::unique_ptr<ExpressionAST> compile() {
        auto ast = parse_comparison();
        skip_spaces();
        if (!is_at_end()) {
            return nullptr;  // 无法完整解析
        }
        return ast;
    }

private:
    std::unique_ptr<ExpressionAST> parse_comparison() {
        auto left = parse_expression();
        if (!left) return nullptr;

        while (true) {
            skip_spaces();
            std::string op;

            if (match_string("==")) op = "==";
            else if (match_string("!=")) op = "!=";
            else if (match_string("<=")) op = "<=";
            else if (match_string(">=")) op = ">=";
            else if (match('<')) op = "<";
            else if (match('>')) op = ">";
            else break;

            auto right = parse_expression();
            if (!right) return nullptr;

            auto node = std::make_unique<ExpressionAST>(ExprKind::kComparison);
            node->comparison_op = op;
            node->children.push_back(std::move(left));
            node->children.push_back(std::move(right));
            left = std::move(node);
        }
        return left;
    }

    std::unique_ptr<ExpressionAST> parse_expression() {
        auto left = parse_term();
        if (!left) return nullptr;

        while (true) {
            skip_spaces();
            char op = '\0';
            if (match('+')) op = '+';
            else if (match('-')) op = '-';
            else break;

            auto right = parse_term();
            if (!right) return nullptr;

            auto node = std::make_unique<ExpressionAST>(ExprKind::kBinaryOp);
            node->op_char = op;
            node->children.push_back(std::move(left));
            node->children.push_back(std::move(right));
            left = std::move(node);
        }
        return left;
    }

    std::unique_ptr<ExpressionAST> parse_term() {
        auto left = parse_unary();
        if (!left) return nullptr;

        while (true) {
            skip_spaces();
            char op = '\0';
            if (match('*')) op = '*';
            else if (match('/')) op = '/';
            else break;

            auto right = parse_unary();
            if (!right) return nullptr;

            auto node = std::make_unique<ExpressionAST>(ExprKind::kBinaryOp);
            node->op_char = op;
            node->children.push_back(std::move(left));
            node->children.push_back(std::move(right));
            left = std::move(node);
        }
        return left;
    }

    std::unique_ptr<ExpressionAST> parse_power() {
        auto base = parse_primary();
        if (!base) return nullptr;

        skip_spaces();
        if (match('^')) {
            auto exponent = parse_unary();
            if (!exponent) return nullptr;

            auto node = std::make_unique<ExpressionAST>(ExprKind::kBinaryOp);
            node->op_char = '^';
            node->children.push_back(std::move(base));
            node->children.push_back(std::move(exponent));
            return node;
        }
        return base;
    }

    std::unique_ptr<ExpressionAST> parse_unary() {
        skip_spaces();
        if (match('+')) {
            return parse_unary();
        }
        if (match('-')) {
            auto operand = parse_unary();
            if (!operand) return nullptr;

            auto node = std::make_unique<ExpressionAST>(ExprKind::kUnaryOp);
            node->op_char = '-';
            node->children.push_back(std::move(operand));
            return node;
        }
        return parse_power();
    }

    std::unique_ptr<ExpressionAST> parse_primary() {
        skip_spaces();

        // 括号表达式
        if (match('(')) {
            auto expr = parse_expression();
            if (!expr) return nullptr;
            skip_spaces();
            if (!match(')')) return nullptr;
            return expr;
        }

        // 标识符或函数调用
        if (peek_is_alpha()) {
            std::string name(parse_identifier());
            skip_spaces();

            // 函数调用
            if (peek('(')) {
                return parse_function_call(name);
            }

            // 变量引用
            auto node = std::make_unique<ExpressionAST>(ExprKind::kVariable);
            node->identifier = name;
            return node;
        }

        // 数字字面量
        return parse_number_ast();
    }

    std::unique_ptr<ExpressionAST> parse_function_call(const std::string& name) {
        expect('(');
        std::vector<std::unique_ptr<ExpressionAST>> args;

        skip_spaces();
        if (!peek(')')) {
            while (true) {
                auto arg = parse_expression();
                if (!arg) return nullptr;
                args.push_back(std::move(arg));

                skip_spaces();
                if (!match(',')) break;
            }
        }
        skip_spaces();
        if (!match(')')) return nullptr;

        auto node = std::make_unique<ExpressionAST>(ExprKind::kFunctionCall);
        node->identifier = name;
        node->children = std::move(args);
        return node;
    }

    std::unique_ptr<ExpressionAST> parse_number_ast() {
        skip_spaces();

        // 十六进制/二进制/八进制
        if (!is_at_end() && source_[pos_] == '0' && pos_ + 1 < source_.size()) {
            int base = 10;
            if (prefixed_base(source_[pos_ + 1], &base)) {
                pos_ += 2;
                while (!is_at_end()) {
                    int digit = digit_value(source_[pos_]);
                    if (digit < 0 || digit >= base) break;
                    ++pos_;
                }
                // 简化处理：直接返回 0，实际使用时回退到 DecimalParser
                return nullptr;
            }
        }

        const std::size_t start = pos_;
        bool has_digit = false;
        bool seen_dot = false;

        while (!is_at_end()) {
            char ch = source_[pos_];
            if (std::isdigit(static_cast<unsigned char>(ch))) {
                has_digit = true;
                ++pos_;
            } else if (ch == '.' && !seen_dot) {
                seen_dot = true;
                ++pos_;
            } else {
                break;
            }
        }

        // 科学计数法
        if (!is_at_end() && (source_[pos_] == 'e' || source_[pos_] == 'E')) {
            ++pos_;
            if (!is_at_end() && (source_[pos_] == '+' || source_[pos_] == '-')) {
                ++pos_;
            }
            while (!is_at_end() && std::isdigit(static_cast<unsigned char>(source_[pos_]))) {
                ++pos_;
            }
        }

        if (!has_digit) {
            return nullptr;
        }

        try {
            double value = std::stod(std::string(source_.substr(start, pos_ - start)));
            auto node = std::make_unique<ExpressionAST>(ExprKind::kNumber);
            node->number_value = value;
            return node;
        } catch (...) {
            return nullptr;
        }
    }

    int prefixed_base(char prefix, int* base) {
        switch (prefix) {
            case 'x': case 'X': *base = 16; return 1;
            case 'b': case 'B': *base = 2; return 1;
            case 'o': case 'O': *base = 8; return 1;
            default: return 0;
        }
    }

    int digit_value(char ch) {
        if (ch >= '0' && ch <= '9') return ch - '0';
        if (ch >= 'a' && ch <= 'f') return ch - 'a' + 10;
        if (ch >= 'A' && ch <= 'F') return ch - 'A' + 10;
        return -1;
    }
};

// ============================================================================
// AST 求值器
// ============================================================================

/**
 * @brief 求值编译后的 AST
 */
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

        case ExprKind::kVariable: {
            // 超快路径：已绑定到槽位索引
            if (ast->is_slot_bound && ast->variable_slot_index >= 0) {
                const StoredValue* found = variables.lookup_by_slot(ast->variable_slot_index);
                if (found) {
                    if (found->is_matrix || found->is_complex ||
                        found->is_string || found->has_symbolic_text) {
                        throw MathError("unsupported variable type: " + ast->identifier);
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
                    throw MathError("unsupported variable type: " + ast->identifier);
                }
                return found->exact ? rational_to_double(found->rational) : found->decimal;
            }
            throw UndefinedError("unknown variable: " + ast->identifier);
        }

        case ExprKind::kBinaryOp: {
            if (ast->children.size() != 2) {
                throw MathError("invalid binary operation");
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
                    if (right == 0.0) throw MathError("division by zero");
                    return left / right;
                case '^': return mymath::pow(left, right);
                default:
                    throw MathError("unknown operator");
            }
        }

        case ExprKind::kUnaryOp: {
            if (ast->children.size() != 1) {
                throw MathError("invalid unary operation");
            }
            double operand = evaluate_ast(ast->children[0].get(), variables,
                                          functions, scalar_functions,
                                          has_script_function, invoke_script_function);
            switch (ast->op_char) {
                case '-': return -operand;
                case '+': return operand;
                default:
                    throw MathError("unknown unary operator");
            }
        }

        case ExprKind::kComparison: {
            if (ast->children.size() != 2) {
                throw MathError("invalid comparison");
            }
            double left = evaluate_ast(ast->children[0].get(), variables,
                                       functions, scalar_functions,
                                       has_script_function, invoke_script_function);
            double right = evaluate_ast(ast->children[1].get(), variables,
                                        functions, scalar_functions,
                                        has_script_function, invoke_script_function);

            if (ast->comparison_op == "==") return mymath::is_near_zero(left - right, 1e-10) ? 1.0 : 0.0;
            if (ast->comparison_op == "!=") return mymath::is_near_zero(left - right, 1e-10) ? 0.0 : 1.0;
            if (ast->comparison_op == "<") return left < right ? 1.0 : 0.0;
            if (ast->comparison_op == ">") return left > right ? 1.0 : 0.0;
            if (ast->comparison_op == "<=") return left <= right ? 1.0 : 0.0;
            if (ast->comparison_op == ">=") return left >= right ? 1.0 : 0.0;

            throw MathError("unknown comparison operator: " + ast->comparison_op);
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
                    if (args.size() != 1) {
                        throw MathError("custom function " + ast->identifier + " expects 1 argument");
                    }
                    std::map<std::string, StoredValue> snapshot = variables.snapshot();
                    StoredValue arg_value;
                    arg_value.decimal = args[0];
                    snapshot[it->second.parameter_name] = arg_value;

                    // 递归求值自定义函数（这里需要完整解析器）
                    // 简化处理：返回 0，实际使用时回退到 DecimalParser
                    throw MathError("custom function in AST not supported");
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

            throw UndefinedError("unknown function: " + ast->identifier);
        }

        case ExprKind::kConditional: {
            if (ast->children.size() != 3) {
                throw MathError("invalid conditional");
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
            throw MathError("unknown AST node kind");
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
    // 排除无法编译的表达式类型
    if (expression.empty()) return false;

    // 包含字符串字面量
    if (expression.find('"') != std::string::npos) return false;

    // 包含矩阵
    if (expression.find('[') != std::string::npos) return false;

    // 包含三元运算符（复杂，暂不支持）
    if (expression.find('?') != std::string::npos) return false;

    // 包含 rat() 调用（需要特殊处理）
    if (expression.find("rat(") != std::string::npos) return false;

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
        return nullptr;
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
