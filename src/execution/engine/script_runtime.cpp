// ============================================================================
// 脚本运行时实现文件
// ============================================================================
//
// 本文件实现了脚本执行的核心功能，包括：
// - 表达式求值（标量、矩阵、字符串、列表、字典等）
// - 命令 AST 执行
// - 脚本语句执行（if、while、for、for-in、match 等）
// - 脚本函数调用
// - 列表推导式和字典构造
// - 索引和切片操作
//
// 是计算器脚本功能的核心实现文件。
// ============================================================================

#include "script_runtime.h"
#include "calculator_module.h"
#include "parser/ast/expression_ast.h"
#include "parser/ast/expression_compiler.h"
#include "parser/grammars/unified_parser_factory.h"
#include "parser/grammars/unified_expression_parser.h"
#include "parser/grammars/symbolic_render_parser.h"
#include "parser/grammars/command_parser.h"
#include "core/services/string_utils.h"
#include "core/services/format_utils.h"
#include "parser/infra/parser_utils.h"
#include "core/services/calculator_service_factory.h"
#include "execution/engine/inline_expander.h"
#include "math/helpers/integer_helpers.h"
#include "mymath.h"
#include "script_parser.h"

#include <iomanip>
#include <algorithm>
#include <cctype>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <utility>

// ============================================================================
// 变量访问辅助函数实现
// ============================================================================

/**
 * @brief 获取当前可见变量的解析器
 * @param impl 计算器实现指针
 * @return 包含所有可见变量的 VariableResolver
 */
VariableResolver visible_variables(const Calculator::Impl* impl) {
    return VariableResolver(&impl->variables, &impl->flat_scopes);
}

/**
 * @brief 检查是否存在可见的脚本函数或原生函数
 * @param impl 计算器实现指针
 * @param name 函数名
 * @return 如果存在返回 true
 */
bool has_visible_script_function(const Calculator::Impl* impl, const std::string& name) {
    if (impl->native_functions.find(name) != impl->native_functions.end()) return true;
    return impl->script_functions.find(name) != impl->script_functions.end();
}

/**
 * @brief 分配可见变量
 * @param impl 计算器实现指针
 * @param name 变量名
 * @param value 变量值
 *
 * 按作用域优先级分配：先尝试更新现有变量，再创建新变量。
 * 优先级：局部作用域 > 全局变量表
 */
void assign_visible_variable(Calculator::Impl* impl,
                             const std::string& name,
                             const StoredValue& value) {
    // 先在局部作用域查找
    if (VariableSlot* existing = impl->flat_scopes.find(name)) {
        existing->value = value;
        return;
    }
    // 再在全局变量表查找
    if (impl->variables.find(name) != impl->variables.end()) {
        impl->variables[name] = value;
        return;
    }
    // 如果有局部作用域，在局部创建新变量
    if (impl->flat_scopes.scope_depth() > 0) {
        impl->flat_scopes.set(name, value);
        return;
    }
    // 否则在全局创建
    impl->variables[name] = value;
}

// ============================================================================
// 匿名命名空间 - 内部辅助函数
// ============================================================================

namespace {

/**
 * @brief 生成缩进文本
 * @param indent 缩进层级
 * @return 缩进字符串
 */
std::string indent_text(int indent) {
    return std::string(static_cast<std::size_t>(indent) * 2, ' ');
}

/**
 * @brief 判断值是否为真值
 * @param value 要判断的值
 * @return 如果为真返回 true
 *
 * 规则：
 * - 矩阵和复数不能用作条件
 * - 字符串：非空为真
 * - 列表/字典：非空为真
 * - 标量：非零为真
 */
bool truthy_value(const StoredValue& value) {
    if (value.is_matrix) {
        throw std::runtime_error("matrix values cannot be used as script conditions");
    }
    if (value.is_complex) {
        throw std::runtime_error("complex values cannot be used as script conditions");
    }
    if (value.is_string) {
        return !value.string_value.empty();
    }
    if (value.is_list) {
        return value.list_value && !value.list_value->empty();
    }
    if (value.is_dict) {
        return value.dict_value && !value.dict_value->empty();
    }
    return !mymath::is_near_zero(value.exact
                                     ? rational_to_double(value.rational)
                                     : value.decimal,
                                 1e-10);
}

/**
 * @brief 求值 AST 并转换为标量值
 * @param calculator 计算器实例
 * @param impl 计算器实现指针
 * @param ast 命令 AST 节点
 * @param context 上下文描述（用于错误信息）
 * @return 标量值
 * @throws 如果值不是标量则抛出异常
 */
Scalar evaluate_scalar(Calculator* calculator, Calculator::Impl* impl, const CommandASTNode& ast, const char* context) {
    StoredValue val = evaluate_command_ast_to_value(calculator, impl, ast, false);
    if (val.is_matrix || val.is_complex || val.is_string) {
        throw std::runtime_error(std::string(context) + " must be a scalar");
    }
    return val.exact ? rational_to_double(val.rational) : val.decimal;
}

/**
 * @brief 检查文本是否被指定字符包裹
 */
bool is_wrapped_by(std::string_view text, char open, char close) {
    return parser_utils::is_wrapped_by(text, open, close);
}

/**
 * @brief 在顶层分割文本（忽略嵌套结构）
 * @param text 输入文本
 * @param delimiter 分隔符
 * @return 分割后的字符串列表
 */
std::vector<std::string> split_script_top_level(std::string_view text, char delimiter) {
    auto parts = parser_utils::split_top_level(text, delimiter);
    for (auto& p : parts) p = trim_copy(p);
    return parts;
}

/**
 * @brief 在顶层查找单词（忽略嵌套结构和字符串）
 * @param text 输入文本
 * @param word 要查找的单词
 * @param start_at 起始位置
 * @return 单词位置，未找到返回 npos
 */
std::size_t find_top_level_word(std::string_view text, std::string_view word, std::size_t start_at = 0) {
    int paren = 0, bracket = 0, brace = 0;
    bool in_string = false, escaping = false;
    for (std::size_t i = start_at; i + word.size() <= text.size(); ++i) {
        const char ch = text[i];
        if (in_string) {
            if (escaping) escaping = false;
            else if (ch == '\\') escaping = true;
            else if (ch == '"') in_string = false;
            continue;
        }
        if (ch == '"') { in_string = true; continue; }
        if (ch == '(') ++paren;
        else if (ch == ')' && paren > 0) --paren;
        else if (ch == '[') ++bracket;
        else if (ch == ']' && bracket > 0) --bracket;
        else if (ch == '{') ++brace;
        else if (ch == '}' && brace > 0) --brace;
        if (paren == 0 && bracket == 0 && brace == 0 &&
            text.substr(i, word.size()) == word) {
            const bool left_ok = i == 0 ||
                !(std::isalnum(static_cast<unsigned char>(text[i - 1])) || text[i - 1] == '_');
            const std::size_t after = i + word.size();
            const bool right_ok = after >= text.size() ||
                !(std::isalnum(static_cast<unsigned char>(text[after])) || text[after] == '_');
            if (left_ok && right_ok) return i;
        }
    }
    return std::string::npos;
}

/**
 * @brief 在顶层查找赋值符号（忽略嵌套结构）
 * @param text 输入文本
 * @return 赋值符号位置，未找到返回 npos
 *
 * 排除比较运算符（==, !=, <=, >=）。
 */
std::size_t find_top_level_assignment(std::string_view text) {
    int paren = 0, bracket = 0, brace = 0;
    bool in_string = false, escaping = false;
    for (std::size_t i = 0; i < text.size(); ++i) {
        const char ch = text[i];
        if (in_string) {
            if (escaping) escaping = false;
            else if (ch == '\\') escaping = true;
            else if (ch == '"') in_string = false;
            continue;
        }
        if (ch == '"') { in_string = true; continue; }
        if (ch == '(') ++paren;
        else if (ch == ')' && paren > 0) --paren;
        else if (ch == '[') ++bracket;
        else if (ch == ']' && bracket > 0) --bracket;
        else if (ch == '{') ++brace;
        else if (ch == '}' && brace > 0) --brace;
        else if (ch == '=' && paren == 0 && bracket == 0 && brace == 0) {
            const char prev = i > 0 ? text[i - 1] : '\0';
            const char next = i + 1 < text.size() ? text[i + 1] : '\0';
            if (prev != '=' && prev != '!' && prev != '<' && prev != '>' && next != '=') return i;
        }
    }
    return std::string::npos;
}

/**
 * @brief 将 StoredValue 转换为整数索引
 * @param value 存储值
 * @param context 上下文描述（用于错误信息）
 * @return 整数索引
 * @throws 如果值不是整数则抛出异常
 */
long long stored_to_index(const StoredValue& value, const char* context) {
    if (value.is_matrix || value.is_complex || value.is_string || value.is_list || value.is_dict) {
        throw std::runtime_error(std::string(context) + " must be an integer");
    }
    const Scalar scalar = value.exact ? rational_to_double(value.rational) : value.decimal;
    if (!is_integer_double(scalar)) {
        throw std::runtime_error(std::string(context) + " must be an integer");
    }
    return round_to_long_long(scalar);
}

/**
 * @brief 将 StoredValue 转换为字典键
 * @param value 存储值
 * @return 字符串形式的键
 * @throws 如果值不能用作键则抛出异常
 */
std::string stored_to_key(const StoredValue& value) {
    if (value.is_string) return value.string_value;
    if (value.is_matrix || value.is_complex || value.is_list || value.is_dict) {
        throw std::runtime_error("dictionary key must be a string or scalar");
    }
    return value.exact ? value.rational.to_string() : format_stored_value(value, false);
}

/**
 * @brief 解析索引表达式
 * @param expression 表达式字符串
 * @param base 输出：基础表达式
 * @param index 输出：索引表达式
 * @return 如果成功解析返回 true
 *
 * 例如："arr[0]" -> base="arr", index="0"
 */
bool parse_index_expression(std::string_view expression, std::string* base, std::string* index) {
    const std::string text = trim_copy(expression);
    if (text.empty() || text.back() != ']') return false;
    int paren = 0, bracket = 0, brace = 0;
    bool in_string = false, escaping = false;
    std::size_t open_pos = std::string::npos;
    for (std::size_t i = 0; i < text.size(); ++i) {
        const char ch = text[i];
        if (in_string) {
            if (escaping) escaping = false;
            else if (ch == '\\') escaping = true;
            else if (ch == '"') in_string = false;
            continue;
        }
        if (ch == '"') { in_string = true; continue; }
        if (ch == '(') ++paren;
        else if (ch == ')' && paren > 0) --paren;
        else if (ch == '{') ++brace;
        else if (ch == '}' && brace > 0) --brace;
        else if (ch == '[') {
            if (paren == 0 && bracket == 0 && brace == 0) {
                open_pos = i;
            }
            ++bracket;
        } else if (ch == ']' && bracket > 0) {
            --bracket;
            if (bracket == 0 && i + 1 != text.size()) {
                open_pos = std::string::npos;
            }
        }
    }
    if (open_pos == std::string::npos) return false;
    *base = trim_copy(text.substr(0, open_pos));
    *index = trim_copy(text.substr(open_pos + 1, text.size() - open_pos - 2));
    return !base->empty();
}

/**
 * @brief 检查文本是否有顶层分号
 */
bool has_top_level_semicolon(std::string_view text) {
    return split_script_top_level(text, ';').size() > 1;
}

// 前向声明
StoredValue evaluate_script_value_expression(Calculator* calculator,
                                             Calculator::Impl* impl,
                                             const std::string& expression,
                                             bool exact_mode);

// ============================================================================
// 列表和字典构造辅助函数
// ============================================================================

/**
 * @brief 创建列表值
 * @param values 元素列表
 * @return 列表类型的 StoredValue
 */
StoredValue make_list_value(std::vector<StoredValue> values) {
    StoredValue stored;
    stored.is_list = true;
    stored.list_value = std::make_shared<std::vector<StoredValue>>(std::move(values));
    return stored;
}

/**
 * @brief 创建字典值
 * @param values 键值对映射
 * @return 字典类型的 StoredValue
 */
StoredValue make_dict_value(std::map<std::string, StoredValue> values) {
    StoredValue stored;
    stored.is_dict = true;
    stored.dict_value = std::make_shared<std::map<std::string, StoredValue>>(std::move(values));
    return stored;
}

/**
 * @brief 求值 range() 函数调用
 * @param calculator 计算器实例
 * @param impl 计算器实现指针
 * @param args_text 参数文本
 * @param exact_mode 是否精确模式
 * @return 包含范围内数值的列表
 *
 * 支持 1-3 个参数：range(stop)、range(start, stop)、range(start, stop, step)
 */
StoredValue evaluate_range_list(Calculator* calculator,
                                Calculator::Impl* impl,
                                const std::string& args_text,
                                bool exact_mode) {
    std::vector<std::string> args;
    if (!trim_copy(args_text).empty()) args = split_script_top_level(args_text, ',');
    if (args.empty() || args.size() > 3) throw std::runtime_error("range expects 1-3 arguments");
    auto eval_scalar = [&](const std::string& arg) {
        StoredValue value = evaluate_script_value_expression(calculator, impl, arg, exact_mode);
        return (stored_to_index(value, "range argument"));
    };
    Scalar start = 0.0L, stop = 0.0L, step = 1.0L;
    if (args.size() == 1) {
        stop = eval_scalar(args[0]);
    } else if (args.size() == 2) {
        start = eval_scalar(args[0]);
        stop = eval_scalar(args[1]);
    } else {
        start = eval_scalar(args[0]);
        stop = eval_scalar(args[1]);
        step = eval_scalar(args[2]);
    }
    if (step == 0.0L) throw std::runtime_error("range step cannot be zero");
    std::vector<StoredValue> values;
    for (Scalar current = start; step > 0 ? current < stop : current > stop; current += step) {
        StoredValue item;
        item.decimal = current;
        values.push_back(item);
    }
    return make_list_value(std::move(values));
}

/**
 * @brief 求值列表推导式
 * @param calculator 计算器实例
 * @param impl 计算器实现指针
 * @param body 推导式体（如 "x*2 for x in list if x > 0"）
 * @param exact_mode 是否精确模式
 * @return 结果列表
 *
 * 语法：[expr for var in iterable if condition]
 * 其中 if condition 是可选的。
 */
StoredValue evaluate_list_comprehension(Calculator* calculator,
                                        Calculator::Impl* impl,
                                        const std::string& body,
                                        bool exact_mode) {
    const std::size_t for_pos = find_top_level_word(body, "for");
    if (for_pos == std::string::npos) throw std::runtime_error("invalid list comprehension");
    const std::string element_expr = trim_copy(body.substr(0, for_pos));
    const std::string rest = trim_copy(body.substr(for_pos + 3));
    const std::size_t in_pos = find_top_level_word(rest, "in");
    if (in_pos == std::string::npos) throw std::runtime_error("invalid list comprehension");
    const std::string var = trim_copy(rest.substr(0, in_pos));
    if (var.empty()) throw std::runtime_error("invalid list comprehension variable");
    std::string iterable_expr = trim_copy(rest.substr(in_pos + 2));
    std::string filter_expr;
    const std::size_t if_pos = find_top_level_word(iterable_expr, "if");
    if (if_pos != std::string::npos) {
        filter_expr = trim_copy(iterable_expr.substr(if_pos + 2));
        iterable_expr = trim_copy(iterable_expr.substr(0, if_pos));
    }

    StoredValue iterable = evaluate_script_value_expression(calculator, impl, iterable_expr, exact_mode);
    if (!iterable.is_list || !iterable.list_value) {
        throw std::runtime_error("list comprehension requires a list iterable");
    }

    std::vector<StoredValue> result;
    impl->flat_scopes.push_scope();
    try {
        for (const StoredValue& item : *iterable.list_value) {
            impl->flat_scopes.set(var, item);
            if (!filter_expr.empty() &&
                !truthy_value(evaluate_script_value_expression(calculator, impl, filter_expr, exact_mode))) {
                continue;
            }
            result.push_back(evaluate_script_value_expression(calculator, impl, element_expr, exact_mode));
        }
        impl->flat_scopes.pop_scope();
    } catch (...) {
        impl->flat_scopes.pop_scope();
        throw;
    }
    return make_list_value(std::move(result));
}

/**
 * @brief 求值列表字面量
 * @param calculator 计算器实例
 * @param impl 计算器实现指针
 * @param expression 列表表达式（如 "[1, 2, 3]"）
 * @param exact_mode 是否精确模式
 * @return 列表值
 *
 * 支持列表推导式语法。
 */
StoredValue evaluate_list_literal(Calculator* calculator,
                                  Calculator::Impl* impl,
                                  const std::string& expression,
                                  bool exact_mode) {
    const std::string inner = trim_copy(expression.substr(1, expression.size() - 2));
    if (inner.empty()) return make_list_value({});
    if (has_top_level_semicolon(inner)) {
        throw std::runtime_error("not a script list literal");
    }
    if (find_top_level_word(inner, "for") != std::string::npos) {
        return evaluate_list_comprehension(calculator, impl, inner, exact_mode);
    }
    std::vector<StoredValue> values;
    for (const std::string& part : split_script_top_level(inner, ',')) {
        if (!part.empty()) values.push_back(evaluate_script_value_expression(calculator, impl, part, exact_mode));
    }
    return make_list_value(std::move(values));
}

/**
 * @brief 求值字典字面量
 * @param calculator 计算器实例
 * @param impl 计算器实现指针
 * @param expression 字典表达式（如 "{a: 1, b: 2}"）
 * @param exact_mode 是否精确模式
 * @return 字典值
 */
StoredValue evaluate_dict_literal(Calculator* calculator,
                                  Calculator::Impl* impl,
                                  const std::string& expression,
                                  bool exact_mode) {
    const std::string inner = trim_copy(expression.substr(1, expression.size() - 2));
    std::map<std::string, StoredValue> values;
    if (inner.empty()) return make_dict_value(std::move(values));
    for (const std::string& item : split_script_top_level(inner, ',')) {
        const std::vector<std::string> kv = split_script_top_level(item, ':');
        if (kv.size() != 2) throw std::runtime_error("invalid dictionary literal");
        StoredValue key = evaluate_script_value_expression(calculator, impl, kv[0], exact_mode);
        values[stored_to_key(key)] = evaluate_script_value_expression(calculator, impl, kv[1], exact_mode);
    }
    return make_dict_value(std::move(values));
}

/**
 * @brief 求值索引或切片操作
 * @param calculator 计算器实例
 * @param impl 计算器实现指针
 * @param expression 索引表达式（如 "arr[0]" 或 "list[1:3]"）
 * @param exact_mode 是否精确模式
 * @return 索引或切片结果
 *
 * 支持列表、字典、矩阵的索引访问，以及列表的切片操作。
 */
StoredValue evaluate_index_or_slice(Calculator* calculator,
                                    Calculator::Impl* impl,
                                    const std::string& expression,
                                    bool exact_mode) {
    std::string base_expr;
    std::string index_expr;
    if (!parse_index_expression(expression, &base_expr, &index_expr)) {
        throw std::runtime_error("invalid index expression");
    }
    StoredValue base = evaluate_script_value_expression(calculator, impl, base_expr, exact_mode);
    if (base.is_list) {
        if (!base.list_value) throw std::runtime_error("invalid list value");
        const auto& list = *base.list_value;
        const std::vector<std::string> slice_parts = split_script_top_level(index_expr, ':');
        if (slice_parts.size() > 1) {
            if (slice_parts.size() > 3) throw std::runtime_error("slice expects start:stop[:step]");
            long long start = 0;
            long long stop = static_cast<long long>(list.size());
            long long step = 1;
            if (!trim_copy(slice_parts[0]).empty()) start = stored_to_index(evaluate_script_value_expression(calculator, impl, slice_parts[0], exact_mode), "slice start");
            if (slice_parts.size() > 1 && !trim_copy(slice_parts[1]).empty()) stop = stored_to_index(evaluate_script_value_expression(calculator, impl, slice_parts[1], exact_mode), "slice stop");
            if (slice_parts.size() > 2 && !trim_copy(slice_parts[2]).empty()) step = stored_to_index(evaluate_script_value_expression(calculator, impl, slice_parts[2], exact_mode), "slice step");
            if (step == 0) throw std::runtime_error("slice step cannot be zero");
            auto normalize = [&](long long idx) {
                if (idx < 0) idx += static_cast<long long>(list.size());
                return std::max(0LL, std::min(idx, static_cast<long long>(list.size())));
            };
            start = normalize(start);
            stop = normalize(stop);
            std::vector<StoredValue> result;
            if (step > 0) {
                for (long long i = start; i < stop; i += step) result.push_back(list[static_cast<std::size_t>(i)]);
            } else {
                for (long long i = start; i > stop; i += step) result.push_back(list[static_cast<std::size_t>(i)]);
            }
            return make_list_value(std::move(result));
        }
        long long index = stored_to_index(evaluate_script_value_expression(calculator, impl, index_expr, exact_mode), "list index");
        if (index < 0) index += static_cast<long long>(list.size());
        if (index < 0 || index >= static_cast<long long>(list.size())) throw std::runtime_error("list index out of range");
        return list[static_cast<std::size_t>(index)];
    }
    if (base.is_dict) {
        if (!base.dict_value) throw std::runtime_error("invalid dictionary value");
        const std::string key = stored_to_key(evaluate_script_value_expression(calculator, impl, index_expr, exact_mode));
        auto it = base.dict_value->find(key);
        if (it == base.dict_value->end()) throw std::runtime_error("dictionary key not found: " + key);
        return it->second;
    }
    if (base.is_matrix) {
        const std::vector<std::string> parts = split_script_top_level(index_expr, ',');
        if (parts.size() == 1) {
            long long index = stored_to_index(evaluate_script_value_expression(calculator, impl, parts[0], exact_mode), "matrix index");
            if (index < 0) index += static_cast<long long>(base.matrix.data.size());
            if (index < 0 || index >= static_cast<long long>(base.matrix.data.size())) throw std::runtime_error("matrix index out of range");
            StoredValue result;
            result.decimal = base.matrix.data[static_cast<std::size_t>(index)];
            return result;
        } else if (parts.size() == 2) {
            long long row = stored_to_index(evaluate_script_value_expression(calculator, impl, parts[0], exact_mode), "matrix row");
            long long col = stored_to_index(evaluate_script_value_expression(calculator, impl, parts[1], exact_mode), "matrix col");
            if (row < 0) row += static_cast<long long>(base.matrix.rows);
            if (col < 0) col += static_cast<long long>(base.matrix.cols);
            if (row < 0 || row >= static_cast<long long>(base.matrix.rows) || col < 0 || col >= static_cast<long long>(base.matrix.cols)) {
                throw std::runtime_error("matrix index out of range");
            }
            StoredValue result;
            result.decimal = base.matrix.at(static_cast<std::size_t>(row), static_cast<std::size_t>(col));
            return result;
        } else {
            throw std::runtime_error("matrix indexing requires 1 or 2 indices");
        }
    }
    throw std::runtime_error("indexing requires a list, dictionary, or matrix");
}

/**
 * @brief 求值脚本值表达式
 * @param calculator 计算器实例
 * @param impl 计算器实现指针
 * @param expression 表达式字符串
 * @param exact_mode 是否精确模式
 * @return 求值结果
 *
 * 支持多种表达式类型：
 * - 内置函数：len(), range(), append()
 * - 列表字面量：[1, 2, 3]
 * - 字典字面量：{a: 1, b: 2}
 * - 索引/切片：arr[0], list[1:3]
 * - 普通数学表达式
 */
StoredValue evaluate_script_value_expression(Calculator* calculator,
                                             Calculator::Impl* impl,
                                             const std::string& expression,
                                             bool exact_mode) {
    const std::string trimmed = trim_copy(expression);

    // len() 函数 - 返回列表、字典或字符串的长度
    if (trimmed.rfind("len(", 0) == 0 && trimmed.back() == ')') {
        std::string arg = trim_copy(trimmed.substr(4, trimmed.size() - 5));
        if (arg.empty()) {
            throw std::runtime_error("len() requires one argument");
        }
        StoredValue value = evaluate_script_value_expression(calculator, impl, arg, exact_mode);
        StoredValue result;
        result.decimal = Scalar(0.0L);
        if (value.is_list && value.list_value) {
            result.decimal = Scalar(static_cast<long long>(value.list_value->size()));
        } else if (value.is_dict && value.dict_value) {
            result.decimal = Scalar(static_cast<long long>(value.dict_value->size()));
        } else if (value.is_string) {
            result.decimal = Scalar(static_cast<long long>(value.string_value.size()));
        } else {
            throw std::runtime_error("len() requires a list, dict, or string");
        }
        return result;
    }

    // range() 函数 - 生成数值范围列表
    if (trimmed.rfind("range(", 0) == 0 && is_wrapped_by(trimmed.substr(5), '(', ')')) {
        return evaluate_range_list(calculator, impl, trimmed.substr(6, trimmed.size() - 7), exact_mode);
    }

    // append() 函数 - 返回追加元素后的新列表
    if (trimmed.rfind("append(", 0) == 0 && trimmed.back() == ')') {
        std::string args_text = trim_copy(trimmed.substr(7, trimmed.size() - 8));
        std::vector<std::string> args = split_script_top_level(args_text, ',');
        if (args.size() != 2) {
            throw std::runtime_error("append() requires two arguments: list and element");
        }
        StoredValue list_val = evaluate_script_value_expression(calculator, impl, args[0], exact_mode);
        if (!list_val.is_list || !list_val.list_value) {
            throw std::runtime_error("append() first argument must be a list");
        }
        StoredValue elem = evaluate_script_value_expression(calculator, impl, args[1], exact_mode);
        std::vector<StoredValue> new_list = *list_val.list_value;
        new_list.push_back(elem);
        return make_list_value(std::move(new_list));
    }

    if (is_wrapped_by(trimmed, '[', ']')) {
        const std::string inner = trim_copy(trimmed.substr(1, trimmed.size() - 2));
        if (!has_top_level_semicolon(inner)) {
            return evaluate_list_literal(calculator, impl, trimmed, exact_mode);
        }
    }
    if (is_wrapped_by(trimmed, '{', '}')) {
        return evaluate_dict_literal(calculator, impl, trimmed, exact_mode);
    }
    std::string base;
    std::string index;
    if (parse_index_expression(trimmed, &base, &index)) {
        return evaluate_index_or_slice(calculator, impl, trimmed, exact_mode);
    }
    return evaluate_expression_value(calculator, impl, trimmed, exact_mode);
}

bool try_execute_index_assignment(Calculator* calculator,
                                  Calculator::Impl* impl,
                                  const std::string& text,
                                  bool exact_mode,
                                  std::string* output) {
    const std::size_t eq = find_top_level_assignment(text);
    if (eq == std::string::npos) return false;
    std::string target = trim_copy(text.substr(0, eq));
    std::string value_expr = trim_copy(text.substr(eq + 1));
    std::string base_name;
    std::string index_expr;
    if (!parse_index_expression(target, &base_name, &index_expr)) return false;
    if (base_name.empty() || (!std::isalpha(static_cast<unsigned char>(base_name[0])) && base_name[0] != '_')) return false;

    VariableSlot* slot = impl->flat_scopes.find(base_name);
    StoredValue* base_value = slot ? &slot->value : nullptr;
    if (!base_value) {
        auto it = impl->variables.find(base_name);
        if (it != impl->variables.end()) base_value = &it->second;
    }
    if (!base_value) throw std::runtime_error("unknown variable: " + base_name);

    StoredValue new_value = evaluate_script_value_expression(calculator, impl, value_expr, exact_mode);
    if (base_value->is_list) {
        if (!base_value->list_value) base_value->list_value = std::make_shared<std::vector<StoredValue>>();
        long long index = stored_to_index(evaluate_script_value_expression(calculator, impl, index_expr, exact_mode), "list index");
        auto& list = *base_value->list_value;
        if (index < 0) index += static_cast<long long>(list.size());
        if (index < 0 || index >= static_cast<long long>(list.size())) throw std::runtime_error("list index out of range");
        list[static_cast<std::size_t>(index)] = new_value;
        *output = base_name + "[" + index_expr + "] = " + format_stored_value(new_value, impl->symbolic_constants_mode);
        return true;
    }
    if (base_value->is_dict) {
        if (!base_value->dict_value) base_value->dict_value = std::make_shared<std::map<std::string, StoredValue>>();
        const std::string key = stored_to_key(evaluate_script_value_expression(calculator, impl, index_expr, exact_mode));
        (*base_value->dict_value)[key] = new_value;
        *output = base_name + "[" + index_expr + "] = " + format_stored_value(new_value, impl->symbolic_constants_mode);
        return true;
    }
    if (base_value->is_matrix) {
        if (new_value.is_matrix || new_value.is_list || new_value.is_dict || new_value.is_string) {
            throw std::runtime_error("matrix element assignment requires a scalar value");
        }
        Scalar val = new_value.exact ? rational_to_double(new_value.rational) : new_value.decimal;
        const std::vector<std::string> parts = split_script_top_level(index_expr, ',');
        if (parts.size() == 1) {
            long long index = stored_to_index(evaluate_script_value_expression(calculator, impl, parts[0], exact_mode), "matrix index");
            if (index < 0) index += static_cast<long long>(base_value->matrix.data.size());
            if (index < 0 || index >= static_cast<long long>(base_value->matrix.data.size())) throw std::runtime_error("matrix index out of range");
            base_value->matrix.data[static_cast<std::size_t>(index)] = val;
            *output = base_name + "[" + index_expr + "] = " + format_stored_value(new_value, impl->symbolic_constants_mode);
            return true;
        } else if (parts.size() == 2) {
            long long row = stored_to_index(evaluate_script_value_expression(calculator, impl, parts[0], exact_mode), "matrix row");
            long long col = stored_to_index(evaluate_script_value_expression(calculator, impl, parts[1], exact_mode), "matrix col");
            if (row < 0) row += static_cast<long long>(base_value->matrix.rows);
            if (col < 0) col += static_cast<long long>(base_value->matrix.cols);
            if (row < 0 || row >= static_cast<long long>(base_value->matrix.rows) || col < 0 || col >= static_cast<long long>(base_value->matrix.cols)) {
                throw std::runtime_error("matrix index out of range");
            }
            base_value->matrix.at(static_cast<std::size_t>(row), static_cast<std::size_t>(col)) = val;
            *output = base_name + "[" + index_expr + "] = " + format_stored_value(new_value, impl->symbolic_constants_mode);
            return true;
        } else {
            throw std::runtime_error("matrix indexing requires 1 or 2 indices");
        }
    }
    throw std::runtime_error("indexed assignment requires a list, dictionary, or matrix");
}

} // namespace

// ============================================================================
// 表达式求值
// ============================================================================

StoredValue evaluate_expression_value(Calculator* calculator,
                                      Calculator::Impl* impl,
                                      const std::string& expression,
                                      bool exact_mode,
                                      std::shared_ptr<ExpressionCache>* cache) {
    const std::string special_expr = trim_copy(expression);
    std::string index_base;
    std::string index_expr;
    if (parse_index_expression(special_expr, &index_base, &index_expr)) {
        return evaluate_index_or_slice(calculator, impl, special_expr, exact_mode);
    }

    const VariableResolver variables = visible_variables(impl);
    
    std::shared_ptr<ExpressionCache> expr_cache;
    if (cache) {
        if (!*cache) {
            *cache = std::make_shared<ExpressionCache>(expression);
            std::string expanded = expand_inline_function_commands(calculator, expression);
            if (expanded != expression) (*cache)->set_expanded(std::move(expanded));
        } else if (!(*cache)->has_expanded) {
            std::string expanded = expand_inline_function_commands(calculator, std::string((*cache)->original_text));
            if (expanded != (*cache)->original_text) (*cache)->set_expanded(std::move(expanded));
        }
        expr_cache = *cache;
    }

    const std::string target_expr =
        expr_cache ? std::string(expr_cache->effective_text())
                   : expand_inline_function_commands(calculator, expression);

    std::string trimmed_expr = trim_copy(target_expr);
    if (is_string_literal(trimmed_expr)) {
        StoredValue res; res.is_string = true; res.string_value = parse_string_literal_value(trimmed_expr);
        return res;
    }

    // 检查是否是原生函数调用（顶级函数调用）
    // 避免与 UnifiedExpressionParser 冲突
    auto is_cmd = [impl](std::string_view name) {
        return impl->command_registry.has_command(std::string(name)) || impl->native_functions.count(std::string(name)) > 0;
    };
    try {
        CommandASTNode ast = parse_command(trimmed_expr, is_cmd);
        if (ast.kind == CommandKind::kFunctionCall) {
            const auto* call = ast.as_function_call();
            if (impl->native_functions.count(std::string(call->name)) > 0) {
                return evaluate_command_ast_to_value(calculator, impl, ast, exact_mode);
            }
        }
    } catch (...) {
        // 解析失败则回退到统一表达式解析器
    }

    // 隐式求值（模块特定逻辑）
    if (!exact_mode) {
        StoredValue implicit_value;
        if (calculator->try_evaluate_implicit(target_expr, &implicit_value, variables.snapshot())) {
            return implicit_value;
        }
    }

    // 进制转换检测
    std::string converted;
    if (try_base_conversion_expression(target_expr, variables, &impl->functions, {impl->hex_prefix_mode, impl->hex_uppercase_mode}, &converted)) {
        StoredValue res; res.is_string = true; res.string_value = converted;
        return res;
    }

    const HasScriptFunctionCallback has_script_function = [impl](const std::string& name) {
        return has_visible_script_function(impl, name);
    };
    const InvokeScriptFunctionCallback invoke_script_function = [calculator, impl](const std::string& name, const std::vector<Scalar>& arguments) {
        return invoke_script_function_decimal(calculator, impl, name, arguments);
    };

    UnifiedExpressionParser parser(variables, &impl->functions, &impl->scalar_functions, &impl->matrix_functions, &impl->value_functions, has_script_function, invoke_script_function);
    
    StoredValue result = parser.evaluate_stored(target_expr, exact_mode, impl->symbolic_constants_mode);
    
    // 符号常量模式下的额外处理（如果解析器没处理完）
    if (impl->symbolic_constants_mode && !result.has_symbolic_text && !result.is_string && !result.is_matrix) {
        std::string symbolic_text;
        if (try_symbolic_constant_expression(target_expr, variables, &impl->functions, &symbolic_text)) {
            bool symbolic_text_is_plain_decimal = false;
            try {
                const Scalar parsed = mymath::from_string(symbolic_text);
                symbolic_text_is_plain_decimal =
                    mymath::precise128::abs(parsed - result.decimal) < Scalar(1e-12L);
            } catch (...) {
                symbolic_text_is_plain_decimal = false;
            }
            if (!symbolic_text_is_plain_decimal) {
                result.has_symbolic_text = true;
                result.symbolic_text = symbolic_text;
            }
        }
    }
    
    return result;
}

// ============================================================================
// 命令 AST 执行与求值
// ============================================================================

std::string execute_command_ast(Calculator* calculator,
                                Calculator::Impl* impl,
                                const CommandASTNode& ast,
                                bool exact_mode) {
    if (ast.kind == CommandKind::kEmpty) return "";
    
    if (ast.kind == CommandKind::kSequence) {
        const auto* nodes = ast.as_sequence();
        std::ostringstream oss;
        for (std::size_t i = 0; i < nodes->size(); ++i) {
            std::string out = execute_command_ast(calculator, impl, (*nodes)[i], exact_mode);
            if (!out.empty()) {
                oss << out;
                if (i + 1 < nodes->size()) oss << "\n";
            }
        }
        return oss.str();
    }

    if (ast.kind == CommandKind::kFunctionDefinition) {
        const FunctionDefinitionInfo* def = ast.as_function_definition();
        std::string name(def->name);
        if (is_reserved_user_function_name(impl, name)) throw std::runtime_error("function name is reserved: " + name);
        std::vector<std::string> params;
        std::string params_display;
        for (std::size_t i = 0; i < def->parameters.size(); ++i) {
            params.emplace_back(def->parameters[i]);
            params_display += def->parameters[i];
            if (i + 1 < def->parameters.size()) params_display += ", ";
        }
        CustomFunction function;
        function.parameter_names = std::move(params);
        function.expression = std::string(def->body.text);
        impl->functions[name] = std::move(function);
        return name + "(" + params_display + ") = " + std::string(def->body.text);
    }

    if (ast.kind == CommandKind::kMetaCommand || ast.kind == CommandKind::kFunctionCall) {
        std::string_view command_name;
        std::vector<std::string_view> arg_views;
        if (ast.kind == CommandKind::kMetaCommand) {
            const auto* meta = ast.as_meta_command();
            command_name = meta->command;
            arg_views = meta->arguments;
        } else {
            const auto* call = ast.as_function_call();
            command_name = call->name;
            for (const auto& arg : call->arguments) arg_views.push_back(arg.text);
        }

        std::string cmd_name = (ast.kind == CommandKind::kMetaCommand) ? ":" + std::string(command_name) : std::string(command_name);
        if (ast.kind == CommandKind::kFunctionCall && cmd_name == "print") {
            std::ostringstream out;
            for (std::size_t i = 0; i < arg_views.size(); ++i) {
                if (i != 0) out << ' ';
                const StoredValue value =
                    evaluate_script_value_expression(calculator, impl, std::string(arg_views[i]), exact_mode);
                out << format_print_value(value, impl->symbolic_constants_mode);
            }
            return out.str();
        }

        const CoreServices& svc = calculator->get_core_services();
        std::string output;
        if (impl->command_registry.has_command(cmd_name)) {
            if (impl->command_registry.try_process(cmd_name, arg_views, &output, exact_mode, svc)) return output;
        }
        if (ast.kind == CommandKind::kFunctionCall) {
            return format_stored_value(evaluate_command_ast_to_value(calculator, impl, ast, exact_mode), impl->symbolic_constants_mode);
        }
        throw std::runtime_error("unknown command: " + cmd_name);
    }

    if (ast.kind == CommandKind::kAssignment) {
        const auto* assign = ast.as_assignment();
        StoredValue val = evaluate_script_value_expression(calculator, impl, std::string(assign->expression.text), exact_mode);
        assign_visible_variable(impl, std::string(assign->variable), val);
        return std::string(assign->variable) + " = " + format_stored_value(val, impl->symbolic_constants_mode);
    }

    if (ast.kind == CommandKind::kExpression) {
        const auto* expr = ast.as_expression();
        return format_stored_value(evaluate_script_value_expression(calculator, impl, std::string(expr->text), exact_mode), impl->symbolic_constants_mode);
    }

    if (ast.kind == CommandKind::kStringLiteral) return "\"" + *ast.as_string_literal() + "\"";
    return "";
}

StoredValue evaluate_command_ast_to_value(Calculator* calculator,
                                          Calculator::Impl* impl,
                                          const CommandASTNode& ast,
                                          bool exact_mode) {
    if (ast.kind == CommandKind::kExpression) {
        const auto* expr = ast.as_expression();
        return evaluate_script_value_expression(calculator, impl, std::string(expr->text), exact_mode);
    }
    if (ast.kind == CommandKind::kAssignment) {
        const auto* assign = ast.as_assignment();
        StoredValue val = evaluate_script_value_expression(calculator, impl, std::string(assign->expression.text), exact_mode);
        assign_visible_variable(impl, std::string(assign->variable), val);
        return val;
    }
    if (ast.kind == CommandKind::kFunctionCall) {
        const auto* call = ast.as_function_call();
        std::vector<StoredValue> args;
        for (const auto& arg : call->arguments) args.push_back(evaluate_script_value_expression(calculator, impl, std::string(arg.text), exact_mode));
        
        auto it = impl->native_functions.find(std::string(call->name));
        if (it != impl->native_functions.end()) {
            return it->second(args);
        }
        
        return invoke_script_function(calculator, impl, std::string(call->name), args);
    }
    if (ast.kind == CommandKind::kStringLiteral) {
        StoredValue value;
        value.is_string = true;
        value.string_value = *ast.as_string_literal();
        return value;
    }
    std::string out = execute_command_ast(calculator, impl, ast, exact_mode);
    StoredValue res; res.is_string = true; res.string_value = out;
    return res;
}

// ============================================================================
// 脚本执行
// ============================================================================

ScriptSignal execute_script_statement(Calculator* calculator,
                                      Calculator::Impl* impl,
                                      const script::Statement& statement,
                                      bool exact_mode,
                                      std::string* last_output,
                                      bool create_scope) {
    switch (statement.kind) {
        case script::Statement::Kind::kBlock:
            return execute_script_block(calculator, impl, static_cast<const script::BlockStatement&>(statement), exact_mode, last_output, create_scope);
        
        case script::Statement::Kind::kSimple: {
            const auto& simple = static_cast<const script::SimpleStatement&>(statement);
            try {
                if (try_execute_index_assignment(calculator, impl, simple.text, exact_mode, last_output)) {
                    return {};
                }
                if (simple.command_ast.kind != CommandKind::kEmpty) {
                    *last_output = execute_command_ast(calculator, impl, simple.command_ast, exact_mode);
                }
            } catch (const std::exception& e) {
                throw std::runtime_error("Line " + std::to_string(simple.line) + ": " + e.what());
            }
            return {};
        }

        case script::Statement::Kind::kIf: {
            const auto& if_stmt = static_cast<const script::IfStatement&>(statement);
            try {
                if (truthy_value(evaluate_command_ast_to_value(calculator, impl, if_stmt.condition_ast, false))) {
                    return execute_script_statement(calculator, impl, *if_stmt.then_branch, exact_mode, last_output, true);
                } else if (if_stmt.else_branch) {
                    return execute_script_statement(calculator, impl, *if_stmt.else_branch, exact_mode, last_output, true);
                }
            } catch (const std::exception& e) {
                throw std::runtime_error("Line " + std::to_string(if_stmt.line) + ": " + e.what());
            }
            return {};
        }

        case script::Statement::Kind::kWhile: {
            const auto& while_stmt = static_cast<const script::WhileStatement&>(statement);
            try {
                while (truthy_value(evaluate_command_ast_to_value(calculator, impl, while_stmt.condition_ast, false))) {
                    const ScriptSignal signal = execute_script_statement(calculator, impl, *while_stmt.body, exact_mode, last_output, true);
                    if (signal.kind == ScriptSignal::Kind::kReturn) return signal;
                    if (signal.kind == ScriptSignal::Kind::kBreak) break;
                    if (signal.kind == ScriptSignal::Kind::kContinue) continue;
                }
            } catch (const std::exception& e) {
                throw std::runtime_error("Line " + std::to_string(while_stmt.line) + ": " + e.what());
            }
            return {};
        }

        case script::Statement::Kind::kFor: {
            const auto& for_stmt = static_cast<const script::ForStatement&>(statement);
            impl->flat_scopes.push_scope();
            try {
                if (for_stmt.init_ast.kind != CommandKind::kEmpty) (void)evaluate_command_ast_to_value(calculator, impl, for_stmt.init_ast, exact_mode);
                while (for_stmt.cond_ast.kind == CommandKind::kEmpty || truthy_value(evaluate_command_ast_to_value(calculator, impl, for_stmt.cond_ast, false))) {
                    const ScriptSignal signal = execute_script_statement(calculator, impl, *for_stmt.body, exact_mode, last_output, true);
                    if (signal.kind == ScriptSignal::Kind::kReturn) { impl->flat_scopes.pop_scope(); return signal; }
                    if (signal.kind == ScriptSignal::Kind::kBreak) break;
                    if (for_stmt.step_ast.kind != CommandKind::kEmpty) (void)evaluate_command_ast_to_value(calculator, impl, for_stmt.step_ast, exact_mode);
                    if (signal.kind == ScriptSignal::Kind::kContinue) continue;
                }
                impl->flat_scopes.pop_scope();
                return {};
            } catch (const std::exception& e) {
                impl->flat_scopes.pop_scope();
                throw std::runtime_error("Line " + std::to_string(for_stmt.line) + ": " + e.what());
            }
        }

        case script::Statement::Kind::kForRange: {
            const auto& for_range = static_cast<const script::ForRangeStatement&>(statement);
            impl->flat_scopes.push_scope();
            try {
                const Scalar start = evaluate_scalar(calculator, impl, for_range.start_ast, "range start");
                const Scalar stop = evaluate_scalar(calculator, impl, for_range.stop_ast, "range stop");
                const Scalar step = evaluate_scalar(calculator, impl, for_range.step_ast, "range step");
                if (step == 0.0L) throw std::runtime_error("range step cannot be zero");
                bool ascending = step > 0.0L;
                for (Scalar current = start; (ascending ? current < stop : current > stop); current += step) {
                    StoredValue loop_val; loop_val.decimal = current; loop_val.exact = false;
                    impl->flat_scopes.set(for_range.variable, loop_val);
                    const ScriptSignal signal = execute_script_statement(calculator, impl, *for_range.body, exact_mode, last_output, true);
                    if (signal.kind == ScriptSignal::Kind::kReturn) { impl->flat_scopes.pop_scope(); return signal; }
                    if (signal.kind == ScriptSignal::Kind::kBreak) break;
                    if (signal.kind == ScriptSignal::Kind::kContinue) continue;
                }
                impl->flat_scopes.pop_scope();
                return {};
            } catch (const std::exception& e) {
                impl->flat_scopes.pop_scope();
                throw std::runtime_error("Line " + std::to_string(for_range.line) + ": " + e.what());
            }
        }

        case script::Statement::Kind::kForIn: {
            const auto& for_in = static_cast<const script::ForInStatement&>(statement);
            impl->flat_scopes.push_scope();
            try {
                StoredValue iterable = evaluate_command_ast_to_value(calculator, impl, for_in.iterable_ast, exact_mode);

                // 支持列表迭代
                if (iterable.is_list && iterable.list_value) {
                    for (const StoredValue& item : *iterable.list_value) {
                        impl->flat_scopes.set(for_in.variable, item);
                        const ScriptSignal signal = execute_script_statement(calculator, impl, *for_in.body, exact_mode, last_output, true);
                        if (signal.kind == ScriptSignal::Kind::kReturn) { impl->flat_scopes.pop_scope(); return signal; }
                        if (signal.kind == ScriptSignal::Kind::kBreak) break;
                        if (signal.kind == ScriptSignal::Kind::kContinue) continue;
                    }
                }
                // 支持矩阵行迭代
                else if (iterable.is_matrix) {
                    for (std::size_t row = 0; row < iterable.matrix.rows; ++row) {
                        StoredValue row_value;
                        row_value.is_matrix = true;
                        row_value.matrix = matrix::Matrix(1, iterable.matrix.cols);
                        for (std::size_t col = 0; col < iterable.matrix.cols; ++col) {
                            row_value.matrix.at(0, col) = iterable.matrix.at(row, col);
                        }
                        impl->flat_scopes.set(for_in.variable, row_value);
                        const ScriptSignal signal = execute_script_statement(calculator, impl, *for_in.body, exact_mode, last_output, true);
                        if (signal.kind == ScriptSignal::Kind::kReturn) { impl->flat_scopes.pop_scope(); return signal; }
                        if (signal.kind == ScriptSignal::Kind::kBreak) break;
                        if (signal.kind == ScriptSignal::Kind::kContinue) continue;
                    }
                }
                // 支持字符串字符迭代
                else if (iterable.is_string) {
                    for (char ch : iterable.string_value) {
                        StoredValue char_value;
                        char_value.is_string = true;
                        char_value.string_value = std::string(1, ch);
                        impl->flat_scopes.set(for_in.variable, char_value);
                        const ScriptSignal signal = execute_script_statement(calculator, impl, *for_in.body, exact_mode, last_output, true);
                        if (signal.kind == ScriptSignal::Kind::kReturn) { impl->flat_scopes.pop_scope(); return signal; }
                        if (signal.kind == ScriptSignal::Kind::kBreak) break;
                        if (signal.kind == ScriptSignal::Kind::kContinue) continue;
                    }
                }
                else {
                    throw std::runtime_error("for-in loop requires a list, matrix, or string iterable");
                }
                impl->flat_scopes.pop_scope();
                return {};
            } catch (const std::exception& e) {
                impl->flat_scopes.pop_scope();
                throw std::runtime_error("Line " + std::to_string(for_in.line) + ": " + e.what());
            }
        }

        case script::Statement::Kind::kReturn: {
            const auto& ret = static_cast<const script::ReturnStatement&>(statement);
            if (!ret.has_expression) return ScriptSignal::make_return({});
            try { return ScriptSignal::make_return(evaluate_command_ast_to_value(calculator, impl, ret.expr_ast, exact_mode)); }
            catch (const std::exception& e) { throw std::runtime_error("Line " + std::to_string(ret.line) + ": " + e.what()); }
        }

        case script::Statement::Kind::kBreak: return ScriptSignal::make_break();
        case script::Statement::Kind::kContinue: return ScriptSignal::make_continue();
        case script::Statement::Kind::kImport: {
            const auto& import_stmt = static_cast<const script::ImportStatement&>(statement);
            try {
                StoredValue path_value =
                    evaluate_command_ast_to_value(calculator, impl, import_stmt.path_ast, exact_mode);
                if (!path_value.is_string) {
                    throw std::runtime_error("import path must be a string");
                }
                *last_output = calculator->execute_script_file(path_value.string_value, exact_mode, true);
                return {};
            } catch (const std::exception& e) {
                throw std::runtime_error("Line " + std::to_string(import_stmt.line) + ": " + e.what());
            }
        }
        case script::Statement::Kind::kFunction: {
            const auto& fs = static_cast<const script::FunctionStatement&>(statement);
            ScriptFunction function; function.parameter_names = fs.parameters; function.body = fs.body;
            impl->script_functions[fs.name] = function;
            *last_output = fs.name + "(...)";
            return {};
        }

        case script::Statement::Kind::kMatch: {
            const auto& match_stmt = static_cast<const script::MatchStatement&>(statement);
            try {
                // 求值匹配主体
                StoredValue subject = evaluate_command_ast_to_value(calculator, impl, match_stmt.subject_ast, exact_mode);

                // 遍历 case 分支
                for (const auto& clause : match_stmt.cases) {
                    bool matches = false;

                    if (clause.is_default) {
                        // 默认分支总是匹配
                        matches = true;
                    } else {
                        // 求值模式并比较
                        StoredValue pattern = evaluate_command_ast_to_value(calculator, impl, clause.pattern_ast, exact_mode);

                        // 比较值
                        if (subject.is_matrix || pattern.is_matrix) {
                            if (subject.is_matrix && pattern.is_matrix) {
                                if (subject.matrix.rows == pattern.matrix.rows && subject.matrix.cols == pattern.matrix.cols) {
                                    matches = true;
                                    for (std::size_t i = 0; i < subject.matrix.data.size(); ++i) {
                                        if (!mymath::is_near_zero(subject.matrix.data[i] - pattern.matrix.data[i], 1e-10)) {
                                            matches = false;
                                            break;
                                        }
                                    }
                                } else {
                                    matches = false;
                                }
                            } else {
                                matches = false;
                            }
                        } else if (subject.is_complex || pattern.is_complex) {
                            // 复数比较
                            if (subject.is_complex && pattern.is_complex) {
                                matches = mymath::is_near_zero(subject.complex.real - pattern.complex.real, 1e-10) &&
                                          mymath::is_near_zero(subject.complex.imag - pattern.complex.imag, 1e-10);
                            } else {
                                matches = false;
                            }
                        } else if (subject.is_string || pattern.is_string) {
                            // 字符串比较
                            if (subject.is_string && pattern.is_string) {
                                matches = subject.string_value == pattern.string_value;
                            } else {
                                matches = false;
                            }
                        } else {
                            // 标量比较
                            Scalar subj_val = subject.exact ? rational_to_double(subject.rational) : subject.decimal;
                            Scalar pat_val = pattern.exact ? rational_to_double(pattern.rational) : pattern.decimal;
                            matches = mymath::is_near_zero(subj_val - pat_val, 1e-10);
                        }

                    }

                    // 检查守卫条件。默认分支也允许 guard，例如 case _ if x > 0:
                    if (matches && clause.is_guarded) {
                        StoredValue guard_val = evaluate_command_ast_to_value(calculator, impl, clause.guard_ast, exact_mode);
                        matches = truthy_value(guard_val);
                    }

                    if (matches) {
                        // 执行匹配的 case 体
                        return execute_script_statement(calculator, impl, *clause.body, exact_mode, last_output, true);
                    }
                }
                // 没有匹配的分支，什么都不做
            } catch (const std::exception& e) {
                throw std::runtime_error("Line " + std::to_string(match_stmt.line) + ": " + e.what());
            }
            return {};
        }
    }
    return {};
}

ScriptSignal execute_script_block(Calculator* calculator,
                                  Calculator::Impl* impl,
                                  const script::BlockStatement& block,
                                  bool exact_mode,
                                  std::string* last_output,
                                  bool create_scope) {
    if (create_scope) impl->flat_scopes.push_scope();
    try {
        for (const auto& stmt : block.statements) {
            const ScriptSignal signal = execute_script_statement(calculator, impl, *stmt, exact_mode, last_output, true);
            if (signal.kind != ScriptSignal::Kind::kNone) {
                if (create_scope) impl->flat_scopes.pop_scope();
                return signal;
            }
        }
        if (create_scope) impl->flat_scopes.pop_scope();
        return {};
    } catch (...) {
        if (create_scope) impl->flat_scopes.pop_scope();
        throw;
    }
}

std::string execute_simple_script_line(Calculator* calculator,
                                       Calculator::Impl* impl,
                                       const std::string& text,
                                       bool exact_mode) {
    const std::string trimmed = trim_copy(text);
    if (trimmed.empty()) return "";
    auto is_command = [impl](std::string_view name) {
        return impl->command_registry.has_command(std::string(name)) || impl->command_registry.has_command(":" + std::string(name));
    };
    CommandASTNode ast = parse_command(trimmed, is_command);
    return execute_command_ast(calculator, impl, ast, exact_mode);
}

Scalar invoke_script_function_decimal(Calculator* calculator,
                                      Calculator::Impl* impl,
                                      const std::string& name,
                                      const std::vector<Scalar>& arguments) {
    std::vector<StoredValue> stored_args;
    for (Scalar arg : arguments) { StoredValue sv; sv.decimal = arg; sv.exact = false; stored_args.push_back(sv); }
    StoredValue result = invoke_script_function(calculator, impl, name, stored_args);
    return result.exact ? rational_to_double(result.rational) : result.decimal;
}

StoredValue invoke_script_function(Calculator* calculator,
                                   Calculator::Impl* impl,
                                   const std::string& name,
                                   const std::vector<StoredValue>& arguments) {
    auto nit = impl->native_functions.find(name);
    if (nit != impl->native_functions.end()) {
        return nit->second(arguments);
    }
    
    auto it = impl->script_functions.find(name);
    if (it == impl->script_functions.end()) throw std::runtime_error("unknown function: " + name);
    const ScriptFunction& function = it->second;
    impl->flat_scopes.push_scope();
    for (std::size_t i = 0; i < arguments.size(); ++i) impl->flat_scopes.set(function.parameter_names[i], arguments[i]);
    std::string ignored;
    const ScriptSignal signal = execute_script_block(calculator, impl, *function.body, false, &ignored, false);
    impl->flat_scopes.pop_scope();
    if (signal.kind != ScriptSignal::Kind::kReturn || !signal.has_value) throw std::runtime_error("script function " + name + " must return a value");
    return signal.value;
}

// 克隆与渲染逻辑 (省略或保持最简，脚本中目前主要使用 shared_ptr)
script::StatementPtr clone_statement(const script::Statement& statement) {
    (void)statement;
    return nullptr;
}

namespace {

std::string command_ast_to_source(const CommandASTNode& ast) {
    switch (ast.kind) {
        case CommandKind::kEmpty:
            return "";
        case CommandKind::kMetaCommand: {
            const auto& meta = std::get<MetaCommandInfo>(ast.data);
            std::string text = ":" + std::string(meta.command);
            for (std::string_view arg : meta.arguments) {
                text += " ";
                text += std::string(arg);
            }
            return text;
        }
        case CommandKind::kFunctionDefinition: {
            const auto& def = std::get<FunctionDefinitionInfo>(ast.data);
            std::string text = std::string(def.name) + "(";
            for (std::size_t i = 0; i < def.parameters.size(); ++i) {
                if (i != 0) text += ", ";
                text += std::string(def.parameters[i]);
            }
            text += ") = ";
            text += std::string(def.body.text);
            return text;
        }
        case CommandKind::kFunctionCall: {
            const auto& call = std::get<FunctionCallInfo>(ast.data);
            std::string text = std::string(call.name) + "(";
            for (std::size_t i = 0; i < call.arguments.size(); ++i) {
                if (i != 0) text += ", ";
                text += std::string(call.arguments[i].text);
            }
            text += ")";
            return text;
        }
        case CommandKind::kAssignment: {
            const auto& assignment = std::get<AssignmentInfo>(ast.data);
            return std::string(assignment.variable) + " = " +
                   std::string(assignment.expression.text);
        }
        case CommandKind::kIndexAssignment: {
            const auto& idx_assign = std::get<IndexAssignmentInfo>(ast.data);
            std::string text = std::string(idx_assign.variable) + "[";
            for (std::size_t i = 0; i < idx_assign.indices.size(); ++i) {
                if (i != 0) text += ", ";
                text += std::string(idx_assign.indices[i].text);
            }
            text += "] = " + std::string(idx_assign.value.text);
            return text;
        }
        case CommandKind::kExpression:
            return std::string(std::get<ExpressionInfo>(ast.data).text);
        case CommandKind::kStringLiteral:
            return "\"" + std::get<std::string>(ast.data) + "\"";
        case CommandKind::kSequence: {
            const auto& nodes = std::get<std::vector<CommandASTNode>>(ast.data);
            std::string text;
            for (std::size_t i = 0; i < nodes.size(); ++i) {
                if (i != 0) text += "; ";
                text += command_ast_to_source(nodes[i]);
            }
            return text;
        }
    }
    return "";
}

}  // namespace

std::string render_script_statement(const script::Statement& statement, int indent) {
    const std::string pad = indent_text(indent);
    switch (statement.kind) {
        case script::Statement::Kind::kBlock:
            return render_script_block(static_cast<const script::BlockStatement&>(statement), indent);
        case script::Statement::Kind::kSimple: {
            const auto& simple = static_cast<const script::SimpleStatement&>(statement);
            const std::string text = simple.text.empty()
                                         ? command_ast_to_source(simple.command_ast)
                                         : simple.text;
            return pad + (text.empty() ? "pass" : text);
        }
        case script::Statement::Kind::kIf: {
            const auto& if_stmt = static_cast<const script::IfStatement&>(statement);
            std::string text = pad + "if " + command_ast_to_source(if_stmt.condition_ast) +
                               " " + render_script_statement(*if_stmt.then_branch, indent);
            if (if_stmt.else_branch) {
                text += "\n" + pad + "else " +
                        render_script_statement(*if_stmt.else_branch, indent);
            }
            return text;
        }
        case script::Statement::Kind::kWhile: {
            const auto& while_stmt = static_cast<const script::WhileStatement&>(statement);
            return pad + "while " + command_ast_to_source(while_stmt.condition_ast) +
                   " " + render_script_statement(*while_stmt.body, indent);
        }
        case script::Statement::Kind::kFor: {
            const auto& for_stmt = static_cast<const script::ForStatement&>(statement);
            return pad + "for (" + command_ast_to_source(for_stmt.init_ast) + "; " +
                   command_ast_to_source(for_stmt.cond_ast) + "; " +
                   command_ast_to_source(for_stmt.step_ast) + ") " +
                   render_script_statement(*for_stmt.body, indent);
        }
        case script::Statement::Kind::kForRange: {
            const auto& for_stmt = static_cast<const script::ForRangeStatement&>(statement);
            return pad + "for " + for_stmt.variable + " in range(" +
                   command_ast_to_source(for_stmt.start_ast) + ", " +
                   command_ast_to_source(for_stmt.stop_ast) + ", " +
                   command_ast_to_source(for_stmt.step_ast) + ") " +
                   render_script_statement(*for_stmt.body, indent);
        }
        case script::Statement::Kind::kForIn: {
            const auto& for_in = static_cast<const script::ForInStatement&>(statement);
            return pad + "for " + for_in.variable + " in " +
                   command_ast_to_source(for_in.iterable_ast) + ": " +
                   render_script_statement(*for_in.body, indent);
        }
        case script::Statement::Kind::kFunction: {
            const auto& function = static_cast<const script::FunctionStatement&>(statement);
            std::string text = pad + "fn " + function.name + "(";
            for (std::size_t i = 0; i < function.parameters.size(); ++i) {
                if (i != 0) text += ", ";
                text += function.parameters[i];
            }
            text += ") " + render_script_block(*function.body, indent);
            return text;
        }
        case script::Statement::Kind::kReturn: {
            const auto& ret = static_cast<const script::ReturnStatement&>(statement);
            return pad + "return" +
                   (ret.has_expression ? " " + command_ast_to_source(ret.expr_ast) : "");
        }
        case script::Statement::Kind::kBreak:
            return pad + "break";
        case script::Statement::Kind::kContinue:
            return pad + "continue";
        case script::Statement::Kind::kImport: {
            const auto& import_stmt = static_cast<const script::ImportStatement&>(statement);
            return pad + "import " +
                   (import_stmt.path_text.empty()
                        ? command_ast_to_source(import_stmt.path_ast)
                        : import_stmt.path_text);
        }
        case script::Statement::Kind::kMatch: {
            const auto& match_stmt = static_cast<const script::MatchStatement&>(statement);
            std::ostringstream out;
            out << pad << "match " << command_ast_to_source(match_stmt.subject_ast) << ":\n";
            for (const auto& clause : match_stmt.cases) {
                out << pad << "  case ";
                if (clause.is_default) {
                    out << "_";
                } else {
                    out << command_ast_to_source(clause.pattern_ast);
                    if (clause.is_guarded) {
                        out << " if " << command_ast_to_source(clause.guard_ast);
                    }
                }
                out << ": " << render_script_statement(*clause.body, indent + 1) << "\n";
            }
            return out.str();
        }
    }
    return pad + "pass";
}

std::string render_script_block(const script::BlockStatement& block, int indent) {
    std::ostringstream out;
    out << "{";
    if (!block.statements.empty()) out << "\n";
    for (const auto& stmt : block.statements) {
        out << render_script_statement(*stmt, indent + 1) << "\n";
    }
    out << indent_text(indent) << "}";
    return out.str();
}

// ============================================================================
// 统一执行接口实现
// ============================================================================

#include "core/api/executable_node.h"

namespace script {

StoredValue Statement::execute(Calculator* calculator, bool exact_mode) const {
    std::string dummy_output;
    ScriptSignal signal = execute_script_statement(calculator, calculator->get_impl_internal(), *this, exact_mode, &dummy_output, false);
    return signal.has_value ? signal.value : StoredValue();
}

} // namespace script

namespace execution {

/**
 * @class CommandExecutable
 * @brief 包装 CommandASTNode 的执行器
 */
class CommandExecutable : public ExecutableNode {
public:
    explicit CommandExecutable(CommandASTNode node) : node_(std::move(node)) {}

    StoredValue execute(Calculator* calculator, bool exact_mode) const override {
        return evaluate_command_ast_to_value(calculator, calculator->get_impl_internal(), node_, exact_mode);
    }

private:
    CommandASTNode node_;
};

/**
 * @brief 创建命令执行器
 */
std::unique_ptr<ExecutableNode> create_command_executable(CommandASTNode node) {
    return std::make_unique<CommandExecutable>(std::move(node));
}

} // namespace execution
