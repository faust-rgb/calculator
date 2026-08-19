// ============================================================================
// 内联函数展开器实现文件
// ============================================================================
//
// 本文件实现了内联函数命令（如 diff, integral, taylor 等）的展开功能。
// 这些命令可以在表达式中内联执行，将结果替换回表达式。
//
// 主要功能：
// - 识别表达式中的内联函数命令
// - 递归展开嵌套的内联函数
// - 将展开结果替换回原始表达式
//
// 支持的内联命令：
// - 微积分：diff, integral, double_integral, triple_integral
// - 数值计算：limit, ode, ode_system
// - 级数展开：taylor
// - 多项式运算：poly_add, poly_sub, poly_mul
// ============================================================================

#include "inline_expander.h"
#include "core/services/core_manager_interfaces.h"

#include <cctype>
#include <stdexcept>

// ============================================================================
// 匿名命名空间 - 内部辅助函数
// ============================================================================

namespace {

/**
 * @brief 检查名称是否为内联函数命令
 * @param ctx 执行上下文
 * @param name 要检查的名称
 * @return 如果是内联函数命令返回 true
 */
bool is_inline_function_command(IExecutionContext* ctx, std::string_view name) {
    if (ctx && ctx->commands().is_inlineable(std::string(name))) {
        return true;
    }
    // 兼容内置数学命令的回退检查
    return name == "diff" ||
           name == "double_integral" ||
           name == "double_integral_cyl" ||
           name == "double_integral_polar" ||
           name == "integral" ||
           name == "integrate" ||
           name == "limit" ||
           name == "ode" ||
           name == "ode_system" ||
           name == "ode_table" ||
           name == "ode_system_table" ||
           name == "taylor" ||
           name == "triple_integral" ||
           name == "triple_integral_cyl" ||
           name == "triple_integral_sph" ||
           name == "poly_add" ||
           name == "poly_sub" ||
           name == "poly_mul";
}

/**
 * @brief 查找匹配的右括号
 * @param text 输入文本
 * @param open_pos 左括号位置
 * @return 右括号位置，如果未找到返回 npos
 *
 * 正确处理嵌套括号、字符串字面量和转义字符。
 */
std::size_t find_matching_paren(std::string_view text, std::size_t open_pos) {
    if (open_pos >= text.size() || text[open_pos] != '(') {
        return std::string::npos;
    }

    int depth = 0;
    bool in_string = false;
    bool escaping = false;
    for (std::size_t i = open_pos; i < text.size(); ++i) {
        const char ch = text[i];
        // 处理字符串字面量内部
        if (in_string) {
            if (escaping) {
                escaping = false;
            } else if (ch == '\\') {
                escaping = true;
            } else if (ch == '"') {
                in_string = false;
            }
            continue;
        }
        if (ch == '"') {
            in_string = true;
            continue;
        }
        // 计算括号深度
        if (ch == '(') {
            ++depth;
        } else if (ch == ')') {
            --depth;
            if (depth == 0) {
                return i;
            }
        }
    }

    return std::string::npos;
}

} // namespace

// ============================================================================
// 主展开函数
// ============================================================================

/**
 * @brief 展开表达式中的内联函数命令
 * @param ctx 计算器实例
 * @param expression 输入表达式
 * @return 展开后的表达式
 *
 * 递归地扫描表达式，识别并执行内联函数命令，将结果替换回表达式。
 * 设置最大递归深度限制，防止无限递归。
 */
std::string expand_inline_function_commands(
                                            IExecutionContext* ctx,
                                            std::string_view expression) {
    // 递归深度限制，防止无限递归
    static thread_local int depth = 0;
    static constexpr int kMaxDepth = 32;
    if (++depth > kMaxDepth) {
        --depth;
        throw std::runtime_error("inline function expansion too deep");
    }

    // RAII 守卫，确保退出时减少深度计数
    struct DepthGuard {
        int* d;
        ~DepthGuard() { if (d) (*d)--; }
    } guard{&depth};

    std::string expanded;
    expanded.reserve(expression.size());

    // 扫描表达式
    for (std::size_t i = 0; i < expression.size();) {
        const char ch = expression[i];
        // 非字母字符直接追加
        if (!std::isalpha(static_cast<unsigned char>(ch))) {
            expanded.push_back(ch);
            ++i;
            continue;
        }

        // 提取标识符名称
        std::size_t name_end = i;
        while (name_end < expression.size()) {
            const char name_ch = expression[name_end];
            if (std::isalnum(static_cast<unsigned char>(name_ch)) || name_ch == '_') {
                ++name_end;
            } else {
                break;
            }
        }

        const std::string_view name = expression.substr(i, name_end - i);
        // 如果不是内联函数命令，直接追加
        if (!is_inline_function_command(ctx, name)) {
            expanded.append(expression.substr(i, name_end - i));
            i = name_end;
            continue;
        }

        // 跳过空白字符
        std::size_t open_pos = name_end;
        while (open_pos < expression.size() &&
               std::isspace(static_cast<unsigned char>(expression[open_pos]))) {
            ++open_pos;
        }
        // 检查是否有括号
        if (open_pos >= expression.size() || expression[open_pos] != '(') {
            expanded.append(expression.substr(i, name_end - i));
            i = name_end;
            continue;
        }

        // 查找匹配的右括号
        const std::size_t close_pos = find_matching_paren(expression, open_pos);
        if (close_pos == std::string::npos) {
            expanded.append(expression.substr(i, name_end - i));
            i = name_end;
            continue;
        }

        // 递归展开括号内的内容
        const std::string inner =
            expand_inline_function_commands(ctx,
                                            expression.substr(open_pos + 1,
                                                              close_pos - open_pos - 1));
        // 重建命令并执行
        const std::string rebuilt = std::string(name) + "(" + inner + ")";
        std::string command_output;
        if (ctx->try_process_function_command(rebuilt, &command_output)) {
            expanded += "(" + command_output + ")";
        } else {
            expanded += rebuilt;
        }
        i = close_pos + 1;
    }

    return expanded;
}
