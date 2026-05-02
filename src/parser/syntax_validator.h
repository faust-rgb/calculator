// ============================================================================
// 语法验证器 - 增强的语法错误检测
// ============================================================================
//
// 设计目标：
// 1. 在解析前检测明显的语法错误
// 2. 提供精确的错误位置和消息
// 3. 支持多种错误类型检测
//
// 检测规则：
// - 括号平衡
// - 运算符序列合法性
// - 操作数上下文
// - 函数调用语法
// - 字符串终止
// - 无效字符序列
// ============================================================================

#ifndef CORE_SYNTAX_VALIDATOR_H
#define CORE_SYNTAX_VALIDATOR_H

#include <string>
#include <string_view>
#include <vector>

/**
 * @enum Severity
 * @brief 错误严重程度
 */
enum class Severity {
    kError,     ///< 错误（阻止解析）
    kWarning,   ///< 警告（可能有问题）
    kInfo,      ///< 信息（建议）
};

/**
 * @struct SyntaxError
 * @brief 语法错误信息
 */
struct SyntaxErrorInfo {
    std::string message;        ///< 错误消息
    std::size_t position;       ///< 错误位置
    Severity severity;          ///< 严重程度
    std::string context;        ///< 上下文片段

    SyntaxErrorInfo()
        : position(0), severity(Severity::kError) {}

    SyntaxErrorInfo(const std::string& msg, std::size_t pos, Severity sev)
        : message(msg), position(pos), severity(sev) {}
};

/**
 * @class SyntaxValidator
 * @brief 语法验证器，检测表达式中的语法错误
 *
 * 在解析前进行语法预检，提供精确的错误报告。
 */
class SyntaxValidator {
public:
    /**
     * @brief 验证表达式语法
     * @param expression 表达式字符串
     * @return 错误列表
     */
    std::vector<SyntaxErrorInfo> validate(std::string_view expression);

    /**
     * @brief 检查表达式是否有错误
     * @param expression 表达式字符串
     * @return 如果有错误返回 true
     */
    bool has_errors(std::string_view expression);

    /**
     * @brief 获取第一个错误消息
     * @param expression 表达式字符串
     * @return 错误消息，如果没有错误返回空字符串
     */
    std::string get_first_error(std::string_view expression);

    /**
     * @brief 格式化错误报告
     * @param expression 表达式字符串
     * @param errors 错误列表
     * @return 格式化的错误报告
     */
    static std::string format_errors(std::string_view expression,
                                      const std::vector<SyntaxErrorInfo>& errors);

private:
    // ========================================================================
    // 检测规则
    // ========================================================================

    /**
     * @brief 检查括号平衡
     */
    bool check_bracket_balance(std::string_view expr, std::vector<SyntaxErrorInfo>& errors);

    /**
     * @brief 检查运算符序列
     */
    bool check_operator_sequences(std::string_view expr, std::vector<SyntaxErrorInfo>& errors);

    /**
     * @brief 检查操作数上下文
     */
    bool check_operand_context(std::string_view expr, std::vector<SyntaxErrorInfo>& errors);

    /**
     * @brief 检查函数调用语法
     */
    bool check_function_syntax(std::string_view expr, std::vector<SyntaxErrorInfo>& errors);

    /**
     * @brief 检查字符串终止
     */
    bool check_string_termination(std::string_view expr, std::vector<SyntaxErrorInfo>& errors);

    /**
     * @brief 检查无效字符序列
     */
    bool check_invalid_chars(std::string_view expr, std::vector<SyntaxErrorInfo>& errors);

    /**
     * @brief 检查表达式结尾
     */
    bool check_expression_end(std::string_view expr, std::vector<SyntaxErrorInfo>& errors);

    // ========================================================================
    // 辅助方法
    // ========================================================================

    /**
     * @brief 添加错误
     */
    void add_error(std::vector<SyntaxErrorInfo>& errors,
                   const std::string& message,
                   std::size_t position,
                   Severity severity = Severity::kError);

    /**
     * @brief 检查字符是否是运算符
     */
    bool is_operator_char(char ch) const;

    /**
     * @brief 检查字符是否是括号
     */
    bool is_bracket_char(char ch) const;

    /**
     * @brief 获取括号的配对
     */
    char get_matching_bracket(char ch) const;

    /**
     * @brief 检查是否是有效的双运算符序列
     */
    bool is_valid_double_operator(char first, char second) const;
};

#endif // CORE_SYNTAX_VALIDATOR_H