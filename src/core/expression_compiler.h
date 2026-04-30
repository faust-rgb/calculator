// ============================================================================
// 表达式预编译系统
// ============================================================================
//
// 提供表达式类型分析、词法分析和预编译功能，解决以下问题：
// 1. 避免循环中重复解析表达式
// 2. 在解析阶段完成函数展开
// 3. 统一表达式词法分析
// ============================================================================

#ifndef EXPRESSION_COMPILER_H
#define EXPRESSION_COMPILER_H

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

// 前向声明
class Calculator;
struct StoredValue;
struct Rational;

// ============================================================================
// 表达式类型标记
// ============================================================================

/**
 * @enum ExpressionHint
 * @brief 表达式类型提示，用于快速分发到正确的求值路径
 */
enum class ExpressionHint {
    kStringLiteral,     ///< 字符串字面量 "..."
    kIdentifier,        ///< 单个标识符（变量名）
    kComplexCandidate,  ///< 可能是复数表达式（包含独立的 i）
    kMatrixCandidate,   ///< 可能是矩阵表达式（包含 [...]）
    kRatCall,           ///< rat(...) 函数调用
    kAssignment,        ///< 赋值表达式 x = ...
    kScalar,            ///< 纯标量表达式
    kUnknown,           ///< 需要动态判断
};

/**
 * @enum ExpressionFeature
 * @brief 表达式特征位掩码
 */
enum class ExpressionFeature : uint32_t {
    kNone           = 0,
    kHasI           = 1 << 0,   ///< 包含独立的 i
    kHasBracket     = 1 << 1,   ///< 包含 [...]
    kHasString      = 1 << 2,   ///< 包含字符串字面量
    kHasRatCall     = 1 << 3,   ///< 包含 rat() 调用
    kHasAssignment  = 1 << 4,   ///< 包含赋值
    kHasFunction    = 1 << 5,   ///< 包含函数调用
    kHasOperator    = 1 << 6,   ///< 包含运算符
    kHasNumber      = 1 << 7,   ///< 包含数字
    kHasIdentifier  = 1 << 8,   ///< 包含标识符
    kHasComparison  = 1 << 9,   ///< 包含比较运算符
    kHasMatrixFunc  = 1 << 10,  ///< 包含矩阵函数 (mat, vec, zeros, etc.)
};

inline ExpressionFeature operator|(ExpressionFeature a, ExpressionFeature b) {
    return static_cast<ExpressionFeature>(static_cast<uint32_t>(a) | static_cast<uint32_t>(b));
}

inline ExpressionFeature operator&(ExpressionFeature a, ExpressionFeature b) {
    return static_cast<ExpressionFeature>(static_cast<uint32_t>(a) & static_cast<uint32_t>(b));
}

inline bool has_feature(ExpressionFeature features, ExpressionFeature feature) {
    return (features & feature) == feature;
}

// ============================================================================
// 表达式 Token
// ============================================================================

/**
 * @enum ExpressionTokenKind
 * @brief 表达式 Token 类型
 */
enum class ExpressionTokenKind {
    kNumber,        ///< 数字字面量
    kIdentifier,    ///< 标识符（变量名或函数名）
    kString,        ///< 字符串字面量
    kOperator,      ///< 运算符
    kLParen,        ///< 左圆括号
    kRParen,        ///< 右圆括号
    kLBracket,      ///< 左方括号
    kRBracket,      ///< 右方括号
    kComma,         ///< 逗号
    kSemicolon,     ///< 分号
    kColon,         ///< 冒号
    kQuestion,      ///< 问号（三元运算符）
    kEqual,         ///< 等号（赋值或比较）
    kEOF,           ///< 结束标记
};

/**
 * @struct ExpressionToken
 * @brief 表达式 Token 结构
 */
struct ExpressionToken {
    ExpressionTokenKind kind;
    std::string text;
    double number_value = 0.0;  ///< 当 kind == kNumber 时的数值
    std::size_t position = 0;   ///< 在源字符串中的位置
};

// ============================================================================
// 字节码指令
// ============================================================================

/**
 * @enum OpCode
 * @brief 表达式字节码操作码
 */
enum class OpCode : uint8_t {
    // 数据操作
    kPushNumber,    ///< 推入数值常量
    kPushVar,       ///< 推入变量值
    kPushString,    ///< 推入字符串常量
    kPushI,         ///< 推入虚数单位 i

    // 算术运算
    kAdd,           ///< 加法
    kSub,           ///< 减法
    kMul,           ///< 乘法
    kDiv,           ///< 除法
    kPow,           ///< 幂运算
    kNeg,           ///< 负号
    kPos,           ///< 正号

    // 比较运算
    kEq,            ///< ==
    kNe,            ///< !=
    kLt,            ///< <
    kLe,            ///< <=
    kGt,            ///< >
    kGe,            ///< >=

    // 函数调用
    kCallFunc,      ///< 调用内置函数
    kCallScript,    ///< 调用脚本函数
    kCallCustom,    ///< 调用自定义函数

    // 矩阵操作
    kBuildMatrix,   ///< 构建矩阵
    kIndexMatrix,   ///< 矩阵索引

    // 特殊操作
    kRat,           ///< rat() 特殊处理
    kConditional,   ///< 三元条件运算

    // 结束
    kHalt,          ///< 停止执行
};

/**
 * @struct Instruction
 * @brief 字节码指令
 */
struct Instruction {
    OpCode op;
    uint32_t operand;  ///< 操作数索引（指向常量池或变量池）
};

// ============================================================================
// 预编译表达式
// ============================================================================

/**
 * @class CompiledExpression
 * @brief 预编译的表达式，包含字节码和常量池
 */
class CompiledExpression {
public:
    CompiledExpression() = default;

    // 分析并编译表达式
    bool compile(const std::string& expression);

    // 获取类型提示
    ExpressionHint hint() const { return hint_; }

    // 获取原始表达式
    const std::string& source() const { return source_; }

    // 获取展开后的表达式
    const std::string& expanded() const { return expanded_; }

    // 设置展开后的表达式
    void set_expanded(const std::string& expanded) { expanded_ = expanded; }

    // 获取特征
    ExpressionFeature features() const { return features_; }

    // 获取 Token 序列
    const std::vector<ExpressionToken>& tokens() const { return tokens_; }

    // 获取字节码
    const std::vector<Instruction>& code() const { return code_; }

    // 获取数值常量池
    const std::vector<double>& number_pool() const { return number_pool_; }

    // 获取字符串/标识符池
    const std::vector<std::string>& string_pool() const { return string_pool_; }

    // 是否有效
    bool is_valid() const { return valid_; }

    // 是否已编译
    bool is_compiled() const { return compiled_; }

private:
    // 词法分析
    bool tokenize();

    // 类型分析
    void analyze_features();

    // 生成字节码
    bool generate_code();

    std::string source_;
    std::string expanded_;
    std::vector<ExpressionToken> tokens_;
    std::vector<Instruction> code_;
    std::vector<double> number_pool_;
    std::vector<std::string> string_pool_;
    ExpressionHint hint_ = ExpressionHint::kUnknown;
    ExpressionFeature features_ = ExpressionFeature::kNone;
    bool valid_ = false;
    bool compiled_ = false;
};

// ============================================================================
// 表达式分析器
// ============================================================================

/**
 * @brief 分析表达式类型（不完整编译）
 * @param expression 表达式字符串
 * @return 类型提示
 */
ExpressionHint analyze_expression_hint(const std::string& expression);

/**
 * @brief 分析表达式特征
 * @param expression 表达式字符串
 * @return 特征位掩码
 */
ExpressionFeature analyze_expression_features(const std::string& expression);

/**
 * @brief 对表达式进行词法分析
 * @param expression 表达式字符串
 * @param tokens 输出的 Token 序列
 * @return 是否成功
 */
bool tokenize_expression(const std::string& expression,
                         std::vector<ExpressionToken>* tokens);

// ============================================================================
// 表达式缓存
// ============================================================================

/**
 * @struct ExpressionCache
 * @brief 表达式缓存结构，用于存储预编译结果
 */
struct ExpressionCache {
    std::string expanded;               ///< 展开后的表达式
    ExpressionHint hint;                ///< 类型提示
    ExpressionFeature features;         ///< 特征
    std::vector<ExpressionToken> tokens; ///< Token 序列（可选）

    // 预编译结果（可选，用于高性能场景）
    std::shared_ptr<CompiledExpression> compiled;

    ExpressionCache() : hint(ExpressionHint::kUnknown),
                        features(ExpressionFeature::kNone) {}

    explicit ExpressionCache(const std::string& expr)
        : expanded(expr), hint(ExpressionHint::kUnknown),
          features(ExpressionFeature::kNone) {}
};

// ============================================================================
// 快速路径检测
// ============================================================================

/**
 * @struct QuickPathResult
 * @brief 快速路径检测结果
 */
struct QuickPathResult {
    enum class Type {
        kNone,          ///< 不是快速路径
        kStringLiteral, ///< 字符串字面量
        kIdentifier,    ///< 单个标识符
        kNumber,        ///< 纯数字
        kRatCall,       ///< rat() 调用
    };

    Type type = Type::kNone;
    std::string value;  ///< 提取的值（标识符名、字符串内容等）
    double number = 0.0; ///< 数值（当 type == kNumber）
};

/**
 * @brief 检测表达式是否可以走快速路径
 * @param expression 表达式字符串
 * @return 检测结果
 */
QuickPathResult detect_quick_path(const std::string& expression);

#endif // EXPRESSION_COMPILER_H
