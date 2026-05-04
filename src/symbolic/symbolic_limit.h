#ifndef SYMBOLIC_LIMIT_H
#define SYMBOLIC_LIMIT_H

#include "symbolic/symbolic_expression.h"
#include "symbolic/symbolic_expression_internal.h"

#include <string>
#include <optional>

namespace symbolic_limit {

// ============================================================================
// 不定型枚举
// ============================================================================

enum class IndeterminateForm {
    kNone,          // 定型
    kZeroOverZero,  // 0/0
    kInfOverInf,    // inf/inf
    kZeroTimesInf,  // 0*inf
    kInfMinusInf,   // inf-inf
    kZeroToZero,    // 0^0
    kInfToZero,     // inf^0
    kOneToInf       // 1^inf
};

// ============================================================================
// 极限结果结构
// ============================================================================

struct LimitResult {
    SymbolicExpression value;
    bool is_elementary = true;      // 是否为初等表达式
    bool is_definite = true;        // 极限是否存在且确定
    bool is_infinite = false;       // 是否为无穷大
    bool is_oscillating = false;    // 是否振荡（极限不存在）
    std::string method_used;

    static LimitResult elementary(const SymbolicExpression& val, const std::string& method) {
        LimitResult r;
        r.value = val;
        r.method_used = method;
        return r;
    }

    static LimitResult infinite(bool positive, const std::string& method) {
        LimitResult r;
        r.value = positive ? SymbolicExpression::parse("inf") : SymbolicExpression::parse("-inf");
        r.is_infinite = true;
        r.method_used = method;
        return r;
    }

    static LimitResult does_not_exist(const std::string& reason) {
        LimitResult r;
        r.is_definite = false;
        r.is_oscillating = true;
        r.method_used = "does_not_exist: " + reason;
        return r;
    }

    static LimitResult unknown() {
        LimitResult r;
        r.is_definite = false;
        r.method_used = "unknown";
        return r;
    }
};

// ============================================================================
// 符号极限引擎
// ============================================================================

class SymbolicLimitEngine {
public:
    /**
     * @brief 计算符号极限
     * @param expr 表达式
     * @param var 变量名
     * @param point 极限点
     * @param direction 方向: 0=双侧, -1=左极限, 1=右极限
     * @return 极限结果
     */
    LimitResult compute_limit(
        const SymbolicExpression& expr,
        const std::string& var,
        const BoundArgument& point,
        int direction = 0);

    /**
     * @brief 计算符号极限（字符串接口）
     */
    static std::optional<SymbolicExpression> limit(
        const std::string& expr_str,
        const std::string& var,
        const BoundArgument& point,
        int direction = 0);

private:
    // 检测不定型
    IndeterminateForm detect_indeterminate_form(
        const SymbolicExpression& expr,
        const std::string& var,
        const BoundArgument& point);

    // 应用 L'Hôpital 法则
    bool apply_lhopital(
        const SymbolicExpression& numerator,
        const SymbolicExpression& denominator,
        const std::string& var,
        const BoundArgument& point,
        int direction,
        LimitResult* result);

    // 尝试已知模式
    bool try_known_pattern(
        const SymbolicExpression& expr,
        const std::string& var,
        const BoundArgument& point,
        int direction,
        LimitResult* result);

    // 直接代入
    std::optional<double> try_direct_substitution(
        const SymbolicExpression& expr,
        const std::string& var,
        double point);

    // 处理无穷远点极限
    LimitResult limit_at_infinity(
        const SymbolicExpression& expr,
        const std::string& var,
        int direction);

    // 处理 1^inf 型极限
    bool handle_one_to_infinity(
        const SymbolicExpression& base,
        const SymbolicExpression& exponent,
        const std::string& var,
        const BoundArgument& point,
        LimitResult* result);

    // 检查表达式在极限点是否为无穷
    bool is_infinite_at_point(
        const SymbolicExpression& expr,
        const std::string& var,
        const BoundArgument& point);

    // 检查表达式在极限点是否为零
    bool is_zero_at_point(
        const SymbolicExpression& expr,
        const std::string& var,
        const BoundArgument& point);
};

// ============================================================================
// 辅助函数
// ============================================================================

/**
 * @brief 解析极限命令参数
 * @param args 参数列表 [expr, var, point, direction?]
 * @param expr 输出表达式
 * @param var 输出变量名
 * @param point 输出极限点
 * @param direction 输出方向
 * @return 是否解析成功
 */
bool parse_limit_arguments(
    const std::vector<std::string>& args,
    SymbolicExpression* expr,
    std::string* var,
    BoundArgument* point,
    int* direction);

}  // namespace symbolic_limit

#endif  // SYMBOLIC_LIMIT_H
