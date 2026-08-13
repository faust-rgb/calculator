/**
 * @file symbolic_node_types.h
 * @brief 符号表达式节点类型公共接口
 *
 * 提供符号表达式节点类型的公共定义，供跨模块安全使用。
 * 包括：
 * - NodeType 枚举：表达式节点类型
 * - BoundKind 枚举：边界值类型
 * - BoundArgument 结构体：边界参数解析结果
 *
 * 此文件不依赖任何内部实现头文件，模块间可以安全包含。
 */

#ifndef SYMBOLIC_PUBLIC_NODE_TYPES_H
#define SYMBOLIC_PUBLIC_NODE_TYPES_H

#include "types/scalar_type.h"
#include <string>

// ============================================================================
// 节点类型枚举
// ============================================================================

/**
 * @enum NodeType
 * @brief 表达式树节点类型
 *
 * 每种类型对应一种表达式构造：
 * - kNumber: 数值常量
 * - kVariable: 变量
 * - kAdd/kSubtract/kMultiply/kDivide/kPower: 二元运算
 * - kNegate: 一元取负
 * - kFunction: 函数调用
 * - kVector: 向量表达式
 * - kTensor: 张量（矩阵）表达式
 * - kDifferentialOp: 微分算子
 * - kRootOf: 代数数
 */
enum class NodeType {
    kNumber,         ///< 数值常量节点
    kVariable,       ///< 变量节点
    kPi,             ///< 精确常数 pi
    kE,              ///< 精确常数 e
    kInfinity,       ///< 无穷大节点 (+inf 或 -inf)
    kAdd,            ///< 加法节点
    kSubtract,       ///< 减法节点
    kMultiply,       ///< 乘法节点
    kDivide,         ///< 除法节点
    kPower,          ///< 幂运算节点
    kNegate,         ///< 取负节点
    kFunction,       ///< 函数调用节点
    kVector,         ///< 向量节点
    kTensor,         ///< 张量节点
    kDifferentialOp, ///< 微分算子节点
    kRootOf,         ///< 代数数节点
};

// ============================================================================
// 边界参数类型
// ============================================================================

/**
 * @enum BoundKind
 * @brief 边界值类型
 */
enum class BoundKind {
    kFinite,   ///< 有限数值
    kPosInf,   ///< 正无穷大 (+inf)
    kNegInf,   ///< 负无穷大 (-inf)
};

/**
 * @struct BoundArgument
 * @brief 边界参数解析结果
 *
 * 统一表示有限值和无穷大，用于 integral、limit 等命令的边界解析。
 */
struct BoundArgument {
    BoundKind kind = BoundKind::kFinite;
    Scalar value = Scalar(0);  ///< 仅当 kind == kFinite 时有效

    bool is_finite() const { return kind == BoundKind::kFinite; }
    bool is_infinite() const { return kind != BoundKind::kFinite; }
    bool is_pos_inf() const { return kind == BoundKind::kPosInf; }
    bool is_neg_inf() const { return kind == BoundKind::kNegInf; }
    Scalar to_scalar() const;

    static BoundArgument finite(Scalar v);
    static BoundArgument pos_inf();
    static BoundArgument neg_inf();
};

/**
 * @brief 解析边界参数字符串
 */
BoundArgument parse_bound_argument(const std::string& text);

/**
 * @brief 检查字符串是否为无穷大字面量
 */
bool is_infinity_literal(const std::string& text);

#endif  // SYMBOLIC_PUBLIC_NODE_TYPES_H