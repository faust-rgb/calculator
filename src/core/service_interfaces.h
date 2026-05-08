// ============================================================================
// 核心服务接口定义
// ============================================================================
//
// 定义计算器核心提供的各种服务接口，采用依赖注入的方式解耦模块间的依赖。
// 主要接口包括：
// - IEvaluationService: 表达式求值服务
// - ISymbolicService: 符号计算服务
// - IEnvironmentService: 环境管理服务
//
// 设计原则：
// - 接口分离：每个服务职责单一
// - 依赖注入：通过函数对象实现服务绑定
// - 易于测试：可替换为 mock 实现
// ============================================================================

#ifndef CORE_SERVICE_INTERFACES_H
#define CORE_SERVICE_INTERFACES_H

#include "matrix/matrix.h"
#include "symbolic/symbolic_expression.h"
#include "types/stored_value.h"
#include <string>
#include <vector>
#include <functional>
#include <map>

// 前向声明
template <typename T> class TFunctionAnalysis;
using FunctionAnalysis = TFunctionAnalysis<Scalar>;

/**
 * @struct IEvaluationService
 * @brief 核心提供的求值服务接口
 *
 * 提供表达式解析和求值功能，包括：
 * - 小数表达式解析
 * - 存储值求值
 * - 作用域求值器构造
 */
struct IEvaluationService {
    std::function<Scalar(const std::string&)> parse_decimal;           ///< 解析小数表达式
    std::function<StoredValue(const std::string&, bool)> evaluate_value;    ///< 求值为存储值
    std::function<Scalar(Scalar)> normalize_result;                    ///< 规范化结果

    // 作用域求值器构造
    std::function<std::function<Scalar(const std::vector<std::pair<std::string, Scalar>>&)>(const std::string&)> build_decimal_evaluator;   ///< 构建小数求值器
    std::function<std::function<Scalar(const std::vector<std::pair<std::string, StoredValue>>&)>(const std::string&)> build_scalar_evaluator;   ///< 构建标量求值器
    std::function<std::function<matrix::Matrix(const std::vector<std::pair<std::string, StoredValue>>&)>(const std::string&)> build_matrix_evaluator; ///< 构建矩阵求值器
};

/**
 * @struct ISymbolicService
 * @brief 符号计算服务接口
 *
 * 提供符号表达式的解析、简化和求值功能，包括：
 * - 符号表达式解析
 * - 内联函数展开
 * - 符号简化
 * - 函数分析构建
 */
struct ISymbolicService {
    std::function<void(const std::string&, bool, std::string*, SymbolicExpression*)> resolve_symbolic;    ///< 解析符号表达式
    std::function<std::string(const std::string&)> expand_inline;                                         ///< 展开内联函数
    std::function<std::string(const std::string&)> simplify_symbolic;                                     ///< 简化符号表达式
    std::function<Scalar(const SymbolicExpression&, const std::string&, Scalar)> evaluate_symbolic_at; ///< 在指定点求值
    std::function<std::vector<SymbolicExpression>(const std::string&)> parse_symbolic_expr_list;          ///< 解析符号表达式列表
    std::function<FunctionAnalysis(const std::string&)> build_analysis;                                   ///< 构建函数分析对象
};

/**
 * @struct IEnvironmentService
 * @brief 环境管理服务接口
 *
 * 提供变量、函数和状态的管理功能，包括：
 * - 变量/函数的查询、列出、清除
 * - 状态的保存和加载
 * - 配置模式的设置
 */
struct IEnvironmentService {
    std::function<bool(const std::string&)> has_variable;          ///< 检查变量是否存在
    std::function<bool(const std::string&)> has_function;          ///< 检查函数是否存在
    std::function<std::string()> list_variables;                   ///< 列出所有变量
    std::function<std::string()> list_functions;                   ///< 列出所有函数
    std::function<std::string(const std::string&)> clear_variable; ///< 清除指定变量
    std::function<std::string(const std::string&)> clear_function; ///< 清除指定函数
    std::function<std::string()> clear_all_variables;              ///< 清除所有变量
    std::function<std::string()> clear_all_functions;              ///< 清除所有函数

    // 状态管理
    std::function<std::string(const std::string&)> save_state;           ///< 保存状态到文件
    std::function<std::string(const std::string&)> load_state;           ///< 从文件加载状态
    std::function<std::string(const std::string&)> export_variable;      ///< 导出变量到文件
    std::function<std::string(const std::string&, bool)> execute_script; ///< 执行脚本代码
    std::function<std::string(const std::string&, bool)> execute_script_file; ///< 执行脚本文件

    // 模式与配置管理
    std::function<std::string(bool)> set_exact_mode;      ///< 设置精确模式
    std::function<std::string(bool)> set_symbolic_mode;   ///< 设置符号常量模式
    std::function<std::string(int)> set_precision;        ///< 设置显示精度
    std::function<std::string(bool)> set_hex_prefix;      ///< 设置十六进制前缀模式
    std::function<std::string(bool)> set_hex_uppercase;   ///< 设置十六进制大小写模式
};

/**
 * @struct CoreServices
 * @brief 汇总后的核心服务，作为各模块访问核心功能的唯一入口
 *
 * 包含三个主要服务接口：
 * - evaluation: 表达式求值服务
 * - symbolic: 符号计算服务
 * - env: 环境管理服务
 *
 * 还提供一些通用工具函数，如参数解析、矩阵处理等。
 */
struct CoreServices {
    IEvaluationService evaluation;   ///< 求值服务
    ISymbolicService symbolic;       ///< 符号计算服务
    IEnvironmentService env;         ///< 环境管理服务

    // 参数解析与辅助（保留作为通用工具）
    std::function<std::vector<std::string>(const std::vector<std::string>&, std::size_t, const std::vector<std::string>&)> parse_symbolic_vars; ///< 解析符号变量
    std::function<bool(const std::string&)> is_matrix_argument;     ///< 检查参数是否为矩阵
    std::function<matrix::Matrix(const std::string&, const std::string&)> parse_matrix_argument; ///< 解析矩阵参数
    std::function<std::string(const std::vector<std::string>&, bool)> render_plot; ///< 渲染图形
    std::function<bool(Scalar, Scalar)> is_integer_double;          ///< 检查是否为整数
    std::function<long long(Scalar)> round_to_long_long;            ///< 四舍五入为长整型
};


#endif
