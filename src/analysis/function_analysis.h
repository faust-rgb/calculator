#ifndef FUNCTION_ANALYSIS_H
#define FUNCTION_ANALYSIS_H

#include <string>
#include <vector>

class Calculator;

/**
 * @file function_analysis.h
 * @brief 函数分析库，提供微积分运算功能
 *
 * 支持单变量函数的数值分析，包括求值、微分、积分和极值查找。
 * 所有运算基于数值方法，不需要符号微分。
 */

/**
 * @struct ExtremumPoint
 * @brief 极值点数据结构
 */
struct ExtremumPoint {
    double x = 0.0;           ///< 极值点位置
    double value = 0.0;       ///< 函数在该点的值
    bool is_maximum = false;  ///< true 表示极大值，false 表示极小值
};

/**
 * @class FunctionAnalysis
 * @brief 函数分析器
 *
 * 用于分析单变量数学函数，支持：
 * - 函数求值
 * - 数值微分（一阶和二阶导数）
 * - 极限计算
 * - 定积分和不定积分
 * - 极值点查找
 */
class FunctionAnalysis {
public:
    /**
     * @brief 构造函数
     * @param variable_name 变量名，默认为 "x"
     * @throw std::runtime_error 当变量名不合法时抛出
     */
    explicit FunctionAnalysis(std::string variable_name = "x");

    /**
     * @brief 定义要分析的函数
     * @param expression 数学表达式字符串
     * @throw std::runtime_error 当表达式为空时抛出
     */
    void define(const std::string& expression);

    /**
     * @brief 计算函数值
     * @param x 自变量值
     * @return f(x)
     */
    double evaluate(double x) const;

    /**
     * @brief 计算导数（数值微分）
     * @param x 求导点
     * @return f'(x)
     *
     * 使用自适应步长的中心差分和 4 层 Richardson 外推。
     */
    double derivative(double x) const;

    /**
     * @brief 计算极限
     * @param x 极限点
     * @param direction 方向：-1 左极限，1 右极限，0 双侧极限
     * @return 极限值
     * @throw std::runtime_error 当极限不存在或不收敛时抛出
     */
    double limit(double x, int direction = 0) const;

    /**
     * @brief 计算定积分
     * @param lower_bound 下限
     * @param upper_bound 上限
     * @return ∫[lower, upper] f(x) dx
     *
     * 使用自适应 Gauss-Kronrod G7-K15 积分法。
     */
    double definite_integral(double lower_bound, double upper_bound) const;

    /**
     * @brief 计算不定积分在指定点的值
     * @param x 计算点
     * @param anchor 积分锚点（积分常数确定的点）
     * @param constant 积分常数
     * @return F(x) = C + ∫[anchor, x] f(t) dt
     */
    double indefinite_integral_at(double x,
                                  double anchor = 0.0,
                                  double constant = 0.0) const;

    /**
     * @brief 查找区间内的极值点
     * @param left_bound 左边界
     * @param right_bound 右边界
     * @param scan_segments 扫描分段数，默认 128
     * @return 极值点列表
     * @throw std::runtime_error 当参数无效时抛出
     *
     * 算法：
     * 1. 扫描区间寻找导数变号点
     * 2. 使用二分法精确定位驻点
     * 3. 通过二阶导数判断极大/极小值
     */
    std::vector<ExtremumPoint> solve_extrema(double left_bound,
                                             double right_bound,
                                             int scan_segments = 128) const;

    /** @brief 获取当前函数表达式 */
    const std::string& expression() const;

    /** @brief 获取变量名 */
    const std::string& variable_name() const;

private:
    /**
     * @brief 带变量赋值的求值
     * @param x 变量值
     * @return f(x)
     *
     * 内部使用 Calculator 实例，先将变量赋值，再求值表达式。
     */
    double evaluate_with_variable(double x) const;

    /**
     * @brief 计算二阶导数
     * @param x 求导点
     * @return f''(x)
     *
     * 使用中心差分法。
     */
    double second_derivative(double x) const;

    /**
     * @brief 二分法查找驻点（导数为零的点）
     * @param left 左端点（导数与右端点异号）
     * @param right 右端点
     * @return 驻点位置
     */
    double bisect_stationary_point(double left, double right) const;

    /**
     * @brief 自适应 Gauss-Kronrod G7-K15 积分
     * @param left 左端点
     * @param right 右端点
     * @param eps 精度要求
     * @param max_depth 最大递归深度
     * @return 积分值
     */
    double adaptive_gauss_kronrod(double left,
                                  double right,
                                  double eps,
                                  int max_depth) const;

    /**
     * @brief 自适应 Gauss-Kronrod 积分的递归实现
     */
    double adaptive_gauss_kronrod_recursive(double left,
                                            double right,
                                            double eps,
                                            double whole,
                                            double error,
                                            int depth) const;

    /**
     * @brief G7-K15 单区间积分，并输出误差估计
     */
    double gauss_kronrod_15(double left,
                            double right,
                            double* error_estimate) const;

    std::string expression_;      ///< 函数表达式
    std::string variable_name_;   ///< 变量名
};

#endif
