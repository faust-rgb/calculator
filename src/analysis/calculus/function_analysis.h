#ifndef FUNCTION_ANALYSIS_H
#define FUNCTION_ANALYSIS_H

#include <functional>
#include <list>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "app/scalar_type.h"

class Calculator;

/**
 * @file function_analysis.h
 * @brief 函数分析库，提供微积分运算功能
 */

/**
 * @struct ExtremumPoint
 * @brief 极值点数据结构
 */
struct ExtremumPoint {
    Scalar x = Scalar(0);               ///< 极值点位置
    Scalar value = Scalar(0);           ///< 函数在该点的值
    bool is_maximum = false;  ///< true 表示极大值，false 表示极小值
};

/**
 * @class FunctionAnalysis
 * @brief 函数分析器
 */
class FunctionAnalysis {
public:
    explicit FunctionAnalysis(std::string variable_name = "x");
    FunctionAnalysis(const FunctionAnalysis& other);
    FunctionAnalysis& operator=(const FunctionAnalysis& other);
    FunctionAnalysis(FunctionAnalysis&& other) noexcept;
    FunctionAnalysis& operator=(FunctionAnalysis&& other) noexcept;
    ~FunctionAnalysis();

    void define(const std::string& expression);

    void set_evaluator(std::function<Scalar(const std::vector<std::pair<std::string, Scalar>>&)> evaluator);

    /**
     * @brief 设置外部变量查找函数
     * 用于极限计算时查询非极限变量的值
     * @param lookup 变量查找函数，返回变量值或默认值
     */
    void set_variable_lookup(std::function<Scalar(const std::string&)> lookup);

    Scalar evaluate(Scalar x) const;

    /**
     * @brief 计算导数（数值微分）
     * 使用自适应步长的中心差分和 4 层 Richardson 外推。
     */
    Scalar derivative(Scalar x) const;

    /**
     * @brief 计算极限
     * @param x 极限点
     * @param direction 方向：-1 左极限，1 右极限，0 双侧极限
     */
    Scalar limit(Scalar x, int direction = 0) const;

    /**
     * @brief 计算定积分
     * 使用自适应 Gauss-Kronrod G7-K15 积分法。
     */
    Scalar definite_integral(Scalar lower_bound, Scalar upper_bound) const;

    Scalar indefinite_integral_at(Scalar x, Scalar anchor = Scalar(0), Scalar constant = Scalar(0)) const;

    std::vector<ExtremumPoint> solve_extrema(Scalar left_bound,
                                                  Scalar right_bound,
                                                  int scan_segments = 128) const;

    const std::string& expression() const;
    const std::string& variable_name() const;

private:
    Scalar evaluate_with_variable(Scalar x) const;
    void ensure_evaluator_initialized() const;
    Scalar second_derivative(Scalar x) const;
    Scalar bisect_stationary_point(Scalar left, Scalar right) const;

    Scalar adaptive_gauss_kronrod(Scalar left, Scalar right, Scalar eps, int max_depth) const;

    Scalar adaptive_gauss_kronrod_recursive(Scalar left, Scalar right, Scalar eps, Scalar whole, Scalar error, int depth) const;

    Scalar gauss_kronrod_15(Scalar left, Scalar right, Scalar* error_estimate) const;

    // 高精度积分辅助函数
    Scalar simpson_rule(Scalar a, Scalar b) const;
    Scalar adaptive_simpson_precise(Scalar a, Scalar b, Scalar eps, int max_depth) const;
    Scalar adaptive_simpson_recursive(Scalar a, Scalar b, Scalar whole, Scalar left, Scalar right, Scalar eps, int depth) const;

    Scalar compute_numerical_limit(Scalar x, int direction) const;

    std::string expression_;
    std::string variable_name_;
    mutable std::function<Scalar(const std::vector<std::pair<std::string, Scalar>>&)> evaluator_;
    mutable std::shared_ptr<Calculator> fallback_calculator_;
    mutable std::function<Scalar(const std::string&)> variable_lookup_;  ///< 外部变量查找函数
    mutable std::list<std::pair<std::string, Scalar>> evaluation_cache_entries_;
    mutable std::unordered_map<std::string, std::list<std::pair<std::string, Scalar>>::iterator> evaluation_cache_index_;
};

#endif
