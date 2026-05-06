#ifndef FUNCTION_ANALYSIS_H
#define FUNCTION_ANALYSIS_H

#include <functional>
#include <list>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

/**
 * @file function_analysis.h
 * @brief 函数分析库，提供微积分运算功能 (泛型版)
 */

/**
 * @struct TExtremumPoint
 * @brief 泛型极值点数据结构
 */
template <typename T>
struct TExtremumPoint {
    T x = T(0);               ///< 极值点位置
    T value = T(0);           ///< 函数在该点的值
    bool is_maximum = false;  ///< true 表示极大值，false 表示极小值
};

using ExtremumPoint = TExtremumPoint<double>;

/**
 * @class TFunctionAnalysis
 * @brief 泛型函数分析器
 */
template <typename T>
class TFunctionAnalysis {
public:
    explicit TFunctionAnalysis(std::string variable_name = "x");
    TFunctionAnalysis(const TFunctionAnalysis& other);
    TFunctionAnalysis& operator=(const TFunctionAnalysis& other);
    TFunctionAnalysis(TFunctionAnalysis&& other) noexcept;
    TFunctionAnalysis& operator=(TFunctionAnalysis&& other) noexcept;
    ~TFunctionAnalysis();

    void define(const std::string& expression);

    void set_evaluator(std::function<T(const std::vector<std::pair<std::string, T>>&)> evaluator);

    T evaluate(T x) const;

    /**
     * @brief 计算导数（数值微分）
     * 使用自适应步长的中心差分和 4 层 Richardson 外推。
     */
    T derivative(T x) const;

    /**
     * @brief 计算极限
     * @param x 极限点
     * @param direction 方向：-1 左极限，1 右极限，0 双侧极限
     */
    T limit(T x, int direction = 0) const;

    /**
     * @brief 计算定积分
     * 使用自适应 Gauss-Kronrod G7-K15 积分法。
     */
    T definite_integral(T lower_bound, T upper_bound) const;

    T indefinite_integral_at(T x, T anchor = T(0), T constant = T(0)) const;

    std::vector<TExtremumPoint<T>> solve_extrema(T left_bound,
                                                  T right_bound,
                                                  int scan_segments = 128) const;

    const std::string& expression() const;
    const std::string& variable_name() const;

private:
    T evaluate_with_variable(T x) const;
    T second_derivative(T x) const;
    T bisect_stationary_point(T left, T right) const;

    T adaptive_gauss_kronrod(T left, T right, T eps, int max_depth) const;

    T adaptive_gauss_kronrod_recursive(T left, T right, T eps, T whole, T error, int depth) const;

    T gauss_kronrod_15(T left, T right, T* error_estimate) const;

    // 高精度积分辅助函数（用于 PreciseDecimal）
    T simpson_rule(T a, T b) const;
    T adaptive_simpson_precise(T a, T b, T eps, int max_depth) const;
    T adaptive_simpson_recursive(T a, T b, T whole, T left, T right, T eps, int depth) const;

    T compute_numerical_limit(T x, int direction) const;

    std::string expression_;
    std::string variable_name_;
    std::function<T(const std::vector<std::pair<std::string, T>>&)> evaluator_;
    mutable std::list<std::pair<std::string, T>> evaluation_cache_entries_;
    mutable std::unordered_map<std::string, typename std::list<std::pair<std::string, T>>::iterator> evaluation_cache_index_;
};

using FunctionAnalysis = TFunctionAnalysis<double>;

#endif
