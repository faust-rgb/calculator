// ============================================================================
// 多重积分器
// ============================================================================
//
// 本文件实现了多重积分的数值计算功能。
// 使用 Simpson 法则进行递归积分，支持任意维数和变边界。

#ifndef MULTIVARIABLE_INTEGRATOR_H
#define MULTIVARIABLE_INTEGRATOR_H

#include <functional>
#include <utility>
#include <vector>

/**
 * @class MultivariableIntegrator
 * @brief 多重积分器
 *
 * 用于计算多维区域上的数值积分。
 * 支持变边界（即边界可以是其他变量的函数）。
 *
 * 积分形式：
 *   ∫ ... ∫ f(x1, ..., xn) dx1 ... dxn
 *
 * 使用 Simpson 法则进行递归积分。
 */
class MultivariableIntegrator {
public:
    /// 被积函数类型：接受点坐标，返回函数值
    using Integrand = std::function<double(const std::vector<double>&)>;

    /// 边界函数类型：接受外层变量值，返回该维度的积分上下界
    using BoundFunc = std::function<std::pair<double, double>(const std::vector<double>&)>;

    /**
     * @brief 构造函数
     * @param integrand 被积函数
     * @throw std::runtime_error 当被积函数为空时抛出
     */
    explicit MultivariableIntegrator(Integrand integrand);

    /**
     * @brief 计算多重积分
     *
     * @param bounds 各维度的边界函数列表
     * @param subdivisions 各维度的分割数列表
     * @return 积分值
     * @throw std::runtime_error 当参数无效时抛出
     */
    double integrate(const std::vector<BoundFunc>& bounds,
                     const std::vector<int>& subdivisions) const;

private:
    /**
     * @brief 计算 Simpson 法则的权重
     * @param index 当前点索引
     * @param subdivisions 总分割数
     * @return Simpson 权重（1, 4, 2, 4, ..., 2, 4, 1）
     */
    static double simpson_weight(int index, int subdivisions);

    /**
     * @brief 规范化分割数为偶数
     * @param subdivisions 输入分割数
     * @return 偶数分割数（如果不为偶数则加 1）
     * @throw std::runtime_error 当分割数非正时抛出
     */
    static int normalize_subdivision_count(int subdivisions);

    /**
     * @brief 递归积分实现
     *
     * @param bounds 边界函数列表
     * @param subdivisions 分割数列表
     * @param point 当前积分点（用于传递给被积函数）
     * @param dimension 当前积分维度
     * @param accumulated_weight 累积的 Simpson 权重
     * @return 该维度及内层维度的积分值
     */
    double integrate_recursive(const std::vector<BoundFunc>& bounds,
                               const std::vector<int>& subdivisions,
                               std::vector<double>* point,
                               std::size_t dimension,
                               double accumulated_weight) const;

    Integrand integrand_;  ///< 被积函数
};

#endif
