// ============================================================================
// Gröbner 基计算模块
// ============================================================================
//
// 本文件实现了 Gröbner 基算法，用于多元多项式方程组的求解和简化。
// 主要功能包括：
// - Buchberger 算法计算 Gröbner 基
// - S-多项式计算和约简
// - 多项式理想的基变换
//
// Gröbner 基是多项式理想的特殊生成集，具有以下性质：
// - 可以判断多项式是否属于理想
// - 可以求解多项式方程组
// - 可以消除变量（消元定理）
//
// 相关文件：
// - groebner_basis.h: Gröbner 基接口定义
// - symbolic_solver.cpp: 方程求解器（使用 Gröbner 基）
// ============================================================================

#include "symbolic/algebra/groebner/groebner_basis.h"
#include "symbolic/core/symbolic_expression_internal.h"
#include "math/mymath.h"
#include <algorithm>
#include <cmath>
#include <map>
#include <vector>

namespace symbolic_groebner {

using namespace symbolic_expression_internal;

// 表示一个单项式: x^a * y^b * ...
typedef std::map<std::string, int> Monomial;

// 比较单项式 (词典序 lex)
struct MonomialComparator {
    std::vector<std::string> vars;
    bool operator()(const Monomial& a, const Monomial& b) const {
        for (const auto& var : vars) {
            int deg_a = a.count(var) ? a.at(var) : 0;
            int deg_b = b.count(var) ? b.at(var) : 0;
            if (deg_a != deg_b) return deg_a > deg_b;
        }
        return false;
    }
};

// 表示一个多元多项式: sum(coeff * monomial)
struct MultivariatePoly {
    std::vector<std::string> vars;
    std::map<Monomial, SymbolicExpression, MonomialComparator> terms;

    MultivariatePoly(const std::vector<std::string>& v) : vars(v), terms(MonomialComparator{v}) {}

    bool is_zero() const {
        for (const auto& [mon, coeff] : terms) {
            if (!expr_is_zero(coeff)) return false;
        }
        return true;
    }

    void simplify() {
        for (auto it = terms.begin(); it != terms.end(); ) {
            it->second = it->second.simplify();
            if (expr_is_zero(it->second)) it = terms.erase(it);
            else ++it;
        }
    }

    std::pair<Monomial, SymbolicExpression> leading_term() const {
        if (terms.empty()) return {{}, SymbolicExpression::number(0.0L)};
        return *terms.begin();
    }
};

// 将 SymbolicExpression 转换为 MultivariatePoly
MultivariatePoly from_expression(const SymbolicExpression& expr, const std::vector<std::string>& vars) {
    MultivariatePoly res(vars);

    if (expr.node_ == nullptr) {
        return res;
    }

    // 处理常数
    Scalar val = Scalar(0.0L);
    if (expr.is_number(&val)) {
        res.terms[{}] = expr;
        return res;
    }

    // 处理变量
    if (expr.node_->type == NodeType::kVariable) {
        Monomial m;
        m[expr.node_->text] = 1;
        res.terms[m] = SymbolicExpression::number(Scalar(1.0L));
        return res;
    }

    // 处理加法: 递归处理每个子项
    if (expr.node_->type == NodeType::kAdd) {
        auto left = from_expression(SymbolicExpression(expr.node_->left), vars);
        auto right = from_expression(SymbolicExpression(expr.node_->right), vars);
        // 合并两个多项式
        for (const auto& [mon, coeff] : left.terms) {
            res.terms[mon] = (res.terms.count(mon) ? res.terms[mon] : SymbolicExpression::number(0.0L)) + coeff;
        }
        for (const auto& [mon, coeff] : right.terms) {
            res.terms[mon] = (res.terms.count(mon) ? res.terms[mon] : SymbolicExpression::number(0.0L)) + coeff;
        }
        res.simplify();
        return res;
    }

    // 处理减法
    if (expr.node_->type == NodeType::kSubtract) {
        auto left = from_expression(SymbolicExpression(expr.node_->left), vars);
        auto right = from_expression(SymbolicExpression(expr.node_->right), vars);
        for (const auto& [mon, coeff] : left.terms) {
            res.terms[mon] = (res.terms.count(mon) ? res.terms[mon] : SymbolicExpression::number(0.0L)) + coeff;
        }
        for (const auto& [mon, coeff] : right.terms) {
            res.terms[mon] = (res.terms.count(mon) ? res.terms[mon] : SymbolicExpression::number(0.0L)) - coeff;
        }
        res.simplify();
        return res;
    }

    // 处理乘法: 多项式乘法
    if (expr.node_->type == NodeType::kMultiply) {
        auto left = from_expression(SymbolicExpression(expr.node_->left), vars);
        auto right = from_expression(SymbolicExpression(expr.node_->right), vars);

        // 多项式乘法
        for (const auto& [mon1, coeff1] : left.terms) {
            for (const auto& [mon2, coeff2] : right.terms) {
                // 合并单项式
                Monomial prod_mon = mon1;
                for (const auto& [v, d] : mon2) {
                    prod_mon[v] += d;
                }
                res.terms[prod_mon] = (res.terms.count(prod_mon) ? res.terms[prod_mon] : SymbolicExpression::number(0.0L)) + coeff1 * coeff2;
            }
        }
        res.simplify();
        return res;
    }

    // 处理幂次: x^n
    if (expr.node_->type == NodeType::kPower) {
        SymbolicExpression base(expr.node_->left);
        SymbolicExpression exp_expr(expr.node_->right);

        Scalar exp_val = Scalar(0.0L);
        if (exp_expr.is_number(&exp_val) && exp_val >= Scalar(0.0L)) {
            // 检查是否为整数
            Scalar exp_floor = mymath::precise128::floor(exp_val);
            if (mymath::precise128::abs(exp_val - exp_floor) < Scalar(1e-15L)) {
                int n = static_cast<int>(static_cast<long double>(exp_floor));
                auto base_poly = from_expression(base, vars);

                // 计算 base_poly^n
                res.terms[{}] = SymbolicExpression::number(Scalar(1.0L));  // 初始化为 1
                for (int i = 0; i < n; ++i) {
                    MultivariatePoly temp(vars);
                    for (const auto& [mon1, coeff1] : res.terms) {
                        for (const auto& [mon2, coeff2] : base_poly.terms) {
                            Monomial prod_mon = mon1;
                            for (const auto& [v, d] : mon2) {
                                prod_mon[v] += d;
                            }
                            temp.terms[prod_mon] = (temp.terms.count(prod_mon) ? temp.terms[prod_mon] : SymbolicExpression::number(0.0L)) + coeff1 * coeff2;
                        }
                    }
                    res = temp;
                }
                res.simplify();
                return res;
            }
        }
    }

    // 处理负号
    if (expr.node_->type == NodeType::kNegate) {
        auto inner = from_expression(SymbolicExpression(expr.node_->left), vars);
        for (const auto& [mon, coeff] : inner.terms) {
            res.terms[mon] = -coeff;
        }
        return res;
    }

    // 处理函数节点：将函数视为新的变量
    // 例如 sin(x) 被视为关于 vars 的多项式中的一个"系数"
    if (expr.node_->type == NodeType::kFunction) {
        // 函数表达式作为常数项（系数为函数本身）
        res.terms[{}] = expr;
        return res;
    }

    // 处理除法：如果分母不含变量，作为系数处理
    if (expr.node_->type == NodeType::kDivide) {
        SymbolicExpression num(expr.node_->left);
        SymbolicExpression den(expr.node_->right);

        // 检查分母是否不含任何变量
        std::string den_str = den.to_string();
        bool den_is_const = true;
        for (const auto& v : vars) {
            if (den_str.find(v) != std::string::npos) {
                den_is_const = false;
                break;
            }
        }

        if (den_is_const) {
            // 分母是常数，分子可能是多项式
            auto num_poly = from_expression(num, vars);
            for (const auto& [mon, coeff] : num_poly.terms) {
                res.terms[mon] = make_divide(coeff, den).simplify();
            }
            return res;
        }

        // 分母含变量，作为整体处理（视为常数项）
        res.terms[{}] = expr;
        return res;
    }

    // 其他无法处理的表达式，作为常数项处理
    // 这确保不会丢失信息，虽然可能不是最优的多项式表示
    res.terms[{}] = expr;
    return res;
}

SymbolicExpression to_expression(const MultivariatePoly& poly) {
    SymbolicExpression res = SymbolicExpression::number(0.0L);
    for (const auto& [mon, coeff] : poly.terms) {
        SymbolicExpression term = coeff;
        for (const auto& [var, deg] : mon) {
            if (deg > 0) {
                term = term * make_power(SymbolicExpression::variable(var), SymbolicExpression::number(deg));
            }
        }
        res = res + term;
    }
    return res.simplify();
}

// 多元多项式除法 (Reduction)
MultivariatePoly reduce(MultivariatePoly p, const std::vector<MultivariatePoly>& basis) {
    MultivariatePoly res(p.vars);
    bool changed = true;
    while (changed) {
        changed = false;
        if (p.is_zero()) break;
        auto [lt_p, c_p] = p.leading_term();
        
        for (const auto& f : basis) {
            auto [lt_f, c_f] = f.leading_term();
            // 检查 lt_f 是否整除 lt_p
            bool divides = true;
            Monomial quot_mon;
            for (const auto& v : p.vars) {
                int d_p = lt_p.count(v) ? lt_p.at(v) : 0;
                int d_f = lt_f.count(v) ? lt_f.at(v) : 0;
                if (d_p < d_f) { divides = false; break; }
                if (d_p > d_f) quot_mon[v] = d_p - d_f;
            }
            
            if (divides) {
                SymbolicExpression quot_coeff = (c_p / c_f).simplify();
                // p = p - (quot_coeff * quot_mon) * f
                for (const auto& [f_mon, f_coeff] : f.terms) {
                    Monomial new_mon = quot_mon;
                    for (const auto& [v, d] : f_mon) new_mon[v] += d;
                    p.terms[new_mon] = (p.terms[new_mon] - quot_coeff * f_coeff).simplify();
                }
                p.simplify();
                changed = true;
                break;
            }
        }
        
        if (!changed && !p.is_zero()) {
            // 将领先项移入结果，继续处理剩余项
            auto [lt, c] = p.leading_term();
            res.terms[lt] = c;
            p.terms.erase(p.terms.begin());
            changed = true;
        }
    }
    return res;
}

// 计算 S-多项式
MultivariatePoly s_polynomial(const MultivariatePoly& f, const MultivariatePoly& g) {
    auto [lt_f, c_f] = f.leading_term();
    auto [lt_g, c_g] = g.leading_term();

    // 计算 LCM(leading monomial of f, leading monomial of g)
    Monomial lcm_mon;
    for (const auto& v : f.vars) {
        lcm_mon[v] = std::max(lt_f.count(v) ? lt_f.at(v) : 0, lt_g.count(v) ? lt_g.at(v) : 0);
    }

    // 计算 x^a / LT(f) 和 x^b / LT(g)
    Monomial quot_f, quot_g;
    for (const auto& v : f.vars) {
        int d_lcm = lcm_mon.count(v) ? lcm_mon.at(v) : 0;
        int d_f = lt_f.count(v) ? lt_f.at(v) : 0;
        int d_g = lt_g.count(v) ? lt_g.at(v) : 0;
        if (d_lcm > d_f) quot_f[v] = d_lcm - d_f;
        if (d_lcm > d_g) quot_g[v] = d_lcm - d_g;
    }

    // S(f, g) = (LCM / LT(f)) * f / c_f - (LCM / LT(g)) * g / c_g
    //         = (x^quot_f * f) / c_f - (x^quot_g * g) / c_g
    MultivariatePoly res(f.vars);

    // 第一项: (x^quot_f * f) / c_f
    for (const auto& [mon, coeff] : f.terms) {
        Monomial new_mon = quot_f;
        for (const auto& [v, d] : mon) new_mon[v] += d;
        res.terms[new_mon] = (res.terms.count(new_mon) ? res.terms[new_mon] : SymbolicExpression::number(0.0L)) + coeff / c_f;
    }

    // 第二项: -(x^quot_g * g) / c_g
    for (const auto& [mon, coeff] : g.terms) {
        Monomial new_mon = quot_g;
        for (const auto& [v, d] : mon) new_mon[v] += d;
        res.terms[new_mon] = (res.terms.count(new_mon) ? res.terms[new_mon] : SymbolicExpression::number(0.0L)) - coeff / c_g;
    }

    res.simplify();
    return res;
}

std::vector<SymbolicExpression> compute_groebner_basis(
    const std::vector<SymbolicExpression>& polynomials,
    const std::vector<std::string>& variables) {

    std::vector<MultivariatePoly> F;
    for (const auto& p : polynomials) F.push_back(from_expression(p, variables));

    // 添加终止条件：最大迭代次数和最大基大小
    constexpr int kMaxIterations = 1000;
    constexpr size_t kMaxBasisSize = 100;
    int iterations = 0;

    bool changed = true;
    while (changed && iterations < kMaxIterations && F.size() < kMaxBasisSize) {
        changed = false;
        ++iterations;

        std::vector<std::pair<size_t, size_t>> pairs;
        for (size_t i = 0; i < F.size(); ++i) {
            for (size_t j = i + 1; j < F.size(); ++j) pairs.push_back({i, j});
        }

        // 使用简单的选择策略：优先处理低索引的对
        for (const auto& pair : pairs) {
            if (F.size() >= kMaxBasisSize) break;

            MultivariatePoly s = s_polynomial(F[pair.first], F[pair.second]);
            MultivariatePoly r = reduce(s, F);
            if (!r.is_zero()) {
                F.push_back(r);
                changed = true;
                break;
            }
        }
    }

    std::vector<SymbolicExpression> res;
    for (const auto& p : F) res.push_back(to_expression(p));
    return res;
}

SymbolicExpression groebner_reduce(
    const SymbolicExpression& p,
    const std::vector<SymbolicExpression>& basis,
    const std::vector<std::string>& variables) {
    std::vector<MultivariatePoly> B;
    for (const auto& b : basis) B.push_back(from_expression(b, variables));
    return to_expression(reduce(from_expression(p, variables), B));
}

} // namespace symbolic_groebner
