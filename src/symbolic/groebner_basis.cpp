#include "symbolic/groebner_basis.h"
#include "symbolic/symbolic_expression_internal.h"
#include "math/mymath.h"
#include <algorithm>
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
    std::map<Monomial, SymbolicExpression, MonomialComparator> terms;
    std::vector<std::string> vars;

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
    // 简化处理：目前仅支持已经展开的多项式
    // 实际实现需要递归遍历表达式树
    // 这里做个非常简单的提取
    long double val = 0.0L;
    if (expr.is_number(&val)) {
        res.terms[{}] = expr;
    } else if (expr.node_->type == NodeType::kVariable) {
        Monomial m; m[expr.node_->text] = 1;
        res.terms[m] = SymbolicExpression::number(1.0L);
    }
    // TODO: 递归解析 Add, Multiply, Power
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
    
    Monomial lcm_mon;
    for (const auto& v : f.vars) {
        lcm_mon[v] = std::max(lt_f.count(v) ? lt_f.at(v) : 0, lt_g.count(v) ? lt_g.at(v) : 0);
    }
    
    MultivariatePoly res(f.vars);
    // TODO: 实现 S-多项式具体合并逻辑
    return res;
}

std::vector<SymbolicExpression> compute_groebner_basis(
    const std::vector<SymbolicExpression>& polynomials,
    const std::vector<std::string>& variables) {
    
    std::vector<MultivariatePoly> F;
    for (const auto& p : polynomials) F.push_back(from_expression(p, variables));
    
    bool changed = true;
    while (changed) {
        changed = false;
        std::vector<std::pair<size_t, size_t>> pairs;
        for (size_t i = 0; i < F.size(); ++i) {
            for (size_t j = i + 1; j < F.size(); ++j) pairs.push_back({i, j});
        }
        
        for (const auto& pair : pairs) {
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
