// ============================================================================
// 变换规则框架实现
// ============================================================================

#include "symbolic/transformation/rules/transform_rules.h"
#include "symbolic/transformation/transform_common.h"
#include "symbolic/core/symbolic_expression_internal.h"
#include "math/mymath.h"

#include <algorithm>

namespace transform_rules {

using namespace symbolic_expression_internal;
using Scalar = mymath::Scalar;

namespace dsl {

SymbolicExpression _any(const std::string& name) {
    return SymbolicExpression::variable("__wild_any__" + name);
}

SymbolicExpression _c(const std::string& name) {
    return SymbolicExpression::variable("__wild_const__" + name);
}

SymbolicExpression _var(const std::string& name) {
    return SymbolicExpression::variable("__wild_var__" + name);
}

SymbolicExpression _num(const std::string& name) {
    return SymbolicExpression::variable("__wild_num__" + name);
}

SymbolicExpression _step(const SymbolicExpression& arg) {
    return make_function("step", arg);
}

SymbolicExpression _delta(const SymbolicExpression& arg) {
    return make_function("delta", arg);
}

SymbolicExpression _exp(const SymbolicExpression& arg) {
    return make_function("exp", arg);
}

SymbolicExpression _sin(const SymbolicExpression& arg) {
    return make_function("sin", arg);
}

SymbolicExpression _cos(const SymbolicExpression& arg) {
    return make_function("cos", arg);
}

SymbolicExpression _sinc(const SymbolicExpression& arg) {
    return make_function("sinc", arg);
}

SymbolicExpression _sinh(const SymbolicExpression& arg) {
    return make_function("sinh", arg);
}

SymbolicExpression _cosh(const SymbolicExpression& arg) {
    return make_function("cosh", arg);
}

SymbolicExpression _sgn(const SymbolicExpression& arg) {
    return make_function("sgn", arg);
}

SymbolicExpression _gamma(const SymbolicExpression& arg) {
    return make_function("gamma", arg);
}

SymbolicExpression _rect(const SymbolicExpression& arg) {
    return make_function("rect", arg);
}

SymbolicExpression _tri(const SymbolicExpression& arg) {
    return make_function("tri", arg);
}

SymbolicExpression _abs(const SymbolicExpression& arg) {
    return make_function("abs", arg);
}

SymbolicExpression _sqrt(const SymbolicExpression& arg) {
    return make_function("sqrt", arg);
}

}  // namespace dsl

namespace {

void flatten_ast_terms(NodeType type, const SymbolicExpression& expr, std::vector<SymbolicExpression>* list) {
    if (expr.node_ && expr.node_->type == type) {
        flatten_ast_terms(type, SymbolicExpression(expr.node_->left), list);
        flatten_ast_terms(type, SymbolicExpression(expr.node_->right), list);
    } else {
        list->push_back(expr);
    }
}

std::vector<SymbolicExpression> get_node_args(const SymbolicExpression::Node* node) {
    std::vector<SymbolicExpression> args;
    if (!node) return args;
    if (!node->children.empty()) {
        for (const auto& child : node->children) {
            args.push_back(SymbolicExpression(child));
        }
    } else if (node->left) {
        args.push_back(SymbolicExpression(node->left));
    }
    return args;
}

} // namespace

bool match_pattern(const SymbolicExpression& target,
                   const SymbolicExpression& pattern,
                   const std::string& input_var,
                   BindingsMap* bindings) {
    if (pattern.node_ == nullptr || target.node_ == nullptr || bindings == nullptr) {
        return false;
    }

    if (pattern.node_->type == NodeType::kVariable) {
        const std::string& name = pattern.node_->text;
        if (name.rfind("__wild_any__", 0) == 0) {
            std::string var_id = name.substr(12);
            auto it = bindings->find(var_id);
            if (it != bindings->end()) {
                return match_zero_expression((target - it->second).simplify());
            }
            (*bindings)[var_id] = target;
            return true;
        }
        if (name.rfind("__wild_const__", 0) == 0) {
            if (!is_constant_expression(target, input_var)) return false;
            std::string var_id = name.substr(14);
            auto it = bindings->find(var_id);
            if (it != bindings->end()) {
                return match_zero_expression((target - it->second).simplify());
            }
            (*bindings)[var_id] = target;
            return true;
        }
        if (name.rfind("__wild_var__", 0) == 0) {
            if (!target.is_variable_named(input_var)) return false;
            std::string var_id = name.substr(12);
            if (!var_id.empty()) {
                (*bindings)[var_id] = target;
            }
            return true;
        }
        if (name.rfind("__wild_num__", 0) == 0) {
            Scalar s = Scalar(0.0L);
            if (!target.is_number(&s)) return false;
            std::string var_id = name.substr(12);
            if (!var_id.empty()) {
                (*bindings)[var_id] = target;
            }
            return true;
        }
        return target.is_variable_named(name);
    }

    Scalar p_num = Scalar(0.0L), t_num = Scalar(0.0L);
    if (pattern.is_number(&p_num) && target.is_number(&t_num)) {
        return mymath::is_near_zero(p_num - t_num, Scalar(1e-9L));
    }

    // 隐式单位元退化匹配：_c("c") * pattern 形式匹配 1 * target
    if (pattern.node_->type == NodeType::kMultiply) {
        SymbolicExpression p_left(pattern.node_->left);
        SymbolicExpression p_right(pattern.node_->right);
        if (p_left.node_ && p_left.node_->type == NodeType::kVariable &&
            p_left.node_->text.rfind("__wild_const__", 0) == 0) {
            std::string bid = p_left.node_->text.substr(14);
            BindingsMap backup = *bindings;
            if (match_pattern(target, p_right, input_var, bindings)) {
                if (bindings->find(bid) == bindings->end()) {
                    (*bindings)[bid] = SymbolicExpression::number(Scalar(1.0L));
                    return true;
                }
            }
            *bindings = backup;
        }
        if (p_right.node_ && p_right.node_->type == NodeType::kVariable &&
            p_right.node_->text.rfind("__wild_const__", 0) == 0) {
            std::string bid = p_right.node_->text.substr(14);
            BindingsMap backup = *bindings;
            if (match_pattern(target, p_left, input_var, bindings)) {
                if (bindings->find(bid) == bindings->end()) {
                    (*bindings)[bid] = SymbolicExpression::number(Scalar(1.0L));
                    return true;
                }
            }
            *bindings = backup;
        }
    }

    // 隐式零加法偏置退化匹配：pattern + _c("c") 形式匹配 target + 0
    if (pattern.node_->type == NodeType::kAdd) {
        SymbolicExpression p_left(pattern.node_->left);
        SymbolicExpression p_right(pattern.node_->right);
        if (p_left.node_ && p_left.node_->type == NodeType::kVariable &&
            p_left.node_->text.rfind("__wild_const__", 0) == 0) {
            std::string bid = p_left.node_->text.substr(14);
            BindingsMap backup = *bindings;
            if (match_pattern(target, p_right, input_var, bindings)) {
                if (bindings->find(bid) == bindings->end()) {
                    (*bindings)[bid] = SymbolicExpression::number(Scalar(0.0L));
                    return true;
                }
            }
            *bindings = backup;
        }
        if (p_right.node_ && p_right.node_->type == NodeType::kVariable &&
            p_right.node_->text.rfind("__wild_const__", 0) == 0) {
            std::string bid = p_right.node_->text.substr(14);
            BindingsMap backup = *bindings;
            if (match_pattern(target, p_left, input_var, bindings)) {
                if (bindings->find(bid) == bindings->end()) {
                    (*bindings)[bid] = SymbolicExpression::number(Scalar(0.0L));
                    return true;
                }
            }
            *bindings = backup;
        }
    }

    // 隐式零减法偏置退化匹配：pattern - _c("c") 匹配 target
    if (pattern.node_->type == NodeType::kSubtract) {
        SymbolicExpression p_right(pattern.node_->right);
        if (p_right.node_ && p_right.node_->type == NodeType::kVariable &&
            p_right.node_->text.rfind("__wild_const__", 0) == 0) {
            std::string bid = p_right.node_->text.substr(14);
            BindingsMap backup = *bindings;
            if (match_pattern(target, SymbolicExpression(pattern.node_->left), input_var, bindings)) {
                if (bindings->find(bid) == bindings->end()) {
                    (*bindings)[bid] = SymbolicExpression::number(Scalar(0.0L));
                    return true;
                }
            }
            *bindings = backup;
        }
    }

    if (pattern.node_->type != target.node_->type) {
        return false;
    }

    // 多参数/单参数函数完整匹配
    if (pattern.node_->type == NodeType::kFunction) {
        if (pattern.node_->text != target.node_->text) return false;
        std::vector<SymbolicExpression> p_args = get_node_args(pattern.node_.get());
        std::vector<SymbolicExpression> t_args = get_node_args(target.node_.get());
        if (p_args.size() != t_args.size()) return false;
        for (size_t i = 0; i < p_args.size(); ++i) {
            if (!match_pattern(t_args[i], p_args[i], input_var, bindings)) {
                return false;
            }
        }
        return true;
    }

    if (pattern.node_->type == NodeType::kNegate) {
        return match_pattern(SymbolicExpression(target.node_->left),
                             SymbolicExpression(pattern.node_->left),
                             input_var, bindings);
    }

    // AC 结合与对易多项匹配
    if (pattern.node_->type == NodeType::kAdd || pattern.node_->type == NodeType::kMultiply) {
        NodeType op_type = pattern.node_->type;
        std::vector<SymbolicExpression> p_terms;
        std::vector<SymbolicExpression> t_terms;
        flatten_ast_terms(op_type, pattern, &p_terms);
        flatten_ast_terms(op_type, target, &t_terms);

        if (p_terms.size() == t_terms.size()) {
            std::vector<bool> used(t_terms.size(), false);
            std::function<bool(size_t)> match_perm = [&](size_t p_idx) -> bool {
                if (p_idx == p_terms.size()) return true;
                for (size_t t_idx = 0; t_idx < t_terms.size(); ++t_idx) {
                    if (!used[t_idx]) {
                        BindingsMap backup = *bindings;
                        if (match_pattern(t_terms[t_idx], p_terms[p_idx], input_var, bindings)) {
                            used[t_idx] = true;
                            if (match_perm(p_idx + 1)) return true;
                            used[t_idx] = false;
                        }
                        *bindings = backup;
                    }
                }
                return false;
            };

            BindingsMap backup = *bindings;
            if (match_perm(0)) return true;
            *bindings = backup;
        }

        // 回退到基本二叉对易尝试
        BindingsMap backup = *bindings;
        if (match_pattern(SymbolicExpression(target.node_->left), SymbolicExpression(pattern.node_->left), input_var, bindings) &&
            match_pattern(SymbolicExpression(target.node_->right), SymbolicExpression(pattern.node_->right), input_var, bindings)) {
            return true;
        }
        *bindings = backup;
        if (match_pattern(SymbolicExpression(target.node_->left), SymbolicExpression(pattern.node_->right), input_var, bindings) &&
            match_pattern(SymbolicExpression(target.node_->right), SymbolicExpression(pattern.node_->left), input_var, bindings)) {
            return true;
        }
        *bindings = backup;
        return false;
    }

    if (pattern.node_->type == NodeType::kSubtract ||
        pattern.node_->type == NodeType::kDivide ||
        pattern.node_->type == NodeType::kPower) {
        return match_pattern(SymbolicExpression(target.node_->left), SymbolicExpression(pattern.node_->left), input_var, bindings) &&
               match_pattern(SymbolicExpression(target.node_->right), SymbolicExpression(pattern.node_->right), input_var, bindings);
    }

    return false;
}

// ============================================================================
// TransformRuleRegistry 实现
// ============================================================================

TransformRuleRegistry& TransformRuleRegistry::instance() {
    static TransformRuleRegistry instance_;
    return instance_;
}

TransformRuleRegistry::TransformRuleRegistry() {
    initialize_builtin_rules();
}

void TransformRuleRegistry::register_rule(const TransformRule& rule) {
    rules_[rule.transform_type].push_back(rule);
    std::sort(rules_[rule.transform_type].begin(),
              rules_[rule.transform_type].end(),
              [](const TransformRule& a, const TransformRule& b) {
                  return a.priority > b.priority;
              });
}

void TransformRuleRegistry::register_declarative_rule(
    const std::string& name,
    const std::string& transform_type,
    const SymbolicExpression& pattern,
    std::function<SymbolicExpression(const BindingsMap& bindings,
                                     const std::string& output_var)> template_generator,
    int priority,
    const std::string& description,
    const std::string& formula) {
    register_rule({
        .name = name,
        .transform_type = transform_type,
        .matcher = [pattern](const SymbolicExpression& expr, const std::string& var) {
            BindingsMap bindings;
            return match_pattern(expr, pattern, var, &bindings);
        },
        .transformer = [pattern, template_generator](const SymbolicExpression& expr,
                                                    const std::string& input_var,
                                                    const std::string& output_var) {
            BindingsMap bindings;
            if (match_pattern(expr, pattern, input_var, &bindings)) {
                return template_generator(bindings, output_var).simplify();
            }
            return expr;
        },
        .priority = priority,
        .description = description,
        .formula = formula
    });
}

std::vector<TransformRule> TransformRuleRegistry::get_rules(
    const std::string& transform_type) const {

    std::vector<TransformRule> result;

    auto it = rules_.find(transform_type);
    if (it != rules_.end()) {
        result = it->second;
    }

    for (const auto& rule : user_rules_) {
        if (rule.transform_type == transform_type) {
            result.push_back(rule);
        }
    }

    std::sort(result.begin(), result.end(),
              [](const TransformRule& a, const TransformRule& b) {
                  return a.priority > b.priority;
              });

    return result;
}

void TransformRuleRegistry::add_user_rule(const TransformRule& rule) {
    user_rules_.push_back(rule);
}

void TransformRuleRegistry::clear_user_rules() {
    user_rules_.clear();
}

void TransformRuleRegistry::initialize_builtin_rules() {
    if (initialized_) return;

    using namespace dsl;

    // ========================================================================
    // Laplace 变换声明式规则
    // ========================================================================

    // L{c} = c / s
    register_declarative_rule(
        "laplace_constant",
        "laplace",
        _c("c"),
        [](const BindingsMap& b, const std::string& s) {
            return make_divide(b.at("c"), SymbolicExpression::variable(s));
        },
        100,
        "L{c} = c / s",
        "c / s"
    );

    // L{t} = 1 / s^2
    register_declarative_rule(
        "laplace_t",
        "laplace",
        _var("t"),
        [](const BindingsMap&, const std::string& s) {
            return make_power(SymbolicExpression::variable(s),
                              SymbolicExpression::number(Scalar(-2.0L)));
        },
        90,
        "L{t} = 1 / s^2",
        "1 / s^2"
    );

    // L{t^k} = k! / s^(k+1)
    register_declarative_rule(
        "laplace_power",
        "laplace",
        make_power(_var("t"), _c("k")),
        [](const BindingsMap& b, const std::string& s) {
            Scalar exp_val;
            if (b.at("k").is_number(&exp_val)) {
                int power = static_cast<int>(mymath::round(exp_val.to_long_double()));
                if (mymath::abs(exp_val - Scalar(static_cast<long long>(power))) < Scalar(1e-9L) && power >= 1) {
                    Scalar factorial = Scalar(1.0L);
                    for (int i = 2; i <= power; ++i) factorial *= Scalar(static_cast<long long>(i));
                    return make_divide(
                        SymbolicExpression::number(factorial),
                        make_power(SymbolicExpression::variable(s),
                                   SymbolicExpression::number(Scalar(static_cast<long long>(power + 1)))));
                }
            }
            throw std::runtime_error("defer to algebraic engine for fractional powers");
        },
        80,
        "L{t^k} = k! / s^(k+1)",
        "k! / s^(k+1)"
    );

    // L{exp(a*t)} = 1 / (s - a)
    register_declarative_rule(
        "laplace_exp",
        "laplace",
        _exp(_c("a") * _var("t")),
        [](const BindingsMap& b, const std::string& s) {
            return make_divide(
                SymbolicExpression::number(Scalar(1.0L)),
                make_subtract(SymbolicExpression::variable(s), b.at("a")));
        },
        70,
        "L{exp(a*t)} = 1 / (s - a)",
        "1 / (s - a)"
    );

    // L{exp(a*t) * step(t)} = 1 / (s - a)
    register_declarative_rule(
        "laplace_exp_step",
        "laplace",
        _exp(_c("a") * _var("t")) * _step(_var("t")),
        [](const BindingsMap& b, const std::string& s) {
            return make_divide(
                SymbolicExpression::number(Scalar(1.0L)),
                make_subtract(SymbolicExpression::variable(s), b.at("a")));
        },
        110,
        "L{exp(a*t) * step(t)} = 1 / (s - a)",
        "1 / (s - a)"
    );

    // L{delta(t)} = 1
    register_declarative_rule(
        "laplace_delta",
        "laplace",
        _delta(_var("t")),
        [](const BindingsMap&, const std::string&) {
            return SymbolicExpression::number(Scalar(1.0L));
        },
        100,
        "L{delta(t)} = 1",
        "1"
    );

    // L{delta(t - tau)} = exp(-s * tau)
    register_declarative_rule(
        "laplace_delta_shift",
        "laplace",
        _delta(_var("t") - _c("tau")),
        [](const BindingsMap& b, const std::string& s) {
            return make_function("exp",
                make_multiply(make_negate(SymbolicExpression::variable(s)), b.at("tau")));
        },
        95,
        "L{delta(t - tau)} = exp(-s * tau)",
        "exp(-s * tau)"
    );

    // L{step(t)} = 1 / s
    register_declarative_rule(
        "laplace_step",
        "laplace",
        _step(_var("t")),
        [](const BindingsMap&, const std::string& s) {
            return make_divide(SymbolicExpression::number(Scalar(1.0L)),
                               SymbolicExpression::variable(s));
        },
        100,
        "L{step(t)} = 1 / s",
        "1 / s"
    );

    // L{sin(w*t)} = w / (s^2 + w^2)
    register_declarative_rule(
        "laplace_sin",
        "laplace",
        _sin(_c("w") * _var("t")),
        [](const BindingsMap& b, const std::string& s) {
            SymbolicExpression s_var = SymbolicExpression::variable(s);
            SymbolicExpression w = b.at("w");
            return make_divide(
                w,
                make_add(make_power(s_var, SymbolicExpression::number(Scalar(2.0L))),
                         make_power(w, SymbolicExpression::number(Scalar(2.0L)))));
        },
        60,
        "L{sin(w*t)} = w / (s^2 + w^2)",
        "w / (s^2 + w^2)"
    );

    // L{sin(w*t) * step(t)} = w / (s^2 + w^2)
    register_declarative_rule(
        "laplace_sin_step",
        "laplace",
        _sin(_c("w") * _var("t")) * _step(_var("t")),
        [](const BindingsMap& b, const std::string& s) {
            SymbolicExpression s_var = SymbolicExpression::variable(s);
            SymbolicExpression w = b.at("w");
            return make_divide(
                w,
                make_add(make_power(s_var, SymbolicExpression::number(Scalar(2.0L))),
                         make_power(w, SymbolicExpression::number(Scalar(2.0L)))));
        },
        110,
        "L{sin(w*t) * step(t)} = w / (s^2 + w^2)",
        "w / (s^2 + w^2)"
    );

    // L{cos(w*t)} = s / (s^2 + w^2)
    register_declarative_rule(
        "laplace_cos",
        "laplace",
        _cos(_c("w") * _var("t")),
        [](const BindingsMap& b, const std::string& s) {
            SymbolicExpression s_var = SymbolicExpression::variable(s);
            SymbolicExpression w = b.at("w");
            return make_divide(
                s_var,
                make_add(make_power(s_var, SymbolicExpression::number(Scalar(2.0L))),
                         make_power(w, SymbolicExpression::number(Scalar(2.0L)))));
        },
        60,
        "L{cos(w*t)} = s / (s^2 + w^2)",
        "s / (s^2 + w^2)"
    );

    // L{cos(w*t) * step(t)} = s / (s^2 + w^2)
    register_declarative_rule(
        "laplace_cos_step",
        "laplace",
        _cos(_c("w") * _var("t")) * _step(_var("t")),
        [](const BindingsMap& b, const std::string& s) {
            SymbolicExpression s_var = SymbolicExpression::variable(s);
            SymbolicExpression w = b.at("w");
            return make_divide(
                s_var,
                make_add(make_power(s_var, SymbolicExpression::number(Scalar(2.0L))),
                         make_power(w, SymbolicExpression::number(Scalar(2.0L)))));
        },
        110,
        "L{cos(w*t) * step(t)} = s / (s^2 + w^2)",
        "s / (s^2 + w^2)"
    );

    // L{sinh(a*t)} = a / (s^2 - a^2)
    register_declarative_rule(
        "laplace_sinh",
        "laplace",
        _sinh(_c("a") * _var("t")),
        [](const BindingsMap& b, const std::string& s) {
            SymbolicExpression s_var = SymbolicExpression::variable(s);
            SymbolicExpression a = b.at("a");
            return make_divide(
                a,
                make_subtract(make_power(s_var, SymbolicExpression::number(Scalar(2.0L))),
                              make_power(a, SymbolicExpression::number(Scalar(2.0L)))));
        },
        60,
        "L{sinh(a*t)} = a / (s^2 - a^2)",
        "a / (s^2 - a^2)"
    );

    // L{sinh(a*t) * step(t)} = a / (s^2 - a^2)
    register_declarative_rule(
        "laplace_sinh_step",
        "laplace",
        _sinh(_c("a") * _var("t")) * _step(_var("t")),
        [](const BindingsMap& b, const std::string& s) {
            SymbolicExpression s_var = SymbolicExpression::variable(s);
            SymbolicExpression a = b.at("a");
            return make_divide(
                a,
                make_subtract(make_power(s_var, SymbolicExpression::number(Scalar(2.0L))),
                              make_power(a, SymbolicExpression::number(Scalar(2.0L)))));
        },
        110,
        "L{sinh(a*t) * step(t)} = a / (s^2 - a^2)",
        "a / (s^2 - a^2)"
    );

    // L{cosh(a*t)} = s / (s^2 - a^2)
    register_declarative_rule(
        "laplace_cosh",
        "laplace",
        _cosh(_c("a") * _var("t")),
        [](const BindingsMap& b, const std::string& s) {
            SymbolicExpression s_var = SymbolicExpression::variable(s);
            SymbolicExpression a = b.at("a");
            return make_divide(
                s_var,
                make_subtract(make_power(s_var, SymbolicExpression::number(Scalar(2.0L))),
                              make_power(a, SymbolicExpression::number(Scalar(2.0L)))));
        },
        60,
        "L{cosh(a*t)} = s / (s^2 - a^2)",
        "s / (s^2 - a^2)"
    );

    // L{cosh(a*t) * step(t)} = s / (s^2 - a^2)
    register_declarative_rule(
        "laplace_cosh_step",
        "laplace",
        _cosh(_c("a") * _var("t")) * _step(_var("t")),
        [](const BindingsMap& b, const std::string& s) {
            SymbolicExpression s_var = SymbolicExpression::variable(s);
            SymbolicExpression a = b.at("a");
            return make_divide(
                s_var,
                make_subtract(make_power(s_var, SymbolicExpression::number(Scalar(2.0L))),
                              make_power(a, SymbolicExpression::number(Scalar(2.0L)))));
        },
        110,
        "L{cosh(a*t) * step(t)} = s / (s^2 - a^2)",
        "s / (s^2 - a^2)"
    );

    // L{t * exp(a*t)} = 1 / (s - a)^2
    register_declarative_rule(
        "laplace_t_exp",
        "laplace",
        _var("t") * _exp(_c("a") * _var("t")),
        [](const BindingsMap& b, const std::string& s) {
            SymbolicExpression s_var = SymbolicExpression::variable(s);
            SymbolicExpression a = b.at("a");
            return make_divide(
                SymbolicExpression::number(Scalar(1.0L)),
                make_power(make_subtract(s_var, a), SymbolicExpression::number(Scalar(2.0L))));
        },
        80,
        "L{t * exp(a*t)} = 1 / (s - a)^2",
        "1 / (s - a)^2"
    );

    // L{t * exp(a*t) * step(t)} = 1 / (s - a)^2
    register_declarative_rule(
        "laplace_t_exp_step",
        "laplace",
        _var("t") * _exp(_c("a") * _var("t")) * _step(_var("t")),
        [](const BindingsMap& b, const std::string& s) {
            SymbolicExpression s_var = SymbolicExpression::variable(s);
            SymbolicExpression a = b.at("a");
            return make_divide(
                SymbolicExpression::number(Scalar(1.0L)),
                make_power(make_subtract(s_var, a), SymbolicExpression::number(Scalar(2.0L))));
        },
        120,
        "L{t * exp(a*t) * step(t)} = 1 / (s - a)^2",
        "1 / (s - a)^2"
    );

    // ========================================================================
    // 逆 Laplace 变换声明式规则
    // ========================================================================

    // L^-1{1 / s} = step(t)
    register_declarative_rule(
        "ilaplace_one_over_s",
        "ilaplace",
        make_divide(SymbolicExpression::number(Scalar(1.0L)), _var("s")),
        [](const BindingsMap&, const std::string& t) {
            return make_function("step", SymbolicExpression::variable(t));
        },
        100,
        "L^-1{1 / s} = step(t)",
        "step(t)"
    );

    // L^-1{1 / (s - a)} = exp(a*t) * step(t)
    register_declarative_rule(
        "ilaplace_exp",
        "ilaplace",
        make_divide(SymbolicExpression::number(Scalar(1.0L)), make_subtract(_var("s"), _c("a"))),
        [](const BindingsMap& b, const std::string& t) {
            SymbolicExpression t_var = SymbolicExpression::variable(t);
            return make_multiply(
                make_function("exp", make_multiply(b.at("a"), t_var)),
                make_function("step", t_var));
        },
        90,
        "L^-1{1 / (s - a)} = exp(a*t) * step(t)",
        "exp(a*t) * step(t)"
    );

    // L^-1{1 / (s + a)} = exp(-a*t) * step(t)
    register_declarative_rule(
        "ilaplace_exp_plus",
        "ilaplace",
        make_divide(SymbolicExpression::number(Scalar(1.0L)), make_add(_var("s"), _c("a"))),
        [](const BindingsMap& b, const std::string& t) {
            SymbolicExpression t_var = SymbolicExpression::variable(t);
            return make_multiply(
                make_function("exp", make_multiply(make_negate(b.at("a")), t_var)),
                make_function("step", t_var));
        },
        90,
        "L^-1{1 / (s + a)} = exp(-a*t) * step(t)",
        "exp(-a*t) * step(t)"
    );

    // L^-1{1 / (s - a)^2} = t * exp(a*t) * step(t)
    register_declarative_rule(
        "ilaplace_t_exp",
        "ilaplace",
        make_divide(SymbolicExpression::number(Scalar(1.0L)),
                    make_power(make_subtract(_var("s"), _c("a")), SymbolicExpression::number(Scalar(2.0L)))),
        [](const BindingsMap& b, const std::string& t) {
            SymbolicExpression t_var = SymbolicExpression::variable(t);
            return make_multiply(
                t_var,
                make_multiply(
                    make_function("exp", make_multiply(b.at("a"), t_var)),
                    make_function("step", t_var)));
        },
        95,
        "L^-1{1 / (s - a)^2} = t * exp(a*t) * step(t)",
        "t * exp(a*t) * step(t)"
    );

    // ========================================================================
    // Fourier 变换声明式规则
    // ========================================================================

    // F{c} = 2*pi*c*delta(w)
    register_declarative_rule(
        "fourier_constant",
        "fourier",
        _c("c"),
        [](const BindingsMap& b, const std::string& w) {
            return make_multiply(
                b.at("c"),
                make_multiply(
                    make_multiply(SymbolicExpression::number(Scalar(2.0L)),
                                  SymbolicExpression::variable("pi")),
                    make_function("delta", SymbolicExpression::variable(w))));
        },
        100,
        "F{c} = 2*pi*c*delta(w)",
        "2*pi*c*delta(w)"
    );

    // F{delta(t)} = 1
    register_declarative_rule(
        "fourier_delta",
        "fourier",
        _delta(_var("t")),
        [](const BindingsMap&, const std::string&) {
            return SymbolicExpression::number(Scalar(1.0L));
        },
        100,
        "F{delta(t)} = 1",
        "1"
    );

    // F{delta(t - tau)} = exp(-i * w * tau)
    register_declarative_rule(
        "fourier_delta_shift",
        "fourier",
        _delta(_var("t") - _c("tau")),
        [](const BindingsMap& b, const std::string& w) {
            return make_function(
                "exp",
                make_multiply(
                    make_multiply(SymbolicExpression::variable("i"),
                                  make_negate(SymbolicExpression::variable(w))),
                    b.at("tau")));
        },
        95,
        "F{delta(t - tau)} = exp(-i * w * tau)",
        "exp(-i * w * tau)"
    );

    // F{sgn(t)} = 2 / (i * w)
    register_declarative_rule(
        "fourier_sgn",
        "fourier",
        _sgn(_var("t")),
        [](const BindingsMap&, const std::string& w) {
            return make_divide(
                SymbolicExpression::number(Scalar(2.0L)),
                make_multiply(SymbolicExpression::variable("i"),
                              SymbolicExpression::variable(w)));
        },
        90,
        "F{sgn(t)} = 2 / (i * w)",
        "2 / (i * w)"
    );

    // F{step(t)} = pi * delta(w) + 1 / (i * w)
    register_declarative_rule(
        "fourier_step",
        "fourier",
        _step(_var("t")),
        [](const BindingsMap&, const std::string& w) {
            SymbolicExpression w_var = SymbolicExpression::variable(w);
            return make_add(
                make_multiply(SymbolicExpression::variable("pi"), make_function("delta", w_var)),
                make_divide(SymbolicExpression::number(Scalar(1.0L)),
                            make_multiply(SymbolicExpression::variable("i"), w_var)));
        },
        80,
        "F{step(t)} = pi * delta(w) + 1 / (i * w)",
        "pi * delta(w) + 1 / (i * w)"
    );

    // ========================================================================
    // 逆 Fourier 变换声明式规则
    // ========================================================================

    // IF{delta(w)} = 1 / (2 * pi)
    register_declarative_rule(
        "ifourier_delta",
        "ifourier",
        _delta(_var("w")),
        [](const BindingsMap&, const std::string&) {
            return make_divide(
                SymbolicExpression::number(Scalar(1.0L)),
                make_multiply(SymbolicExpression::number(Scalar(2.0L)), SymbolicExpression::variable("pi")));
        },
        100,
        "IF{delta(w)} = 1 / (2 * pi)",
        "1 / (2 * pi)"
    );

    // IF{delta(w - w0)} = 1 / (2 * pi) * exp(i * w0 * t)
    register_declarative_rule(
        "ifourier_delta_shift",
        "ifourier",
        _delta(_var("w") - _c("w0")),
        [](const BindingsMap& b, const std::string& t) {
            SymbolicExpression t_var = SymbolicExpression::variable(t);
            SymbolicExpression w0 = b.at("w0");
            SymbolicExpression one_over_2pi = make_divide(
                SymbolicExpression::number(Scalar(1.0L)),
                make_multiply(SymbolicExpression::number(Scalar(2.0L)), SymbolicExpression::variable("pi")));
            SymbolicExpression exp_term = make_function(
                "exp",
                make_multiply(SymbolicExpression::variable("i"), make_multiply(w0, t_var)));
            return make_multiply(one_over_2pi, exp_term);
        },
        95,
        "IF{delta(w - w0)} = 1 / (2 * pi) * exp(i * w0 * t)",
        "1 / (2 * pi) * exp(i * w0 * t)"
    );

    // ========================================================================
    // Z 变换声明式规则
    // ========================================================================

    // Z{c} = c*z / (z - 1)
    register_declarative_rule(
        "z_constant",
        "z",
        _c("c"),
        [](const BindingsMap& b, const std::string& z) {
            SymbolicExpression z_var = SymbolicExpression::variable(z);
            return make_divide(
                make_multiply(b.at("c"), z_var),
                make_subtract(z_var, SymbolicExpression::number(Scalar(1.0L))));
        },
        100,
        "Z{c} = c*z / (z - 1)",
        "c*z / (z - 1)"
    );

    // Z{delta(n)} = 1
    register_declarative_rule(
        "z_delta",
        "z",
        _delta(_var("n")),
        [](const BindingsMap&, const std::string&) {
            return SymbolicExpression::number(Scalar(1.0L));
        },
        100,
        "Z{delta(n)} = 1",
        "1"
    );

    // Z{delta(n - k)} = z^(-k)
    register_declarative_rule(
        "z_delta_shift",
        "z",
        _delta(_var("n") - _c("k")),
        [](const BindingsMap& b, const std::string& z) {
            return make_power(SymbolicExpression::variable(z), make_negate(b.at("k")));
        },
        110,
        "Z{delta(n - k)} = z^(-k)",
        "z^(-k)"
    );

    // Z{step(n)} = z / (z - 1)
    register_declarative_rule(
        "z_step",
        "z",
        _step(_var("n")),
        [](const BindingsMap&, const std::string& z) {
            SymbolicExpression z_var = SymbolicExpression::variable(z);
            return make_divide(
                z_var,
                make_subtract(z_var, SymbolicExpression::number(Scalar(1.0L))));
        },
        90,
        "Z{step(n)} = z / (z - 1)",
        "z / (z - 1)"
    );

    // Z{a^n * step(n)} = z / (z - a)
    register_declarative_rule(
        "z_geometric_step",
        "z",
        make_power(_c("a"), _var("n")) * _step(_var("n")),
        [](const BindingsMap& b, const std::string& z) {
            SymbolicExpression z_var = SymbolicExpression::variable(z);
            return make_divide(z_var, make_subtract(z_var, b.at("a")));
        },
        110,
        "Z{a^n * step(n)} = z / (z - a)",
        "z / (z - a)"
    );

    // Z{n * step(n)} = z / (z - 1)^2
    register_declarative_rule(
        "z_ramp_step",
        "z",
        _var("n") * _step(_var("n")),
        [](const BindingsMap&, const std::string& z) {
            SymbolicExpression z_var = SymbolicExpression::variable(z);
            return make_divide(
                z_var,
                make_power(make_subtract(z_var, SymbolicExpression::number(Scalar(1.0L))),
                           SymbolicExpression::number(Scalar(2.0L))));
        },
        105,
        "Z{n * step(n)} = z / (z - 1)^2",
        "z / (z - 1)^2"
    );

    // Z{n * a^n * step(n)} = a*z / (z - a)^2
    register_declarative_rule(
        "z_geom_ramp_step",
        "z",
        _var("n") * make_power(_c("a"), _var("n")) * _step(_var("n")),
        [](const BindingsMap& b, const std::string& z) {
            SymbolicExpression z_var = SymbolicExpression::variable(z);
            SymbolicExpression a = b.at("a");
            return make_divide(
                make_multiply(a, z_var),
                make_power(make_subtract(z_var, a), SymbolicExpression::number(Scalar(2.0L))));
        },
        115,
        "Z{n * a^n * step(n)} = a*z / (z - a)^2",
        "a*z / (z - a)^2"
    );

    // ========================================================================
    // 逆 Z 变换声明式规则
    // ========================================================================

    // IZ{1} = delta(n)
    register_declarative_rule(
        "iz_one",
        "iz",
        SymbolicExpression::number(Scalar(1.0L)),
        [](const BindingsMap&, const std::string& n) {
            return make_function("delta", SymbolicExpression::variable(n));
        },
        100,
        "IZ{1} = delta(n)",
        "delta(n)"
    );

    // IZ{z / (z - 1)} = step(n)
    register_declarative_rule(
        "iz_step",
        "iz",
        make_divide(_var("z"), make_subtract(_var("z"), SymbolicExpression::number(Scalar(1.0L)))),
        [](const BindingsMap&, const std::string& n) {
            return make_function("step", SymbolicExpression::variable(n));
        },
        100,
        "IZ{z / (z - 1)} = step(n)",
        "step(n)"
    );

    // IZ{z / (z - a)} = a^n * step(n)
    register_declarative_rule(
        "iz_geom_step",
        "iz",
        make_divide(_var("z"), make_subtract(_var("z"), _c("a"))),
        [](const BindingsMap& b, const std::string& n) {
            SymbolicExpression n_var = SymbolicExpression::variable(n);
            return make_multiply(
                make_power(b.at("a"), n_var),
                make_function("step", n_var));
        },
        90,
        "IZ{z / (z - a)} = a^n * step(n)",
        "a^n * step(n)"
    );

    // IZ{z / (z - 1)^2} = n * step(n)
    register_declarative_rule(
        "iz_ramp_step",
        "iz",
        make_divide(_var("z"),
                    make_power(make_subtract(_var("z"), SymbolicExpression::number(Scalar(1.0L))),
                               SymbolicExpression::number(Scalar(2.0L)))),
        [](const BindingsMap&, const std::string& n) {
            SymbolicExpression n_var = SymbolicExpression::variable(n);
            return make_multiply(n_var, make_function("step", n_var));
        },
        95,
        "IZ{z / (z - 1)^2} = n * step(n)",
        "n * step(n)"
    );

    initialized_ = true;
}

// ============================================================================
// 辅助函数实现
// ============================================================================

bool is_constant_expression(const SymbolicExpression& expr, const std::string& var) {
    // Check if expression contains the variable
    std::set<std::string> vars;
    // Simple check: if it's a number or doesn't contain the variable
    Scalar val;
    if (expr.is_number(&val)) return true;

    // Check if it's a variable different from var
    if (expr.node_->type == NodeType::kVariable && expr.node_->text != var) return true;
    if (expr.node_->type == NodeType::kVariable && expr.node_->text == var) return false;

    // For other cases, check recursively
    if (expr.node_->type == NodeType::kFunction) {
        return is_constant_expression(SymbolicExpression(expr.node_->left), var);
    }
    if (expr.node_->type == NodeType::kNegate) {
        return is_constant_expression(SymbolicExpression(expr.node_->left), var);
    }
    if (expr.node_->type == NodeType::kAdd || expr.node_->type == NodeType::kSubtract ||
        expr.node_->type == NodeType::kMultiply || expr.node_->type == NodeType::kDivide ||
        expr.node_->type == NodeType::kPower) {
        return is_constant_expression(SymbolicExpression(expr.node_->left), var) &&
               is_constant_expression(SymbolicExpression(expr.node_->right), var);
    }

    return true;
}

bool is_variable_expression(const SymbolicExpression& expr, const std::string& var) {
    if (expr.node_->type == NodeType::kVariable && expr.node_->text == var) {
        return true;
    }
    return false;
}

bool is_power_of_variable(const SymbolicExpression& expr, const std::string& var, int* power) {
    if (power) *power = 0;

    if (is_variable_expression(expr, var)) {
        if (power) *power = 1;
        return true;
    }

    if (expr.node_->type == NodeType::kPower) {
        SymbolicExpression base(expr.node_->left);
        SymbolicExpression exp(expr.node_->right);

        if (is_variable_expression(base, var)) {
            Scalar exp_val;
            if (exp.is_number(&exp_val)) {
                int n = static_cast<int>(mymath::round(exp_val.to_long_double()));
                if (mymath::abs(exp_val - Scalar(static_cast<long long>(n))) < Scalar(1e-9L) && n > 0) {
                    if (power) *power = n;
                    return true;
                }
            }
        }
    }

    return false;
}

bool is_exponential_form(const SymbolicExpression& expr, const std::string& var,
                         SymbolicExpression* coefficient, SymbolicExpression* constant) {
    if (expr.node_->type != NodeType::kFunction || expr.node_->text != "exp") {
        return false;
    }

    SymbolicExpression arg = SymbolicExpression(expr.node_->left);

    // Check for a*var + b form
    if (arg.node_->type == NodeType::kMultiply) {
        SymbolicExpression left(arg.node_->left);
        SymbolicExpression right(arg.node_->right);

        if (is_variable_expression(left, var) && is_constant_expression(right, var)) {
            if (coefficient) *coefficient = right;
            if (constant) *constant = SymbolicExpression::number(Scalar(0));
            return true;
        }
        if (is_variable_expression(right, var) && is_constant_expression(left, var)) {
            if (coefficient) *coefficient = left;
            if (constant) *constant = SymbolicExpression::number(Scalar(0));
            return true;
        }
    }

    if (is_variable_expression(arg, var)) {
        if (coefficient) *coefficient = SymbolicExpression::number(Scalar(1));
        if (constant) *constant = SymbolicExpression::number(Scalar(0));
        return true;
    }

    return false;
}

bool is_trig_form(const SymbolicExpression& expr, const std::string& func_name,
                  const std::string& var, SymbolicExpression* coefficient,
                  SymbolicExpression* constant) {
    if (expr.node_->type != NodeType::kFunction || expr.node_->text != func_name) {
        return false;
    }

    SymbolicExpression arg = SymbolicExpression(expr.node_->left);

    // Check for a*var form
    if (arg.node_->type == NodeType::kMultiply) {
        SymbolicExpression left(arg.node_->left);
        SymbolicExpression right(arg.node_->right);

        if (is_variable_expression(left, var) && is_constant_expression(right, var)) {
            if (coefficient) *coefficient = right;
            if (constant) *constant = SymbolicExpression::number(Scalar(0));
            return true;
        }
        if (is_variable_expression(right, var) && is_constant_expression(left, var)) {
            if (coefficient) *coefficient = left;
            if (constant) *constant = SymbolicExpression::number(Scalar(0));
            return true;
        }
    }

    if (is_variable_expression(arg, var)) {
        if (coefficient) *coefficient = SymbolicExpression::number(Scalar(1));
        if (constant) *constant = SymbolicExpression::number(Scalar(0));
        return true;
    }

    return false;
}

} // namespace transform_rules
