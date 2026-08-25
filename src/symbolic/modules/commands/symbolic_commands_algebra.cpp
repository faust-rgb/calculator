/**
 * @file symbolic_commands_algebra.cpp
 * @brief 代数运算命令实现
 *
 * 本文件实现了符号代数相关的命令处理：
 * - simplify: 表达式简化
 * - expand: 表达式展开
 * - factor: 因式分解
 * - latex: 转换为 LaTeX 格式
 * - groebner: Gröbner 基计算
 *
 * 这些命令用于符号表达式的代数变换。
 */

#include "symbolic/modules/commands/symbolic_commands_internal.h"
#include "symbolic/base/assumptions.h"
#include "symbolic/algebra/groebner/groebner_basis.h"
#include "symbolic/solver/symbolic_solver.h"
#include "core/services/string_utils.h"
#include <vector>
#include <set>
#include <sstream>

namespace symbolic_commands {

bool handle_algebra_commands(const SymbolicCommandContext& ctx,
                            const std::string& command,
                            const std::string& inside,
                            const std::vector<std::string>& arguments,
                            std::string* output) {
    if (command == "groebner") {
        if (arguments.empty()) throw std::runtime_error("groebner expects polynomials and variables");
        std::vector<SymbolicExpression> polys;
        std::vector<std::string> vars;

        if (arguments.size() >= 2 && !arguments[0].empty() && arguments[0].front() == '[' &&
            !arguments.back().empty() && arguments.back().front() == '[') {
            symbolic_solver::parse_equation_system(arguments[0], &polys);
            std::vector<SymbolicExpression> var_exprs;
            symbolic_solver::parse_equation_system(arguments.back(), &var_exprs);
            for (const auto& ve : var_exprs) {
                if (ve.node_type() == NodeType::kVariable) {
                    vars.push_back(ve.node_text());
                }
            }
        } else {
            auto parts = ctx.parse_symbolic_expression_list(inside);
            for (const auto& p : parts) {
                if (p.node_type() == NodeType::kVariable) {
                    vars.push_back(p.node_text());
                } else {
                    polys.push_back(p);
                }
            }
            if (vars.empty()) {
                std::set<std::string> var_set;
                for (const auto& p : polys) {
                    for (const auto& iv : p.identifier_variables()) var_set.insert(iv);
                }
                vars.assign(var_set.begin(), var_set.end());
            }
        }

        auto basis = symbolic_groebner::compute_groebner_basis(polys, vars);
        std::ostringstream oss;
        oss << "[";
        for (size_t i = 0; i < basis.size(); ++i) {
            if (i > 0) oss << ", ";
            oss << basis[i].simplify().to_string();
        }
        oss << "]";
        *output = oss.str();
        return true;
    }

    if (command == "simplify" || command == "latex" || command == "to_latex") {
        std::string var;
        SymbolicExpression expr;
        ctx.resolve_symbolic(inside, false, &var, &expr);
        if (command == "latex" || command == "to_latex") {
            *output = expr.to_latex();
        } else {
            *output = expr.simplify().to_string();
        }
        return true;
    }

    if (command == "expand" || command == "cse") {
        std::string var;
        SymbolicExpression expr;
        ctx.resolve_symbolic(arguments[0], false, &var, &expr);
        if (command == "expand") {
            *output = expr.expand().simplify().to_string();
        } else if (command == "cse") {
            std::vector<std::pair<std::string, SymbolicExpression>> bindings;
            SymbolicExpression rewritten = expr.cse_rewrite(bindings);
            std::ostringstream out;
            for (const auto& binding : bindings) {
                out << binding.first << " = " << binding.second.simplify().to_string() << "\n";
            }
            out << rewritten.simplify().to_string();
            *output = out.str();
        }
        return true;
    }

    return false;
}

} // namespace symbolic_commands
