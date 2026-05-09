#include "symbolic/modules/commands/symbolic_commands_internal.h"
#include "symbolic/base/assumptions.h"
#include "symbolic/algebra/groebner/groebner_basis.h"
#include "core/string_utils.h"
#include <vector>

namespace symbolic_commands {

bool handle_algebra_commands(const SymbolicCommandContext& ctx,
                            const std::string& command,
                            const std::string& inside,
                            const std::vector<std::string>& arguments,
                            std::string* output) {
    if (command == "groebner") {
        auto parts = ctx.parse_symbolic_expression_list(inside);
        if (parts.size() < 2) throw std::runtime_error("groebner expects list of polynomials and list of variables");
        // Simplified handling for now as in original code
        return false;
    }

    if (command == "simplify" || command == "latex") {
        std::string var;
        SymbolicExpression expr;
        ctx.resolve_symbolic(inside, false, &var, &expr);
        if (command == "latex") {
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
