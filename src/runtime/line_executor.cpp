#include "line_executor.h"

#include "evaluator.h"
#include "parser.h"

#include <cctype>
#include <sstream>
#include <stdexcept>

namespace runtime {
namespace {

std::string trim(const std::string& text) {
    std::size_t begin = 0;
    while (begin < text.size() &&
           std::isspace(static_cast<unsigned char>(text[begin]))) {
        ++begin;
    }
    std::size_t end = text.size();
    while (end > begin &&
           std::isspace(static_cast<unsigned char>(text[end - 1]))) {
        --end;
    }
    return text.substr(begin, end - begin);
}

bool is_identifier(const std::string& text) {
    if (text.empty() ||
        !(std::isalpha(static_cast<unsigned char>(text[0])) || text[0] == '_')) {
        return false;
    }
    for (char ch : text) {
        if (!(std::isalnum(static_cast<unsigned char>(ch)) || ch == '_')) {
            return false;
        }
    }
    return true;
}

bool split_top_level_assignment(const std::string& line,
                                std::string* lhs,
                                std::string* rhs) {
    int depth = 0;
    for (std::size_t i = 0; i < line.size(); ++i) {
        const char ch = line[i];
        if (ch == '(' || ch == '[' || ch == '{') {
            ++depth;
        } else if (ch == ')' || ch == ']' || ch == '}') {
            --depth;
            if (depth < 0) {
                throw std::runtime_error("unbalanced closing delimiter");
            }
        } else if (ch == '=' && depth == 0) {
            const bool comparison =
                (i > 0 && (line[i - 1] == '<' || line[i - 1] == '>' ||
                           line[i - 1] == '!' || line[i - 1] == '=')) ||
                (i + 1 < line.size() && line[i + 1] == '=');
            if (comparison) {
                continue;
            }
            *lhs = trim(line.substr(0, i));
            *rhs = trim(line.substr(i + 1));
            return true;
        }
    }
    if (depth != 0) {
        throw std::runtime_error("unbalanced opening delimiter");
    }
    return false;
}

}  // namespace

std::string LineResult::to_string() const {
    if (!assigned) {
        return value.to_string();
    }
    return name + " = " + value.to_string();
}

LineResult evaluate_line(const std::string& line, Environment& env) {
    const std::string trimmed = trim(line);
    if (trimmed.empty()) {
        throw std::runtime_error("empty input");
    }

    std::string lhs;
    std::string rhs;
    if (split_top_level_assignment(trimmed, &lhs, &rhs)) {
        if (!is_identifier(lhs)) {
            throw std::runtime_error("invalid assignment target: " + lhs);
        }
        if (rhs.empty()) {
            throw std::runtime_error("assignment requires a value");
        }
        const expression::Expr expr = expression::parse_expression(rhs);
        const Value value = expression::evaluate(expr, env);
        env.set(lhs, value, expr);
        return {true, lhs, value};
    }

    const expression::Expr expr = expression::parse_expression(trimmed);
    return {false, "", expression::evaluate(expr, env)};
}

std::string list_variables(const Environment& env) {
    if (env.variables().empty()) {
        return "No variables defined.";
    }
    std::ostringstream out;
    bool first = true;
    for (const auto& [name, binding] : env.variables()) {
        if (!first) {
            out << '\n';
        }
        first = false;
        out << name << " = " << binding.value.to_string();
    }
    return out.str();
}

}  // namespace runtime
