// ============================================================================
// script_utils.cpp - 脚本运行时辅助函数实现
// ============================================================================

#include "script_runtime_internal.h"
#include "execution/resolver/variable_resolver.h"
#include "parser/infra/parser_utils.h"
#include "core/services/string_utils.h"
#include "core/services/format_utils.h"
#include "core/services/core_manager_interfaces.h"
#include "math/helpers/integer_helpers.h"
#include "mymath.h"
#include <algorithm>
#include <cctype>

VariableResolver visible_variables(IExecutionContext* ctx) {
    return ctx->variables().create_resolver();
}

bool has_visible_script_function(IExecutionContext* ctx, const std::string& name) {
    if (ctx->functions().get_native_functions()->count(name) > 0) return true;
    return ctx->functions().get_script(name) != nullptr;
}

void assign_visible_variable(IExecutionContext* ctx,
                             const std::string& name,
                             const StoredValue& value) {
    ctx->variables().assign_visible(name, value);
}

bool truthy_value(const StoredValue& value) {
    if (value.is_matrix) {
        throw std::runtime_error("matrix values cannot be used as script conditions");
    }
    if (value.is_complex) {
        throw std::runtime_error("complex values cannot be used as script conditions");
    }
    if (value.is_string) {
        return !value.string_value.empty();
    }
    if (value.is_list) {
        return value.list_value && !value.list_value->empty();
    }
    if (value.is_dict) {
        return value.dict_value && !value.dict_value->empty();
    }
    return !mymath::is_near_zero(value.exact
                                     ? rational_to_double(value.rational)
                                     : value.decimal,
                                 1e-10);
}

Scalar evaluate_scalar( IExecutionContext* ctx, const CommandASTNode& ast, const char* context) {
    StoredValue val = evaluate_command_ast_to_value(ctx, ast, false);
    if (val.is_matrix || val.is_complex || val.is_string) {
        throw std::runtime_error(std::string(context) + " must be a scalar");
    }
    return val.exact ? rational_to_double(val.rational) : val.decimal;
}

bool is_wrapped_by(std::string_view text, char open, char close) {
    return parser_utils::is_wrapped_by(text, open, close);
}

std::vector<std::string> split_script_top_level(std::string_view text, char delimiter) {
    auto parts = parser_utils::split_top_level(text, delimiter);
    for (auto& p : parts) p = utils::trim_copy(p);
    return parts;
}

std::size_t find_top_level_word(std::string_view text, std::string_view word, std::size_t start_at) {
    int paren = 0, bracket = 0, brace = 0;
    bool in_string = false, escaping = false;
    for (std::size_t i = start_at; i + word.size() <= text.size(); ++i) {
        const char ch = text[i];
        if (in_string) {
            if (escaping) escaping = false;
            else if (ch == '\\') escaping = true;
            else if (ch == '"') in_string = false;
            continue;
        }
        if (ch == '"') { in_string = true; continue; }
        if (ch == '(') ++paren;
        else if (ch == ')' && paren > 0) --paren;
        else if (ch == '[') ++bracket;
        else if (ch == ']' && bracket > 0) --bracket;
        else if (ch == '{') ++brace;
        else if (ch == '}' && brace > 0) --brace;
        if (paren == 0 && bracket == 0 && brace == 0 &&
            text.substr(i, word.size()) == word) {
            const bool left_ok = i == 0 ||
                !(std::isalnum(static_cast<unsigned char>(text[i - 1])) || text[i - 1] == '_');
            const std::size_t after = i + word.size();
            const bool right_ok = after >= text.size() ||
                !(std::isalnum(static_cast<unsigned char>(text[after])) || text[after] == '_');
            if (left_ok && right_ok) return i;
        }
    }
    return std::string::npos;
}

std::size_t find_top_level_assignment(std::string_view text) {
    int paren = 0, bracket = 0, brace = 0;
    bool in_string = false, escaping = false;
    for (std::size_t i = 0; i < text.size(); ++i) {
        const char ch = text[i];
        if (in_string) {
            if (escaping) escaping = false;
            else if (ch == '\\') escaping = true;
            else if (ch == '"') in_string = false;
            continue;
        }
        if (ch == '"') { in_string = true; continue; }
        if (ch == '(') ++paren;
        else if (ch == ')' && paren > 0) --paren;
        else if (ch == '[') ++bracket;
        else if (ch == ']' && bracket > 0) --bracket;
        else if (ch == '{') ++brace;
        else if (ch == '}' && brace > 0) --brace;
        else if (ch == '=' && paren == 0 && bracket == 0 && brace == 0) {
            const char prev = i > 0 ? text[i - 1] : '\0';
            const char next = i + 1 < text.size() ? text[i + 1] : '\0';
            if (prev != '=' && prev != '!' && prev != '<' && prev != '>' && next != '=') return i;
        }
    }
    return std::string::npos;
}

long long stored_to_index(const StoredValue& value, const char* context) {
    if (value.is_matrix || value.is_complex || value.is_string || value.is_list || value.is_dict) {
        throw std::runtime_error(std::string(context) + " must be an integer");
    }
    const Scalar scalar = value.exact ? rational_to_double(value.rational) : value.decimal;
    if (!is_integer_double(scalar)) {
        throw std::runtime_error(std::string(context) + " must be an integer");
    }
    return round_to_long_long(scalar);
}

std::string stored_to_key(const StoredValue& value) {
    if (value.is_string) return value.string_value;
    if (value.is_matrix || value.is_complex || value.is_list || value.is_dict) {
        throw std::runtime_error("dictionary key must be a string or scalar");
    }
    return value.exact ? value.rational.to_string() : format_stored_value(value, false);
}

bool parse_index_expression(std::string_view expression, std::string* base, std::string* index) {
    const std::string text = utils::trim_copy(expression);
    if (text.empty() || text.back() != ']') return false;
    int paren = 0, bracket = 0, brace = 0;
    bool in_string = false, escaping = false;
    std::size_t open_pos = std::string::npos;
    for (std::size_t i = 0; i < text.size(); ++i) {
        const char ch = text[i];
        if (in_string) {
            if (escaping) escaping = false;
            else if (ch == '\\') escaping = true;
            else if (ch == '"') in_string = false;
            continue;
        }
        if (ch == '"') { in_string = true; continue; }
        if (ch == '(') ++paren;
        else if (ch == ')' && paren > 0) --paren;
        else if (ch == '{') ++brace;
        else if (ch == '}' && brace > 0) --brace;
        else if (ch == '[') {
            if (paren == 0 && bracket == 0 && brace == 0) {
                open_pos = i;
            }
            ++bracket;
        } else if (ch == ']' && bracket > 0) {
            --bracket;
            if (bracket == 0 && i + 1 != text.size()) {
                open_pos = std::string::npos;
            }
        }
    }
    if (open_pos == std::string::npos) return false;
    *base = utils::trim_copy(text.substr(0, open_pos));
    *index = utils::trim_copy(text.substr(open_pos + 1, text.size() - open_pos - 2));
    return !base->empty();
}

bool has_top_level_semicolon(std::string_view text) {
    return split_script_top_level(text, ';').size() > 1;
}

StoredValue make_list_value(std::vector<StoredValue> values) {
    StoredValue stored;
    stored.is_list = true;
    stored.list_value = std::make_shared<std::vector<StoredValue>>(std::move(values));
    return stored;
}

StoredValue make_dict_value(std::map<std::string, StoredValue> values) {
    StoredValue stored;
    stored.is_dict = true;
    stored.dict_value = std::make_shared<std::map<std::string, StoredValue>>(std::move(values));
    return stored;
}
