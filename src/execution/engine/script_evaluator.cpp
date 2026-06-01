#include "core/services/core_manager_interfaces.h"
#include "types/function.h"
#include "core/services/core_manager_interfaces.h"
// ============================================================================
// script_evaluator.cpp - 脚本值求值与表达式处理实现
// ============================================================================

#include "script_runtime_internal.h"
#include "execution/resolver/variable_resolver.h"
#include "parser/grammars/unified_expression_parser.h"
#include "parser/grammars/symbolic_render_parser.h"
#include "core/services/string_utils.h"
#include "core/services/format_utils.h"
#include "execution/engine/inline_expander.h"
#include "math/helpers/integer_helpers.h"
#include "mymath.h"
#include "matrix/matrix.h"
#include <map>
#include <stdexcept>

/**
 * @brief 求值 range() 函数调用
 */
StoredValue evaluate_range_list(
                                IExecutionContext* ctx,
                                const std::string& args_text,
                                bool exact_mode) {
    std::vector<std::string> args;
    if (!utils::trim_copy(args_text).empty()) args = split_script_top_level(args_text, ',');
    if (args.empty() || args.size() > 3) throw std::runtime_error("range expects 1-3 arguments");
    auto eval_scalar = [&](const std::string& arg) {
        StoredValue value = evaluate_script_value_expression(ctx, arg, exact_mode);
        return (stored_to_index(value, "range argument"));
    };
    Scalar start = 0.0L, stop = 0.0L, step = 1.0L;
    if (args.size() == 1) {
        stop = eval_scalar(args[0]);
    } else if (args.size() == 2) {
        start = eval_scalar(args[0]);
        stop = eval_scalar(args[1]);
    } else {
        start = eval_scalar(args[0]);
        stop = eval_scalar(args[1]);
        step = eval_scalar(args[2]);
    }
    if (step == 0.0L) throw std::runtime_error("range step cannot be zero");
    std::vector<StoredValue> values;
    for (Scalar current = start; step > 0 ? current < stop : current > stop; current += step) {
        StoredValue item;
        item.decimal = current;
        values.push_back(item);
    }
    return make_list_value(std::move(values));
}

/**
 * @brief 求值列表推导式
 */
StoredValue evaluate_list_comprehension(
                                        IExecutionContext* ctx,
                                        const std::string& body,
                                        bool exact_mode) {
    const std::size_t for_pos = find_top_level_word(body, "for");
    if (for_pos == std::string::npos) throw std::runtime_error("invalid list comprehension");
    const std::string element_expr = utils::trim_copy(body.substr(0, for_pos));
    const std::string rest = utils::trim_copy(body.substr(for_pos + 3));
    const std::size_t in_pos = find_top_level_word(rest, "in");
    if (in_pos == std::string::npos) throw std::runtime_error("invalid list comprehension");
    const std::string var = utils::trim_copy(rest.substr(0, in_pos));
    if (var.empty()) throw std::runtime_error("invalid list comprehension variable");
    std::string iterable_expr = utils::trim_copy(rest.substr(in_pos + 2));
    std::string filter_expr;
    const std::size_t if_pos = find_top_level_word(iterable_expr, "if");
    if (if_pos != std::string::npos) {
        filter_expr = utils::trim_copy(iterable_expr.substr(if_pos + 2));
        iterable_expr = utils::trim_copy(iterable_expr.substr(0, if_pos));
    }

    StoredValue iterable = evaluate_script_value_expression(ctx, iterable_expr, exact_mode);
    if (!iterable.is_list || !iterable.list_value) {
        throw std::runtime_error("list comprehension requires a list iterable");
    }

    std::vector<StoredValue> result;
    ctx->variables().push_scope();
    try {
        for (const StoredValue& item : *iterable.list_value) {
            ctx->variables().set_local(var, item);
            if (!filter_expr.empty() &&
                !truthy_value(evaluate_script_value_expression(ctx, filter_expr, exact_mode))) {
                continue;
            }
            result.push_back(evaluate_script_value_expression(ctx, element_expr, exact_mode));
        }
        ctx->variables().pop_scope();
    } catch (...) {
        ctx->variables().pop_scope();
        throw;
    }
    return make_list_value(std::move(result));
}

/**
 * @brief 求值列表字面量
 */
StoredValue evaluate_list_literal(
                                  IExecutionContext* ctx,
                                  const std::string& expression,
                                  bool exact_mode) {
    const std::string inner = utils::trim_copy(expression.substr(1, expression.size() - 2));
    if (inner.empty()) return make_list_value({});
    if (has_top_level_semicolon(inner)) {
        throw std::runtime_error("not a script list literal");
    }
    if (find_top_level_word(inner, "for") != std::string::npos) {
        return evaluate_list_comprehension(ctx, inner, exact_mode);
    }
    std::vector<StoredValue> values;
    for (const std::string& part : split_script_top_level(inner, ',')) {
        if (!part.empty()) values.push_back(evaluate_script_value_expression(ctx, part, exact_mode));
    }
    return make_list_value(std::move(values));
}

/**
 * @brief 求值字典字面量
 */
StoredValue evaluate_dict_literal(
                                  IExecutionContext* ctx,
                                  const std::string& expression,
                                  bool exact_mode) {
    const std::string inner = utils::trim_copy(expression.substr(1, expression.size() - 2));
    std::map<std::string, StoredValue> values;
    if (inner.empty()) return make_dict_value(std::move(values));
    for (const std::string& item : split_script_top_level(inner, ',')) {
        const std::vector<std::string> kv = split_script_top_level(item, ':');
        if (kv.size() != 2) throw std::runtime_error("invalid dictionary literal");
        StoredValue key = evaluate_script_value_expression(ctx, kv[0], exact_mode);
        values[stored_to_key(key)] = evaluate_script_value_expression(ctx, kv[1], exact_mode);
    }
    return make_dict_value(std::move(values));
}

/**
 * @brief 求值索引或切片操作
 */
StoredValue evaluate_index_or_slice(
                                    IExecutionContext* ctx,
                                    const std::string& expression,
                                    bool exact_mode) {
    std::string base_expr;
    std::string index_expr;
    if (!parse_index_expression(expression, &base_expr, &index_expr)) {
        throw std::runtime_error("invalid index expression");
    }
    StoredValue base = evaluate_script_value_expression(ctx, base_expr, exact_mode);
    if (base.is_list) {
        if (!base.list_value) throw std::runtime_error("invalid list value");
        const auto& list = *base.list_value;
        const std::vector<std::string> slice_parts = split_script_top_level(index_expr, ':');
        if (slice_parts.size() > 1) {
            if (slice_parts.size() > 3) throw std::runtime_error("slice expects start:stop[:step]");
            long long start = 0;
            long long stop = static_cast<long long>(list.size());
            long long step = 1;
            if (!utils::trim_copy(slice_parts[0]).empty()) start = stored_to_index(evaluate_script_value_expression(ctx, slice_parts[0], exact_mode), "slice start");
            if (slice_parts.size() > 1 && !utils::trim_copy(slice_parts[1]).empty()) stop = stored_to_index(evaluate_script_value_expression(ctx, slice_parts[1], exact_mode), "slice stop");
            if (slice_parts.size() > 2 && !utils::trim_copy(slice_parts[2]).empty()) step = stored_to_index(evaluate_script_value_expression(ctx, slice_parts[2], exact_mode), "slice step");
            if (step == 0) throw std::runtime_error("slice step cannot be zero");
            auto normalize = [&](long long idx) {
                if (idx < 0) idx += static_cast<long long>(list.size());
                return std::max(0LL, std::min(idx, static_cast<long long>(list.size())));
            };
            start = normalize(start);
            stop = normalize(stop);
            std::vector<StoredValue> result;
            if (step > 0) {
                for (long long i = start; i < stop; i += step) result.push_back(list[static_cast<std::size_t>(i)]);
            } else {
                for (long long i = start; i > stop; i += step) result.push_back(list[static_cast<std::size_t>(i)]);
            }
            return make_list_value(std::move(result));
        }
        long long index = stored_to_index(evaluate_script_value_expression(ctx, index_expr, exact_mode), "list index");
        if (index < 0) index += static_cast<long long>(list.size());
        if (index < 0 || index >= static_cast<long long>(list.size())) throw std::runtime_error("list index out of range");
        return list[static_cast<std::size_t>(index)];
    }
    if (base.is_dict) {
        if (!base.dict_value) throw std::runtime_error("invalid dictionary value");
        const std::string key = stored_to_key(evaluate_script_value_expression(ctx, index_expr, exact_mode));
        auto it = base.dict_value->find(key);
        if (it == base.dict_value->end()) throw std::runtime_error("dictionary key not found: " + key);
        return it->second;
    }
    if (base.is_matrix) {
        if (!base.matrix_ptr) throw std::runtime_error("invalid matrix value");
        const std::vector<std::string> parts = split_script_top_level(index_expr, ',');
        if (parts.size() == 1) {
            long long index = stored_to_index(evaluate_script_value_expression(ctx, parts[0], exact_mode), "matrix index");
            if (index < 0) index += static_cast<long long>(base.matrix_ptr->data.size());
            if (index < 0 || index >= static_cast<long long>(base.matrix_ptr->data.size())) throw std::runtime_error("matrix index out of range");
            StoredValue result;
            result.decimal = base.matrix_ptr->data[static_cast<std::size_t>(index)];
            return result;
        } else if (parts.size() == 2) {
            long long row = stored_to_index(evaluate_script_value_expression(ctx, parts[0], exact_mode), "matrix row");
            long long col = stored_to_index(evaluate_script_value_expression(ctx, parts[1], exact_mode), "matrix col");
            if (row < 0) row += static_cast<long long>(base.matrix_ptr->rows);
            if (col < 0) col += static_cast<long long>(base.matrix_ptr->cols);
            if (row < 0 || row >= static_cast<long long>(base.matrix_ptr->rows) || col < 0 || col >= static_cast<long long>(base.matrix_ptr->cols)) {
                throw std::runtime_error("matrix index out of range");
            }
            StoredValue result;
            result.decimal = base.matrix_ptr->at(static_cast<std::size_t>(row), static_cast<std::size_t>(col));
            return result;
        } else {
            throw std::runtime_error("matrix indexing requires 1 or 2 indices");
        }
    }
    throw std::runtime_error("indexing requires a list, dictionary, or matrix");
}

StoredValue evaluate_script_value_expression(
                                             IExecutionContext* ctx,
                                             const std::string& expression,
                                             bool exact_mode) {
    const std::string trimmed = utils::trim_copy(expression);

    if (trimmed.rfind("len(", 0) == 0 && trimmed.back() == ')') {
        std::string arg = utils::trim_copy(trimmed.substr(4, trimmed.size() - 5));
        if (arg.empty()) throw std::runtime_error("len() requires one argument");
        StoredValue value = evaluate_script_value_expression(ctx, arg, exact_mode);
        StoredValue result;
        result.decimal = Scalar(0.0L);
        if (value.is_list && value.list_value) {
            result.decimal = Scalar(static_cast<long long>(value.list_value->size()));
        } else if (value.is_dict && value.dict_value) {
            result.decimal = Scalar(static_cast<long long>(value.dict_value->size()));
        } else if (value.is_string) {
            result.decimal = Scalar(static_cast<long long>(value.string_value.size()));
        } else {
            throw std::runtime_error("len() requires a list, dict, or string");
        }
        return result;
    }

    if (trimmed.rfind("range(", 0) == 0 && is_wrapped_by(trimmed.substr(5), '(', ')')) {
        return evaluate_range_list(ctx, trimmed.substr(6, trimmed.size() - 7), exact_mode);
    }

    if (trimmed.rfind("append(", 0) == 0 && trimmed.back() == ')') {
        std::string args_text = utils::trim_copy(trimmed.substr(7, trimmed.size() - 8));
        std::vector<std::string> args = split_script_top_level(args_text, ',');
        if (args.size() != 2) throw std::runtime_error("append() requires two arguments: list and element");
        StoredValue list_val = evaluate_script_value_expression(ctx, args[0], exact_mode);
        if (!list_val.is_list || !list_val.list_value) throw std::runtime_error("append() first argument must be a list");
        StoredValue elem = evaluate_script_value_expression(ctx, args[1], exact_mode);
        std::vector<StoredValue> new_list = *list_val.list_value;
        new_list.push_back(elem);
        return make_list_value(std::move(new_list));
    }

    if (is_wrapped_by(trimmed, '[', ']')) {
        const std::string inner = utils::trim_copy(trimmed.substr(1, trimmed.size() - 2));
        if (!has_top_level_semicolon(inner)) return evaluate_list_literal(ctx, trimmed, exact_mode);
    }
    if (is_wrapped_by(trimmed, '{', '}')) return evaluate_dict_literal(ctx, trimmed, exact_mode);
    
    std::string base;
    std::string index;
    if (parse_index_expression(trimmed, &base, &index)) return evaluate_index_or_slice(ctx, trimmed, exact_mode);
    
    return evaluate_expression_value(ctx, trimmed, exact_mode);
}

bool try_execute_index_assignment(
                                  IExecutionContext* ctx,
                                  const std::string& text,
                                  bool exact_mode,
                                  std::string* output) {
    const std::size_t eq = find_top_level_assignment(text);
    if (eq == std::string::npos) return false;
    std::string target = utils::trim_copy(text.substr(0, eq));
    std::string value_expr = utils::trim_copy(text.substr(eq + 1));
    std::string base_name;
    std::string index_expr;
    if (!parse_index_expression(target, &base_name, &index_expr)) return false;
    if (base_name.empty() || (!std::isalpha(static_cast<unsigned char>(base_name[0])) && base_name[0] != '_')) return false;

    auto opt_val = ctx->variables().get(base_name);
    if (!opt_val) throw std::runtime_error("unknown variable: " + base_name);
    StoredValue base_value = *opt_val;

    StoredValue new_value = evaluate_script_value_expression(ctx, value_expr, exact_mode);
    if (base_value.is_list) {
        if (!base_value.list_value) base_value.list_value = std::make_shared<std::vector<StoredValue>>();
        long long index = stored_to_index(evaluate_script_value_expression(ctx, index_expr, exact_mode), "list index");
        auto& list = *base_value.list_value;
        if (index < 0) index += static_cast<long long>(list.size());
        if (index < 0 || index >= static_cast<long long>(list.size())) throw std::runtime_error("list index out of range");
        list[static_cast<std::size_t>(index)] = new_value;
        *output = base_name + "[" + index_expr + "] = " + format_stored_value(new_value, ctx->config().is_symbolic_constants_mode());
        ctx->variables().assign_visible(base_name, base_value);
        return true;
    }
    if (base_value.is_dict) {
        if (!base_value.dict_value) base_value.dict_value = std::make_shared<std::map<std::string, StoredValue>>();
        const std::string key = stored_to_key(evaluate_script_value_expression(ctx, index_expr, exact_mode));
        (*base_value.dict_value)[key] = new_value;
        *output = base_name + "[" + index_expr + "] = " + format_stored_value(new_value, ctx->config().is_symbolic_constants_mode());
        ctx->variables().assign_visible(base_name, base_value);
        return true;
    }
    if (base_value.is_matrix) {
        if (!base_value.matrix_ptr) throw std::runtime_error("invalid matrix value");
        if (new_value.is_matrix || new_value.is_list || new_value.is_dict || new_value.is_string) {
            throw std::runtime_error("matrix element assignment requires a scalar value");
        }
        Scalar val = new_value.exact ? rational_to_double(new_value.rational) : new_value.decimal;
        const std::vector<std::string> parts = split_script_top_level(index_expr, ',');
        if (parts.size() == 1) {
            long long index = stored_to_index(evaluate_script_value_expression(ctx, parts[0], exact_mode), "matrix index");
            if (index < 0) index += static_cast<long long>(base_value.matrix_ptr->data.size());
            if (index < 0 || index >= static_cast<long long>(base_value.matrix_ptr->data.size())) throw std::runtime_error("matrix index out of range");
            base_value.matrix_ptr->data[static_cast<std::size_t>(index)] = val;
            *output = base_name + "[" + index_expr + "] = " + format_stored_value(new_value, ctx->config().is_symbolic_constants_mode());
            ctx->variables().assign_visible(base_name, base_value);
            return true;
        } else if (parts.size() == 2) {
            long long row = stored_to_index(evaluate_script_value_expression(ctx, parts[0], exact_mode), "matrix row");
            long long col = stored_to_index(evaluate_script_value_expression(ctx, parts[1], exact_mode), "matrix col");
            if (row < 0) row += static_cast<long long>(base_value.matrix_ptr->rows);
            if (col < 0) col += static_cast<long long>(base_value.matrix_ptr->cols);
            if (row < 0 || row >= static_cast<long long>(base_value.matrix_ptr->rows) || col < 0 || col >= static_cast<long long>(base_value.matrix_ptr->cols)) {
                throw std::runtime_error("matrix index out of range");
            }
            base_value.matrix_ptr->at(static_cast<std::size_t>(row), static_cast<std::size_t>(col)) = val;
            *output = base_name + "[" + index_expr + "] = " + format_stored_value(new_value, ctx->config().is_symbolic_constants_mode());
            ctx->variables().assign_visible(base_name, base_value);
            return true;
        }
    }
    return false;
}

StoredValue evaluate_expression_value(
                                      IExecutionContext* ctx,
                                      const std::string& expression,
                                      bool exact_mode,
                                      std::shared_ptr<ExpressionCache>* cache) {
    const std::string special_expr = utils::trim_copy(expression);
    std::string index_base;
    std::string index_expr;
    if (parse_index_expression(special_expr, &index_base, &index_expr)) {
        return evaluate_index_or_slice(ctx, special_expr, exact_mode);
    }

    const VariableResolver variables = visible_variables(ctx);
    
    std::shared_ptr<ExpressionCache> expr_cache;
    if (cache) {
        if (!*cache) {
            *cache = std::make_shared<ExpressionCache>(expression);
            std::string expanded = ctx->expand_inline(expression);
            if (expanded != expression) (*cache)->set_expanded(std::move(expanded));
        } else if (!(*cache)->has_expanded) {
            std::string expanded = ctx->expand_inline(std::string((*cache)->original_text));
            if (expanded != (*cache)->original_text) (*cache)->set_expanded(std::move(expanded));
        }
        expr_cache = *cache;
    }

    const std::string target_expr =
        expr_cache ? std::string(expr_cache->effective_text())
                   : ctx->expand_inline(expression);

    std::string trimmed_expr = utils::trim_copy(target_expr);
    if (is_string_literal(trimmed_expr)) {
        StoredValue res; res.is_string = true; res.string_value = parse_string_literal_value(trimmed_expr);
        return res;
    }

    CommandParser::IsCommandCallback is_cmd = [ctx](std::string_view name) {
        return ctx->commands().has_command(std::string(name)) || ctx->functions().get_native_functions()->count(std::string(name)) > 0;
    };
    try {
        CommandASTNode ast = parse_command(trimmed_expr, is_cmd);
        if (ast.kind == CommandKind::kFunctionCall) {
            const auto* call = ast.as_function_call();
            if (ctx->functions().get_native_functions()->count(std::string(call->name)) > 0) {
                return evaluate_command_ast_to_value(ctx, ast, exact_mode);
            }
        }
    } catch (...) {}

    if (!exact_mode) {
        StoredValue implicit_value;
        if (ctx->try_evaluate_implicit(target_expr, &implicit_value, variables.snapshot())) return implicit_value;
    }

    std::string converted;
    if (try_base_conversion_expression(target_expr, variables, ctx->functions().get_custom_functions_map(), {ctx->config().is_hex_prefix_mode(), ctx->config().is_hex_uppercase_mode()}, &converted)) {
        StoredValue res; res.is_string = true; res.string_value = converted;
        return res;
    }

    const HasScriptFunctionCallback has_script_function = [ctx](const std::string& name) {
        return has_visible_script_function(ctx, name);
    };
    const InvokeScriptFunctionCallback invoke_script_function = [ctx](const std::string& name, const std::vector<Scalar>& arguments) {
        return invoke_script_function_decimal(ctx, name, arguments);
    };

    UnifiedExpressionParser parser(variables, ctx->functions().get_custom_functions_map(),
                                   nullptr,
                                   nullptr,
                                   nullptr,
                                   ctx->functions().get_native_functions(),
                                   has_script_function, invoke_script_function);
    StoredValue result = parser.evaluate_stored(target_expr, exact_mode, ctx->config().is_symbolic_constants_mode());

    if (ctx->config().is_symbolic_constants_mode() && !result.has_symbolic_text && !result.is_string && !result.is_matrix) {
        std::string symbolic_text;
        if (try_symbolic_constant_expression(target_expr, variables, ctx->functions().get_custom_functions_map(), &symbolic_text)) {
            bool symbolic_text_is_plain_decimal = false;
            try {
                const Scalar parsed = mymath::scalar_from_string(symbolic_text);
                symbolic_text_is_plain_decimal = mymath::abs(parsed - result.decimal) < Scalar(1e-12L);
            } catch (...) { symbolic_text_is_plain_decimal = false; }
            if (!symbolic_text_is_plain_decimal) {
                result.has_symbolic_text = true;
                result.symbolic_text = symbolic_text;
            }
        }
    }
    
    return result;
}
