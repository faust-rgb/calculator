#include "core/services/core_manager_interfaces.h"
#include "execution/functions/user_function.h"
#include "core/services/core_manager_interfaces.h"
// ============================================================================
// script_executor.cpp - 脚本语句与块执行实现
// ============================================================================

#include "script_runtime_internal.h"
#include "core/services/string_utils.h"
#include "math/functions/integer/integer_helpers.h"
#include "mymath.h"
#include "matrix/matrix.h"
#include <sstream>
#include <stdexcept>

ScriptSignal execute_script_statement(
                                      IExecutionContext* ctx,
                                      const script::Statement& statement,
                                      bool exact_mode,
                                      std::string* last_output,
                                      bool create_scope) {
    switch (statement.kind) {
        case script::Statement::Kind::kBlock:
            return execute_script_block(ctx, static_cast<const script::BlockStatement&>(statement), exact_mode, last_output, create_scope);
        
        case script::Statement::Kind::kSimple: {
            const auto& simple = static_cast<const script::SimpleStatement&>(statement);
            try {
                if (simple.command_ast.kind == CommandKind::kIndexAssignment) {
                    const auto* assign = simple.command_ast.as_index_assignment();
                    std::string output;
                    if (execute_index_assignment_direct(ctx, assign->variable, assign->indices, assign->value, exact_mode, &output)) {
                        *last_output = output;
                        return {};
                    }
                }
                if (try_execute_index_assignment(ctx, simple.text, exact_mode, last_output)) {
                    return {};
                }
                if (simple.command_ast.kind != CommandKind::kEmpty) {
                    *last_output = execute_command_ast(ctx, simple.command_ast, exact_mode);
                }
            } catch (const std::exception& e) {
                throw std::runtime_error(format_script_error(simple.line, e));
            }
            return {};
        }

        case script::Statement::Kind::kIf: {
            const auto& if_stmt = static_cast<const script::IfStatement&>(statement);
            try {
                if (truthy_value(evaluate_command_ast_to_value(ctx, if_stmt.condition_ast, false))) {
                    return execute_script_statement(ctx, *if_stmt.then_branch, exact_mode, last_output, false);
                } else if (if_stmt.else_branch) {
                    return execute_script_statement(ctx, *if_stmt.else_branch, exact_mode, last_output, false);
                }
            } catch (const std::exception& e) {
                throw std::runtime_error(format_script_error(if_stmt.line, e));
            }
            return {};
        }

        case script::Statement::Kind::kWhile: {
            const auto& while_stmt = static_cast<const script::WhileStatement&>(statement);
            try {
                while (truthy_value(evaluate_command_ast_to_value(ctx, while_stmt.condition_ast, false))) {
                    const ScriptSignal signal = execute_script_statement(ctx, *while_stmt.body, exact_mode, last_output, false);
                    if (signal.kind == ScriptSignal::Kind::kReturn) return signal;
                    if (signal.kind == ScriptSignal::Kind::kBreak) break;
                    if (signal.kind == ScriptSignal::Kind::kContinue) continue;
                }
            } catch (const std::exception& e) {
                throw std::runtime_error(format_script_error(while_stmt.line, e));
            }
            return {};
        }

        case script::Statement::Kind::kFor: {
            const auto& for_stmt = static_cast<const script::ForStatement&>(statement);
            try {
                if (for_stmt.init_ast.kind != CommandKind::kEmpty) (void)evaluate_command_ast_to_value(ctx, for_stmt.init_ast, exact_mode);
                while (for_stmt.cond_ast.kind == CommandKind::kEmpty || truthy_value(evaluate_command_ast_to_value(ctx, for_stmt.cond_ast, false))) {
                    const ScriptSignal signal = execute_script_statement(ctx, *for_stmt.body, exact_mode, last_output, false);
                    if (signal.kind == ScriptSignal::Kind::kReturn) return signal;
                    if (signal.kind == ScriptSignal::Kind::kBreak) break;
                    if (for_stmt.step_ast.kind != CommandKind::kEmpty) (void)evaluate_command_ast_to_value(ctx, for_stmt.step_ast, exact_mode);
                    if (signal.kind == ScriptSignal::Kind::kContinue) continue;
                }
                return {};
            } catch (const std::exception& e) {
                throw std::runtime_error(format_script_error(for_stmt.line, e));
            }
        }

        case script::Statement::Kind::kForRange: {
            const auto& for_range = static_cast<const script::ForRangeStatement&>(statement);
            try {
                const Scalar start = evaluate_scalar(ctx, for_range.start_ast, "range start");
                const Scalar stop = evaluate_scalar(ctx, for_range.stop_ast, "range stop");
                const Scalar step = evaluate_scalar(ctx, for_range.step_ast, "range step");
                if (step == 0.0L) throw std::runtime_error("range step cannot be zero");
                bool ascending = step > 0.0L;
                for (Scalar current = start; (ascending ? current < stop : current > stop); current += step) {
                    StoredValue loop_val(current);
                    ctx->variables().set_local(for_range.variable, loop_val);
                    const ScriptSignal signal = execute_script_statement(ctx, *for_range.body, exact_mode, last_output, false);
                    if (signal.kind == ScriptSignal::Kind::kReturn) return signal;
                    if (signal.kind == ScriptSignal::Kind::kBreak) break;
                    if (signal.kind == ScriptSignal::Kind::kContinue) continue;
                }
                return {};
            } catch (const std::exception& e) {
                throw std::runtime_error(format_script_error(for_range.line, e));
            }
        }

        case script::Statement::Kind::kForIn: {
            const auto& for_in = static_cast<const script::ForInStatement&>(statement);
            try {
                StoredValue iterable = evaluate_command_ast_to_value(ctx, for_in.iterable_ast, exact_mode);

                if (iterable.is_list && iterable.list_value) {
                    for (const StoredValue& item : *iterable.list_value) {
                        ctx->variables().set_local(for_in.variable, item);
                        const ScriptSignal signal = execute_script_statement(ctx, *for_in.body, exact_mode, last_output, false);
                        if (signal.kind == ScriptSignal::Kind::kReturn) return signal;
                        if (signal.kind == ScriptSignal::Kind::kBreak) break;
                        if (signal.kind == ScriptSignal::Kind::kContinue) continue;
                    }
                }
                else if (iterable.is_matrix) {
                    if (!iterable.matrix_ptr) throw std::runtime_error("invalid matrix value");
                    for (std::size_t row = 0; row < iterable.matrix_ptr->rows; ++row) {
                        StoredValue row_value;
                        row_value.is_matrix = true;
                        row_value.matrix_ptr = std::make_shared<matrix::Matrix>(1, iterable.matrix_ptr->cols);
                        for (std::size_t col = 0; col < iterable.matrix_ptr->cols; ++col) {
                            row_value.matrix_ptr->at(0, col) = iterable.matrix_ptr->at(row, col);
                        }
                        ctx->variables().set_local(for_in.variable, row_value);
                        const ScriptSignal signal = execute_script_statement(ctx, *for_in.body, exact_mode, last_output, false);
                        if (signal.kind == ScriptSignal::Kind::kReturn) return signal;
                        if (signal.kind == ScriptSignal::Kind::kBreak) break;
                        if (signal.kind == ScriptSignal::Kind::kContinue) continue;
                    }
                }
                else if (iterable.is_string) {
                    for (char ch : iterable.get_string_value()) {
                        StoredValue char_value;
                        char_value.is_string = true;
                        char_value.string_value = std::string(1, ch);
                        ctx->variables().set_local(for_in.variable, char_value);
                        const ScriptSignal signal = execute_script_statement(ctx, *for_in.body, exact_mode, last_output, false);
                        if (signal.kind == ScriptSignal::Kind::kReturn) return signal;
                        if (signal.kind == ScriptSignal::Kind::kBreak) break;
                        if (signal.kind == ScriptSignal::Kind::kContinue) continue;
                    }
                }
                else {
                    throw std::runtime_error("for-in loop requires a list, matrix, or string iterable");
                }
                return {};
            } catch (const std::exception& e) {
                throw std::runtime_error(format_script_error(for_in.line, e));
            }
        }

        case script::Statement::Kind::kReturn: {
            const auto& ret = static_cast<const script::ReturnStatement&>(statement);
            if (!ret.has_expression) return ScriptSignal::make_return({});
            try { return ScriptSignal::make_return(evaluate_command_ast_to_value(ctx, ret.expr_ast, exact_mode)); }
            catch (const std::exception& e) { throw std::runtime_error(format_script_error(ret.line, e)); }
        }

        case script::Statement::Kind::kBreak: return ScriptSignal::make_break();
        case script::Statement::Kind::kContinue: return ScriptSignal::make_continue();
        case script::Statement::Kind::kImport: {
            const auto& import_stmt = static_cast<const script::ImportStatement&>(statement);
            try {
                StoredValue path_value = evaluate_command_ast_to_value(ctx, import_stmt.path_ast, exact_mode);
                if (!path_value.is_string) throw std::runtime_error("import path must be a string");
                *last_output = ctx->execute_script_file(path_value.string_value, exact_mode, true);
                return {};
            } catch (const std::exception& e) {
                throw std::runtime_error(format_script_error(import_stmt.line, e));
            }
        }
        case script::Statement::Kind::kFunction: {
            const auto& fs = static_cast<const script::FunctionStatement&>(statement);
            ScriptFunction function;
            function.parameter_names = fs.parameters;
            function.default_values = fs.default_values;
            function.body = fs.body;
            ctx->functions().add_script_function(fs.name, function);
            *last_output = fs.name + "(...)";
            return {};
        }

        case script::Statement::Kind::kMatch: {
            const auto& match_stmt = static_cast<const script::MatchStatement&>(statement);
            try {
                StoredValue subject = evaluate_command_ast_to_value(ctx, match_stmt.subject_ast, exact_mode);
                for (const auto& clause : match_stmt.cases) {
                    bool matches = false;
                    if (clause.is_default) {
                        matches = true;
                    } else {
                        StoredValue pattern = evaluate_command_ast_to_value(ctx, clause.pattern_ast, exact_mode);
                        if (subject.is_matrix || pattern.is_matrix) {
                            if (subject.is_matrix && pattern.is_matrix) {
                                if (subject.matrix_ptr && pattern.matrix_ptr &&
                                    subject.matrix_ptr->rows == pattern.matrix_ptr->rows && 
                                    subject.matrix_ptr->cols == pattern.matrix_ptr->cols) {
                                    matches = true;
                                    for (std::size_t i = 0; i < subject.matrix_ptr->data.size(); ++i) {
                                        if (!mymath::is_near_zero(subject.matrix_ptr->data[i] - pattern.matrix_ptr->data[i], 1e-10)) {
                                            matches = false;
                                            break;
                                        }
                                    }
                                } else matches = false;
                            } else matches = false;
                        } else if (subject.is_complex || pattern.is_complex) {
                            if (subject.is_complex && pattern.is_complex) {
                                matches = mymath::is_near_zero(subject.complex.real() - pattern.complex.real(), 1e-10) &&
                                          mymath::is_near_zero(subject.complex.imag() - pattern.complex.imag(), 1e-10);
                            } else matches = false;
                        } else if (subject.is_string || pattern.is_string) {
                            if (subject.is_string && pattern.is_string) matches = subject.string_value == pattern.string_value;
                            else matches = false;
                        } else {
                            Scalar subj_val = subject.get_decimal();
                            Scalar pat_val = pattern.get_decimal();
                            matches = mymath::is_near_zero(subj_val - pat_val, 1e-10);
                        }
                    }
                    if (matches && clause.is_guarded) {
                        StoredValue guard_val = evaluate_command_ast_to_value(ctx, clause.guard_ast, exact_mode);
                        matches = truthy_value(guard_val);
                    }
                    if (matches) return execute_script_statement(ctx, *clause.body, exact_mode, last_output, false);
                }
            } catch (const std::exception& e) {
                throw std::runtime_error(format_script_error(match_stmt.line, e));
            }
            return {};
        }
    }
    return {};
}

#include "parser/ast/unified_ast.h"

ScriptSignal execute_script_block(
                                  IExecutionContext* ctx,
                                  const script::BlockStatement& block,
                                  bool exact_mode,
                                  std::string* last_output,
                                  bool create_scope) {
    std::unique_ptr<core::ExecutionContext::ScopeGuard> guard;
    if (create_scope && ctx) {
        guard = std::make_unique<core::ExecutionContext::ScopeGuard>(ctx->core_context());
        ctx->variables().push_scope();
    }
    try {
        for (const auto& stmt : block.statements) {
            const ScriptSignal signal = execute_script_statement(ctx, *stmt, exact_mode, last_output, false);
            if (signal.kind != ScriptSignal::Kind::kNone) {
                if (create_scope) ctx->variables().pop_scope();
                return signal;
            }
        }
        if (create_scope) ctx->variables().pop_scope();
        return {};
    } catch (const core::ScriptRuntimeError&) {
        if (create_scope) ctx->variables().pop_scope();
        throw;
    } catch (const std::exception& e) {
        if (create_scope) ctx->variables().pop_scope();
        core::SourceLocation loc;
        loc.line = block.line;
        throw core::ScriptRuntimeError(e.what(), loc);
    }
}

std::string execute_simple_script_line(
                                       IExecutionContext* ctx,
                                       const std::string& text,
                                       bool exact_mode) {
    const std::string trimmed = utils::trim_copy(text);
    if (trimmed.empty()) return "";
    auto is_command = [ctx](std::string_view name) {
        return ctx->commands().has_command(std::string(name)) || ctx->commands().has_command(":" + std::string(name));
    };
    CommandASTNode ast = parse_command(trimmed, is_command);
    return execute_command_ast(ctx, ast, exact_mode);
}

namespace {

std::string command_ast_to_source(const CommandASTNode& ast) {
    switch (ast.kind) {
        case CommandKind::kEmpty: return "";
        case CommandKind::kMetaCommand: {
            const auto& meta = std::get<MetaCommandInfo>(ast.data);
            std::string text = ":" + std::string(meta.command);
            for (std::string_view arg : meta.arguments) { text += " "; text += std::string(arg); }
            return text;
        }
        case CommandKind::kFunctionDefinition: {
            const auto& def = std::get<FunctionDefinitionInfo>(ast.data);
            std::string text = std::string(def.name) + "(";
            for (std::size_t i = 0; i < def.parameters.size(); ++i) {
                if (i != 0) text += ", ";
                text += std::string(def.parameters[i]);
            }
            text += ") = "; text += std::string(def.body.text);
            return text;
        }
        case CommandKind::kFunctionCall: {
            const auto& call = std::get<FunctionCallInfo>(ast.data);
            std::string text = std::string(call.name) + "(";
            for (std::size_t i = 0; i < call.arguments.size(); ++i) {
                if (i != 0) text += ", ";
                text += std::string(call.arguments[i].text);
            }
            text += ")"; return text;
        }
        case CommandKind::kAssignment: {
            const auto& assignment = std::get<AssignmentInfo>(ast.data);
            return std::string(assignment.variable) + " = " + std::string(assignment.expression.text);
        }
        case CommandKind::kIndexAssignment: {
            const auto& idx_assign = std::get<IndexAssignmentInfo>(ast.data);
            std::string text = std::string(idx_assign.variable) + "[";
            for (std::size_t i = 0; i < idx_assign.indices.size(); ++i) {
                if (i != 0) text += ", ";
                text += std::string(idx_assign.indices[i].text);
            }
            text += "] = " + std::string(idx_assign.value.text);
            return text;
        }
        case CommandKind::kExpression: return std::string(std::get<ExpressionInfo>(ast.data).text);
        case CommandKind::kStringLiteral: return "\"" + std::get<std::string>(ast.data) + "\"";
        case CommandKind::kSequence: {
            const auto& nodes = std::get<std::vector<CommandASTNode>>(ast.data);
            std::string text;
            for (std::size_t i = 0; i < nodes.size(); ++i) {
                if (i != 0) text += "; ";
                text += command_ast_to_source(nodes[i]);
            }
            return text;
        }
    }
    return "";
}

} // namespace

std::string render_script_statement(const script::Statement& statement, int indent) {
    const std::string pad = indent_text(indent);
    switch (statement.kind) {
        case script::Statement::Kind::kBlock:
            return render_script_block(static_cast<const script::BlockStatement&>(statement), indent);
        case script::Statement::Kind::kSimple: {
            const auto& simple = static_cast<const script::SimpleStatement&>(statement);
            const std::string text = simple.text.empty() ? command_ast_to_source(simple.command_ast) : simple.text;
            return pad + (text.empty() ? "pass" : text);
        }
        case script::Statement::Kind::kIf: {
            const auto& if_stmt = static_cast<const script::IfStatement&>(statement);
            std::string text = pad + "if " + command_ast_to_source(if_stmt.condition_ast) + " " + render_script_statement(*if_stmt.then_branch, indent);
            if (if_stmt.else_branch) text += "\n" + pad + "else " + render_script_statement(*if_stmt.else_branch, indent);
            return text;
        }
        case script::Statement::Kind::kWhile: {
            const auto& while_stmt = static_cast<const script::WhileStatement&>(statement);
            return pad + "while " + command_ast_to_source(while_stmt.condition_ast) + " " + render_script_statement(*while_stmt.body, indent);
        }
        case script::Statement::Kind::kFor: {
            const auto& for_stmt = static_cast<const script::ForStatement&>(statement);
            return pad + "for (" + command_ast_to_source(for_stmt.init_ast) + "; " +
                   command_ast_to_source(for_stmt.cond_ast) + "; " +
                   command_ast_to_source(for_stmt.step_ast) + ") " +
                   render_script_statement(*for_stmt.body, indent);
        }
        case script::Statement::Kind::kForRange: {
            const auto& for_stmt = static_cast<const script::ForRangeStatement&>(statement);
            return pad + "for " + for_stmt.variable + " in range(" +
                   command_ast_to_source(for_stmt.start_ast) + ", " +
                   command_ast_to_source(for_stmt.stop_ast) + ", " +
                   command_ast_to_source(for_stmt.step_ast) + ") " +
                   render_script_statement(*for_stmt.body, indent);
        }
        case script::Statement::Kind::kForIn: {
            const auto& for_in = static_cast<const script::ForInStatement&>(statement);
            return pad + "for " + for_in.variable + " in " + command_ast_to_source(for_in.iterable_ast) + ": " + render_script_statement(*for_in.body, indent);
        }
        case script::Statement::Kind::kFunction: {
            const auto& function = static_cast<const script::FunctionStatement&>(statement);
            std::string text = pad + "fn " + function.name + "(";
            for (std::size_t i = 0; i < function.parameters.size(); ++i) {
                if (i != 0) text += ", ";
                text += function.parameters[i];
            }
            text += ") " + render_script_block(*function.body, indent);
            return text;
        }
        case script::Statement::Kind::kReturn: {
            const auto& ret = static_cast<const script::ReturnStatement&>(statement);
            return pad + "return" + (ret.has_expression ? " " + command_ast_to_source(ret.expr_ast) : "");
        }
        case script::Statement::Kind::kBreak: return pad + "break";
        case script::Statement::Kind::kContinue: return pad + "continue";
        case script::Statement::Kind::kImport: {
            const auto& import_stmt = static_cast<const script::ImportStatement&>(statement);
            return pad + "import " + (import_stmt.path_text.empty() ? command_ast_to_source(import_stmt.path_ast) : import_stmt.path_text);
        }
        case script::Statement::Kind::kMatch: {
            const auto& match_stmt = static_cast<const script::MatchStatement&>(statement);
            std::ostringstream out;
            out << pad << "match " << command_ast_to_source(match_stmt.subject_ast) << ":\n";
            for (const auto& clause : match_stmt.cases) {
                out << pad << "  case ";
                if (clause.is_default) out << "_";
                else {
                    out << command_ast_to_source(clause.pattern_ast);
                    if (clause.is_guarded) out << " if " << command_ast_to_source(clause.guard_ast);
                }
                out << ": " << render_script_statement(*clause.body, indent + 1) << "\n";
            }
            return out.str();
        }
    }
    return pad + "pass";
}

std::string render_script_block(const script::BlockStatement& block, int indent) {
    std::ostringstream out;
    out << "{";
    if (!block.statements.empty()) out << "\n";
    for (const auto& stmt : block.statements) out << render_script_statement(*stmt, indent + 1) << "\n";
    out << indent_text(indent) << "}";
    return out.str();
}

namespace script {

StoredValue Statement::execute(IExecutionContext* ctx, bool exact_mode) const {
    std::string dummy_output;
    ScriptSignal signal = execute_script_statement(ctx, *this, exact_mode, &dummy_output, false);
    return signal.has_value ? signal.value : StoredValue();
}

} // namespace script
