// ============================================================================
// state_persistence.cpp - 状态持久化服务实现
// ============================================================================
//
// 从 calculator_core.cpp 中提取的状态保存/加载逻辑。
// ============================================================================

#include "state_persistence.h"
#include "parser/grammars/command_parser.h"
#include "parser/grammars/script_parser.h"
#include "execution/engine/script_runtime.h"
#include "matrix/matrix.h"
#include "core/services/string_utils.h"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

StatePersistenceService::StatePersistenceService(
    std::shared_ptr<IVariableManager> variables,
    std::shared_ptr<IFunctionManager> functions)
    : variables_(std::move(variables))
    , functions_(std::move(functions)) {}

void StatePersistenceService::set_script_executor(std::function<void(const std::string&)> executor) {
    script_executor_ = std::move(executor);
}

std::string StatePersistenceService::encode_field(const std::string& text) {
    std::string result;
    result.reserve(text.size());
    for (char c : text) {
        switch (c) {
            case '\n': result += "\\n"; break;
            case '\r': result += "\\r"; break;
            case '\t': result += "\\t"; break;
            case '\\': result += "\\\\"; break;
            default: result += c; break;
        }
    }
    return result;
}

std::string StatePersistenceService::decode_field(const std::string& text) {
    std::string result;
    result.reserve(text.size());
    for (std::size_t i = 0; i < text.size(); ++i) {
        if (text[i] == '\\' && i + 1 < text.size()) {
            switch (text[i + 1]) {
                case 'n': result += '\n'; ++i; break;
                case 'r': result += '\r'; ++i; break;
                case 't': result += '\t'; ++i; break;
                case '\\': result += '\\'; ++i; break;
                default: result += text[i]; break;
            }
        } else {
            result += text[i];
        }
    }
    return result;
}

std::vector<std::string> StatePersistenceService::split_tab_fields(const std::string& row_text) {
    std::vector<std::string> parts;
    std::size_t start = 0;
    for (std::size_t i = 0; i <= row_text.size(); ++i) {
        if (i == row_text.size() || row_text[i] == '\t') {
            parts.push_back(row_text.substr(start, i - start));
            start = i + 1;
        }
    }
    return parts;
}

std::string StatePersistenceService::save_state(const std::string& path) const {
    const std::filesystem::path target_path(path);
    const std::filesystem::path temp_path =
        target_path.parent_path() /
        (target_path.filename().string() + ".tmp-save");
    std::ofstream out(temp_path);
    if (!out) {
        throw std::runtime_error("unable to open file for writing: " + path);
    }

    auto all_vars = variables_->get_all_globals();
    auto custom_names = functions_->get_custom_names();
    auto script_names = functions_->get_script_names();

    out << "STATE_V5\n";

    // 保存变量
    for (const auto& [name, value] : all_vars) {
        if (value.is_matrix) {
            out << "VAR\t" << encode_field(name)
                << "\tMATRIX\t" << value.matrix_ptr->rows
                << '\t' << value.matrix_ptr->cols;
            for (Scalar element : value.matrix_ptr->data) {
                out << '\t' << std::setprecision(17) << element;
            }
            out << '\n';
        } else if (value.is_complex) {
            out << "VAR\t" << encode_field(name)
                << "\tCOMPLEX\t" << std::setprecision(17) << value.complex.real()
                << '\t' << std::setprecision(17) << value.complex.imag() << '\n';
        } else if (value.is_string) {
            out << "VAR\t" << encode_field(name)
                << "\tSTRING\t" << encode_field(value.string_value) << '\n';
        } else if (value.exact) {
            out << "VAR\t" << encode_field(name)
                << "\tEXACT\t" << value.rational.numerator
                << '\t' << value.rational.denominator
                << '\t' << std::setprecision(17) << value.decimal << '\n';
        } else {
            out << "VAR\t" << encode_field(name)
                << "\tDECIMAL\t" << std::setprecision(17) << value.decimal << '\n';
        }
        if (value.has_precise_decimal_text) {
            out << "PRECISE\t" << encode_field(name)
                << '\t' << encode_field(value.precise_decimal_text) << '\n';
        }
        if (value.has_symbolic_text) {
            out << "SYMBOLIC\t" << encode_field(name)
                << '\t' << encode_field(value.symbolic_text) << '\n';
        }
    }

    // 保存自定义函数
    for (const auto& name : custom_names) {
        const CustomFunction* func = functions_->get_custom(name);
        if (!func) continue;
        std::string params_str;
        for (std::size_t i = 0; i < func->parameter_names.size(); ++i) {
            if (i != 0) params_str += ", ";
            params_str += func->parameter_names[i];
        }
        out << "EXPRFUNC\t"
            << encode_field(name + "(" + params_str + ") = " + func->expression)
            << '\n';
    }

    // 保存脚本函数
    for (const auto& name : script_names) {
        const ScriptFunction* func = functions_->get_script(name);
        if (!func) continue;
        std::ostringstream source;
        source << "fn " << name << "(";
        for (std::size_t i = 0; i < func->parameter_names.size(); ++i) {
            if (i != 0) {
                source << ", ";
            }
            source << func->parameter_names[i];
        }
        source << ") " << render_script_block(*func->body, 0);
        out << "SCRIPT\t" << encode_field(source.str()) << '\n';
    }

    out.close();
    if (!out) {
        std::error_code remove_error;
        std::filesystem::remove(temp_path, remove_error);
        throw std::runtime_error("unable to finish writing state file: " + path);
    }

    std::error_code rename_error;
    std::filesystem::rename(temp_path, target_path, rename_error);
    if (rename_error) {
        std::error_code remove_existing_error;
        std::filesystem::remove(target_path, remove_existing_error);
        rename_error.clear();
        std::filesystem::rename(temp_path, target_path, rename_error);
    }
    if (rename_error) {
        std::error_code remove_error;
        std::filesystem::remove(temp_path, remove_error);
        throw std::runtime_error("unable to replace state file: " + path);
    }

    return "Saved variables to: " + path;
}

std::string StatePersistenceService::load_state(const std::string& path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("unable to open file for reading: " + path);
    }

    std::map<std::string, StoredValue> loaded;
    std::map<std::string, CustomFunction> loaded_functions;
    std::vector<std::string> loaded_script_sources;
    std::string line;
    int state_version = 1;

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        if (line == "STATE_V2") {
            state_version = 2;
            continue;
        }
        if (line == "STATE_V3") {
            state_version = 3;
            continue;
        }
        if (line == "STATE_V4") {
            state_version = 4;
            continue;
        }
        if (line == "STATE_V5") {
            state_version = 5;
            continue;
        }

        const std::vector<std::string> parts = split_tab_fields(line);
        if (state_version >= 2) {
            if (parts.empty()) {
                continue;
            }
            if (parts[0] == "VAR") {
                if (parts.size() < 4) {
                    throw std::runtime_error("invalid save file format");
                }

                StoredValue value;
                const std::string name = decode_field(parts[1]);
                if (parts[2] == "STRING") {
                    if (parts.size() != 4) {
                        throw std::runtime_error("invalid save file format");
                    }
                    value.is_string = true;
                    value.string_value = decode_field(parts[3]);
                } else if (parts[2] == "MATRIX" && state_version >= 4) {
                    if (parts.size() < 5) {
                        throw std::runtime_error("invalid save file format");
                    }
                    const std::size_t rows =
                        static_cast<std::size_t>(std::stoull(parts[3]));
                    const std::size_t cols =
                        static_cast<std::size_t>(std::stoull(parts[4]));
                    if (parts.size() != rows * cols + 5) {
                        throw std::runtime_error("invalid save file format");
                    }
                    value.is_matrix = true;
                    value.matrix_ptr = std::make_shared<matrix::Matrix>(rows, cols, 0.0L);
                    for (std::size_t i = 0; i < value.matrix_ptr->data.size(); ++i) {
                        value.matrix_ptr->data[i] = Scalar(parts[i + 5]);
                    }
                } else if (parts[2] == "COMPLEX" && state_version >= 5) {
                    if (parts.size() != 5) {
                        throw std::runtime_error("invalid save file format");
                    }
                    value.is_complex = true;
                    value.complex.real(Scalar(parts[3]));
                    value.complex.imag(Scalar(parts[4]));
                } else if (parts[2] == "EXACT") {
                    if (parts.size() != 6) {
                        throw std::runtime_error("invalid save file format");
                    }
                    value = StoredValue(Rational(std::stoll(parts[3]), std::stoll(parts[4])));
                } else if (parts[2] == "DECIMAL") {
                    if ((state_version == 2 && parts.size() != 4 && parts.size() != 5) ||
                        (state_version >= 3 && parts.size() != 4)) {
                        throw std::runtime_error("invalid save file format");
                    }
                    value.decimal = Scalar(parts[3]);
                    if (state_version == 2 && parts.size() == 5) {
                        value.has_precise_decimal_text = true;
                        value.precise_decimal_text = decode_field(parts[4]);
                    }
                } else {
                    throw std::runtime_error("invalid save file format");
                }
                loaded[name] = value;
                continue;
            }

            if (state_version >= 3 && parts[0] == "PRECISE") {
                if (parts.size() != 3) {
                    throw std::runtime_error("invalid save file format");
                }
                const std::string name = decode_field(parts[1]);
                auto it = loaded.find(name);
                if (it == loaded.end() || it->second.is_matrix || it->second.is_complex ||
                    it->second.is_string) {
                    throw std::runtime_error("invalid save file format");
                }
                it->second.has_precise_decimal_text = true;
                it->second.precise_decimal_text = decode_field(parts[2]);
                continue;
            }

            if (state_version >= 3 && parts[0] == "SYMBOLIC") {
                if (parts.size() != 3) {
                    throw std::runtime_error("invalid save file format");
                }
                const std::string name = decode_field(parts[1]);
                auto it = loaded.find(name);
                if (it == loaded.end() || it->second.is_matrix || it->second.is_complex ||
                    it->second.is_string) {
                    throw std::runtime_error("invalid save file format");
                }
                it->second.has_symbolic_text = true;
                it->second.symbolic_text = decode_field(parts[2]);
                continue;
            }

            if (parts[0] == "EXPRFUNC") {
                if (parts.size() != 2) {
                    throw std::runtime_error("invalid save file format");
                }
                const std::string definition = decode_field(parts[1]);
                CommandASTNode ast = parse_command(definition);
                if (ast.kind == CommandKind::kFunctionDefinition) {
                    const FunctionDefinitionInfo* def = ast.as_function_definition();
                    if (def) {
                        std::vector<std::string> params;
                        for (auto p : def->parameters) {
                            params.emplace_back(p);
                        }
                        CustomFunction function;
                        function.parameter_names = std::move(params);
                        function.expression = std::string(def->body.text);
                        loaded_functions[std::string(def->name)] = std::move(function);
                    }
                } else {
                    throw std::runtime_error("invalid save file format");
                }
                continue;
            }

            if (parts[0] == "SCRIPT") {
                if (parts.size() != 2) {
                    throw std::runtime_error("invalid save file format");
                }
                loaded_script_sources.push_back(decode_field(parts[1]));
                continue;
            }

            throw std::runtime_error("invalid save file format");
        }

        // V1 格式（已弃用但保持兼容）
        if (parts.size() != 5) {
            throw std::runtime_error("invalid save file format");
        }

        StoredValue value = std::stoi(parts[1]) != 0
            ? StoredValue(Rational(std::stoll(parts[2]), std::stoll(parts[3])))
            : StoredValue(Scalar(parts[4]));
        loaded[parts[0]] = value;
    }

    // 清除并重新加载变量
    variables_->clear_all();
    for (const auto& [name, value] : loaded) {
        variables_->set_global(name, value);
    }

    // 清除并重新加载函数
    functions_->clear_all();
    for (const auto& [name, func] : loaded_functions) {
        functions_->add_custom_function(name, func);
    }
    for (const std::string& script_source : loaded_script_sources) {
        if (script_executor_) {
            script_executor_(script_source);
        }
    }

    return "Loaded variables from: " + path;
}
