#include "calculator_internal_types.h"

#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <utility>

std::string Calculator::save_state(const std::string& path) const {
    for (const auto& [name, value] : impl_->variables) {
        (void)name;
        if (value.is_matrix) {
            throw std::runtime_error("save_state does not yet support matrix variables");
        }
    }

    const std::filesystem::path target_path(path);
    const std::filesystem::path temp_path =
        target_path.parent_path() /
        (target_path.filename().string() + ".tmp-save");
    std::ofstream out(temp_path);
    if (!out) {
        throw std::runtime_error("unable to open file for writing: " + path);
    }

    out << "STATE_V3\n";

    for (const auto& [name, value] : impl_->variables) {
        if (value.is_string) {
            out << "VAR\t" << encode_state_field(name)
                << "\tSTRING\t" << encode_state_field(value.string_value) << '\n';
        } else if (value.exact) {
            out << "VAR\t" << encode_state_field(name)
                << "\tEXACT\t" << value.rational.numerator
                << '\t' << value.rational.denominator
                << '\t' << std::setprecision(17) << value.decimal << '\n';
        } else {
            out << "VAR\t" << encode_state_field(name)
                << "\tDECIMAL\t" << std::setprecision(17) << value.decimal << '\n';
        }
        if (value.has_precise_decimal_text) {
            out << "PRECISE\t" << encode_state_field(name)
                << '\t' << encode_state_field(value.precise_decimal_text) << '\n';
        }
        if (value.has_symbolic_text) {
            out << "SYMBOLIC\t" << encode_state_field(name)
                << '\t' << encode_state_field(value.symbolic_text) << '\n';
        }
    }

    for (const auto& [name, function] : impl_->functions) {
        out << "EXPRFUNC\t"
            << encode_state_field(name + "(" + function.parameter_name + ") = " +
                                  function.expression)
            << '\n';
    }

    for (const auto& [name, function] : impl_->script_functions) {
        std::ostringstream source;
        source << "fn " << name << "(";
        for (std::size_t i = 0; i < function.parameter_names.size(); ++i) {
            if (i != 0) {
                source << ", ";
            }
            source << function.parameter_names[i];
        }
        source << ") " << render_script_block(*function.body, 0);
        out << "SCRIPT\t" << encode_state_field(source.str()) << '\n';
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

std::string Calculator::load_state(const std::string& path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("unable to open file for reading: " + path);
    }

    std::map<std::string, StoredValue> loaded;
    std::map<std::string, CustomFunction> loaded_functions;
    std::map<std::string, ScriptFunction> loaded_script_functions;
    std::string line;
    int state_version = 1;

    auto split_tab_fields = [](const std::string& row_text) {
        std::vector<std::string> parts;
        std::size_t start = 0;
        for (std::size_t i = 0; i <= row_text.size(); ++i) {
            if (i == row_text.size() || row_text[i] == '\t') {
                parts.push_back(row_text.substr(start, i - start));
                start = i + 1;
            }
        }
        return parts;
    };

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
                const std::string name = decode_state_field(parts[1]);
                if (parts[2] == "STRING") {
                    if (parts.size() != 4) {
                        throw std::runtime_error("invalid save file format");
                    }
                    value.is_string = true;
                    value.string_value = decode_state_field(parts[3]);
                } else if (parts[2] == "EXACT") {
                    if (parts.size() != 6) {
                        throw std::runtime_error("invalid save file format");
                    }
                    value.exact = true;
                    value.rational = Rational(std::stoll(parts[3]), std::stoll(parts[4]));
                    value.decimal = std::stod(parts[5]);
                } else if (parts[2] == "DECIMAL") {
                    if ((state_version == 2 && parts.size() != 4 && parts.size() != 5) ||
                        (state_version >= 3 && parts.size() != 4)) {
                        throw std::runtime_error("invalid save file format");
                    }
                    value.decimal = std::stod(parts[3]);
                    if (state_version == 2 && parts.size() == 5) {
                        value.has_precise_decimal_text = true;
                        value.precise_decimal_text = decode_state_field(parts[4]);
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
                const std::string name = decode_state_field(parts[1]);
                auto it = loaded.find(name);
                if (it == loaded.end() || it->second.is_matrix || it->second.is_string) {
                    throw std::runtime_error("invalid save file format");
                }
                it->second.has_precise_decimal_text = true;
                it->second.precise_decimal_text = decode_state_field(parts[2]);
                continue;
            }

            if (state_version >= 3 && parts[0] == "SYMBOLIC") {
                if (parts.size() != 3) {
                    throw std::runtime_error("invalid save file format");
                }
                const std::string name = decode_state_field(parts[1]);
                auto it = loaded.find(name);
                if (it == loaded.end() || it->second.is_matrix || it->second.is_string) {
                    throw std::runtime_error("invalid save file format");
                }
                it->second.has_symbolic_text = true;
                it->second.symbolic_text = decode_state_field(parts[2]);
                continue;
            }

            if (parts[0] == "EXPRFUNC") {
                if (parts.size() != 2) {
                    throw std::runtime_error("invalid save file format");
                }
                std::string function_name;
                std::string parameter_name;
                std::string body;
                const std::string definition = decode_state_field(parts[1]);
                if (!split_function_definition(definition,
                                               &function_name,
                                               &parameter_name,
                                               &body)) {
                    throw std::runtime_error("invalid save file format");
                }
                loaded_functions[function_name] = {parameter_name, body};
                continue;
            }

            if (parts[0] == "SCRIPT") {
                if (parts.size() != 2) {
                    throw std::runtime_error("invalid save file format");
                }
                Calculator temp;
                temp.execute_script(decode_state_field(parts[1]), false);
                for (const auto& [name, function] : temp.impl_->script_functions) {
                    loaded_script_functions[name] = function;
                }
                continue;
            }

            throw std::runtime_error("invalid save file format");
        }

        if (parts.size() != 5) {
            throw std::runtime_error("invalid save file format");
        }

        StoredValue value;
        value.exact = std::stoi(parts[1]) != 0;
        value.rational = Rational(std::stoll(parts[2]), std::stoll(parts[3]));
        value.decimal = std::stod(parts[4]);
        loaded[parts[0]] = value;
    }

    impl_->variables = loaded;
    if (state_version >= 2) {
        impl_->functions = loaded_functions;
        impl_->script_functions = loaded_script_functions;
    }
    return "Loaded variables from: " + path;
}
