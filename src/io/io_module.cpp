#include "io_module.h"
#include "core/calculator_exceptions.h"
#include "core/string_utils.h"
#include "matrix/matrix.h"

#include <filesystem>
#include <sstream>
#include <stdexcept>

namespace {
    double get_scalar(const StoredValue& val, const std::string& ctx) {
        if (val.is_matrix || val.is_complex || val.is_string || val.is_list || val.is_dict) {
            throw std::runtime_error(ctx + " must be a scalar");
        }
        return val.exact ? rational_to_double(val.rational) : val.decimal;
    }

    std::string get_string(const StoredValue& val, const std::string& ctx) {
        if (val.is_string) return val.string_value;
        throw std::runtime_error(ctx + " must be a string");
    }

    matrix::Matrix get_matrix(const StoredValue& val, const std::string& ctx) {
        if (val.is_matrix) return val.matrix;
        throw std::runtime_error(ctx + " must be a matrix");
    }

    // JSON 序列化辅助函数
    std::string value_to_json_impl(const StoredValue& val) {
        if (val.is_string) {
            std::string escaped;
            for (char c : val.string_value) {
                switch (c) {
                    case '"': escaped += "\\\""; break;
                    case '\\': escaped += "\\\\"; break;
                    case '\n': escaped += "\\n"; break;
                    case '\r': escaped += "\\r"; break;
                    case '\t': escaped += "\\t"; break;
                    default: escaped += c;
                }
            }
            return "\"" + escaped + "\"";
        }
        if (val.is_matrix) {
            std::string result = "[";
            for (std::size_t r = 0; r < val.matrix.rows; ++r) {
                if (r > 0) result += ", ";
                if (val.matrix.cols > 1) result += "[";
                for (std::size_t c = 0; c < val.matrix.cols; ++c) {
                    if (c > 0) result += ", ";
                    result += std::to_string(val.matrix.at(r, c));
                }
                if (val.matrix.cols > 1) result += "]";
            }
            result += "]";
            return result;
        }
        if (val.is_list && val.list_value) {
            std::string result = "[";
            for (std::size_t i = 0; i < val.list_value->size(); ++i) {
                if (i > 0) result += ", ";
                result += value_to_json_impl((*val.list_value)[i]);
            }
            result += "]";
            return result;
        }
        if (val.is_dict && val.dict_value) {
            std::string result = "{";
            bool first = true;
            for (const auto& [key, value] : *val.dict_value) {
                if (!first) result += ", ";
                first = false;
                result += "\"" + key + "\": " + value_to_json_impl(value);
            }
            result += "}";
            return result;
        }
        if (val.is_complex) {
            return "[" + std::to_string(val.complex.real) + ", " + std::to_string(val.complex.imag) + "]";
        }
        // 标量
        return std::to_string(val.exact ? rational_to_double(val.rational) : val.decimal);
    }

    // JSON 解析辅助函数
    StoredValue parse_json_value_impl(std::string_view& s) {
        // 跳过空白
        while (!s.empty() && std::isspace(static_cast<unsigned char>(s.front()))) {
            s.remove_prefix(1);
        }
        if (s.empty()) throw std::runtime_error("Unexpected end of JSON");

        StoredValue result;

        // null
        if (s.substr(0, 4) == "null") {
            s.remove_prefix(4);
            result.decimal = 0;
            return result;
        }

        // true
        if (s.substr(0, 4) == "true") {
            s.remove_prefix(4);
            result.decimal = 1;
            return result;
        }

        // false
        if (s.substr(0, 5) == "false") {
            s.remove_prefix(5);
            result.decimal = 0;
            return result;
        }

        // 数字
        if (std::isdigit(static_cast<unsigned char>(s.front())) || s.front() == '-' || s.front() == '+') {
            std::size_t i = 0;
            if (s.front() == '-' || s.front() == '+') i = 1;
            while (i < s.size() && (std::isdigit(static_cast<unsigned char>(s[i])) || s[i] == '.' || s[i] == 'e' || s[i] == 'E' || s[i] == '+' || s[i] == '-')) {
                ++i;
            }
            result.decimal = std::stod(std::string(s.substr(0, i)));
            result.exact = false;
            s.remove_prefix(i);
            return result;
        }

        // 字符串
        if (s.front() == '"') {
            s.remove_prefix(1);
            std::string str;
            while (!s.empty() && s.front() != '"') {
                if (s.front() == '\\' && s.size() > 1) {
                    s.remove_prefix(1);
                    switch (s.front()) {
                        case 'n': str += '\n'; break;
                        case 'r': str += '\r'; break;
                        case 't': str += '\t'; break;
                        case '"': str += '"'; break;
                        case '\\': str += '\\'; break;
                        default: str += s.front();
                    }
                } else {
                    str += s.front();
                }
                s.remove_prefix(1);
            }
            if (!s.empty()) s.remove_prefix(1);  // 跳过结束引号
            result.is_string = true;
            result.string_value = str;
            return result;
        }

        // 数组
        if (s.front() == '[') {
            s.remove_prefix(1);
            std::vector<StoredValue> items;
            while (!s.empty() && s.front() != ']') {
                while (!s.empty() && std::isspace(static_cast<unsigned char>(s.front()))) s.remove_prefix(1);
                if (s.empty() || s.front() == ']') break;
                items.push_back(parse_json_value_impl(s));
                while (!s.empty() && std::isspace(static_cast<unsigned char>(s.front()))) s.remove_prefix(1);
                if (!s.empty() && s.front() == ',') s.remove_prefix(1);
            }
            if (!s.empty()) s.remove_prefix(1);  // 跳过 ]

            // 检查是否为矩阵
            // 情况1：所有元素都是数字 -> 列向量
            // 情况2：所有元素都是相同长度的数字列表 -> 2D矩阵
            bool is_numeric = true;
            for (const auto& item : items) {
                if (item.is_string || item.is_list || item.is_dict || item.is_matrix) {
                    is_numeric = false;
                    break;
                }
            }

            if (is_numeric && !items.empty()) {
                // 检查这是否可能是矩阵的行（延迟决定，先保存为列表）
                result.is_list = true;
                result.list_value = std::make_shared<std::vector<StoredValue>>(std::move(items));
            } else if (!items.empty()) {
                // 检查是否为2D矩阵（所有元素都是相同长度的数字列表）
                bool is_2d_matrix = true;
                std::size_t cols = 0;
                for (const auto& item : items) {
                    if (item.is_list && item.list_value) {
                        // 检查列表元素是否都是数字
                        for (const auto& sub : *item.list_value) {
                            if (sub.is_string || sub.is_list || sub.is_dict || sub.is_matrix) {
                                is_2d_matrix = false;
                                break;
                            }
                        }
                        if (!is_2d_matrix) break;
                        if (cols == 0) cols = item.list_value->size();
                        else if (item.list_value->size() != cols) { is_2d_matrix = false; break; }
                    } else {
                        is_2d_matrix = false;
                        break;
                    }
                }

                if (is_2d_matrix && cols > 0) {
                    // 转换为矩阵
                    matrix::Matrix mat(items.size(), cols, 0.0);
                    for (std::size_t r = 0; r < items.size(); ++r) {
                        const auto& item = items[r];
                        for (std::size_t c = 0; c < item.list_value->size(); ++c) {
                            mat.at(r, c) = (*item.list_value)[c].decimal;
                        }
                    }
                    result.is_matrix = true;
                    result.matrix = mat;
                } else {
                    result.is_list = true;
                    result.list_value = std::make_shared<std::vector<StoredValue>>(std::move(items));
                }
            } else {
                result.is_list = true;
                result.list_value = std::make_shared<std::vector<StoredValue>>(std::move(items));
            }
            return result;
        }

        // 对象
        if (s.front() == '{') {
            s.remove_prefix(1);
            std::map<std::string, StoredValue> obj;
            while (!s.empty() && s.front() != '}') {
                while (!s.empty() && std::isspace(static_cast<unsigned char>(s.front()))) s.remove_prefix(1);
                if (s.empty() || s.front() == '}') break;

                // 解析键
                if (s.front() != '"') throw std::runtime_error("Expected string key in JSON object");
                s.remove_prefix(1);
                std::string key;
                while (!s.empty() && s.front() != '"') {
                    if (s.front() == '\\' && s.size() > 1) {
                        s.remove_prefix(1);
                        key += s.front();
                    } else {
                        key += s.front();
                    }
                    s.remove_prefix(1);
                }
                if (!s.empty()) s.remove_prefix(1);  // 跳过结束引号

                // 跳过冒号
                while (!s.empty() && std::isspace(static_cast<unsigned char>(s.front()))) s.remove_prefix(1);
                if (!s.empty() && s.front() == ':') s.remove_prefix(1);
                while (!s.empty() && std::isspace(static_cast<unsigned char>(s.front()))) s.remove_prefix(1);

                // 解析值
                obj[key] = parse_json_value_impl(s);

                // 跳过逗号
                while (!s.empty() && std::isspace(static_cast<unsigned char>(s.front()))) s.remove_prefix(1);
                if (!s.empty() && s.front() == ',') s.remove_prefix(1);
            }
            if (!s.empty()) s.remove_prefix(1);  // 跳过 }
            result.is_dict = true;
            result.dict_value = std::make_shared<std::map<std::string, StoredValue>>(std::move(obj));
            return result;
        }

        throw std::runtime_error("Unexpected character in JSON");
    }
}

std::map<std::string, std::function<StoredValue(const std::vector<StoredValue>&)>> IoModule::get_native_functions() const {
    std::map<std::string, std::function<StoredValue(const std::vector<StoredValue>&)>> funcs;

    funcs["open"] = [this](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.empty()) throw std::runtime_error("open expects at least 1 argument");
        std::string path;
        if (args[0].is_string) {
            path = args[0].string_value;
        } else {
            throw std::runtime_error("open expects a string path");
        }
        
        std::string mode = "r";
        if (args.size() > 1) {
            if (!args[1].is_string) throw std::runtime_error("open mode must be a string");
            mode = args[1].string_value;
        }
        
        std::ios_base::openmode std_mode = std::ios_base::in;
        if (mode == "w") std_mode = std::ios_base::out | std::ios_base::trunc;
        else if (mode == "a") std_mode = std::ios_base::out | std::ios_base::app;
        else if (mode == "rw" || mode == "r+") std_mode = std::ios_base::in | std::ios_base::out;
        else if (mode != "r") throw std::runtime_error("Invalid open mode: " + mode);

        auto fs = std::make_shared<std::fstream>(path, std_mode);
        if (!fs->is_open()) {
            throw std::runtime_error("Failed to open file: " + path);
        }
        
        int fd = next_fd_++;
        files_[fd] = fs;
        
        StoredValue res;
        res.decimal = fd;
        res.exact = false;
        return res;
    };

    funcs["close"] = [this](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.empty()) throw std::runtime_error("close expects 1 argument");
        int fd = static_cast<int>(get_scalar(args[0], "close fd"));
        auto it = files_.find(fd);
        if (it != files_.end()) {
            it->second->close();
            files_.erase(it);
            StoredValue res; res.decimal = 1; return res;
        }
        throw std::runtime_error("Invalid file descriptor");
    };

    funcs["read"] = [this](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.empty()) throw std::runtime_error("read expects at least 1 argument");
        int fd = static_cast<int>(get_scalar(args[0], "read fd"));
        auto it = files_.find(fd);
        if (it == files_.end()) throw std::runtime_error("Invalid file descriptor");
        
        std::ostringstream ss;
        ss << it->second->rdbuf();
        
        StoredValue res;
        res.is_string = true;
        res.string_value = ss.str();
        return res;
    };

    funcs["write"] = [this](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() < 2) throw std::runtime_error("write expects 2 arguments");
        int fd = static_cast<int>(get_scalar(args[0], "write fd"));
        auto it = files_.find(fd);
        if (it == files_.end()) throw std::runtime_error("Invalid file descriptor");
        
        std::string content;
        if (args[1].is_string) content = args[1].string_value;
        else content = std::to_string(get_scalar(args[1], "write data")); // Basic fallback
        
        *(it->second) << content;
        it->second->flush();
        
        StoredValue res; res.decimal = 1; return res;
    };

    funcs["read_lines"] = [this](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.empty()) throw std::runtime_error("read_lines expects 1 argument");
        int fd = static_cast<int>(get_scalar(args[0], "read_lines fd"));
        auto it = files_.find(fd);
        if (it == files_.end()) throw std::runtime_error("Invalid file descriptor");

        std::vector<StoredValue> lines;
        std::string line;
        while (std::getline(*(it->second), line)) {
            StoredValue lval;
            lval.is_string = true;
            lval.string_value = line;
            lines.push_back(lval);
        }

        StoredValue res;
        res.is_list = true;
        res.list_value = std::make_shared<std::vector<StoredValue>>(std::move(lines));
        return res;
    };

    // readline(fd) - 读取单行
    funcs["readline"] = [this](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.empty()) throw std::runtime_error("readline expects 1 argument");
        int fd = static_cast<int>(get_scalar(args[0], "readline fd"));
        auto it = files_.find(fd);
        if (it == files_.end()) throw std::runtime_error("Invalid file descriptor");

        std::string line;
        if (std::getline(*(it->second), line)) {
            StoredValue res;
            res.is_string = true;
            res.string_value = line;
            return res;
        }

        // EOF reached
        StoredValue res;
        res.is_string = true;
        res.string_value = "";
        return res;
    };

    // seek(fd, pos) - 定位文件指针
    funcs["seek"] = [this](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() < 2) throw std::runtime_error("seek expects 2 arguments");
        int fd = static_cast<int>(get_scalar(args[0], "seek fd"));
        std::streampos pos = static_cast<std::streampos>(get_scalar(args[1], "seek position"));

        auto it = files_.find(fd);
        if (it == files_.end()) throw std::runtime_error("Invalid file descriptor");

        it->second->seekg(pos);
        it->second->seekp(pos);

        StoredValue res;
        res.decimal = 1;
        res.exact = false;
        return res;
    };

    // tell(fd) - 获取当前位置
    funcs["tell"] = [this](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.empty()) throw std::runtime_error("tell expects 1 argument");
        int fd = static_cast<int>(get_scalar(args[0], "tell fd"));

        auto it = files_.find(fd);
        if (it == files_.end()) throw std::runtime_error("Invalid file descriptor");

        std::streampos pos = it->second->tellg();

        StoredValue res;
        res.decimal = static_cast<double>(pos);
        res.exact = false;
        return res;
    };

    // exists(path) - 检查文件是否存在
    funcs["exists"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.empty()) throw std::runtime_error("exists expects 1 argument");
        std::string path = get_string(args[0], "exists path");

        StoredValue res;
        res.decimal = std::filesystem::exists(path) ? 1.0 : 0.0;
        res.exact = false;
        return res;
    };

    // delete(path) - 删除文件
    funcs["delete"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.empty()) throw std::runtime_error("delete expects 1 argument");
        std::string path = get_string(args[0], "delete path");

        bool removed = std::filesystem::remove(path);

        StoredValue res;
        res.decimal = removed ? 1.0 : 0.0;
        res.exact = false;
        return res;
    };

    // read_csv(path) - 从 CSV 文件读取矩阵
    funcs["read_csv"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.empty()) throw std::runtime_error("read_csv expects 1 argument");
        std::string path = get_string(args[0], "read_csv path");

        std::ifstream file(path);
        if (!file.is_open()) {
            throw std::runtime_error("Failed to open file: " + path);
        }

        std::vector<std::vector<double>> rows;
        std::string line;

        while (std::getline(file, line)) {
            if (line.empty()) continue;

            std::vector<double> row;
            std::stringstream ss(line);
            std::string cell;

            while (std::getline(ss, cell, ',')) {
                try {
                    row.push_back(std::stod(trim_copy(cell)));
                } catch (...) {
                    row.push_back(0.0);
                }
            }

            if (!row.empty()) {
                rows.push_back(std::move(row));
            }
        }

        if (rows.empty()) {
            StoredValue res;
            res.is_matrix = true;
            res.matrix = matrix::Matrix(0, 0, 0.0);
            return res;
        }

        std::size_t cols = rows[0].size();
        matrix::Matrix result(rows.size(), cols, 0.0);

        for (std::size_t r = 0; r < rows.size(); ++r) {
            for (std::size_t c = 0; c < std::min(rows[r].size(), cols); ++c) {
                result.at(r, c) = rows[r][c];
            }
        }

        StoredValue res;
        res.is_matrix = true;
        res.matrix = result;
        return res;
    };

    // write_csv(path, matrix) - 将矩阵写入 CSV 文件
    funcs["write_csv"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() < 2) throw std::runtime_error("write_csv expects 2 arguments");
        std::string path = get_string(args[0], "write_csv path");
        matrix::Matrix mat = get_matrix(args[1], "write_csv matrix");

        std::ofstream file(path);
        if (!file.is_open()) {
            throw std::runtime_error("Failed to open file for writing: " + path);
        }

        for (std::size_t r = 0; r < mat.rows; ++r) {
            for (std::size_t c = 0; c < mat.cols; ++c) {
                if (c > 0) file << ",";
                file << mat.at(r, c);
            }
            file << "\n";
        }

        StoredValue res;
        res.decimal = 1;
        res.exact = false;
        return res;
    };

    // read_json(path) - 从 JSON 文件读取数据
    funcs["read_json"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.empty()) throw std::runtime_error("read_json expects 1 argument");
        std::string path = get_string(args[0], "read_json path");

        std::ifstream file(path);
        if (!file.is_open()) {
            throw std::runtime_error("Failed to open file: " + path);
        }

        std::stringstream ss;
        ss << file.rdbuf();
        std::string content = ss.str();

        std::string_view s = content;
        return parse_json_value_impl(s);
    };

    // write_json(path, data) - 将数据写入 JSON 文件
    funcs["write_json"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() < 2) throw std::runtime_error("write_json expects 2 arguments");
        std::string path = get_string(args[0], "write_json path");

        std::ofstream file(path);
        if (!file.is_open()) {
            throw std::runtime_error("Failed to open file for writing: " + path);
        }

        file << value_to_json_impl(args[1]);

        StoredValue res;
        res.decimal = 1;
        res.exact = false;
        return res;
    };

    return funcs;
}

std::string IoModule::execute_args(const std::string& command,
                                   const std::vector<std::string>& args,
                                   const CoreServices& services) {
    // For command-line usage without parentheses
    // e.g. open "file.txt" w
    std::vector<StoredValue> s_args;
    for (const auto& arg : args) {
        // Strip quotes if any
        std::string parsed = trim_copy(arg);
        if (parsed.size() >= 2 && parsed.front() == '"' && parsed.back() == '"') {
            StoredValue sv;
            sv.is_string = true;
            sv.string_value = parsed.substr(1, parsed.size() - 2);
            s_args.push_back(sv);
        } else {
            // 使用 evaluate_value 而不是 parse_decimal，以支持矩阵和复数变量
            StoredValue sv = services.evaluation.evaluate_value(parsed, false);
            s_args.push_back(sv);
        }
    }
    
    auto funcs = get_native_functions();
    auto it = funcs.find(command);
    if (it != funcs.end()) {
        StoredValue res = it->second(s_args);
        if (res.is_string) return res.string_value;
        if (res.is_list) {
            std::string out = "[";
            for (size_t i = 0; i < res.list_value->size(); ++i) {
                if (i > 0) out += ", ";
                out += "\"" + (*res.list_value)[i].string_value + "\"";
            }
            out += "]";
            return out;
        }
        return std::to_string(res.exact ? rational_to_double(res.rational) : res.decimal);
    }
    
    return "";
}

std::string IoModule::get_help_snippet(const std::string& topic) const {
    if (topic == "io" || topic == "file") {
        return "File I/O Functions:\n"
               "  open(path, [mode])    - Open file, returns fd (modes: r, w, a, rw)\n"
               "  close(fd)             - Close file\n"
               "  read(fd)              - Read entire content\n"
               "  read_lines(fd)        - Read all lines into list\n"
               "  readline(fd)          - Read single line\n"
               "  write(fd, text)       - Write text to file\n"
               "  seek(fd, pos)         - Set file position\n"
               "  tell(fd)              - Get current position\n"
               "  exists(path)          - Check if file exists (returns 1/0)\n"
               "  delete(path)          - Delete file\n"
               "  read_csv(path)        - Read matrix from CSV\n"
               "  write_csv(path, mat)  - Write matrix to CSV";
    }
    if (topic == "open") return "open(path, [mode]) - Open a file (modes: r, w, a, rw)";
    if (topic == "close") return "close(fd) - Close a file";
    if (topic == "read") return "read(fd) - Read entire file content";
    if (topic == "write") return "write(fd, text) - Write text to file";
    if (topic == "read_lines") return "read_lines(fd) - Read all lines into a list";
    if (topic == "readline") return "readline(fd) - Read single line";
    if (topic == "seek") return "seek(fd, pos) - Set file position";
    if (topic == "tell") return "tell(fd) - Get current position";
    if (topic == "exists") return "exists(path) - Check if file exists (returns 1 or 0)";
    if (topic == "delete") return "delete(path) - Delete a file";
    if (topic == "read_csv") return "read_csv(path) - Read matrix from CSV file";
    if (topic == "write_csv") return "write_csv(path, matrix) - Write matrix to CSV file";
    return "";
}
