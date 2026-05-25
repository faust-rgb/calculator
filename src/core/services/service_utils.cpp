// ============================================================================
// service_utils.cpp - 服务工具函数实现
// ============================================================================
//
// 实现 service_utils 命名空间中的工具函数。
// ============================================================================

#include "core/services/service_utils.h"
#include "core/services/string_utils.h"
#include <cmath>

namespace service_utils {

std::vector<std::string> parse_symbolic_vars(
    const std::vector<std::string>& arguments,
    std::size_t start_index,
    const std::vector<std::string>& defaults) {

    std::vector<std::string> vars;

    // 如果参数中有显式指定的变量列表
    if (arguments.size() > start_index) {
        for (std::size_t i = start_index; i < arguments.size(); ++i) {
            std::string var = arguments[i];
            // 移除逗号分隔符
            if (!var.empty() && var.back() == ',') {
                var.pop_back();
            }
            if (!var.empty()) {
                vars.push_back(var);
            }
        }
    } else {
        // 使用默认变量列表
        vars = defaults;
    }

    return vars;
}

bool is_matrix_argument(const std::string& arg) {
    // 检查是否以 '[' 开始（矩阵格式）
    if (arg.empty()) return false;
    return arg[0] == '[';
}

// parse_matrix_argument is provided by CoreServices

bool is_integer_double(Scalar value, Scalar epsilon) {
    const Scalar abs_val = value < Scalar(0) ? -value : value;
    const Scalar nearest = Scalar(std::round(static_cast<long double>(abs_val)));
    const Scalar diff = abs_val - nearest;
    return (diff < Scalar(0) ? -diff : diff) < epsilon;
}

long long round_to_long_long(Scalar value) {
    return static_cast<long long>(std::round(static_cast<long double>(value)));
}

bool is_near_zero(Scalar value, Scalar epsilon) {
    return (value < Scalar(0) ? -value : value) < epsilon;
}

} // namespace service_utils