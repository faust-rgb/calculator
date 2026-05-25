/**
 * @file dsp_module.cpp
 * @brief 数字信号处理模块实现
 *
 * 本文件实现了 DSP 模块的命令处理和函数注册：
 * - residue: 计算有理函数在指定点的留数（命令形式）
 * - filter: 数字滤波器（函数形式）
 * - freqz: 频率响应计算（函数形式）
 * - residue: 部分分式展开（函数形式，矩阵版）
 *
 * 该模块将命令路由到具体的信号处理函数实现，
 * 并将矩阵级 DSP 函数注册为全局可用函数。
 */

#include "dsp_module.h"
#include "core/services/service_locator.h"
#include "core/services/core_manager_interfaces.h"
#include "residue.h"
#include "matrix_dsp.h"
#include "core/services/string_utils.h"
#include <stdexcept>

namespace {

using matrix::Matrix;
using matrix::filter;
using matrix::freqz;
using matrix::residue;

Matrix require_matrix(const StoredValue& val, const std::string& func_name) {
    if (!val.is_matrix || !val.matrix_ptr) {
        throw std::runtime_error(func_name + " expects a matrix argument");
    }
    return *val.matrix_ptr;
}

} // namespace

// ============================================================================
// 命令接口实现
// ============================================================================

std::vector<std::string> DspModule::get_commands() const {
    return {"residue"};
}


std::string DspModule::execute_args_view(std::string_view command,
                                    const std::vector<std::string_view>& args,
                                    ServiceLocator& locator) {
    using namespace module_helpers;
    std::vector<std::string> string_args;
    for (const auto& arg : args) string_args.emplace_back(arg);

    return dsp_ops::handle_residue_command(std::string(command), string_args, locator);
}

// ============================================================================
// 函数接口实现 — 矩阵级 DSP 函数注册
// ============================================================================

std::map<std::string, std::function<StoredValue(const std::vector<StoredValue>&)>>
DspModule::get_functions_map() const {
    std::map<std::string, std::function<StoredValue(const std::vector<StoredValue>&)>> funcs;

    // 数字滤波器
    funcs["filter"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() != 3) throw std::runtime_error("filter expects 3 arguments");
        Matrix b = require_matrix(args[0], "filter");
        Matrix a = require_matrix(args[1], "filter");
        Matrix x = require_matrix(args[2], "filter");
        StoredValue res;
        res.is_matrix = true;
        res.matrix_ptr = std::make_shared<Matrix>(filter(b, a, x));
        return res;
    };

    // 频率响应
    funcs["freqz"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() < 2 || args.size() > 3) throw std::runtime_error("freqz expects 2 or 3 arguments");
        Matrix b = require_matrix(args[0], "freqz");
        Matrix a = require_matrix(args[1], "freqz");
        std::size_t n = 512;
        if (args.size() == 3) {
            Matrix nm = require_matrix(args[2], "freqz");
            if (nm.rows > 0 && nm.cols > 0) {
                n = static_cast<std::size_t>(static_cast<long double>(nm.at(0, 0)) + 0.5);
            }
        }
        StoredValue res;
        res.is_matrix = true;
        res.matrix_ptr = std::make_shared<Matrix>(freqz(b, a, n));
        return res;
    };

    // 部分分式展开（矩阵版，与命令版 residue 不同）
    funcs["residue_pf"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() != 2) throw std::runtime_error("residue_pf expects 2 arguments");
        Matrix b = require_matrix(args[0], "residue_pf");
        Matrix a = require_matrix(args[1], "residue_pf");
        StoredValue res;
        res.is_matrix = true;
        res.matrix_ptr = std::make_shared<Matrix>(residue(b, a));
        return res;
    };

    return funcs;
}

std::vector<std::string> DspModule::get_function_names() const {
    return {"filter", "freqz", "residue_pf"};
}

std::string DspModule::get_help_snippet(const std::string& topic) const {
    if (topic == "dsp") {
        return "DSP commands:\n"
               "  residue(expr, var, point) - Compute residue of rational function at point\n"
               "DSP functions:\n"
               "  filter(b, a, x) - Digital filter\n"
               "  freqz(b, a, n) - Frequency response\n"
               "  residue_pf(b, a) - Partial fraction expansion";
    }
    return "";
}

REGISTER_CALCULATOR_MODULE(DspModule)