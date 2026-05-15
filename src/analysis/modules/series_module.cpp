// ============================================================================
// 级数展开模块命令处理
// ============================================================================
//
// 本文件是级数模块的入口，负责：
// - 解析命令参数
// - 调用 series/ 子目录中的核心算法
// - 格式化输出结果
//
// 核心计算已拆分到：
// - psa_engine.cpp: 幂级数代数运算
// - taylor_series.cpp: Taylor 级数展开
// - pade_approximation.cpp: Pade 有理逼近
// - puiseux_series.cpp: Puiseux 级数
// - series_summation.cpp: 级数求和

#include "analysis/modules/series_module.h"
#include "analysis/series/psa_engine.h"
#include "analysis/series/taylor_series.h"
#include "analysis/series/pade_approximation.h"
#include "analysis/series/puiseux_series.h"
#include "analysis/series/series_summation.h"
#include "app/scalar_type.h"
#include "symbolic/core/symbolic_expression_internal.h"
#include "math/mymath.h"
#include "math/helpers/integer_helpers.h"
#include "parser/grammars/unified_expression_parser.h"
#include "expression_utils.h"
#include "string_utils.h"
#include <sstream>
#include <iomanip>

namespace series_ops {

using Scalar = mymath::Scalar;

// ============================================================================
// 级数展开命令处理
// ============================================================================

bool is_series_command(const std::string& command) {
    return command == "taylor" ||
           command == "pade" ||
           command == "puiseux" ||
           command == "series" ||
           command == "series_sum" ||
           command == "summation";
}

bool handle_series_command(const SeriesContext& ctx,
                           const std::string& command,
                           const std::vector<std::string>& arguments,
                           std::string* output) {

    if (command == "taylor") {
        if (arguments.size() != 3) {
            throw std::runtime_error("taylor expects exactly three arguments");
        }

        const Scalar center = ctx.parse_decimal(arguments[1]);
        const Scalar degree_value = ctx.parse_decimal(arguments[2]);
        if (!is_integer_double(degree_value) || degree_value < 0.0L) {
            throw std::runtime_error("taylor degree must be a non-negative integer");
        }
        const int degree = static_cast<int>(round_to_long_long(degree_value));

        *output = taylor::taylor(ctx, arguments[0], center, degree);
        return true;
    }

    if (command == "pade") {
        if (arguments.size() != 3 && arguments.size() != 4) {
            throw std::runtime_error(
                "pade expects expr, m, n or expr, center, m, n");
        }

        const bool explicit_center = arguments.size() == 4;
        const Scalar center = explicit_center
                                  ? ctx.parse_decimal(arguments[1])
                                  : 0.0L;
        const Scalar numerator_degree_value = ctx.parse_decimal(
            arguments[explicit_center ? 2 : 1]);
        const Scalar denominator_degree_value = ctx.parse_decimal(
            arguments[explicit_center ? 3 : 2]);
        if (!is_integer_double(numerator_degree_value) ||
            numerator_degree_value < 0.0L ||
            !is_integer_double(denominator_degree_value) ||
            denominator_degree_value < 0.0L) {
            throw std::runtime_error(
                "pade degrees must be non-negative integers");
        }

        const int numerator_degree =
            static_cast<int>(round_to_long_long(numerator_degree_value));
        const int denominator_degree =
            static_cast<int>(round_to_long_long(denominator_degree_value));

        *output = pade::pade(ctx, arguments[0], center, numerator_degree, denominator_degree);
        return true;
    }

    if (command == "series") {
        if (arguments.size() != 2 && arguments.size() != 3) {
            throw std::runtime_error("series expects expr, degree or expr, center, degree");
        }
        Scalar center = 0.0L;
        int degree = 0;
        if (arguments.size() == 2) {
            degree = static_cast<int>(ctx.parse_decimal(arguments[1]));
        } else {
            center = ctx.parse_decimal(arguments[1]);
            degree = static_cast<int>(ctx.parse_decimal(arguments[2]));
        }
        if (degree < 0) throw std::runtime_error("series degree must be non-negative");
        *output = puiseux::general_series_auto(ctx, arguments[0], center, degree);
        return true;
    }

    if (command == "puiseux") {
        if (arguments.size() != 3 && arguments.size() != 4) {
            throw std::runtime_error(
                "puiseux expects expr, degree, denominator or expr, center, degree, denominator");
        }

        const bool explicit_center = arguments.size() == 4;
        const Scalar center = explicit_center
                                  ? ctx.parse_decimal(arguments[1])
                                  : 0.0L;
        const Scalar degree_value = ctx.parse_decimal(
            arguments[explicit_center ? 2 : 1]);
        const Scalar denominator_value = ctx.parse_decimal(
            arguments[explicit_center ? 3 : 2]);
        if (!is_integer_double(degree_value) || degree_value < 0.0L) {
            throw std::runtime_error(
                "puiseux degree must be a non-negative integer");
        }
        if (!is_integer_double(denominator_value) || denominator_value <= 0.0L) {
            throw std::runtime_error(
                "puiseux denominator must be a positive integer");
        }

        const int degree = static_cast<int>(round_to_long_long(degree_value));
        const int denom = static_cast<int>(round_to_long_long(denominator_value));

        *output = puiseux::puiseux(ctx, arguments[0], center, degree, denom);
        return true;
    }

    if (command == "series_sum" || command == "summation") {
        if (arguments.size() != 4) {
            throw std::runtime_error(
                "series_sum expects expr, index, lower, upper");
        }

        const std::string index_name = trim_copy(arguments[1]);
        if (!is_identifier_text(index_name)) {
            throw std::runtime_error("series_sum index must be an identifier");
        }

        *output = summation::series_sum(ctx,
                             arguments[0],
                             index_name,
                             arguments[2],
                             trim_copy(arguments[3]));
        return true;
    }

    return false;
}


std::string SeriesModule::execute_args(const std::string& command,
                                      const std::vector<std::string>& args,
                                      const CoreServices& services) {
    SeriesContext ctx;
    ctx.resolve_symbolic = services.symbolic.resolve_symbolic;
    ctx.parse_decimal = services.evaluation.parse_decimal;
    ctx.evaluate_at = services.symbolic.evaluate_symbolic_at;
    ctx.simplify_symbolic = services.symbolic.simplify_symbolic;
    ctx.expand_inline = services.symbolic.expand_inline;

    std::string output;
    if (handle_series_command(ctx, command, args, &output)) {
        return output;
    }
    throw std::runtime_error("Series command failed: " + command);
}

std::vector<std::string> SeriesModule::get_commands() const {
    return {"taylor", "pade", "puiseux", "series", "series_sum", "summation"};
}

std::string SeriesModule::get_help_snippet(const std::string& topic) const {
    if (topic == "symbolic") {
        return "Series:\n"
               "  taylor(f, a, n)     Taylor series at x=a up to degree n\n"
               "  pade(f, [a], m, n)  Pade approximation [m/n]\n"
               "  series(f, [a], n)   Auto-detecting series (Laurent/Puiseux)\n"
               "  puiseux(f, n, d)    Puiseux series (fractional powers)\n"
               "  series_sum(f, i, a, b) Finite or infinite summation";
    }
    return "";
}

}  // namespace series_ops
#include "module/module_registration.h"
REGISTER_CALCULATOR_MODULE(series_ops::SeriesModule)
