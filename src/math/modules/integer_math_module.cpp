/**
 * @file integer_math_module.cpp
 * @brief 整数数学模块实现
 *
 * 本文件实现了 IntegerMathModule 类，提供整数数学、数论和进制转换功能。
 * 包括因式分解、进制转换、位运算和数论函数等。
 */

#include "math/modules/integer_math_module.h"
#include "core/services/core_manager_interfaces.h"
#include "core/services/service_locator.h"
#include "math/helpers/integer_helpers.h"
#include "math/helpers/combinatorics.h"
#include "math/helpers/bitwise_helpers.h"
#include "math/helpers/base_conversions.h"
#include "core/common/calculator_exceptions.h"
#include "parser/grammars/command_parser.h"
#include "math/mymath.h"
#include <algorithm>
#include <sstream>
#include <set>

namespace {

/**
 * @brief 确保参数为整数并转换为 long long
 * @param x 待检查的值
 * @param name 参数名称（用于错误消息）
 * @param func 函数名称（用于错误消息）
 * @return 转换后的整数值
 * @throws MathError 如果值不是整数
 */
long long require_integer(Scalar x, const std::string& name, const std::string& func) {
    if (!mymath::is_integer(x)) {
        throw MathError(func + " requires an integer " + name);
    }
    return round_to_long_long(static_cast<long double>(x)); // round_to_long_long only takes long double
}

} // namespace


/**
 * @brief 执行命令式操作（如 factor、bin、oct、hex、base）
 * @param node 命令 AST 节点
 * @param locator 服务定位器
 * @return 执行结果的字符串表示
 * @throws std::runtime_error 如果参数无效
 */
std::string IntegerMathModule::execute_command(const CommandASTNode& node,
                                              ServiceLocator& locator) {
    // 提取命令名和参数
    std::string command;
    std::vector<std::string> args;

    if (node.kind == CommandKind::kMetaCommand) {
        command = ":" + std::string(node.as_meta_command()->command);
        for (const auto& arg : node.as_meta_command()->arguments) {
            if (arg->kind == CommandKind::kExpression && arg->as_expression()) {
                args.push_back(std::string(arg->as_expression()->text));
            }
        }
    } else if (node.kind == CommandKind::kFunctionCall) {
        command = std::string(node.as_function_call()->name);
        for (const auto& arg : node.as_function_call()->arguments) {
            if (arg->kind == CommandKind::kExpression && arg->as_expression()) {
                args.push_back(std::string(arg->as_expression()->text));
            }
        }
    } else {
        throw std::runtime_error("Invalid command node type");
    }

    auto engine = locator.resolve<IEvaluationEngine>();

    if (command == "factor") {
        if (args.empty()) throw std::runtime_error("factor expects 1 argument");
        Scalar val = engine->parse_decimal(args[0]);
        return factor_integer(require_integer(val, "argument", "factor"));
    }

    if (command == "bin" || command == "oct" || command == "hex" || command == "base") {
        if (args.empty()) throw std::runtime_error(command + " expects at least 1 argument");

        Scalar value = engine->parse_decimal(args[0]);
        int base = 10;

        if (command == "bin") base = 2;
        else if (command == "oct") base = 8;
        else if (command == "hex") base = 16;
        else {
            if (args.size() < 2) throw std::runtime_error("base expects 2 arguments: value, base");
            base = static_cast<int>(require_integer(engine->parse_decimal(args[1]), "base", "base"));
        }

        if (base < 2 || base > 36) throw std::runtime_error("base must be in range [2, 36]");

        // 获取 hex 格式配置
        auto config = locator.resolve<IConfigManager>();
        bool uppercase = config->is_hex_uppercase_mode();
        bool prefix = config->is_hex_prefix_mode();

        long long n = require_integer(value, "value", command);
        return convert_to_base(n, base, uppercase, prefix);
    }

    return "";
}

/**
 * @brief 获取模块提供的标量函数映射
 * @return 函数名称到函数实现的映射
 *
 * 包括数论函数（gcd、lcm、factorial、nCr、nPr等）、
 * 位运算函数（and、or、xor、shl、shr、rol、ror等）
 * 以及位统计函数（popcount、bitlen、ctz、clz等）。
 */
std::map<std::string, std::function<Scalar(const std::vector<Scalar>&)>>
IntegerMathModule::get_scalar_functions() const {
    std::map<std::string, std::function<Scalar(const std::vector<Scalar>&)>> funcs;

    // Number Theory
    funcs["gcd"] = [](const std::vector<Scalar>& a) { if(a.size()!=2) throw MathError("gcd expects 2 arguments"); return Scalar(gcd_ll(require_integer(a[0], "a", "gcd"), require_integer(a[1], "b", "gcd"))); };
    funcs["lcm"] = [](const std::vector<Scalar>& a) { if(a.size()!=2) throw MathError("lcm expects 2 arguments"); return Scalar(lcm_ll(require_integer(a[0], "a", "lcm"), require_integer(a[1], "b", "lcm"))); };
    funcs["mod"] = [](const std::vector<Scalar>& a) {
        if(a.size()!=2) throw MathError("mod expects 2 arguments");
        const long long lhs = require_integer(a[0], "lhs", "mod");
        const long long rhs = require_integer(a[1], "rhs", "mod");
        if (rhs == 0) throw MathError("mod by zero");
        return Scalar(lhs % rhs);
    };
    funcs["factorial"] = [](const std::vector<Scalar>& a) { if(a.size()!=1) throw MathError("factorial expects 1 argument"); return factorial_scalar(require_integer(a[0], "argument", "factorial")); };
    funcs["nCr"] = [](const std::vector<Scalar>& a) { if(a.size()!=2) throw MathError("nCr expects 2 arguments"); return combination_scalar(require_integer(a[0], "n", "nCr"), require_integer(a[1], "r", "nCr")); };
    funcs["binom"] = funcs["nCr"];
    funcs["nPr"] = [](const std::vector<Scalar>& a) { if(a.size()!=2) throw MathError("nPr expects 2 arguments"); return permutation_scalar(require_integer(a[0], "n", "nPr"), require_integer(a[1], "r", "nPr")); };
    funcs["fib"] = [](const std::vector<Scalar>& a) { if(a.size()!=1) throw MathError("fib expects 1 argument"); return fibonacci_scalar(require_integer(a[0], "argument", "fib")); };
    funcs["is_prime"] = [](const std::vector<Scalar>& a) { if(a.size()!=1) throw MathError("is_prime expects 1 argument"); return Scalar(is_prime_ll(require_integer(a[0], "argument", "is_prime")) ? 1LL : 0LL); };
    funcs["next_prime"] = [](const std::vector<Scalar>& a) { if(a.size()!=1) throw MathError("next_prime expects 1 argument"); return Scalar(next_prime_ll(require_integer(a[0], "argument", "next_prime"))); };
    funcs["prev_prime"] = [](const std::vector<Scalar>& a) { if(a.size()!=1) throw MathError("prev_prime expects 1 argument"); return Scalar(prev_prime_ll(require_integer(a[0], "argument", "prev_prime"))); };
    funcs["euler_phi"] = [](const std::vector<Scalar>& a) { if(a.size()!=1) throw MathError("euler_phi expects 1 argument"); return Scalar(euler_phi_ll(require_integer(a[0], "argument", "euler_phi"))); };
    funcs["phi"] = funcs["euler_phi"];
    funcs["mobius"] = [](const std::vector<Scalar>& a) { if(a.size()!=1) throw MathError("mobius expects 1 argument"); return Scalar(mobius_ll(require_integer(a[0], "argument", "mobius"))); };
    funcs["prime_pi"] = [](const std::vector<Scalar>& a) { if(a.size()!=1) throw MathError("prime_pi expects 1 argument"); return Scalar(prime_pi_ll(require_integer(a[0], "argument", "prime_pi"))); };

    funcs["egcd"] = [](const std::vector<Scalar>& a) {
        if(a.size()!=2) throw MathError("egcd expects 2 arguments");
        long long x = 0, y = 0;
        return Scalar(extended_gcd_ll(require_integer(a[0], "a", "egcd"), require_integer(a[1], "b", "egcd"), &x, &y));
    };

    // Bitwise
    funcs["and"] = [](const std::vector<Scalar>& a) { if(a.size()!=2) throw MathError("and expects 2 arguments"); return Scalar(require_integer(a[0], "lhs", "and") & require_integer(a[1], "rhs", "and")); };
    funcs["or"] = [](const std::vector<Scalar>& a) { if(a.size()!=2) throw MathError("or expects 2 arguments"); return Scalar(require_integer(a[0], "lhs", "or") | require_integer(a[1], "rhs", "or")); };
    funcs["xor"] = [](const std::vector<Scalar>& a) { if(a.size()!=2) throw MathError("xor expects 2 arguments"); return Scalar(require_integer(a[0], "lhs", "xor") ^ require_integer(a[1], "rhs", "xor")); };
    funcs["not"] = [](const std::vector<Scalar>& a) { if(a.size()!=1) throw MathError("not expects 1 argument"); return Scalar(~require_integer(a[0], "argument", "not")); };
    funcs["shl"] = [](const std::vector<Scalar>& a) {
        if(a.size()!=2) throw MathError("shl expects 2 arguments");
        const long long count = require_integer(a[1], "shift", "shl");
        if (count < 0) throw MathError("shift count cannot be negative");
        return Scalar(require_integer(a[0], "value", "shl") << count);
    };
    funcs["shr"] = [](const std::vector<Scalar>& a) {
        if(a.size()!=2) throw MathError("shr expects 2 arguments");
        const long long count = require_integer(a[1], "shift", "shr");
        if (count < 0) throw MathError("shift count cannot be negative");
        return Scalar(require_integer(a[0], "value", "shr") >> count);
    };
    funcs["rol"] = [](const std::vector<Scalar>& a) {
        if(a.size()!=2) throw MathError("rol expects 2 arguments");
        return Scalar(static_cast<long long>(from_unsigned_bits(rotate_left_bits(
            to_unsigned_bits(require_integer(a[0], "value", "rol")),
            normalize_rotation_count(require_integer(a[1], "shift", "rol"))))));
    };
    funcs["ror"] = [](const std::vector<Scalar>& a) {
        if(a.size()!=2) throw MathError("ror expects 2 arguments");
        return Scalar(static_cast<long long>(from_unsigned_bits(rotate_right_bits(
            to_unsigned_bits(require_integer(a[0], "value", "ror")),
            normalize_rotation_count(require_integer(a[1], "shift", "ror"))))));
    };
    funcs["popcount"] = [](const std::vector<Scalar>& a) { if(a.size()!=1) throw MathError("popcount expects 1 argument"); return Scalar(static_cast<long long>(popcount_bits(to_unsigned_bits(require_integer(a[0], "argument", "popcount"))))); };
    funcs["bitlen"] = [](const std::vector<Scalar>& a) { if(a.size()!=1) throw MathError("bitlen expects 1 argument"); return Scalar(static_cast<long long>(bit_length_bits(to_unsigned_bits(require_integer(a[0], "argument", "bitlen"))))); };
    funcs["ctz"] = [](const std::vector<Scalar>& a) { if(a.size()!=1) throw MathError("ctz expects 1 argument"); return Scalar(static_cast<long long>(trailing_zero_count_bits(to_unsigned_bits(require_integer(a[0], "argument", "ctz"))))); };
    funcs["clz"] = [](const std::vector<Scalar>& a) { if(a.size()!=1) throw MathError("clz expects 1 argument"); return Scalar(static_cast<long long>(leading_zero_count_bits(to_unsigned_bits(require_integer(a[0], "argument", "clz"))))); };
    funcs["parity"] = [](const std::vector<Scalar>& a) { if(a.size()!=1) throw MathError("parity expects 1 argument"); return Scalar(static_cast<long long>(parity_bits(to_unsigned_bits(require_integer(a[0], "argument", "parity"))))); };
    funcs["reverse_bits"] = [](const std::vector<Scalar>& a) { if(a.size()!=1) throw MathError("reverse_bits expects 1 argument"); return Scalar(static_cast<long long>(from_unsigned_bits(reverse_bits(to_unsigned_bits(require_integer(a[0], "argument", "reverse_bits")))))); };

    return funcs;
}

/**
 * @brief 获取模块提供的所有函数名称列表
 * @return 函数名称列表
 */
std::vector<std::string> IntegerMathModule::get_functions() const {
    std::vector<std::string> names;
    auto sfuncs = get_scalar_functions();
    for (const auto& [name, _] : sfuncs) names.push_back(name);
    return names;
}

/**
 * @brief 获取帮助信息片段
 * @param topic 主题名称
 * @return 对应主题的帮助文本
 */
std::string IntegerMathModule::get_help_snippet(const std::string& topic) const {
    if (topic == "programmer") {
        return "Integer & Programmer tools:\n"
               "  factor(n)      Factorize an integer\n"
               "  bin, oct, hex, base(n, b)  Base conversion\n"
               "  Bitwise:       and, or, xor, not, shl, shr, rol, ror\n"
               "  Bit metrics:   popcount, bitlen, ctz, clz, parity, reverse_bits\n"
               "  Number Theory: gcd, lcm, mod, is_prime, next_prime, prev_prime, phi, prime_pi";
    }
    return "";
}

#include "module/module_registration.h"
REGISTER_CALCULATOR_MODULE(IntegerMathModule)
