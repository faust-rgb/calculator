/**
 * @file integer_math_module.cpp
 * @brief 整数数学模块实现
 *
 * 本文件实现了 IntegerMathModule 类，提供整数数学、数论和进制转换功能。
 * 包括因式分解、进制转换、位运算和数论函数等。
 * 所有函数通过统一的 StoredValue 接口注册，
 * 使用 wrap_scalar 辅助函数将标量计算包装为 StoredValue 签名。
 */

#include "integer_math_module.h"
#include "core/services/core_manager_interfaces.h"
#include "core/services/service_locator.h"
#include "math/functions/integer/integer_helpers.h"
#include "math/functions/special/special_functions.h"
#include "math/functions/integer/bitwise_helpers.h"
#include "math/functions/conversion/base_conversions.h"
#include "core/common/calculator_exceptions.h"
#include "math/mymath.h"
#include "matrix/matrix.h"
#include <algorithm>
#include <sstream>
#include <set>
#include <tuple>

namespace {

using Scalar = mymath::Scalar;

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

auto wrap_scalar(std::function<Scalar(const std::vector<Scalar>&)> f,
                 const std::string& name, std::size_t min_args, std::size_t max_args) {
    return [f = std::move(f), name, min_args, max_args]
           (const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() < min_args || args.size() > max_args)
            throw MathError(name + " expects " + std::to_string(min_args) +
                            (min_args == max_args ? "" : " to " + std::to_string(max_args)) +
                            " argument(s)");
        std::vector<Scalar> sa;
        sa.reserve(args.size());
        for (const auto& a : args) sa.push_back(a.decimal);
        StoredValue sv;
        sv.decimal = f(sa);
        return sv;
    };
}

} // namespace


/**
 * @brief 执行命令式操作（如 factor、bin、oct、hex、base）
 * @param command 命令名称
 * @param args 参数列表
 * @param services 核心服务接口
 * @return 执行结果的字符串表示
 * @throws std::runtime_error 如果参数无效
 */
std::string IntegerMathModule::execute_args_view(std::string_view command,
                                          const std::vector<std::string_view>& args,
                                          ServiceLocator& locator) {
    using namespace module_helpers;
    auto& services = *locator.resolve<CoreServices>();

    if (command == "factor") {
        require_args_count(args, 1, 1, command);
        Scalar val = services.evaluation.parse_decimal(std::string(args[0]));
        return factor_integer(require_integer(val, "argument", "factor"));
    }

    if (command == "bin" || command == "oct" || command == "hex" || command == "base") {
        require_args_count(args, 1, 2, command);

        Scalar value = services.evaluation.parse_decimal(std::string(args[0]));
        int base = 10;

        if (command == "bin") base = 2;
        else if (command == "oct") base = 8;
        else if (command == "hex") base = 16;
        else {
            if (args.size() < 2) throw std::runtime_error("base expects 2 arguments: value, base");
            base = static_cast<int>(require_integer(services.evaluation.parse_decimal(std::string(args[1])), "base", "base"));
        }

        if (base < 2 || base > 36) throw std::runtime_error("base must be in range [2, 36]");

        // 获取 hex 格式配置
        auto config = locator.resolve<IConfigManager>();
        bool uppercase = config->is_hex_uppercase_mode();
        bool prefix = config->is_hex_prefix_mode();

        long long n = require_integer(value, "value", std::string(command));
        return convert_to_base(n, base, uppercase, prefix);
    }

    return "";
}

std::vector<std::string> IntegerMathModule::get_help_topics() const {
    return {"programmer"};
}

/**
 * @brief 获取模块提供的函数映射
 * @return 函数名称到函数实现的映射
 *
 * 包括数论函数（gcd、lcm、factorial、nCr、nPr等）、
 * 位运算函数（and、or、xor、shl、shr、rol、ror等）
 * 以及位统计函数（popcount、bitlen、ctz、clz等）。
 */
std::map<std::string, std::function<StoredValue(const std::vector<StoredValue>&)>>
IntegerMathModule::get_functions_map() const {
    std::map<std::string, std::function<StoredValue(const std::vector<StoredValue>&)>> funcs;

    // Number Theory
    funcs["gcd"] = wrap_scalar([](const std::vector<Scalar>& a) {
        return Scalar(gcd_ll(require_integer(a[0], "a", "gcd"), require_integer(a[1], "b", "gcd")));
    }, "gcd", 2, 2);
    funcs["lcm"] = wrap_scalar([](const std::vector<Scalar>& a) {
        return Scalar(lcm_ll(require_integer(a[0], "a", "lcm"), require_integer(a[1], "b", "lcm")));
    }, "lcm", 2, 2);
    funcs["mod"] = wrap_scalar([](const std::vector<Scalar>& a) {
        const long long lhs = require_integer(a[0], "lhs", "mod");
        const long long rhs = require_integer(a[1], "rhs", "mod");
        if (rhs == 0) throw MathError("mod by zero");
        return Scalar(lhs % rhs);
    }, "mod", 2, 2);
    funcs["factorial"] = wrap_scalar([](const std::vector<Scalar>& a) {
        return mymath::factorial_scalar(require_integer(a[0], "argument", "factorial"));
    }, "factorial", 1, 1);
    funcs["nCr"] = wrap_scalar([](const std::vector<Scalar>& a) {
        return mymath::combination_scalar(require_integer(a[0], "n", "nCr"), require_integer(a[1], "r", "nCr"));
    }, "nCr", 2, 2);
    funcs["binom"] = funcs["nCr"];
    funcs["nPr"] = wrap_scalar([](const std::vector<Scalar>& a) {
        return mymath::permutation_scalar(require_integer(a[0], "n", "nPr"), require_integer(a[1], "r", "nPr"));
    }, "nPr", 2, 2);
    funcs["fib"] = wrap_scalar([](const std::vector<Scalar>& a) {
        return mymath::fibonacci_scalar(require_integer(a[0], "argument", "fib"));
    }, "fib", 1, 1);
    funcs["is_prime"] = wrap_scalar([](const std::vector<Scalar>& a) {
        return Scalar(is_prime_ll(require_integer(a[0], "argument", "is_prime")) ? 1LL : 0LL);
    }, "is_prime", 1, 1);
    funcs["next_prime"] = wrap_scalar([](const std::vector<Scalar>& a) {
        return Scalar(next_prime_ll(require_integer(a[0], "argument", "next_prime")));
    }, "next_prime", 1, 1);
    funcs["prev_prime"] = wrap_scalar([](const std::vector<Scalar>& a) {
        return Scalar(prev_prime_ll(require_integer(a[0], "argument", "prev_prime")));
    }, "prev_prime", 1, 1);
    funcs["euler_phi"] = wrap_scalar([](const std::vector<Scalar>& a) {
        return Scalar(euler_phi_ll(require_integer(a[0], "argument", "euler_phi")));
    }, "euler_phi", 1, 1);
    funcs["phi"] = funcs["euler_phi"];
    funcs["mobius"] = wrap_scalar([](const std::vector<Scalar>& a) {
        return Scalar(mobius_ll(require_integer(a[0], "argument", "mobius")));
    }, "mobius", 1, 1);
    funcs["prime_pi"] = wrap_scalar([](const std::vector<Scalar>& a) {
        return Scalar(prime_pi_ll(require_integer(a[0], "argument", "prime_pi")));
    }, "prime_pi", 1, 1);
    funcs["egcd"] = wrap_scalar([](const std::vector<Scalar>& a) {
        long long x = 0, y = 0;
        return Scalar(extended_gcd_ll(require_integer(a[0], "a", "egcd"), require_integer(a[1], "b", "egcd"), &x, &y));
    }, "egcd", 2, 2);

    // Bitwise
    funcs["and"] = wrap_scalar([](const std::vector<Scalar>& a) {
        return Scalar(require_integer(a[0], "lhs", "and") & require_integer(a[1], "rhs", "and"));
    }, "and", 2, 2);
    funcs["or"] = wrap_scalar([](const std::vector<Scalar>& a) {
        return Scalar(require_integer(a[0], "lhs", "or") | require_integer(a[1], "rhs", "or"));
    }, "or", 2, 2);
    funcs["xor"] = wrap_scalar([](const std::vector<Scalar>& a) {
        return Scalar(require_integer(a[0], "lhs", "xor") ^ require_integer(a[1], "rhs", "xor"));
    }, "xor", 2, 2);
    funcs["not"] = wrap_scalar([](const std::vector<Scalar>& a) {
        return Scalar(~require_integer(a[0], "argument", "not"));
    }, "not", 1, 1);
    funcs["shl"] = wrap_scalar([](const std::vector<Scalar>& a) {
        const long long count = require_integer(a[1], "shift", "shl");
        if (count < 0) throw MathError("shift count cannot be negative");
        return Scalar(require_integer(a[0], "value", "shl") << count);
    }, "shl", 2, 2);
    funcs["shr"] = wrap_scalar([](const std::vector<Scalar>& a) {
        const long long count = require_integer(a[1], "shift", "shr");
        if (count < 0) throw MathError("shift count cannot be negative");
        return Scalar(require_integer(a[0], "value", "shr") >> count);
    }, "shr", 2, 2);
    funcs["rol"] = wrap_scalar([](const std::vector<Scalar>& a) {
        return Scalar(static_cast<long long>(from_unsigned_bits(rotate_left_bits(
            to_unsigned_bits(require_integer(a[0], "value", "rol")),
            normalize_rotation_count(require_integer(a[1], "shift", "rol"))))));
    }, "rol", 2, 2);
    funcs["ror"] = wrap_scalar([](const std::vector<Scalar>& a) {
        return Scalar(static_cast<long long>(from_unsigned_bits(rotate_right_bits(
            to_unsigned_bits(require_integer(a[0], "value", "ror")),
            normalize_rotation_count(require_integer(a[1], "shift", "ror"))))));
    }, "ror", 2, 2);
    funcs["popcount"] = wrap_scalar([](const std::vector<Scalar>& a) {
        return Scalar(static_cast<long long>(popcount_bits(to_unsigned_bits(require_integer(a[0], "argument", "popcount")))));
    }, "popcount", 1, 1);
    funcs["bitlen"] = wrap_scalar([](const std::vector<Scalar>& a) {
        return Scalar(static_cast<long long>(bit_length_bits(to_unsigned_bits(require_integer(a[0], "argument", "bitlen")))));
    }, "bitlen", 1, 1);
    funcs["ctz"] = wrap_scalar([](const std::vector<Scalar>& a) {
        return Scalar(static_cast<long long>(trailing_zero_count_bits(to_unsigned_bits(require_integer(a[0], "argument", "ctz")))));
    }, "ctz", 1, 1);
    funcs["clz"] = wrap_scalar([](const std::vector<Scalar>& a) {
        return Scalar(static_cast<long long>(leading_zero_count_bits(to_unsigned_bits(require_integer(a[0], "argument", "clz")))));
    }, "clz", 1, 1);
    funcs["parity"] = wrap_scalar([](const std::vector<Scalar>& a) {
        return Scalar(static_cast<long long>(parity_bits(to_unsigned_bits(require_integer(a[0], "argument", "parity")))));
    }, "parity", 1, 1);
    funcs["reverse_bits"] = wrap_scalar([](const std::vector<Scalar>& a) {
        return Scalar(static_cast<long long>(from_unsigned_bits(reverse_bits(to_unsigned_bits(require_integer(a[0], "argument", "reverse_bits"))))));
    }, "reverse_bits", 1, 1);

    funcs["divisors"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() != 1 || args[0].is_matrix || args[0].is_complex) {
            throw std::runtime_error("divisors expects one integer");
        }
        const Scalar raw = args[0].get_decimal();
        if (!mymath::is_integer(raw) || raw <= Scalar(0)) throw std::runtime_error("divisors expects a positive integer");
        const long long value = static_cast<long long>(raw);
        std::vector<Scalar> values;
        for (long long d = 1; d <= value / d; ++d) {
            if (value % d != 0) continue;
            values.push_back(Scalar(d));
            if (d != value / d) values.push_back(Scalar(value / d));
        }
        std::sort(values.begin(), values.end());
        StoredValue res;
        res.is_matrix = true;
        res.matrix_ptr = std::make_shared<matrix::Matrix>(matrix::Matrix::vector(std::move(values)));
        return res;
    };
    funcs["extended_gcd"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() != 2 || args[0].is_matrix || args[1].is_matrix ||
            args[0].is_complex || args[1].is_complex) {
            throw std::runtime_error("extended_gcd expects two integers");
        }
        if (!mymath::is_integer(args[0].get_decimal()) || !mymath::is_integer(args[1].get_decimal())) {
            throw std::runtime_error("extended_gcd expects two integers");
        }
        long long a = static_cast<long long>(args[0].get_decimal());
        long long b = static_cast<long long>(args[1].get_decimal());
        long long old_r = a, r = b, old_s = 1, s = 0, old_t = 0, t = 1;
        while (r != 0) {
            const long long q = old_r / r;
            std::tie(old_r, r) = std::make_pair(r, old_r - q * r);
            std::tie(old_s, s) = std::make_pair(s, old_s - q * s);
            std::tie(old_t, t) = std::make_pair(t, old_t - q * t);
        }
        StoredValue res;
        res.is_matrix = true;
        res.matrix_ptr = std::make_shared<matrix::Matrix>(matrix::Matrix::vector({Scalar(old_r), Scalar(old_s), Scalar(old_t)}));
        return res;
    };
    funcs["xgcd"] = funcs["extended_gcd"];

    return funcs;
}

/**
 * @brief 获取模块提供的所有函数名称列表
 * @return 函数名称列表
 */
std::vector<std::string> IntegerMathModule::get_function_names() const {
    static const std::vector<std::string> names = {
        "and", "binom", "bitlen", "clz", "ctz", "divisors", "egcd", "euler_phi",
        "extended_gcd", "factorial", "fib", "gcd", "is_prime", "lcm", "mobius",
        "mod", "nCr", "nPr", "next_prime", "not", "or", "parity", "phi",
        "popcount", "prev_prime", "prime_pi", "reverse_bits", "rol", "ror",
        "shl", "shr", "xgcd", "xor"
    };
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
