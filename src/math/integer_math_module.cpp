#include "integer_math_module.h"
#include "helpers/integer_helpers.h"
#include "helpers/combinatorics.h"
#include "helpers/bitwise_helpers.h"
#include "helpers/base_conversions.h"
#include "../core/utils.h"
#include "../core/calculator_exceptions.h"
#include "mymath.h"
#include <cmath>
#include <algorithm>
#include <sstream>
#include <set>

namespace {

using namespace utils;

long long require_integer(double x, const std::string& name, const std::string& func) {
    if (!is_integer_double(x)) {
        throw MathError(func + " requires an integer " + name);
    }
    return round_to_long_long(x);
}

} // namespace

bool IntegerMathModule::can_handle(const std::string& command) const {
    static const std::set<std::string> cmds = {"factor", "bin", "oct", "hex", "base"};
    return cmds.count(command) > 0;
}

std::string IntegerMathModule::execute_args(const std::string& command,
                                          const std::vector<std::string>& args,
                                          const CoreServices& services) {
    if (command == "factor") {
        if (args.empty()) throw std::runtime_error("factor expects 1 argument");
        double val = services.evaluation.parse_decimal(args[0]);
        return factor_integer(require_integer(val, "argument", "factor"));
    }

    if (command == "bin" || command == "oct" || command == "hex" || command == "base") {
        if (args.empty()) throw std::runtime_error(command + " expects at least 1 argument");
        
        double value = services.evaluation.parse_decimal(args[0]);
        int base = 10;
        
        if (command == "bin") base = 2;
        else if (command == "oct") base = 8;
        else if (command == "hex") base = 16;
        else {
            if (args.size() < 2) throw std::runtime_error("base expects 2 arguments: value, base");
            base = static_cast<int>(require_integer(services.evaluation.parse_decimal(args[1]), "base", "base"));
        }

        std::string output;
        // 使用 utils 中的转换逻辑，或者在此重写
        // 这里为了简单，我们调用核心提供的转换服务（如果以后增加了的话）
        // 目前我们直接复用 core/utils.cpp 中的 convert_base_value (如果它是公开的)
        // 由于它是匿名的或非公开的，我们在此实现一个简版或者调用 evaluate_value
        
        // 我们先通过 execute_script 回调回核心，利用现有的 try_base_conversion_expression 逻辑
        // 或者直接在这里重新实现 base conversion 逻辑以彻底解耦。
        // 考虑到解耦目标，我们在这里实现。
        
        if (base < 2 || base > 36) throw std::runtime_error("base must be in range [2, 36]");
        
        long long n = require_integer(value, "value", command);
        bool negative = n < 0;
        if (negative) n = -n;
        
        std::string res;
        if (n == 0) res = "0";
        else {
            while (n > 0) {
                int digit = n % base;
                res += (digit < 10 ? (char)('0' + digit) : (char)('A' + digit - 10));
                n /= base;
            }
            std::reverse(res.begin(), res.end());
        }
        return (negative ? "-" : "") + res;
    }

    return "";
}

std::map<std::string, std::function<double(const std::vector<double>&)>>
IntegerMathModule::get_scalar_functions() const {
    std::map<std::string, std::function<double(const std::vector<double>&)>> funcs;

    // Number Theory
    funcs["gcd"] = [](const std::vector<double>& a) { if(a.size()!=2) throw MathError("gcd expects 2 arguments"); return static_cast<double>(gcd_ll(round_to_long_long(a[0]), round_to_long_long(a[1]))); };
    funcs["lcm"] = [](const std::vector<double>& a) { if(a.size()!=2) throw MathError("lcm expects 2 arguments"); return static_cast<double>(lcm_ll(round_to_long_long(a[0]), round_to_long_long(a[1]))); };
    funcs["mod"] = [](const std::vector<double>& a) { 
        if(a.size()!=2) throw MathError("mod expects 2 arguments");
        const long long lhs = require_integer(a[0], "lhs", "mod");
        const long long rhs = require_integer(a[1], "rhs", "mod");
        if (rhs == 0) throw MathError("mod by zero");
        return static_cast<double>(lhs % rhs);
    };
    funcs["factorial"] = [](const std::vector<double>& a) { if(a.size()!=1) throw MathError("factorial expects 1 argument"); return factorial_value(require_integer(a[0], "argument", "factorial")); };
    funcs["nCr"] = [](const std::vector<double>& a) { if(a.size()!=2) throw MathError("nCr expects 2 arguments"); return combination_value(require_integer(a[0], "n", "nCr"), require_integer(a[1], "r", "nCr")); };
    funcs["binom"] = funcs["nCr"];
    funcs["nPr"] = [](const std::vector<double>& a) { if(a.size()!=2) throw MathError("nPr expects 2 arguments"); return permutation_value(require_integer(a[0], "n", "nPr"), require_integer(a[1], "r", "nPr")); };
    funcs["fib"] = [](const std::vector<double>& a) { if(a.size()!=1) throw MathError("fib expects 1 argument"); return fibonacci_value(round_to_long_long(a[0])); };
    funcs["is_prime"] = [](const std::vector<double>& a) { if(a.size()!=1) throw MathError("is_prime expects 1 argument"); return is_prime_ll(require_integer(a[0], "argument", "is_prime")) ? 1.0 : 0.0; };
    funcs["next_prime"] = [](const std::vector<double>& a) { if(a.size()!=1) throw MathError("next_prime expects 1 argument"); return static_cast<double>(next_prime_ll(require_integer(a[0], "argument", "next_prime"))); };
    funcs["prev_prime"] = [](const std::vector<double>& a) { if(a.size()!=1) throw MathError("prev_prime expects 1 argument"); return static_cast<double>(prev_prime_ll(require_integer(a[0], "argument", "prev_prime"))); };
    funcs["euler_phi"] = [](const std::vector<double>& a) { if(a.size()!=1) throw MathError("euler_phi expects 1 argument"); return static_cast<double>(euler_phi_ll(require_integer(a[0], "argument", "euler_phi"))); };
    funcs["phi"] = funcs["euler_phi"];
    funcs["mobius"] = [](const std::vector<double>& a) { if(a.size()!=1) throw MathError("mobius expects 1 argument"); return static_cast<double>(mobius_ll(require_integer(a[0], "argument", "mobius"))); };
    funcs["prime_pi"] = [](const std::vector<double>& a) { if(a.size()!=1) throw MathError("prime_pi expects 1 argument"); return static_cast<double>(prime_pi_ll(require_integer(a[0], "argument", "prime_pi"))); };
    
    funcs["egcd"] = [](const std::vector<double>& a) { 
        if(a.size()!=2) throw MathError("egcd expects 2 arguments"); 
        long long x = 0, y = 0;
        return static_cast<double>(extended_gcd_ll(require_integer(a[0], "a", "egcd"), require_integer(a[1], "b", "egcd"), &x, &y));
    };

    // Bitwise
    funcs["and"] = [](const std::vector<double>& a) { if(a.size()!=2) throw MathError("and expects 2 arguments"); return static_cast<double>(require_integer(a[0], "lhs", "and") & require_integer(a[1], "rhs", "and")); };
    funcs["or"] = [](const std::vector<double>& a) { if(a.size()!=2) throw MathError("or expects 2 arguments"); return static_cast<double>(require_integer(a[0], "lhs", "or") | require_integer(a[1], "rhs", "or")); };
    funcs["xor"] = [](const std::vector<double>& a) { if(a.size()!=2) throw MathError("xor expects 2 arguments"); return static_cast<double>(require_integer(a[0], "lhs", "xor") ^ require_integer(a[1], "rhs", "xor")); };
    funcs["not"] = [](const std::vector<double>& a) { if(a.size()!=1) throw MathError("not expects 1 argument"); return static_cast<double>(~require_integer(a[0], "argument", "not")); };
    funcs["shl"] = [](const std::vector<double>& a) { 
        if(a.size()!=2) throw MathError("shl expects 2 arguments"); 
        const long long count = require_integer(a[1], "shift", "shl");
        if (count < 0) throw MathError("shift count cannot be negative");
        return static_cast<double>(require_integer(a[0], "value", "shl") << count); 
    };
    funcs["shr"] = [](const std::vector<double>& a) { 
        if(a.size()!=2) throw MathError("shr expects 2 arguments"); 
        const long long count = require_integer(a[1], "shift", "shr");
        if (count < 0) throw MathError("shift count cannot be negative");
        return static_cast<double>(require_integer(a[0], "value", "shr") >> count); 
    };
    funcs["rol"] = [](const std::vector<double>& a) { 
        if(a.size()!=2) throw MathError("rol expects 2 arguments"); 
        return static_cast<double>(from_unsigned_bits(rotate_left_bits(
            to_unsigned_bits(require_integer(a[0], "value", "rol")),
            normalize_rotation_count(require_integer(a[1], "shift", "rol")))));
    };
    funcs["ror"] = [](const std::vector<double>& a) { 
        if(a.size()!=2) throw MathError("ror expects 2 arguments"); 
        return static_cast<double>(from_unsigned_bits(rotate_right_bits(
            to_unsigned_bits(require_integer(a[0], "value", "ror")),
            normalize_rotation_count(require_integer(a[1], "shift", "ror")))));
    };
    funcs["popcount"] = [](const std::vector<double>& a) { if(a.size()!=1) throw MathError("popcount expects 1 argument"); return static_cast<double>(popcount_bits(to_unsigned_bits(require_integer(a[0], "argument", "popcount")))); };
    funcs["bitlen"] = [](const std::vector<double>& a) { if(a.size()!=1) throw MathError("bitlen expects 1 argument"); return static_cast<double>(bit_length_bits(to_unsigned_bits(require_integer(a[0], "argument", "bitlen")))); };
    funcs["ctz"] = [](const std::vector<double>& a) { if(a.size()!=1) throw MathError("ctz expects 1 argument"); return static_cast<double>(trailing_zero_count_bits(to_unsigned_bits(require_integer(a[0], "argument", "ctz")))); };
    funcs["clz"] = [](const std::vector<double>& a) { if(a.size()!=1) throw MathError("clz expects 1 argument"); return static_cast<double>(leading_zero_count_bits(to_unsigned_bits(require_integer(a[0], "argument", "clz")))); };
    funcs["parity"] = [](const std::vector<double>& a) { if(a.size()!=1) throw MathError("parity expects 1 argument"); return static_cast<double>(parity_bits(to_unsigned_bits(require_integer(a[0], "argument", "parity")))); };
    funcs["reverse_bits"] = [](const std::vector<double>& a) { if(a.size()!=1) throw MathError("reverse_bits expects 1 argument"); return static_cast<double>(from_unsigned_bits(reverse_bits(to_unsigned_bits(require_integer(a[0], "argument", "reverse_bits"))))); };

    return funcs;
}

std::vector<std::string> IntegerMathModule::get_functions() const {
    std::vector<std::string> names;
    auto sfuncs = get_scalar_functions();
    for (const auto& [name, _] : sfuncs) names.push_back(name);
    return names;
}

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
