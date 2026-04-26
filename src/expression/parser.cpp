#include "parser.h"

#include <cctype>
#include <stdexcept>
#include <utility>

namespace expression {
namespace {

numeric::BigDecimal decimal_from_scientific(const std::string& token) {
    std::string digits;
    bool negative = false;
    bool seen_dot = false;
    int scale = 0;
    int exponent = 0;
    std::size_t pos = 0;
    bool has_exponent = false;
    if (pos < token.size() && (token[pos] == '+' || token[pos] == '-')) {
        negative = token[pos] == '-';
        ++pos;
    }
    bool has_digit = false;
    for (; pos < token.size(); ++pos) {
        const char ch = token[pos];
        if (ch == 'e' || ch == 'E') {
            has_exponent = true;
            ++pos;
            break;
        }
        if (ch == '.') {
            if (seen_dot) {
                throw std::runtime_error("invalid number literal: " + token);
            }
            seen_dot = true;
            continue;
        }
        if (!std::isdigit(static_cast<unsigned char>(ch))) {
            throw std::runtime_error("invalid number literal: " + token);
        }
        digits.push_back(ch);
        if (seen_dot) {
            ++scale;
        }
        has_digit = true;
    }
    if (!has_digit) {
        throw std::runtime_error("invalid number literal: " + token);
    }
    if (has_exponent) {
        bool negative_exponent = false;
        if (pos < token.size() && (token[pos] == '+' || token[pos] == '-')) {
            negative_exponent = token[pos] == '-';
            ++pos;
        }
        if (pos == token.size()) {
            throw std::runtime_error("invalid exponent in number literal: " + token);
        }
        for (; pos < token.size(); ++pos) {
            if (!std::isdigit(static_cast<unsigned char>(token[pos]))) {
                throw std::runtime_error("invalid exponent in number literal: " + token);
            }
            exponent = exponent * 10 + (token[pos] - '0');
        }
        if (negative_exponent) {
            exponent = -exponent;
        }
    }
    scale -= exponent;
    numeric::BigInt coefficient = numeric::BigInt::from_string(digits);
    if (negative) {
        coefficient = -coefficient;
    }
    return numeric::BigDecimal(coefficient, scale);
}

class Parser {
public:
    explicit Parser(std::string source) : source_(std::move(source)) {}

    Expr parse() {
        Expr expr = parse_expression();
        skip_spaces();
        if (!at_end()) {
            throw std::runtime_error("unexpected token near: " + source_.substr(pos_, 1));
        }
        return expr;
    }

private:
    Expr parse_expression() {
        Expr value = parse_term();
        while (true) {
            skip_spaces();
            if (match('+')) {
                value = Expr::add(value, parse_term());
            } else if (match('-')) {
                value = Expr::add(value, Expr::mul(Expr::number(numeric::Number(-1)),
                                                  parse_term()));
            } else {
                return value;
            }
        }
    }

    Expr parse_term() {
        Expr value = parse_unary();
        while (true) {
            skip_spaces();
            if (match('*')) {
                value = Expr::mul(value, parse_unary());
            } else if (match('/')) {
                value = Expr::mul(value, Expr::pow(parse_unary(),
                                                  Expr::number(numeric::Number(-1))));
            } else {
                return value;
            }
        }
    }

    Expr parse_unary() {
        skip_spaces();
        if (match('+')) {
            return parse_unary();
        }
        if (match('-')) {
            return Expr::mul(Expr::number(numeric::Number(-1)), parse_unary());
        }
        return parse_power();
    }

    Expr parse_power() {
        Expr value = parse_primary();
        skip_spaces();
        if (match('^')) {
            value = Expr::pow(value, parse_unary());
        }
        return value;
    }

    Expr parse_primary() {
        skip_spaces();
        if (match('(')) {
            Expr value = parse_expression();
            skip_spaces();
            expect(')');
            return value;
        }
        if (match('[')) {
            std::vector<Expr> items;
            skip_spaces();
            if (!peek(']')) {
                while (true) {
                    items.push_back(parse_expression());
                    skip_spaces();
                    if (!match(',')) {
                        break;
                    }
                }
            }
            expect(']');
            return Expr::list(items);
        }
        if (!at_end() && (std::isdigit(static_cast<unsigned char>(source_[pos_])) ||
                          source_[pos_] == '.')) {
            return parse_number();
        }
        if (peek_alpha()) {
            const std::string name = parse_identifier();
            skip_spaces();
            if (!match('(')) {
                return Expr::symbol(name);
            }
            std::vector<Expr> args;
            skip_spaces();
            if (!peek(')')) {
                while (true) {
                    args.push_back(parse_expression());
                    skip_spaces();
                    if (!match(',')) {
                        break;
                    }
                }
            }
            expect(')');
            return Expr::function(name, args);
        }
        throw std::runtime_error("expected expression");
    }

    Expr parse_number() {
        const std::size_t start = pos_;
        bool seen_dot = false;
        while (!at_end()) {
            const char ch = source_[pos_];
            if (std::isdigit(static_cast<unsigned char>(ch))) {
                ++pos_;
            } else if (ch == '.' && !seen_dot) {
                seen_dot = true;
                ++pos_;
            } else if (ch == 'e' || ch == 'E') {
                seen_dot = true;
                ++pos_;
                if (!at_end() && (source_[pos_] == '+' || source_[pos_] == '-')) {
                    ++pos_;
                }
                const std::size_t exponent_start = pos_;
                while (!at_end() && std::isdigit(static_cast<unsigned char>(source_[pos_]))) {
                    ++pos_;
                }
                if (exponent_start == pos_) {
                    throw std::runtime_error("invalid exponent in number literal");
                }
                break;
            } else {
                break;
            }
        }
        std::string token = source_.substr(start, pos_ - start);
        if (token == ".") {
            throw std::runtime_error("invalid number literal");
        }
        skip_spaces();
        if (match('i')) {
            const numeric::BigDecimal imag = seen_dot
                ? decimal_from_scientific(token)
                : numeric::BigDecimal(numeric::BigInt::from_string(token));
            return Expr::number(numeric::Number(numeric::Complex(numeric::BigDecimal(), imag)));
        }
        if (seen_dot) {
            return Expr::number(numeric::Number(decimal_from_scientific(token)));
        }
        return Expr::number(numeric::Number(numeric::BigInt::from_string(token)));
    }

    std::string parse_identifier() {
        const std::size_t start = pos_;
        while (!at_end() &&
               (std::isalnum(static_cast<unsigned char>(source_[pos_])) ||
                source_[pos_] == '_')) {
            ++pos_;
        }
        return source_.substr(start, pos_ - start);
    }

    void skip_spaces() {
        while (!at_end() && std::isspace(static_cast<unsigned char>(source_[pos_]))) {
            ++pos_;
        }
    }

    bool match(char ch) {
        skip_spaces();
        if (peek(ch)) {
            ++pos_;
            return true;
        }
        return false;
    }

    void expect(char ch) {
        if (!match(ch)) {
            throw std::runtime_error(std::string("expected '") + ch + "'");
        }
    }

    bool peek(char ch) const {
        return !at_end() && source_[pos_] == ch;
    }

    bool peek_alpha() const {
        return !at_end() &&
               (std::isalpha(static_cast<unsigned char>(source_[pos_])) ||
                source_[pos_] == '_');
    }

    bool at_end() const {
        return pos_ >= source_.size();
    }

    std::string source_;
    std::size_t pos_ = 0;
};

}  // namespace

Expr parse_expression(const std::string& text) {
    return Parser(text).parse();
}

}  // namespace expression
