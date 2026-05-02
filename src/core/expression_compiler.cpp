#include "expression_compiler.h"
#include "base_parser.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <stdexcept>

// ============================================================================
// 表达式词法分析器
// ============================================================================

ExpressionLexer::ExpressionLexer(std::string_view source)
    : source_(source), pos_(0) {}

void ExpressionLexer::skip_spaces() {
    while (pos_ < source_.size() &&
           std::isspace(static_cast<unsigned char>(source_[pos_]))) {
        ++pos_;
    }
}

bool ExpressionLexer::is_at_end() const {
    return pos_ >= source_.size();
}

char ExpressionLexer::peek() const {
    if (is_at_end()) return '\0';
    return source_[pos_];
}

void ExpressionLexer::expect(char expected) {
    if (peek() == expected) {
        ++pos_;
    } else {
        // Simple error handling for lexer
        throw std::runtime_error("Lexer expected character not found.");
    }
}

std::string_view ExpressionLexer::parse_identifier() {
    std::size_t start = pos_;
    while (!is_at_end() && (std::isalnum(static_cast<unsigned char>(source_[pos_])) ||
                            source_[pos_] == '_' || source_[pos_] == '\'')) {
        ++pos_;
    }
    return source_.substr(start, pos_ - start);
}

bool ExpressionLexer::tokenize(std::vector<ExpressionToken>* tokens) {
    tokens->clear();

    while (!is_at_end()) {
        skip_spaces();
        if (is_at_end()) break;

        const char ch = peek();
        const std::size_t pos = pos_;

        // 字符串字面量
        if (ch == '"') {
            if (!parse_string_token(tokens)) return false;
            continue;
        }

        // 数字字面量
        if (std::isdigit(static_cast<unsigned char>(ch)) ||
            (ch == '.' && pos_ + 1 < source_.size() &&
             std::isdigit(static_cast<unsigned char>(source_[pos_ + 1])))) {
            if (!parse_number_token(tokens)) return false;
            continue;
        }

        // 标识符
        if (std::isalpha(static_cast<unsigned char>(ch)) || ch == '_') {
            tokens->push_back({ExpressionTokenKind::kIdentifier, std::string(parse_identifier()), 0.0, pos});
            continue;
        }

        // 运算符和符号
        if (!parse_operator_token(tokens)) return false;
    }

    tokens->push_back({ExpressionTokenKind::kEOF, "", 0.0, pos_});
    return true;
}

bool ExpressionLexer::parse_string_token(std::vector<ExpressionToken>* tokens) {
    const std::size_t start = pos_;
    expect('"');
    bool escaping = false;

    while (!is_at_end()) {
        const char ch = source_[pos_++];
        if (escaping) {
            escaping = false;
        } else if (ch == '\\') {
            escaping = true;
        } else if (ch == '"') {
            std::string text = std::string(source_.substr(start, pos_ - start));
            tokens->push_back({ExpressionTokenKind::kString, text, 0.0, start});
            return true;
        }
    }

    return false; // 未终止的字符串
}

bool ExpressionLexer::parse_number_token(std::vector<ExpressionToken>* tokens) {
    const std::size_t start = pos_;
    
    // 处理进制前缀 (0x, 0b, 0o)
    if (peek() == '0' && pos_ + 1 < source_.size()) {
        char base_char = std::tolower(static_cast<unsigned char>(source_[pos_ + 1]));
        if (base_char == 'x' || base_char == 'b' || base_char == 'o') {
            pos_ += 2;
            while (!is_at_end() && std::isxdigit(static_cast<unsigned char>(peek()))) {
                pos_++;
            }
            goto end_num;
        }
    }

    { // 浮点数/普通数字解析
        bool seen_dot = false;

        // 处理前导小数点
        if (peek() == '.') {
            seen_dot = true;
            ++pos_;
        }

        // 整数部分
        while (!is_at_end() && std::isdigit(static_cast<unsigned char>(peek()))) {
            ++pos_;
        }

        // 小数部分
        if (!seen_dot && !is_at_end() && peek() == '.') {
            seen_dot = true;
            ++pos_;
            while (!is_at_end() && std::isdigit(static_cast<unsigned char>(peek()))) {
                ++pos_;
            }
        }

        // 指数部分
        if (!is_at_end() && (peek() == 'e' || peek() == 'E')) {
            ++pos_;
            if (!is_at_end() && (peek() == '+' || peek() == '-')) {
                ++pos_;
            }
            while (!is_at_end() && std::isdigit(static_cast<unsigned char>(peek()))) {
                ++pos_;
            }
        }
    }

end_num:
    // 虚数单位 i
    if (!is_at_end() && peek() == 'i') {
        // 检查 i 后面是否是标识符字符
        if (pos_ + 1 >= source_.size() ||
            (!std::isalnum(static_cast<unsigned char>(source_[pos_ + 1])) &&
             source_[pos_ + 1] != '_')) {
            // 这是一个虚数，如 3i
            ++pos_;
        }
    }

    std::string text = std::string(source_.substr(start, pos_ - start));
    double value = 0.0;
    
    // 解析数值
    try {
        // 处理进制前缀 (0x, 0b, 0o)
        if (text.size() > 2 && text[0] == '0' && std::isalpha(text[1])) {
            value = static_cast<double>(parse_prefixed_integer_token(text));
        } else {
            // 处理虚数
            std::string num_part = text;
            if (!text.empty() && text.back() == 'i') {
                num_part = text.substr(0, text.size() - 1);
            }
            if (!num_part.empty()) {
                value = std::stod(num_part);
            }
        }
    } catch (...) {
        return false;
    }

    tokens->push_back({ExpressionTokenKind::kNumber, text, value, start});
    return true;
}

bool ExpressionLexer::parse_operator_token(std::vector<ExpressionToken>* tokens) {
    const char ch = peek();
    const std::size_t pos = pos_;

    // 多字符运算符
    if (ch == '=' && pos_ + 1 < source_.size() && source_[pos_ + 1] == '=') {
        pos_ += 2;
        tokens->push_back({ExpressionTokenKind::kOperator, "==", 0.0, pos});
        return true;
    }
    if (ch == '!' && pos_ + 1 < source_.size() && source_[pos_ + 1] == '=') {
        pos_ += 2;
        tokens->push_back({ExpressionTokenKind::kOperator, "!=", 0.0, pos});
        return true;
    }
    if (ch == '<' && pos_ + 1 < source_.size() && source_[pos_ + 1] == '=') {
        pos_ += 2;
        tokens->push_back({ExpressionTokenKind::kOperator, "<=", 0.0, pos});
        return true;
    }
    if (ch == '>' && pos_ + 1 < source_.size() && source_[pos_ + 1] == '=') {
        pos_ += 2;
        tokens->push_back({ExpressionTokenKind::kOperator, ">=", 0.0, pos});
        return true;
    }
    if (ch == '&' && pos_ + 1 < source_.size() && source_[pos_ + 1] == '&') {
        pos_ += 2;
        tokens->push_back({ExpressionTokenKind::kOperator, "&&", 0.0, pos});
        return true;
    }
    if (ch == '|' && pos_ + 1 < source_.size() && source_[pos_ + 1] == '|') {
        pos_ += 2;
        tokens->push_back({ExpressionTokenKind::kOperator, "||", 0.0, pos});
        return true;
    }

    // 单字符符号
    ++pos_;
    switch (ch) {
        case '+':
        case '-':
        case '*':
        case '/':
        case '^':
        case '%':
            tokens->push_back({ExpressionTokenKind::kOperator, std::string(1, ch), 0.0, pos});
            return true;
        case '(':
            tokens->push_back({ExpressionTokenKind::kLParen, "(", 0.0, pos});
            return true;
        case ')':
            tokens->push_back({ExpressionTokenKind::kRParen, ")", 0.0, pos});
            return true;
        case '[':
            tokens->push_back({ExpressionTokenKind::kLBracket, "[", 0.0, pos});
            return true;
        case ']':
            tokens->push_back({ExpressionTokenKind::kRBracket, "]", 0.0, pos});
            return true;
        case ',':
            tokens->push_back({ExpressionTokenKind::kComma, ",", 0.0, pos});
            return true;
        case ';':
            tokens->push_back({ExpressionTokenKind::kSemicolon, ";", 0.0, pos});
            return true;
        case ':':
            tokens->push_back({ExpressionTokenKind::kColon, ":", 0.0, pos});
            return true;
        case '?':
            tokens->push_back({ExpressionTokenKind::kQuestion, "?", 0.0, pos});
            return true;
        case '=':
            tokens->push_back({ExpressionTokenKind::kEqual, "=", 0.0, pos});
            return true;
        case '<':
            tokens->push_back({ExpressionTokenKind::kOperator, "<", 0.0, pos});
            return true;
        case '>':
            tokens->push_back({ExpressionTokenKind::kOperator, ">", 0.0, pos});
            return true;
        case '!':
            tokens->push_back({ExpressionTokenKind::kOperator, "!", 0.0, pos});
            return true;
        default:
            return false;
    }
}

// ============================================================================
// 表达式分析函数
// ============================================================================

ExpressionFeature analyze_expression_features(const std::string& expression) {
    ExpressionFeature features = ExpressionFeature::kNone;
    bool in_string = false;
    bool escaping = false;
    int paren_depth = 0;
    int bracket_depth = 0;

    for (std::size_t i = 0; i < expression.size(); ++i) {
        const char ch = expression[i];

        // 字符串内部
        if (in_string) {
            if (escaping) {
                escaping = false;
            } else if (ch == '\\') {
                escaping = true;
            } else if (ch == '"') {
                in_string = false;
            }
            continue;
        }

        if (ch == '"') {
            in_string = true;
            features = features | ExpressionFeature::kHasString;
            continue;
        }

        // 括号
        if (ch == '(') {
            ++paren_depth;
        } else if (ch == ')') {
            --paren_depth;
        } else if (ch == '[') {
            ++bracket_depth;
            features = features | ExpressionFeature::kHasBracket;
        } else if (ch == ']') {
            --bracket_depth;
        }

        // 赋值（顶层等号）
        if (ch == '=' && paren_depth == 0 && bracket_depth == 0) {
            // 排除比较运算符
            if (i > 0 && expression[i-1] != '!' && expression[i-1] != '<' &&
                expression[i-1] != '>' && expression[i-1] != '=') {
                if (i + 1 >= expression.size() || expression[i+1] != '=') {
                    features = features | ExpressionFeature::kHasAssignment;
                }
            }
        }

        // 数字
        if (std::isdigit(static_cast<unsigned char>(ch)) ||
            (ch == '.' && i + 1 < expression.size() &&
             std::isdigit(static_cast<unsigned char>(expression[i+1])))) {
            features = features | ExpressionFeature::kHasNumber;
        }

        // 标识符和独立的 i
        if (std::isalpha(static_cast<unsigned char>(ch)) || ch == '_') {
            std::size_t start = i;
            while (i < expression.size() &&
                   (std::isalnum(static_cast<unsigned char>(expression[i])) ||
                    expression[i] == '_')) {
                ++i;
            }
            std::string id = expression.substr(start, i - start);
            --i;

            features = features | ExpressionFeature::kHasIdentifier;

            // 检查是否是独立的 i（虚数单位）
            if (id == "i") {
                bool prev_ok = (start == 0 ||
                    (!std::isalnum(static_cast<unsigned char>(expression[start-1])) &&
                     expression[start-1] != '_'));
                bool next_ok = (i + 1 >= expression.size() ||
                    (!std::isalnum(static_cast<unsigned char>(expression[i+1])) &&
                     expression[i+1] != '_'));
                if (prev_ok && next_ok) {
                    features = features | ExpressionFeature::kHasI;
                }
            }

            // 检查是否是 rat 调用
            if (id == "rat") {
                features = features | ExpressionFeature::kHasRatCall;
            }

            // 检查是否是由矩阵/复数表达式引擎处理的函数调用。
            if (id == "mat" || id == "vec" || id == "zeros" || id == "identity" ||
                id == "diag" || id == "eye" || id == "linspace" || id == "logspace" ||
                id == "get" || id == "set" || id == "resize" || id == "append_row" ||
                id == "append_col" || id == "transpose" || id == "inverse" ||
                id == "pinv" || id == "dot" || id == "outer" || id == "kron" ||
                id == "hadamard" || id == "null" || id == "least_squares" ||
                id == "qr_q" || id == "qr_r" || id == "lu_l" || id == "lu_u" ||
                id == "lu_p" || id == "svd_u" || id == "svd_s" || id == "svd_vt" ||
                id == "solve" || id == "norm" || id == "cond" || id == "trace" ||
                id == "det" || id == "rank" || id == "rref" || id == "eigvals" ||
                id == "eigvecs" || id == "reshape" || id == "cholesky" ||
                id == "hessenberg" || id == "schur" || id == "poly_eval" ||
                id == "poly_deriv" || id == "poly_integ" || id == "poly_compose" ||
                id == "poly_gcd" || id == "poly_fit" || id == "polynomial_fit" ||
                id == "lagrange" || id == "linear_regression" || id == "dft" ||
                id == "fft" || id == "idft" || id == "ifft" || id == "convolve" ||
                id == "hann" || id == "hanning" || id == "hamming" ||
                id == "blackman" || id == "divisors" || id == "extended_gcd" ||
                id == "xgcd" || id == "complex" || id == "polar" || id == "real" ||
                id == "imag" || id == "arg" || id == "conj") {
                features = features | ExpressionFeature::kHasMatrixFunc;
            }
        }

        // 运算符
        if (ch == '+' || ch == '-' || ch == '*' || ch == '/' || ch == '^' || ch == '%') {
            features = features | ExpressionFeature::kHasOperator;
        }

        // 比较运算符
        if (ch == '<' || ch == '>' || ch == '!' ||
            (ch == '=' && i > 0 && expression[i-1] == '=')) {
            features = features | ExpressionFeature::kHasComparison;
        }
    }

    return features;
}

ExpressionHint analyze_expression_hint(const std::string& expression) {
    const std::string trimmed = trim_copy(expression);

    // 空表达式
    if (trimmed.empty()) {
        return ExpressionHint::kScalar;
    }

    // 字符串字面量
    if (trimmed.size() >= 2 && trimmed.front() == '"' && trimmed.back() == '"') {
        // 检查是否是合法的字符串字面量
        bool valid = true;
        bool escaping = false;
        for (std::size_t i = 1; i < trimmed.size() - 1; ++i) {
            if (escaping) {
                escaping = false;
            } else if (trimmed[i] == '\\') {
                escaping = true;
            } else if (trimmed[i] == '"') {
                valid = false;
                break;
            }
        }
        if (valid) {
            return ExpressionHint::kStringLiteral;
        }
    }

    // 单个标识符
    if (is_identifier_text(trimmed)) {
        return ExpressionHint::kIdentifier;
    }

    // 纯数字
    bool is_number = true;
    bool has_digit = false;
    for (char ch : trimmed) {
        if (std::isdigit(static_cast<unsigned char>(ch))) {
            has_digit = true;
        } else if (ch != '.' && ch != 'e' && ch != 'E' && ch != '+' && ch != '-') {
            is_number = false;
            break;
        }
    }
    if (is_number && has_digit) {
        return ExpressionHint::kScalar;
    }

    // 分析特征
    ExpressionFeature features = analyze_expression_features(trimmed);

    // rat 调用
    if (has_feature(features, ExpressionFeature::kHasRatCall)) {
        return ExpressionHint::kRatCall;
    }

    // 赋值
    if (has_feature(features, ExpressionFeature::kHasAssignment)) {
        return ExpressionHint::kAssignment;
    }

    // 复数候选
    if (has_feature(features, ExpressionFeature::kHasI)) {
        return ExpressionHint::kComplexCandidate;
    }

    // 矩阵候选（包含方括号或矩阵函数）
    if (has_feature(features, ExpressionFeature::kHasBracket) ||
        has_feature(features, ExpressionFeature::kHasMatrixFunc)) {
        return ExpressionHint::kMatrixCandidate;
    }

    // 默认为标量
    if (has_feature(features, ExpressionFeature::kHasOperator) ||
        has_feature(features, ExpressionFeature::kHasNumber) ||
        has_feature(features, ExpressionFeature::kHasIdentifier)) {
        return ExpressionHint::kScalar;
    }

    return ExpressionHint::kUnknown;
}
