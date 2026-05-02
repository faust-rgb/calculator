// ============================================================================
// 惰性 Token 流实现
// ============================================================================

#include "parser/lazy_token_stream.h"
#include "core/calculator_exceptions.h"
#include "math/helpers/base_conversions.h"
#include <sstream>
#include <algorithm>

// 静态成员定义
Token LazyTokenStream::end_token_;

// ============================================================================
// 构造函数
// ============================================================================

LazyTokenStream::LazyTokenStream(std::string_view source)
    : BaseParser(source, true) { // 启用注释跳过
    cache_.reserve(32);  // 预分配适量空间
}

// ============================================================================
// Token 访问
// ============================================================================

const Token& LazyTokenStream::peek() {
    return peek(0);
}

const Token& LazyTokenStream::peek(std::size_t offset) {
    ensure_cache_size(cache_pos_ + offset + 1);

    if (cache_pos_ + offset >= cache_.size()) {
        return end_token_;
    }

    return cache_[cache_pos_ + offset];
}

Token LazyTokenStream::advance() {
    ensure_cache_size(cache_pos_ + 1);

    if (cache_pos_ >= cache_.size()) {
        return end_token_;
    }

    return cache_[cache_pos_++];
}

bool LazyTokenStream::is_at_end() const {
    if (cache_pos_ < cache_.size()) {
        return cache_[cache_pos_].kind == TokenKind::kEnd;
    }
    return end_reached_;
}

// ============================================================================
// 回溯支持
// ============================================================================

LazyTokenStream::Checkpoint LazyTokenStream::save_checkpoint() const {
    Checkpoint cp;
    cp.cache_pos = cache_pos_;
    cp.source_pos = pos_;
    cp.end_reached = end_reached_;
    return cp;
}

void LazyTokenStream::restore_checkpoint(const Checkpoint& cp) {
    cache_pos_ = cp.cache_pos;
    // Token data beyond cache_pos_ is intentionally kept.  Rewinding the
    // underlying character cursor while keeping that cache makes the lexer
    // regenerate already-cached text once the cache is consumed, which corrupts
    // nested calls such as diff(abs(x)).
}

void LazyTokenStream::reset() {
    cache_.clear();
    cache_pos_ = 0;
    pos_ = 0;
    end_reached_ = false;
}

// ============================================================================
// Token 匹配辅助
// ============================================================================

bool LazyTokenStream::check(TokenKind kind) {
    return peek().kind == kind;
}

bool LazyTokenStream::match(TokenKind kind) {
    if (check(kind)) {
        advance();
        return true;
    }
    return false;
}

Token LazyTokenStream::expect(TokenKind kind, const char* message) {
    if (!check(kind)) {
        const Token& tok = peek();
        std::ostringstream oss;
        oss << message << " at position " << tok.position << "\n";
        oss << "  " << source_ << "\n";
        oss << "  " << std::string(tok.position, ' ') << "^";
        throw SyntaxError(oss.str());
    }
    return advance();
}

// ============================================================================
// 源字符串访问
// ============================================================================

std::string_view LazyTokenStream::source_from(std::size_t start_pos) const {
    if (start_pos >= source_.size()) {
        return {};
    }
    return source_.substr(start_pos);
}

// ============================================================================
// 词法分析核心
// ============================================================================

void LazyTokenStream::ensure_cache_size(std::size_t required_index) {
    while (cache_.size() < required_index && !end_reached_) {
        Token tok = generate_next();
        cache_.push_back(tok);
        if (tok.kind == TokenKind::kEnd) {
            end_reached_ = true;
        }
    }
}

Token LazyTokenStream::generate_next() {
    skip_ignorable();

    if (BaseParser::is_at_end()) {
        Token tok;
        tok.kind = TokenKind::kEnd;
        tok.position = pos_;
        return tok;
    }

    const char ch = BaseParser::peek();
    const std::size_t start_pos = pos_;

    // 字符串字面量
    if (ch == '"') {
        return parse_string_token();
    }

    // 数字字面量
    if (peek_is_digit() || (ch == '.' && std::isdigit(static_cast<unsigned char>(peek_next())))) {
        return parse_number_token();
    }

    // 标识符
    if (peek_is_identifier_start()) {
        return parse_identifier_token();
    }

    // 单字符 Token 和运算符
    switch (ch) {
        case '(': {
            pos_++;
            Token tok;
            tok.kind = TokenKind::kLParen;
            tok.text = source_.substr(start_pos, 1);
            tok.position = start_pos;
            return tok;
        }
        case ')': {
            pos_++;
            Token tok;
            tok.kind = TokenKind::kRParen;
            tok.text = source_.substr(start_pos, 1);
            tok.position = start_pos;
            return tok;
        }
        case '[': {
            pos_++;
            Token tok;
            tok.kind = TokenKind::kLBracket;
            tok.text = source_.substr(start_pos, 1);
            tok.position = start_pos;
            return tok;
        }
        case ']': {
            pos_++;
            Token tok;
            tok.kind = TokenKind::kRBracket;
            tok.text = source_.substr(start_pos, 1);
            tok.position = start_pos;
            return tok;
        }
        case '{': {
            pos_++;
            Token tok;
            tok.kind = TokenKind::kLBrace;
            tok.text = source_.substr(start_pos, 1);
            tok.position = start_pos;
            return tok;
        }
        case '}': {
            pos_++;
            Token tok;
            tok.kind = TokenKind::kRBrace;
            tok.text = source_.substr(start_pos, 1);
            tok.position = start_pos;
            return tok;
        }
        case ',': {
            pos_++;
            Token tok;
            tok.kind = TokenKind::kComma;
            tok.text = source_.substr(start_pos, 1);
            tok.position = start_pos;
            return tok;
        }
        case ';': {
            pos_++;
            Token tok;
            tok.kind = TokenKind::kSemicolon;
            tok.text = source_.substr(start_pos, 1);
            tok.position = start_pos;
            return tok;
        }
        case ':': {
            pos_++;
            Token tok;
            tok.kind = TokenKind::kColon;
            tok.text = source_.substr(start_pos, 1);
            tok.position = start_pos;
            return tok;
        }
        case '?': {
            pos_++;
            Token tok;
            tok.kind = TokenKind::kQuestion;
            tok.text = source_.substr(start_pos, 1);
            tok.position = start_pos;
            return tok;
        }
        case '=': {
            // 检查是否是 ==
            if (peek_next() == '=') {
                return parse_operator_token();
            }
            pos_++;
            Token tok;
            tok.kind = TokenKind::kEqual;
            tok.text = source_.substr(start_pos, 1);
            tok.position = start_pos;
            return tok;
        }
        case '+': case '-': case '*': case '/': case '^': case '%':
        case '<': case '>': case '!': case '&': case '|':
            return parse_operator_token();
        case '\n': {
            pos_++;
            Token tok;
            tok.kind = TokenKind::kNewline;
            tok.text = source_.substr(start_pos, 1);
            tok.position = start_pos;
            return tok;
        }
        default:
            // 未知字符，作为运算符处理
            return parse_operator_token();
    }
}

// ============================================================================
// Token 类型解析
// ============================================================================

Token LazyTokenStream::parse_number_token() {
    const std::size_t start = pos_;
    Token tok;
    tok.kind = TokenKind::kNumber;
    tok.text = parse_number_token_view();
    tok.position = start;
    // 解析数值
    try {
        std::string num_str(tok.text);
        if (num_str.size() > 2 && num_str[0] == '0' &&
            std::isalpha(static_cast<unsigned char>(num_str[1]))) {
            tok.number_value = static_cast<double>(parse_prefixed_integer_token(num_str));
        } else {
            tok.number_value = std::stod(num_str);
        }
    } catch (...) {
        std::ostringstream oss;
        oss << "invalid number literal at position " << tok.position << "\n";
        oss << "  " << source_ << "\n";
        oss << "  " << std::string(tok.position, ' ') << "^";
        throw SyntaxError(oss.str());
    }
    return tok;
}

Token LazyTokenStream::parse_string_token() {
    const std::size_t start = pos_;
    Token tok;
    tok.kind = TokenKind::kString;
    tok.string_value = parse_string_literal();
    tok.text = source_.substr(start, pos_ - start);
    tok.position = start;
    return tok;
}

Token LazyTokenStream::parse_identifier_token() {
    const std::size_t start = pos_;
    Token tok;
    tok.kind = TokenKind::kIdentifier;
    tok.text = parse_identifier();
    tok.position = start;
    return tok;
}

Token LazyTokenStream::parse_operator_token() {
    const std::size_t start = pos_;
    if (BaseParser::is_at_end()) {
        Token tok;
        tok.kind = TokenKind::kEnd;
        tok.position = start;
        return tok;
    }
    const char ch = BaseParser::peek();
    pos_++;  // 消费字符

    // 检查双字符运算符
    if (!BaseParser::is_at_end()) {
        const char next = BaseParser::peek();
        if ((ch == '=' && next == '=') ||
            (ch == '!' && next == '=') ||
            (ch == '<' && next == '=') ||
            (ch == '>' && next == '=') ||
            (ch == '&' && next == '&') ||
            (ch == '|' && next == '|')) {
            pos_++;
        }
    }

    Token tok;
    tok.kind = TokenKind::kOperator;
    tok.text = source_.substr(start, pos_ - start);
    tok.position = start;
    return tok;
}
