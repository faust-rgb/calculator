// ============================================================================
// 统一 Token 类型定义
// ============================================================================
//
// 提供通用的 Token 类型和词法分析基础结构，供多个解析器复用：
// - CommandParser: 命令解析
// - ExpressionLexer: 表达式词法分析
// - 其他需要词法分析的模块
//
// 设计目标：
// - 避免重复定义 Token 类型
// - 统一 Token 分类，便于跨模块使用
// ============================================================================

#ifndef CORE_TOKEN_TYPES_H
#define CORE_TOKEN_TYPES_H

#include <string>
#include <string_view>

// ============================================================================
// 通用 Token 类型枚举
// ============================================================================

/**
 * @enum TokenKind
 * @brief 通用 Token 类型分类
 *
 * 涵盖命令解析和表达式解析所需的所有 Token 类型。
 */
enum class TokenKind {
    // 特殊标记
    kEnd,           ///< 输入结束
    kError,         ///< 词法错误

    // 字面量
    kNumber,        ///< 数字字面量（保留原始文本）
    kString,        ///< 字符串字面量（已解码转义）
    kIdentifier,    ///< 标识符（变量名、函数名）

    // 括号
    kLParen,        ///< 左圆括号 (
    kRParen,        ///< 右圆括号 )
    kLBracket,      ///< 左方括号 [
    kRBracket,      ///< 右方括号 ]
    kLBrace,        ///< 左花括号 {
    kRBrace,        ///< 右花括号 }

    // 分隔符
    kComma,         ///< 逗号 ,
    kSemicolon,     ///< 分号 ;
    kColon,         ///< 冒号 :

    // 运算符
    kOperator,      ///< 运算符 + - * / ^ % 等
    kEqual,         ///< 等号 =（赋值或定义）
    kQuestion,      ///< 问号 ?（三元运算符）

    // 特殊
    kNewline,       ///< 换行符（脚本解析用）
};

// ============================================================================
// 通用 Token 结构
// ============================================================================

/**
 * @struct Token
 * @brief 通用 Token 结构
 *
 * 数字 Token 保留原始文本，不进行数值转换，避免精度丢失。
 * 数值转换延迟到求值阶段。
 */
struct Token {
    TokenKind kind = TokenKind::kEnd;
    std::string_view text;      ///< 原始文本视图
    std::string string_value;   ///< 字符串值（仅 kString 类型使用）
    double number_value = 0.0;  ///< 数值（仅 kNumber 类型，可选）
    std::size_t position = 0;   ///< 在源字符串中的起始位置

    Token() = default;

    static Token make_end(std::size_t pos) {
        Token t;
        t.kind = TokenKind::kEnd;
        t.position = pos;
        return t;
    }

    static Token make_error(std::string_view msg, std::size_t pos) {
        Token t;
        t.kind = TokenKind::kError;
        t.text = msg;
        t.position = pos;
        return t;
    }
};

// ============================================================================
// Token 类型判断辅助函数
// ============================================================================

inline bool is_literal_token(TokenKind kind) {
    return kind == TokenKind::kNumber ||
           kind == TokenKind::kString ||
           kind == TokenKind::kIdentifier;
}

inline bool is_bracket_token(TokenKind kind) {
    return kind == TokenKind::kLParen || kind == TokenKind::kRParen ||
           kind == TokenKind::kLBracket || kind == TokenKind::kRBracket ||
           kind == TokenKind::kLBrace || kind == TokenKind::kRBrace;
}

inline bool is_opening_bracket(TokenKind kind) {
    return kind == TokenKind::kLParen ||
           kind == TokenKind::kLBracket ||
           kind == TokenKind::kLBrace;
}

inline bool is_closing_bracket(TokenKind kind) {
    return kind == TokenKind::kRParen ||
           kind == TokenKind::kRBracket ||
           kind == TokenKind::kRBrace;
}

#endif // CORE_TOKEN_TYPES_H
