// ============================================================================
// 惰性 Token 流 - 按需生成 Token，避免全量词法分析
// ============================================================================
//
// 设计目标：
// 1. Token 按需生成，减少内存分配
// 2. 支持任意前瞻（lookahead）而不触发全量扫描
// 3. 支持检查点（checkpoint）机制用于高效回溯
// 4. 缓存已生成的 Token，避免重复词法分析
//
// 使用场景：
// - CommandParser: 命令解析需要前瞻判断命令类型
// - ExpressionLexer: 表达式解析需要运算符优先级判断
// - ScriptParser: 脚本解析需要块结构识别
// ============================================================================

#ifndef CORE_LAZY_TOKEN_STREAM_H
#define CORE_LAZY_TOKEN_STREAM_H

#include "parser/base_parser.h"
#include "parser/token_types.h"
#include <string_view>
#include <vector>
#include <cctype>
#include <stdexcept>

/**
 * @class LazyTokenStream
 * @brief 惰性 Token 流，按需生成 Token
 */
class LazyTokenStream : public BaseParser {
public:
    enum class WhitespaceMode {
        kCommand,
        kScript
    };

    explicit LazyTokenStream(std::string_view source,
                             WhitespaceMode whitespace_mode = WhitespaceMode::kCommand);

    // ========================================================================
    // Token 访问
    // ========================================================================

    /**
     * @brief 查看当前 Token（不消费）
     * @return 当前 Token 的引用
     *
     * 如果当前 Token 未生成，会触发词法分析。
     */
    const Token& peek();

    /**
     * @brief 查看指定偏移量的 Token（前瞻）
     * @param offset 前瞻偏移量（0 = 当前，1 = 下一个，...）
     * @return Token 的引用
     *
     * 支持任意前瞻深度，会按需生成 Token。
     */
    const Token& peek(std::size_t offset);

    /**
     * @brief 消费并返回当前 Token
     * @return 当前 Token
     *
     * 消费后，当前位置前进一格。
     */
    Token advance();

    /**
     * @brief 检查是否到达输入末尾
     * @return 如果当前 Token 是 kEnd 返回 true
     */
    bool is_at_end() const;

    /**
     * @brief 获取当前位置索引
     * @return 当前 Token 缓存位置
     */
    std::size_t position() const { return cache_pos_; }

    // ========================================================================
    // 回溯支持
    // ========================================================================

    /**
     * @struct Checkpoint
     * @brief 检查点，用于保存和恢复解析状态
     */
    struct Checkpoint {
        std::size_t cache_pos;      ///< Token 缓存位置
        std::size_t source_pos;     ///< 源字符串位置
        bool end_reached;           ///< 是否已到达末尾

        Checkpoint() : cache_pos(0), source_pos(0), end_reached(false) {}
    };

    /**
     * @brief 创建当前状态的检查点
     * @return 检查点对象
     *
     * 用于尝试性解析，失败后可恢复。
     */
    Checkpoint save_checkpoint() const;

    /**
     * @brief 恢复到指定检查点
     * @param cp 检查点对象
     *
     * 恢复后，可以重新从该位置开始解析。
     */
    void restore_checkpoint(const Checkpoint& cp);

    /**
     * @brief 重置到开始位置
     *
     * 清空缓存，重新从源字符串开头开始。
     */
    void reset();

    // ========================================================================
    // Token 匹配辅助
    // ========================================================================

    /**
     * @brief 检查当前 Token 类型
     * @param kind 期望的 Token 类型
     * @return 如果匹配返回 true
     */
    bool check(TokenKind kind);

    /**
     * @brief 匹配并消费指定类型的 Token
     * @param kind 期望的 Token 类型
     * @return 如果匹配并消费成功返回 true
     */
    bool match(TokenKind kind);

    /**
     * @brief 期望指定类型的 Token，失败抛出异常
     * @param kind 期望的 Token 类型
     * @param message 错误消息
     * @return 消费的 Token
     */
    Token expect(TokenKind kind, const char* message);

    // ========================================================================
    // 源字符串访问
    // ========================================================================

    /**
     * @brief 获取源字符串
     * @return 源字符串视图
     */
    std::string_view source() const { return source_; }

    /**
     * @brief 获取从指定位置开始的源字符串片段
     * @param start_pos 起始位置
     * @return 源字符串片段
     */
    std::string_view source_from(std::size_t start_pos) const;

private:
    // ========================================================================
    // 词法分析核心
    // ========================================================================

    /**
     * @brief 生成下一个 Token
     * @return 生成的 Token
     *
     * 从 source_pos_ 位置开始解析，生成一个 Token。
     */
    Token generate_next();

    /**
     * @brief 确保缓存中有足够的 Token
     * @param required_index 需要的最小缓存索引
     *
     * 按需生成 Token 直到缓存满足要求。
     */
    void ensure_cache_size(std::size_t required_index);
    void skip_ignorable_for_mode();

    // ========================================================================
    // Token 类型解析
    // ========================================================================

    Token parse_number_token();
    Token parse_string_token();
    Token parse_identifier_token();
    Token parse_operator_token();

    // ========================================================================
    // 辅助方法
    // ========================================================================

    void skip_spaces();

    bool is_at_source_end() const;

    char peek_char() const;
    char peek_next_char() const;

    bool is_identifier_start(char ch) const;
    bool is_identifier_char(char ch) const;
    bool is_digit_char(char ch) const;

    // ========================================================================
    // 成员变量
    // ========================================================================

    // source_ 继承自 BaseParser
    std::vector<Token> cache_;      ///< Token 缓存
    std::size_t cache_pos_ = 0;     ///< 当前缓存位置（消费进度）
    bool end_reached_ = false;      ///< 是否已生成 kEnd Token
    WhitespaceMode whitespace_mode_ = WhitespaceMode::kCommand;

    // 静态 kEnd Token，用于返回末尾引用
    static Token end_token_;
};

#endif // CORE_LAZY_TOKEN_STREAM_H
