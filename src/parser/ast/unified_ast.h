/**
 * @file unified_ast.h
 * @brief 统一 AST 节点与求值器 (Phase 2 Unified AST & Pipeline)
 */

#ifndef PARSER_UNIFIED_AST_H
#define PARSER_UNIFIED_AST_H

#include "types/stored_value.h"
#include "core/execution_context.h"
#include <memory>
#include <string>
#include <vector>
#include <stdexcept>

namespace core {

/**
 * @struct SourceLocation
 * @brief 源码文件与行列号信息
 */
struct SourceLocation {
    std::string file;
    std::size_t line = 0;
    std::size_t column = 0;

    std::string to_string() const {
        if (file.empty()) {
            return "line " + std::to_string(line) + ", col " + std::to_string(column);
        }
        return file + ":" + std::to_string(line) + ":" + std::to_string(column);
    }
};

/**
 * @class ScriptRuntimeError
 * @brief 带有行列号位置信息的运行时异常
 */
class ScriptRuntimeError : public std::runtime_error {
public:
    ScriptRuntimeError(const std::string& msg, SourceLocation loc)
        : std::runtime_error(loc.to_string() + ": " + msg), location_(std::move(loc)) {}
    const SourceLocation& location() const { return location_; }
private:
    SourceLocation location_;
};

enum class ASTNodeType {
    kScalar,
    kString,
    kVariable,
    kUnaryOp,
    kBinaryOp,
    kFunctionCall,
    kMatrixLiteral,
    kListLiteral,
    kDictLiteral,
    kSlice
};

class ASTNode {
public:
    virtual ~ASTNode() = default;
    virtual ASTNodeType type() const = 0;
    virtual StoredValue evaluate(ExecutionContext& ctx) const = 0;

    void set_location(SourceLocation loc) { location_ = std::move(loc); }
    const SourceLocation& location() const { return location_; }
protected:
    SourceLocation location_;
};

class ScalarNode : public ASTNode {
public:
    explicit ScalarNode(StoredValue val) : val_(std::move(val)) {}
    ASTNodeType type() const override { return ASTNodeType::kScalar; }
    StoredValue evaluate(ExecutionContext& ctx) const override {
        (void)ctx;
        return val_;
    }
private:
    StoredValue val_;
};

class StringNode : public ASTNode {
public:
    explicit StringNode(std::string str) : str_(std::move(str)) {}
    ASTNodeType type() const override { return ASTNodeType::kString; }
    StoredValue evaluate(ExecutionContext& ctx) const override {
        (void)ctx;
        return StoredValue(str_);
    }
private:
    std::string str_;
};

class VariableNode : public ASTNode {
public:
    explicit VariableNode(std::string name) : name_(std::move(name)) {}
    ASTNodeType type() const override { return ASTNodeType::kVariable; }
    StoredValue evaluate(ExecutionContext& ctx) const override;
    const std::string& name() const { return name_; }
private:
    std::string name_;
};

class UnaryOpNode : public ASTNode {
public:
    UnaryOpNode(char op, std::unique_ptr<ASTNode> child)
        : op_(op), child_(std::move(child)) {}
    ASTNodeType type() const override { return ASTNodeType::kUnaryOp; }
    StoredValue evaluate(ExecutionContext& ctx) const override;
private:
    char op_;
    std::unique_ptr<ASTNode> child_;
};

class BinaryOpNode : public ASTNode {
public:
    BinaryOpNode(std::string op, std::unique_ptr<ASTNode> left, std::unique_ptr<ASTNode> right)
        : op_(std::move(op)), left_(std::move(left)), right_(std::move(right)) {}
    ASTNodeType type() const override { return ASTNodeType::kBinaryOp; }
    StoredValue evaluate(ExecutionContext& ctx) const override;
private:
    std::string op_;
    std::unique_ptr<ASTNode> left_;
    std::unique_ptr<ASTNode> right_;
};

class FunctionCallNode : public ASTNode {
public:
    FunctionCallNode(std::string name, std::vector<std::unique_ptr<ASTNode>> args)
        : name_(std::move(name)), args_(std::move(args)) {}
    ASTNodeType type() const override { return ASTNodeType::kFunctionCall; }
    StoredValue evaluate(ExecutionContext& ctx) const override;
private:
    std::string name_;
    std::vector<std::unique_ptr<ASTNode>> args_;
};

class MatrixLiteralNode : public ASTNode {
public:
    explicit MatrixLiteralNode(std::vector<std::vector<std::unique_ptr<ASTNode>>> rows)
        : rows_(std::move(rows)) {}
    ASTNodeType type() const override { return ASTNodeType::kMatrixLiteral; }
    StoredValue evaluate(ExecutionContext& ctx) const override;
private:
    std::vector<std::vector<std::unique_ptr<ASTNode>>> rows_;
};

class ListLiteralNode : public ASTNode {
public:
    explicit ListLiteralNode(std::vector<std::unique_ptr<ASTNode>> elements)
        : elements_(std::move(elements)) {}
    ASTNodeType type() const override { return ASTNodeType::kListLiteral; }
    StoredValue evaluate(ExecutionContext& ctx) const override;
private:
    std::vector<std::unique_ptr<ASTNode>> elements_;
};

class DictLiteralNode : public ASTNode {
public:
    explicit DictLiteralNode(std::vector<std::pair<std::unique_ptr<ASTNode>, std::unique_ptr<ASTNode>>> pairs)
        : pairs_(std::move(pairs)) {}
    ASTNodeType type() const override { return ASTNodeType::kDictLiteral; }
    StoredValue evaluate(ExecutionContext& ctx) const override;
private:
    std::vector<std::pair<std::unique_ptr<ASTNode>, std::unique_ptr<ASTNode>>> pairs_;
};

class SliceNode : public ASTNode {
public:
    SliceNode(std::unique_ptr<ASTNode> start, std::unique_ptr<ASTNode> stop, std::unique_ptr<ASTNode> step)
        : start_(std::move(start)), stop_(std::move(stop)), step_(std::move(step)) {}
    ASTNodeType type() const override { return ASTNodeType::kSlice; }
    StoredValue evaluate(ExecutionContext& ctx) const override;
private:
    std::unique_ptr<ASTNode> start_;
    std::unique_ptr<ASTNode> stop_;
    std::unique_ptr<ASTNode> step_;
};

/**
 * @class ASTEvaluator
 * @brief 统一 AST 求值调度器
 */
class ASTEvaluator {
public:
    static StoredValue evaluate(const ASTNode& node, ExecutionContext& ctx) {
        return node.evaluate(ctx);
    }
};

} // namespace core

struct ExpressionAST;

namespace core {

/**
 * @brief 从 legacy ExpressionAST 转换构建 Unified AST 节点
 */
std::unique_ptr<ASTNode> build_unified_ast(const ExpressionAST* ast);

/**
 * @brief 使用 ASTEvaluator 和 ExecutionContext 评估 AST 节点
 */
StoredValue evaluate_unified_ast(const ExpressionAST* ast, ExecutionContext& ctx);

} // namespace core

#endif // PARSER_UNIFIED_AST_H
