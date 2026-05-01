// ============================================================================
// 命令解析器实现
// ============================================================================

#include "command_parser.h"
#include "calculator_exceptions.h"
#include "utils.h"
#include <algorithm>
#include <cctype>

// ============================================================================
// CommandASTNode 实现
// ============================================================================

CommandASTNode CommandASTNode::make_empty() {
    CommandASTNode node;
    node.kind = CommandKind::kEmpty;
    node.data = std::monostate{};
    return node;
}

CommandASTNode CommandASTNode::make_meta_command(const std::string& cmd,
                                                  const std::vector<std::string>& args) {
    CommandASTNode node;
    node.kind = CommandKind::kMetaCommand;
    node.data = MetaCommandInfo{cmd, args};
    return node;
}

CommandASTNode CommandASTNode::make_function_definition(const std::string& name,
                                                        const std::string& param,
                                                        const std::string& body) {
    CommandASTNode node;
    node.kind = CommandKind::kFunctionDefinition;
    node.data = FunctionDefinitionInfo{name, param, body};
    return node;
}

CommandASTNode CommandASTNode::make_function_call(const std::string& name,
                                                  const std::vector<std::string>& args) {
    CommandASTNode node;
    node.kind = CommandKind::kFunctionCall;
    node.data = FunctionCallInfo{name, args};
    return node;
}

CommandASTNode CommandASTNode::make_assignment(const std::string& var,
                                               const std::string& expr) {
    CommandASTNode node;
    node.kind = CommandKind::kAssignment;
    node.data = AssignmentInfo{var, expr};
    return node;
}

CommandASTNode CommandASTNode::make_expression(const std::string& expr) {
    CommandASTNode node;
    node.kind = CommandKind::kExpression;
    node.data = expr;
    return node;
}

CommandASTNode CommandASTNode::make_string_literal(const std::string& value) {
    CommandASTNode node;
    node.kind = CommandKind::kStringLiteral;
    node.data = value;
    return node;
}

const MetaCommandInfo* CommandASTNode::as_meta_command() const {
    if (kind == CommandKind::kMetaCommand) {
        return &std::get<MetaCommandInfo>(data);
    }
    return nullptr;
}

const FunctionDefinitionInfo* CommandASTNode::as_function_definition() const {
    if (kind == CommandKind::kFunctionDefinition) {
        return &std::get<FunctionDefinitionInfo>(data);
    }
    return nullptr;
}

const FunctionCallInfo* CommandASTNode::as_function_call() const {
    if (kind == CommandKind::kFunctionCall) {
        return &std::get<FunctionCallInfo>(data);
    }
    return nullptr;
}

const AssignmentInfo* CommandASTNode::as_assignment() const {
    if (kind == CommandKind::kAssignment) {
        return &std::get<AssignmentInfo>(data);
    }
    return nullptr;
}

const std::string* CommandASTNode::as_expression() const {
    if (kind == CommandKind::kExpression) {
        return &std::get<std::string>(data);
    }
    return nullptr;
}

const std::string* CommandASTNode::as_string_literal() const {
    if (kind == CommandKind::kStringLiteral) {
        return &std::get<std::string>(data);
    }
    return nullptr;
}

// ============================================================================
// 辅助函数声明（来自 utils.h）
// ============================================================================

// split_function_definition 在 utils.h 中声明
// 我们在这里使用它

// ============================================================================
// CommandParser 实现
// ============================================================================

CommandParser::CommandParser(std::string source)
    : BaseParser(std::move(source)) {}

CommandASTNode CommandParser::parse() {
    skip_spaces();

    // 1. 空输入
    if (is_at_end()) {
        return CommandASTNode::make_empty();
    }

    // 2. 字符串字面量
    if (peek() == '"') {
        std::string content = parse_string_literal_content();
        return CommandASTNode::make_string_literal(content);
    }

    // 3. 元命令（以冒号开头）
    if (peek() == ':') {
        match(':');
        std::string command = parse_identifier();
        skip_spaces();

        std::vector<std::string> arguments;
        if (!is_at_end()) {
            // 元命令参数：剩余全部内容作为一个参数，或按逗号分隔
            std::string remaining = source_.substr(pos_);
            arguments = split_command_arguments(remaining);
        }
        return CommandASTNode::make_meta_command(command, arguments);
    }

    // 4. 检测顶层等号，判断是函数定义还是赋值
    // 使用 utils.cpp 中已有的 split_assignment 逻辑
    std::string lhs;
    std::string rhs;
    if (split_assignment(source_, &lhs, &rhs)) {
        // 检查 lhs 是否为函数定义形式 f(x)
        std::string func_name;
        std::string param_name;
        std::string body;
        if (split_function_definition(source_, &func_name, &param_name, &body)) {
            return CommandASTNode::make_function_definition(func_name, param_name, body);
        }

        // 否则是简单赋值
        return CommandASTNode::make_assignment(lhs, rhs);
    }

    // 5. 解析函数调用或表达式
    return parse_function_call_or_expression();
}


CommandASTNode CommandParser::parse_function_call_or_expression() {
    skip_spaces();

    // 检查是否为标识符开头
    if (!peek_is_identifier_start()) {
        // 不是标识符，直接作为表达式
        std::string expr = trim_copy(source_);
        return CommandASTNode::make_expression(expr);
    }

    // 解析标识符
    std::string name = parse_identifier();
    skip_spaces();

    // 检查是否为函数调用
    if (match('(')) {
        // 函数调用
        std::vector<std::string> arguments = parse_argument_list();
        skip_spaces();
        expect(')');

        // 检查后面是否还有内容
        skip_spaces();
        if (!is_at_end()) {
            // 后面还有内容，整体作为表达式
            // 例如：sin(x) + cos(x)
            std::string expr = trim_copy(source_);
            return CommandASTNode::make_expression(expr);
        }

        return CommandASTNode::make_function_call(name, arguments);
    }

    // 不是函数调用，检查是否为赋值
    skip_spaces();
    if (match('=')) {
        skip_spaces();
        // 检查是否为比较运算符
        if (peek() == '=') {
            // == 比较运算符，整体作为表达式
            std::string expr = trim_copy(source_);
            return CommandASTNode::make_expression(expr);
        }
        // 赋值
        std::string expr = trim_copy(source_.substr(pos_));
        return CommandASTNode::make_assignment(name, expr);
    }

    // 不是赋值，检查后面是否还有内容
    skip_spaces();
    if (!is_at_end()) {
        // 后面还有内容，整体作为表达式
        std::string expr = trim_copy(source_);
        return CommandASTNode::make_expression(expr);
    }

    // 纯变量引用，作为表达式
    return CommandASTNode::make_expression(name);
}

std::vector<std::string> CommandParser::parse_argument_list() {
    std::vector<std::string> arguments;
    skip_spaces();

    if (peek(')')) {
        return arguments;
    }

    while (true) {
        skip_spaces();

        // 解析一个参数（到逗号或右括号）
        std::string arg;
        int paren_depth = 0;
        int bracket_depth = 0;
        bool in_string = false;
        bool escaping = false;

        while (!is_at_end()) {
            char ch = peek();

            if (in_string) {
                arg.push_back(ch);
                if (escaping) {
                    escaping = false;
                } else if (ch == '\\') {
                    escaping = true;
                } else if (ch == '"') {
                    in_string = false;
                }
                ++pos_;
                continue;
            }

            if (ch == '"') {
                in_string = true;
                arg.push_back(ch);
                ++pos_;
                continue;
            }

            if (ch == '(') {
                ++paren_depth;
                arg.push_back(ch);
                ++pos_;
                continue;
            }

            if (ch == '[') {
                ++bracket_depth;
                arg.push_back(ch);
                ++pos_;
                continue;
            }

            if (ch == ')') {
                if (paren_depth == 0 && bracket_depth == 0) {
                    // 参数结束
                    break;
                }
                --paren_depth;
                arg.push_back(ch);
                ++pos_;
                continue;
            }

            if (ch == ']') {
                if (bracket_depth > 0) {
                    --bracket_depth;
                }
                arg.push_back(ch);
                ++pos_;
                continue;
            }

            if (ch == ',' && paren_depth == 0 && bracket_depth == 0) {
                // 参数分隔符
                break;
            }

            arg.push_back(ch);
            ++pos_;
        }

        arguments.push_back(trim_copy(arg));

        skip_spaces();
        if (!match(',')) {
            break;
        }
    }

    return arguments;
}

std::string CommandParser::parse_string_literal_content() {
    expect('"');

    std::string content;
    bool escaping = false;

    while (!is_at_end()) {
        char ch = peek();

        if (escaping) {
            if (ch == 'n') {
                content.push_back('\n');
            } else if (ch == 't') {
                content.push_back('\t');
            } else if (ch == '\\') {
                content.push_back('\\');
            } else if (ch == '"') {
                content.push_back('"');
            } else {
                content.push_back(ch);
            }
            escaping = false;
            ++pos_;
            continue;
        }

        if (ch == '\\') {
            escaping = true;
            ++pos_;
            continue;
        }

        if (ch == '"') {
            ++pos_;
            break;
        }

        content.push_back(ch);
        ++pos_;
    }

    return content;
}

// ============================================================================
// 便捷函数实现
// ============================================================================

CommandASTNode parse_command(const std::string& source) {
    CommandParser parser(source);
    return parser.parse();
}

std::vector<std::string> split_command_arguments(const std::string& text) {
    // 复用 utils.cpp 中的 split_top_level_arguments
    return split_top_level_arguments(text);
}