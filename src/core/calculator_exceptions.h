#ifndef CALCULATOR_EXCEPTIONS_H
#define CALCULATOR_EXCEPTIONS_H

#include <stdexcept>
#include <string>

/**
 * @brief 计算器基础异常类
 */
class CalculatorError : public std::runtime_error {
public:
    explicit CalculatorError(const std::string& message) : std::runtime_error(message) {}
};

/**
 * @brief 语法错误异常
 */
class SyntaxError : public CalculatorError {
public:
    explicit SyntaxError(const std::string& message) : CalculatorError("Syntax error: " + message) {}
};

/**
 * @brief 数学领域错误异常（如除以零、对负数取对数等）
 */
class MathError : public CalculatorError {
public:
    explicit MathError(const std::string& message) : CalculatorError("Math error: " + message) {}
};

/**
 * @brief 变量或函数未定义错误
 */
class UndefinedError : public CalculatorError {
public:
    explicit UndefinedError(const std::string& message) : CalculatorError("Undefined error: " + message) {}
};

/**
 * @brief 矩阵维度不匹配错误
 */
class DimensionError : public CalculatorError {
public:
    explicit DimensionError(const std::string& message) : CalculatorError("Dimension error: " + message) {}
};

/**
 * @brief 命令参数错误
 */
class ArgumentError : public CalculatorError {
public:
    explicit ArgumentError(const std::string& message) : CalculatorError("Argument error: " + message) {}
};

/**
 * @brief 精确模式不支持异常
 */
class ExactModeUnsupported : public CalculatorError {
public:
    explicit ExactModeUnsupported(const std::string& message)
        : CalculatorError("Exact mode unsupported: " + message) {}
};

/**
 * @brief 精确小数不支持异常
 */
class PreciseDecimalUnsupported : public CalculatorError {
public:
    explicit PreciseDecimalUnsupported(const std::string& message)
        : CalculatorError("Precise decimal unsupported: " + message) {}
};

#endif