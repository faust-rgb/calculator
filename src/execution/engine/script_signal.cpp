// ============================================================================
// 脚本信号实现文件
// ============================================================================
//
// 本文件实现了 ScriptSignal 结构体的静态工厂方法。
// 用于创建不同类型的控制流信号。
// ============================================================================

#include "script_signal.h"

// ============================================================================
// 静态工厂方法实现
// ============================================================================

/**
 * @brief 创建带有返回值的 return 信号
 * @param return_value 函数返回值
 * @return 配置好的 return 信号
 */
ScriptSignal ScriptSignal::make_return(const StoredValue& return_value) {
    ScriptSignal signal;
    signal.kind = Kind::kReturn;
    signal.has_value = true;
    signal.value = return_value;
    return signal;
}

/**
 * @brief 创建 break 信号
 * @return 配置好的 break 信号
 */
ScriptSignal ScriptSignal::make_break() {
    ScriptSignal signal;
    signal.kind = Kind::kBreak;
    return signal;
}

/**
 * @brief 创建 continue 信号
 * @return 配置好的 continue 信号
 */
ScriptSignal ScriptSignal::make_continue() {
    ScriptSignal signal;
    signal.kind = Kind::kContinue;
    return signal;
}
