// ============================================================================
// 脚本信号头文件
// ============================================================================
//
// 本文件定义了 ScriptSignal 结构体，用于实现脚本执行过程中的控制流。
//
// 主要功能：
// - 支持 return、break、continue 语句
// - 语句执行返回信号，外层控制结构根据信号类型决定行为
// - 可携带返回值（用于 return 语句）
// ============================================================================

#ifndef SCRIPT_SIGNAL_H
#define SCRIPT_SIGNAL_H

#include "core/value/stored_value.h"

// ============================================================================
// ScriptSignal 结构体定义
// ============================================================================

/**
 * @struct ScriptSignal
 * @brief 脚本执行的控制流信号
 *
 * 在脚本执行过程中，控制流语句（return、break、continue）会产生信号。
 * 外层控制结构（如循环、函数调用）根据信号类型执行相应操作。
 */
struct ScriptSignal {
    /**
     * @brief 信号类型枚举
     */
    enum class Kind {
        kNone,      ///< 无信号，正常执行
        kReturn,    ///< return 语句，需要返回函数值
        kBreak,     ///< break 语句，跳出当前循环
        kContinue,  ///< continue 语句，跳过当前迭代
    };

    Kind kind = Kind::kNone;  ///< 信号类型
    bool has_value = false;   ///< 是否有返回值（仅 return 使用）
    StoredValue value;        ///< 返回值（return 时使用）

    // ========================================================================
    // 静态工厂方法
    // ========================================================================

    /**
     * @brief 创建 return 信号
     * @param return_value 返回值
     * @return 带返回值的 return 信号
     */
    static ScriptSignal make_return(const StoredValue& return_value);

    /**
     * @brief 创建 break 信号
     * @return break 信号
     */
    static ScriptSignal make_break();

    /**
     * @brief 创建 continue 信号
     * @return continue 信号
     */
    static ScriptSignal make_continue();
};

#endif // SCRIPT_SIGNAL_H
