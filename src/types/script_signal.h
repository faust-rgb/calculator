// ============================================================================
// 脚本控制流信号
// ============================================================================
//
// 用于实现 return、break、continue 语句。
// 语句执行返回信号，外层控制结构根据信号类型决定行为。

#ifndef TYPES_SCRIPT_SIGNAL_H
#define TYPES_SCRIPT_SIGNAL_H

#include "stored_value.h"

/**
 * @struct ScriptSignal
 * @brief 脚本执行的控制流信号
 */
struct ScriptSignal {
    enum class Kind {
        kNone,      ///< 无信号，正常执行
        kReturn,    ///< return 语句
        kBreak,     ///< break 语句
        kContinue,  ///< continue 语句
    };

    Kind kind = Kind::kNone;  ///< 信号类型
    bool has_value = false;   ///< 是否有返回值
    StoredValue value;        ///< 返回值（return 时使用）

    /** @brief 创建 return 信号 */
    static ScriptSignal make_return(const StoredValue& return_value);

    /** @brief 创建 break 信号 */
    static ScriptSignal make_break();

    /** @brief 创建 continue 信号 */
    static ScriptSignal make_continue();
};

#endif // TYPES_SCRIPT_SIGNAL_H
