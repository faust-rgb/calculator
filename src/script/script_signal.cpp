#include "script_signal.h"

ScriptSignal ScriptSignal::make_return(const StoredValue& return_value) {
    ScriptSignal signal;
    signal.kind = Kind::kReturn;
    signal.has_value = true;
    signal.value = return_value;
    return signal;
}

ScriptSignal ScriptSignal::make_break() {
    ScriptSignal signal;
    signal.kind = Kind::kBreak;
    return signal;
}

ScriptSignal ScriptSignal::make_continue() {
    ScriptSignal signal;
    signal.kind = Kind::kContinue;
    return signal;
}
