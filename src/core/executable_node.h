#ifndef CORE_EXECUTABLE_NODE_H
#define CORE_EXECUTABLE_NODE_H

#include "types/stored_value.h"

class Calculator;

namespace execution {

/**
 * @class ExecutableNode
 * @brief Unified interface for any component that can be executed.
 */
class ExecutableNode {
public:
    virtual ~ExecutableNode() = default;

    /**
     * @brief Executes the node and returns a result.
     * @param calc The calculator instance providing context and services.
     * @param exact_mode Whether to evaluate in exact rational mode.
     * @return The resulting StoredValue.
     */
    virtual StoredValue execute(Calculator* calc, bool exact_mode) const = 0;
};

} // namespace execution

#endif
