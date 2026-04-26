#ifndef NUMERIC_PRECISION_CONTEXT_H
#define NUMERIC_PRECISION_CONTEXT_H

namespace numeric {

enum class RoundingMode {
    Nearest,
    Floor,
    Ceil,
    Zero,
};

struct PrecisionContext {
    int digits = 50;
    RoundingMode rounding = RoundingMode::Nearest;
    int max_iterations = 1000;
    bool exact_preferred = true;
};

}  // namespace numeric

#endif
