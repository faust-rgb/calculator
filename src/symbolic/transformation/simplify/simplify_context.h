#ifndef SYMBOLIC_SIMPLIFY_CONTEXT_H
#define SYMBOLIC_SIMPLIFY_CONTEXT_H

#include "symbolic/core/symbolic_expression.h"
#include <cstddef>
#include <string>
#include <vector>

enum class SimplifyMode {
    kSafe,
    kExpand,
    kFactor,
    kTogether,
    kRealDomain,
    kComplexPrincipal
};

struct DomainCondition {
    std::vector<SymbolicExpression> nonzero;
    std::vector<SymbolicExpression> positive;
    std::vector<SymbolicExpression> real;

    bool empty() const {
        return nonzero.empty() && positive.empty() && real.empty();
    }
};

struct SimplifyContext {
    SimplifyMode mode = SimplifyMode::kSafe;
    DomainCondition domain;
    std::size_t max_nodes = 10000;
    std::size_t max_steps = 64;
    bool preserve_domain = true;
    bool trace_rules = false;
};

inline const char* simplify_mode_name(SimplifyMode mode) {
    switch (mode) {
        case SimplifyMode::kSafe: return "safe";
        case SimplifyMode::kExpand: return "expand";
        case SimplifyMode::kFactor: return "factor";
        case SimplifyMode::kTogether: return "together";
        case SimplifyMode::kRealDomain: return "real";
        case SimplifyMode::kComplexPrincipal: return "complex-principal";
    }
    return "safe";
}

#endif
