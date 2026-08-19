#ifndef SYMBOLIC_LAPLACE_NORMALIZER_H
#define SYMBOLIC_LAPLACE_NORMALIZER_H

#include "symbolic/core/symbolic_expression.h"

#include <string>
#include <vector>

namespace symbolic_expression_internal {

struct TransformContext {
    std::string time_variable;
    std::string transform_variable;
    int recursion_depth = 0;
    int max_recursion_depth = 256;
};

class TransformContextScope {
public:
    explicit TransformContextScope(const TransformContext& context);
    ~TransformContextScope();

    TransformContextScope(const TransformContextScope&) = delete;
    TransformContextScope& operator=(const TransformContextScope&) = delete;

private:
    int previous_max_depth_;
};

int current_transform_max_depth();

// A deliberately small intermediate form shared by Laplace rules.  The
// numerator and denominator remain symbolic expressions so parameters are not
// forced through Scalar.
struct RationalForm {
    SymbolicExpression numerator;
    SymbolicExpression denominator;
};

bool depends_on(const SymbolicExpression& expression,
                const std::string& variable_name);

bool contains_function(const SymbolicExpression& expression,
                       const std::string& function_name);

void flatten_additive_terms(const SymbolicExpression& expression,
                            std::vector<SymbolicExpression>* terms);

void flatten_multiplicative_factors(const SymbolicExpression& expression,
                                    std::vector<SymbolicExpression>* factors);

bool normalize_rational_form(const SymbolicExpression& expression,
                             const std::string& variable_name,
                             RationalForm* result);

}  // namespace symbolic_expression_internal

#endif
