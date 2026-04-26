#ifndef EXPRESSION_EVALUATOR_H
#define EXPRESSION_EVALUATOR_H

#include "environment.h"
#include "expr.h"
#include "value.h"

namespace expression {

runtime::Value evaluate(const Expr& expr, runtime::Environment& env);

}  // namespace expression

#endif
