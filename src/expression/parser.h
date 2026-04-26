#ifndef EXPRESSION_PARSER_H
#define EXPRESSION_PARSER_H

#include "expr.h"

#include <string>

namespace expression {

Expr parse_expression(const std::string& text);

}  // namespace expression

#endif
