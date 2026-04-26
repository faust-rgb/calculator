#ifndef CAS_DIFFERENTIATE_H
#define CAS_DIFFERENTIATE_H

#include "expr.h"

#include <string>

namespace cas {

expression::Expr differentiate(const expression::Expr& expr, const std::string& variable_name);

}  // namespace cas

#endif
