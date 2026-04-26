#ifndef CAS_INTEGRATE_H
#define CAS_INTEGRATE_H

#include "expr.h"

#include <string>

namespace cas {

expression::Expr integrate(const expression::Expr& expr, const std::string& variable_name);

}  // namespace cas

#endif
