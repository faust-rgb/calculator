#ifndef RUNTIME_ENVIRONMENT_H
#define RUNTIME_ENVIRONMENT_H

#include "precision_context.h"
#include "value.h"

#include <map>
#include <string>

namespace runtime {

struct Binding {
    Value value;
    expression::Expr original_expr;
    bool has_original_expr = false;
};

class Environment {
public:
    numeric::PrecisionContext& precision();
    const numeric::PrecisionContext& precision() const;

    void set(const std::string& name, const Value& value);
    void set(const std::string& name, const Value& value, const expression::Expr& original);
    bool lookup(const std::string& name, Value* value) const;
    const std::map<std::string, Binding>& variables() const;

private:
    std::map<std::string, Binding> variables_;
    numeric::PrecisionContext precision_;
};

}  // namespace runtime

#endif
