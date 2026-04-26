#ifndef RUNTIME_FUNCTION_REGISTRY_H
#define RUNTIME_FUNCTION_REGISTRY_H

#include "number.h"
#include "precision_context.h"

#include <functional>
#include <map>
#include <string>
#include <vector>

namespace runtime {

struct FunctionSpec {
    std::string name;
    std::vector<std::string> aliases;
    std::string category;
    std::string help;
    std::vector<std::string> autocomplete;
    std::string derivative_rule;
    std::string integration_hint;
    std::size_t min_arity = 0;
    std::size_t max_arity = 0;
    std::function<numeric::Number(const std::vector<numeric::Number>&,
                                  const numeric::PrecisionContext&)>
        numeric_eval;
};

class FunctionRegistry {
public:
    static const FunctionRegistry& builtins();

    void register_function(const FunctionSpec& spec);
    const FunctionSpec* find(const std::string& name) const;
    std::vector<FunctionSpec> functions() const;
    numeric::Number evaluate_numeric(const std::string& name,
                                     const std::vector<numeric::Number>& args,
                                     const numeric::PrecisionContext& context) const;

private:
    std::map<std::string, FunctionSpec> functions_;
    std::map<std::string, std::string> aliases_;
};

}  // namespace runtime

#endif
