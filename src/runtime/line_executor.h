#ifndef RUNTIME_LINE_EXECUTOR_H
#define RUNTIME_LINE_EXECUTOR_H

#include "environment.h"
#include "value.h"

#include <string>

namespace runtime {

struct LineResult {
    bool assigned = false;
    std::string name;
    Value value;

    std::string to_string() const;
};

LineResult evaluate_line(const std::string& line, Environment& env);
std::string list_variables(const Environment& env);

}  // namespace runtime

#endif
