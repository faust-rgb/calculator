#include "expression_compiler.h"
#include "core/string_utils.h"
#include "math/helpers/base_conversions.h"
#include "parser/unified_parser_factory.h"
#include "parser/function_categories.h"

#include <algorithm>
#include <cctype>
#include <stdexcept>

// ============================================================================
// 表达式分析函数（委托给 UnifiedParserFactory）
// ============================================================================

namespace {
    UnifiedParserFactory& get_global_factory() {
        static UnifiedParserFactory factory;
        return factory;
    }
}

ExpressionFeature analyze_expression_features(const std::string& expression) {
    return get_global_factory().analyze_features(expression);
}

ExpressionHint analyze_expression_hint(const std::string& expression) {
    return get_global_factory().analyze(expression).hint;
}
