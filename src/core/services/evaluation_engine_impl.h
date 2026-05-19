// ============================================================================
// evaluation_engine_impl.h - IEvaluationEngine 的具体实现
// ============================================================================

#ifndef EVALUATION_ENGINE_IMPL_H
#define EVALUATION_ENGINE_IMPL_H

#include "core/services/core_manager_interfaces.h"
#include "core/api/calculator.h"
#include "core/services/service_locator.h"
#include "analysis/calculus/function_analysis.h"
#include "symbolic/core/symbolic_expression.h"
#include "symbolic/modules/symbolic_module.h"
#include "execution/engine/inline_expander.h"

class EvaluationEngineImpl : public IEvaluationEngine {
public:
    EvaluationEngineImpl(Calculator* calculator, Calculator::Impl* impl)
        : calculator_(calculator), impl_(impl) {}

    Scalar evaluate(const std::string& expression) override {
        return calculator_->evaluate(expression);
    }

    std::string evaluate_for_display(const std::string& expression, bool exact_mode) override {
        return calculator_->evaluate_for_display(expression, exact_mode);
    }

    std::string execute_script(const std::string& source, bool exact_mode) override {
        return calculator_->execute_script(source, exact_mode);
    }

    bool is_matrix_argument(const std::string& arg) override;

    matrix::Matrix parse_matrix_argument(const std::string& arg, const std::string& command) override;

    void resolve_symbolic(const std::string& arg, bool req, std::string* var, SymbolicExpression* expr) override;

    std::string render_plot(const std::vector<std::string>& args, bool gnuplot) override;

    Scalar parse_decimal(const std::string& expr) override;

    Scalar normalize_result(Scalar value) override {
        return Calculator::normalize_result(value);
    }

    std::function<Scalar(const std::vector<std::pair<std::string, Scalar>>&)> build_scoped_evaluator(const std::string& expression) override;

    std::function<Scalar(const std::vector<std::pair<std::string, StoredValue>>&)> build_scoped_scalar_evaluator(const std::string& expression) override;

    std::function<matrix::Matrix(const std::vector<std::pair<std::string, StoredValue>>&)> build_scoped_matrix_evaluator(const std::string& expression) override;

    StoredValue evaluate_expression_value(const std::string& arg, bool exact) override;

    FunctionAnalysis build_analysis(const std::string& expression) override;

    std::vector<std::string> parse_symbolic_vars(const std::vector<std::string>& arguments, std::size_t start_index, const std::vector<std::string>& defaults) override;

    std::string expand_inline(const std::string& expression) override;

private:
    Calculator* calculator_;
    Calculator::Impl* impl_;
};

#endif // EVALUATION_ENGINE_IMPL_H