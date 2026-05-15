// ============================================================================
// evaluation_engine_impl.h - IEvaluationEngine 的具体实现
// ============================================================================

#ifndef EVALUATION_ENGINE_IMPL_H
#define EVALUATION_ENGINE_IMPL_H

#include "core_manager_interfaces.h"
#include "core/api/calculator.h"
#include "analysis/calculus/function_analysis.h"

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

    bool is_matrix_argument(const std::string& arg) override {
        return calculator_->get_core_services().is_matrix_argument(arg);
    }

    matrix::Matrix parse_matrix_argument(const std::string& arg, const std::string& command) override {
        return calculator_->get_core_services().parse_matrix_argument(arg, command);
    }
    
    void resolve_symbolic(const std::string& arg, bool req, std::string* var, SymbolicExpression* expr) override {
        calculator_->get_core_services().symbolic.resolve_symbolic(arg, req, var, expr);
    }

    std::string render_plot(const std::vector<std::string>& args, bool gnuplot) override {
        return calculator_->get_core_services().render_plot(args, gnuplot);
    }

    Scalar parse_decimal(const std::string& expr) override {
        return calculator_->get_core_services().evaluation.parse_decimal(expr);
    }

    Scalar normalize_result(Scalar value) override {
        return calculator_->get_core_services().evaluation.normalize_result(value);
    }

    std::function<Scalar(const std::vector<std::pair<std::string, Scalar>>&)> build_scoped_evaluator(const std::string& expression) override {
        return calculator_->get_core_services().evaluation.build_decimal_evaluator(expression);
    }

    std::function<Scalar(const std::vector<std::pair<std::string, StoredValue>>&)> build_scoped_scalar_evaluator(const std::string& expression) override {
        return calculator_->get_core_services().evaluation.build_scalar_evaluator(expression);
    }

    std::function<matrix::Matrix(const std::vector<std::pair<std::string, StoredValue>>&)> build_scoped_matrix_evaluator(const std::string& expression) override {
        return calculator_->get_core_services().evaluation.build_matrix_evaluator(expression);
    }

    StoredValue evaluate_expression_value(const std::string& arg, bool exact) override {
        return calculator_->get_core_services().evaluation.evaluate_value(arg, exact);
    }

    FunctionAnalysis build_analysis(const std::string& expression) override {
        return calculator_->get_core_services().symbolic.build_analysis(expression);
    }

    std::vector<std::string> parse_symbolic_vars(const std::vector<std::string>& arguments, std::size_t start_index, const std::vector<std::string>& defaults) override {
        return calculator_->get_core_services().parse_symbolic_vars(arguments, start_index, defaults);
    }

    std::string expand_inline(const std::string& expression) override {
        return calculator_->get_core_services().symbolic.expand_inline(expression);
    }

private:
    Calculator* calculator_;
    Calculator::Impl* impl_;
};

#endif // EVALUATION_ENGINE_IMPL_H
