// ============================================================================
// matrix_eval_stats.cpp - 矩阵/向量统计函数求值实现
// ============================================================================

#include "matrix_eval_stats.h"
#include "matrix_internal.h"
#include "mymath.h"
#include <stdexcept>

namespace matrix {
namespace internal {

namespace {

std::vector<Scalar> extract_scalar_vector(
    const std::vector<std::string>& arguments,
    const std::function<bool(const std::string&, Value*)>& eval_value,
    const ScalarEvaluator& scalar_eval,
    const std::string& name) {
    if (arguments.empty()) {
        throw std::runtime_error(name + " expects at least one argument");
    }
    std::vector<Scalar> values;
    if (arguments.size() == 1) {
        Value val;
        if (eval_value(arguments[0], &val) && val.is_matrix) {
            values = as_vector_values(val.matrix, name);
        } else {
            values.push_back(scalar_eval(arguments[0]));
        }
    } else {
        values.reserve(arguments.size());
        for (const std::string& argument : arguments) {
            values.push_back(scalar_eval(argument));
        }
    }
    return values;
}

} // namespace

bool try_evaluate_stats_function(
    const std::string& name,
    const std::vector<std::string>& arguments,
    const std::function<bool(const std::string&, Value*)>& eval_value,
    const ScalarEvaluator& scalar_eval,
    const std::function<Matrix(const std::string&, const std::string&)>& require_matrix,
    Value* result) {
    if (!result) return false;

    if (name == "mean") {
        auto values = extract_scalar_vector(arguments, eval_value, scalar_eval, "mean");
        *result = Value::from_scalar(mean_values(values));
        return true;
    }

    if (name == "median") {
        auto values = extract_scalar_vector(arguments, eval_value, scalar_eval, "median");
        *result = Value::from_scalar(median_values(values));
        return true;
    }

    if (name == "mode") {
        auto values = extract_scalar_vector(arguments, eval_value, scalar_eval, "mode");
        *result = Value::from_scalar(mode_values(values));
        return true;
    }

    if (name == "var") {
        auto values = extract_scalar_vector(arguments, eval_value, scalar_eval, "var");
        *result = Value::from_scalar(variance_values(values));
        return true;
    }

    if (name == "std") {
        auto values = extract_scalar_vector(arguments, eval_value, scalar_eval, "std");
        *result = Value::from_scalar(mymath::sqrt(variance_values(values)));
        return true;
    }

    if (name == "skewness" || name == "skew" || name == "kurtosis") {
        auto values = extract_scalar_vector(arguments, eval_value, scalar_eval, name);
        const Scalar mean = mean_values(values);
        Scalar second_moment = Scalar(0.0L);
        Scalar higher_moment = Scalar(0.0L);
        for (Scalar value : values) {
            const Scalar delta = value - mean;
            const Scalar delta2 = delta * delta;
            second_moment += delta2;
            higher_moment += (name == "kurtosis") ? delta2 * delta2 : delta2 * delta;
        }
        second_moment /= Scalar(static_cast<long long>(values.size()));
        if (mymath::is_near_zero(second_moment, Scalar(1e-12L))) {
            throw std::runtime_error(name + " is undefined for zero variance data");
        }
        higher_moment /= Scalar(static_cast<long long>(values.size()));
        if (name == "kurtosis") {
            *result = Value::from_scalar(
                higher_moment / (second_moment * second_moment) - Scalar(3.0L));
        } else {
            *result = Value::from_scalar(
                higher_moment / mymath::pow(second_moment, Scalar(1.5L)));
        }
        return true;
    }

    if (name == "percentile") {
        if (arguments.size() < 2) {
            throw std::runtime_error("percentile expects vector,p or p,value...");
        }
        if (arguments.size() == 2) {
            Value val;
            if (eval_value(arguments[0], &val) && val.is_matrix) {
                *result = Value::from_scalar(percentile_values(
                    as_vector_values(val.matrix, "percentile"),
                    scalar_eval(arguments[1])));
                return true;
            }
        }
        std::vector<Scalar> values;
        values.reserve(arguments.size() - 1);
        const Scalar p = scalar_eval(arguments[0]);
        for (std::size_t i = 1; i < arguments.size(); ++i) {
            values.push_back(scalar_eval(arguments[i]));
        }
        *result = Value::from_scalar(percentile_values(values, p));
        return true;
    }

    if (name == "quartile") {
        if (arguments.size() < 2) {
            throw std::runtime_error("quartile expects vector,q or q,value...");
        }
        if (arguments.size() == 2) {
            Value val;
            if (eval_value(arguments[0], &val) && val.is_matrix) {
                *result = Value::from_scalar(quartile_values(
                    as_vector_values(val.matrix, "quartile"),
                    scalar_eval(arguments[1])));
                return true;
            }
        }
        std::vector<Scalar> values;
        values.reserve(arguments.size() - 1);
        const Scalar q = scalar_eval(arguments[0]);
        for (std::size_t i = 1; i < arguments.size(); ++i) {
            values.push_back(scalar_eval(arguments[i]));
        }
        *result = Value::from_scalar(quartile_values(values, q));
        return true;
    }

    if (name == "cov") {
        if (arguments.size() != 2) {
            throw std::runtime_error("cov expects exactly two vector arguments");
        }
        *result = Value::from_scalar(covariance_values(
            as_vector_values(require_matrix(arguments[0], "cov"), "cov"),
            as_vector_values(require_matrix(arguments[1], "cov"), "cov")));
        return true;
    }

    if (name == "corr") {
        if (arguments.size() != 2) {
            throw std::runtime_error("corr expects exactly two vector arguments");
        }
        *result = Value::from_scalar(correlation_values(
            as_vector_values(require_matrix(arguments[0], "corr"), "corr"),
            as_vector_values(require_matrix(arguments[1], "corr"), "corr")));
        return true;
    }

    return false;
}

} // namespace internal
} // namespace matrix
