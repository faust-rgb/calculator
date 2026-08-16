/**
 * @file vector_calculus.cpp
 * @brief 矢量场与微积分算子实现
 */

#include "vector_calculus.h"
#include <stdexcept>

namespace analysis::pde {

std::vector<SymbolicExpression> gradient(const SymbolicExpression& f,
                                         const std::vector<std::string>& variables) {
    if (variables.empty()) {
        throw std::runtime_error("gradient requires at least one variable");
    }
    std::vector<SymbolicExpression> grad;
    grad.reserve(variables.size());
    for (const auto& var : variables) {
        grad.push_back(f.derivative(var).simplify());
    }
    return grad;
}

SymbolicExpression divergence(const std::vector<SymbolicExpression>& F,
                             const std::vector<std::string>& variables) {
    if (F.size() != variables.size()) {
        throw std::runtime_error("divergence: field dimension (" + std::to_string(F.size()) +
                                 ") must match number of variables (" + std::to_string(variables.size()) + ")");
    }
    if (F.empty()) {
        return SymbolicExpression::number(0.0L);
    }

    SymbolicExpression div = F[0].derivative(variables[0]);
    for (std::size_t i = 1; i < F.size(); ++i) {
        div = (div + F[i].derivative(variables[i])).simplify();
    }
    return div.simplify();
}

std::vector<SymbolicExpression> curl(const std::vector<SymbolicExpression>& F,
                                     const std::vector<std::string>& variables) {
    if (variables.size() != 2 && variables.size() != 3) {
        throw std::runtime_error("curl is defined for 2D or 3D vector fields");
    }
    if (F.size() != variables.size()) {
        throw std::runtime_error("curl: field dimension must match number of variables");
    }

    if (variables.size() == 2) {
        // 2D curl: dF2/dx - dF1/dy
        SymbolicExpression c2d = (F[1].derivative(variables[0]) - F[0].derivative(variables[1])).simplify();
        return {c2d};
    }

    // 3D curl: [dF3/dy - dF2/dz, dF1/dz - dF3/dx, dF2/dx - dF1/dy]
    const std::string& x = variables[0];
    const std::string& y = variables[1];
    const std::string& z = variables[2];

    SymbolicExpression cx = (F[2].derivative(y) - F[1].derivative(z)).simplify();
    SymbolicExpression cy = (F[0].derivative(z) - F[2].derivative(x)).simplify();
    SymbolicExpression cz = (F[1].derivative(x) - F[0].derivative(y)).simplify();

    return {cx, cy, cz};
}

SymbolicExpression laplacian(const SymbolicExpression& f,
                            const std::vector<std::string>& variables) {
    if (variables.empty()) {
        throw std::runtime_error("laplacian requires at least one variable");
    }
    SymbolicExpression lap = f.derivative(variables[0]).derivative(variables[0]);
    for (std::size_t i = 1; i < variables.size(); ++i) {
        lap = (lap + f.derivative(variables[i]).derivative(variables[i])).simplify();
    }
    return lap.simplify();
}

} // namespace analysis::pde
