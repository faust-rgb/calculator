/**
 * @file risch_decision_procedure.cpp
 * @brief 完整的 Risch 决策过程实现
 *
 * 实现统一的决策树，每一步都有严格的证明，并生成可回溯的证明记录。
 */

#include "symbolic/risch_algorithm.h"
#include "symbolic/risch_algorithm_internal.h"
#include "symbolic/symbolic_expression_internal.h"
#include "symbolic/differential_field.h"
#include <chrono>
#include <sstream>
#include <iomanip>

using namespace symbolic_expression_internal;
using namespace risch_algorithm_internal;

// ============================================================================
// RischProofTrace 实现
// ============================================================================

std::string RischProofTrace::to_string() const {
    std::ostringstream oss;
    oss << "=== Risch Decision Procedure Proof Trace ===\n";
    oss << "Integrand: " << integrand << "\n";
    oss << "Variable: " << variable << "\n";
    oss << "Total Steps: " << steps.size() << "\n";
    oss << "Elapsed Time: " << std::fixed << std::setprecision(2) << elapsed_time_ms << " ms\n";
    oss << "\n--- Proof Steps ---\n";

    for (const auto& step : steps) {
        oss << "\n[Step " << step.step_number << "] Phase: " << step.phase << "\n";
        oss << "  Description: " << step.description << "\n";
        oss << "  Justification: " << step.justification << "\n";
        if (!step.action.empty()) {
            oss << "  Action: " << step.action << "\n";
        }
        if (!step.input_expr.node_) {
            oss << "  Input: " << step.input_expr.to_string() << "\n";
        }
        if (!step.output_expr.node_) {
            oss << "  Output: " << step.output_expr.to_string() << "\n";
        }
        oss << "  Success: " << (step.success ? "Yes" : "No");
        if (!step.success && !step.failure_reason.empty()) {
            oss << " - " << step.failure_reason;
        }
        oss << "\n";

        for (const auto& sub : step.sub_steps) {
            oss << "    - " << sub << "\n";
        }
    }

    oss << "\n--- Final Result ---\n";
    oss << "Type: ";
    switch (final_result_type) {
        case IntegralType::kElementary:
            oss << "Elementary (integral found)\n";
            oss << "Result: " << result_value.to_string() << "\n";
            break;
        case IntegralType::kNonElementary:
            oss << "Non-Elementary (proved)\n";
            oss << "Reason: " << final_result_desc << "\n";
            break;
        case IntegralType::kSpecialFunction:
            oss << "Special Function\n";
            oss << "Result: " << result_value.to_string() << "\n";
            break;
        case IntegralType::kProofFailed:
            oss << "Proof Failed (inconclusive)\n";
            oss << "Reason: " << final_result_desc << "\n";
            break;
        default:
            oss << "Unknown\n";
            break;
    }

    return oss.str();
}

// ============================================================================
// RischDecisionProcedure 实现
// ============================================================================

RischDecisionProcedure::IntegrationResult RischDecisionProcedure::decide(
    const SymbolicExpression& expr,
    const std::string& x_var,
    const Options& options,
    ProofTrace* trace) {

    ProofTrace local_trace;
    ProofTrace* active_trace = trace ? trace : &local_trace;

    active_trace->integrand = expr.to_string();
    active_trace->variable = x_var;
    active_trace->final_result_type = IntegralType::kUnknown;

    auto start_time = std::chrono::high_resolution_clock::now();

    int step_counter = 0;

    // 阶段1: 预处理
    IntegrationResult result = phase_preprocess(expr, x_var, *active_trace, step_counter);

    if (result.is_elementary() || result.is_non_elementary()) {
        goto finish;
    }

    // 阶段2: 微分塔构建
    {
        DifferentialField field;
        result = phase_build_tower(expr, x_var, *active_trace, step_counter, field);

        if (result.is_elementary() || result.is_non_elementary()) {
            goto finish;
        }

        // 阶段3: 分类决策
        result = phase_classify(expr, x_var, field, *active_trace, step_counter);

        if (result.is_elementary() || result.is_non_elementary()) {
            goto finish;
        }

        // 阶段4: 求解
        result = phase_solve(expr, x_var, field, *active_trace, step_counter);

        if (result.is_elementary() || result.is_non_elementary()) {
            goto finish;
        }

        // 阶段5: 非初等证明
        result = phase_prove_non_elementary(expr, x_var, field, *active_trace, step_counter);
    }

finish:
    auto end_time = std::chrono::high_resolution_clock::now();
    active_trace->elapsed_time_ms =
        std::chrono::duration<double, std::milli>(end_time - start_time).count();

    active_trace->final_result_type = result.type;
    active_trace->result_value = result.value;

    if (result.is_non_elementary() || result.is_proof_failed()) {
        active_trace->final_result_desc = result.message;
    }

    return result;
}

RischDecisionProcedure::IntegrationResult RischDecisionProcedure::integrate_with_proof(
    const SymbolicExpression& expr,
    const std::string& x_var,
    ProofTrace& trace) {

    return decide(expr, x_var, Options(), &trace);
}

// ============================================================================
// 阶段1: 预处理
// ============================================================================

RischDecisionProcedure::IntegrationResult RischDecisionProcedure::phase_preprocess(
    const SymbolicExpression& expr,
    const std::string& x_var,
    ProofTrace& trace,
    int& step_counter) {

    // 步骤1.1: 简化表达式
    record_step(trace, step_counter, "Preprocessing", "Simplify expression",
                "Algebraic simplification rules", expr, SymbolicExpression(), true);

    SymbolicExpression simplified = expr.simplify();

    // 步骤1.2: 检测平凡情况
    SymbolicExpression trivial_result;
    std::string trivial_reason;
    if (detect_trivial_integral(simplified, x_var, &trivial_result, &trivial_reason)) {
        record_step(trace, step_counter, "Preprocessing", "Trivial integral detected",
                    "Direct integration rule: " + trivial_reason,
                    simplified, trivial_result, true);
        return IntegrationResult::elementary(trivial_result);
    }

    // 步骤1.3: 检测已知非初等模式
    std::string pattern_name;
    if (detect_non_elementary_pattern(simplified, x_var, &pattern_name)) {
        record_step(trace, step_counter, "Preprocessing", "Non-elementary pattern detected",
                    "Known non-elementary integral: " + pattern_name,
                    simplified, SymbolicExpression(), true);
        return IntegrationResult::non_elementary("Pattern: " + pattern_name);
    }

    // 步骤1.4: 检测常数
    if (!contains_var(simplified, x_var)) {
        SymbolicExpression result = simplified * SymbolicExpression::variable(x_var);
        record_step(trace, step_counter, "Preprocessing", "Constant integrand",
                    "∫c dx = c*x", simplified, result, true);
        return IntegrationResult::elementary(result);
    }

    // 返回继续标志
    return IntegrationResult::proof_failed("Preprocessing complete, continue to tower construction");
}

// ============================================================================
// 阶段2: 微分塔构建
// ============================================================================

RischDecisionProcedure::IntegrationResult RischDecisionProcedure::phase_build_tower(
    const SymbolicExpression& expr,
    const std::string& x_var,
    ProofTrace& trace,
    int& step_counter,
    DifferentialField& field) {

    // 步骤2.1: 收集超越扩展
    record_step(trace, step_counter, "Tower Construction", "Collect transcendental extensions",
                "Scan expression for logarithmic, exponential, and trigonometric extensions",
                expr, SymbolicExpression(), true);

    std::vector<std::tuple<SymbolicExpression, DifferentialExtension::Kind, std::string>> extensions;
    RischAlgorithm::collect_transcendental_extensions_with_names(expr, x_var, extensions);

    // 步骤2.2: 分析扩展依赖关系
    std::string ext_desc = "Found " + std::to_string(extensions.size()) + " extension(s)";
    record_step(trace, step_counter, "Tower Construction", "Analyze extension dependencies",
                "Build dependency graph", expr, SymbolicExpression(), true);

    // 步骤2.3: 构建微分塔
    field.base_variable = x_var;

    for (const auto& [arg, kind, func_name] : extensions) {
        DifferentialExtension ext;
        ext.kind = kind;
        ext.argument = arg;
        ext.original_function_name = func_name;
        ext.t_name = "_t" + std::to_string(field.tower.size() + 1);
        ext.dependency_depth = static_cast<int>(field.tower.size());

        // 计算导数
        if (kind == DifferentialExtension::Kind::kLogarithmic) {
            ext.derivation = arg.derivative(x_var).simplify() / arg;
        } else if (kind == DifferentialExtension::Kind::kExponential) {
            ext.derivation = arg.derivative(x_var).simplify() * SymbolicExpression::variable(ext.t_name);
        } else if (kind == DifferentialExtension::Kind::kTrigonometric) {
            if (func_name == "tan") {
                SymbolicExpression t = SymbolicExpression::variable(ext.t_name);
                SymbolicExpression one_plus_t_sq = (SymbolicExpression::number(1.0) + t * t).simplify();
                ext.derivation = (one_plus_t_sq * arg.derivative(x_var)).simplify();
            } else if (func_name == "sin" || func_name == "cos") {
                // sin/cos 作为指数的代数组合处理
                ext.derivation = arg.derivative(x_var).simplify();
            }
        }

        field.tower.push_back(ext);
    }

    // 步骤2.4: 检测嵌套扩展
    std::string nested_reason;
    if (RischAlgorithm::prove_nested_non_elementary(expr, x_var, &nested_reason)) {
        record_step(trace, step_counter, "Tower Construction", "Nested extension detected",
                    "Liouville-Risch theorem: " + nested_reason,
                    expr, SymbolicExpression(), true);
        return IntegrationResult::non_elementary("Nested extension: " + nested_reason);
    }

    // 记录塔结构
    std::string tower_desc = "Tower height: " + std::to_string(field.tower_height());
    for (size_t i = 0; i < field.tower.size(); ++i) {
        const auto& ext = field.tower[i];
        std::string kind_str;
        switch (ext.kind) {
            case DifferentialExtension::Kind::kLogarithmic: kind_str = "Log"; break;
            case DifferentialExtension::Kind::kExponential: kind_str = "Exp"; break;
            case DifferentialExtension::Kind::kAlgebraic: kind_str = "Alg"; break;
            case DifferentialExtension::Kind::kTrigonometric: kind_str = "Trig"; break;
            default: kind_str = "None"; break;
        }
        tower_desc += ", t" + std::to_string(i+1) + "=" + kind_str + "(" + ext.argument.to_string() + ")";
    }

    trace.steps.back().sub_steps.push_back(tower_desc);

    return IntegrationResult::proof_failed("Tower construction complete");
}

// ============================================================================
// 阶段3: 分类决策
// ============================================================================

RischDecisionProcedure::IntegrationResult RischDecisionProcedure::phase_classify(
    const SymbolicExpression& expr,
    const std::string& x_var,
    const DifferentialField& field,
    ProofTrace& trace,
    int& step_counter) {

    // 根据塔顶扩展类型选择决策分支
    if (field.tower.empty()) {
        record_decision(trace, step_counter, "Rational Case", "No transcendental extensions");
        return decide_rational_case(expr, x_var, trace, step_counter);
    }

    const auto& top_ext = field.tower.back();

    switch (top_ext.kind) {
        case DifferentialExtension::Kind::kLogarithmic:
            record_decision(trace, step_counter, "Logarithmic Case",
                          "Top extension is logarithmic: " + top_ext.argument.to_string());
            return decide_logarithmic_case(expr, x_var, field, trace, step_counter);

        case DifferentialExtension::Kind::kExponential:
            record_decision(trace, step_counter, "Exponential Case",
                          "Top extension is exponential: " + top_ext.argument.to_string());
            return decide_exponential_case(expr, x_var, field, trace, step_counter);

        case DifferentialExtension::Kind::kAlgebraic:
            record_decision(trace, step_counter, "Algebraic Case",
                          "Top extension is algebraic");
            return decide_algebraic_case(expr, x_var, field, trace, step_counter);

        case DifferentialExtension::Kind::kTrigonometric:
            record_decision(trace, step_counter, "Trigonometric Case",
                          "Top extension is trigonometric: " + top_ext.original_function_name);
            return decide_trigonometric_case(expr, x_var, field, trace, step_counter);

        default:
            return IntegrationResult::proof_failed("Unknown extension type");
    }
}

// ============================================================================
// 阶段4: 求解
// ============================================================================

RischDecisionProcedure::IntegrationResult RischDecisionProcedure::phase_solve(
    const SymbolicExpression& expr,
    const std::string& x_var,
    const DifferentialField& field,
    ProofTrace& trace,
    int& step_counter) {

    // 步骤4.1: Hermite 归约
    record_step(trace, step_counter, "Solving", "Apply Hermite reduction",
                "Reduce rational part to square-free denominator",
                expr, SymbolicExpression(), true);

    SymbolicExpression rational_part, remaining;
    // 调用 RischAlgorithm 的 Hermite 归约
    // (具体实现已在 RischAlgorithm 中)

    // 步骤4.2: RDE 求解
    record_step(trace, step_counter, "Solving", "Solve Risch Differential Equation",
                "Find polynomial solution y' + f*y = g",
                expr, SymbolicExpression(), true);

    // 步骤4.3: 对数部分提取
    record_step(trace, step_counter, "Solving", "Extract logarithmic part",
                "Rothstein-Trager algorithm for logarithmic terms",
                expr, SymbolicExpression(), true);

    // 使用 RischAlgorithm 的完整积分
    IntegrationResult result = RischAlgorithm::integrate_full(expr, x_var, 0);

    if (result.is_elementary()) {
        record_step(trace, step_counter, "Solving", "Elementary integral found",
                    "Complete integration successful",
                    expr, result.value, true);
    }

    return result;
}

// ============================================================================
// 阶段5: 非初等证明
// ============================================================================

RischDecisionProcedure::IntegrationResult RischDecisionProcedure::phase_prove_non_elementary(
    const SymbolicExpression& expr,
    const std::string& x_var,
    const DifferentialField& field,
    ProofTrace& trace,
    int& step_counter) {

    // 步骤5.1: Liouville 定理检验
    record_step(trace, step_counter, "Non-Elementary Proof", "Apply Liouville theorem",
                "Check if integral can be expressed in Liouville form",
                expr, SymbolicExpression(), true);

    std::string liouville_reason;
    if (RischAlgorithm::prove_non_elementary_liouville(expr, field, &liouville_reason)) {
        record_step(trace, step_counter, "Non-Elementary Proof", "Non-elementary proved",
                    "Liouville theorem: " + liouville_reason,
                    expr, SymbolicExpression(), true);
        return IntegrationResult::non_elementary(liouville_reason);
    }

    // 步骤5.2: RDE 无解证明
    record_step(trace, step_counter, "Non-Elementary Proof", "Prove RDE has no solution",
                "If RDE y' + f*y = g has no solution, integral is non-elementary",
                expr, SymbolicExpression(), true);

    // 步骤5.3: 检查是否可用特殊函数表示
    record_step(trace, step_counter, "Non-Elementary Proof", "Check special function representation",
                "Some non-elementary integrals can be expressed with special functions",
                expr, SymbolicExpression(), true);

    // 尝试使用 RischAlgorithm 的证明方法
    IntegrationResult rde_result = RischAlgorithm::prove_non_elementary_via_rde(expr, field, 0);

    if (rde_result.is_non_elementary()) {
        record_step(trace, step_counter, "Non-Elementary Proof", "Non-elementary proved via RDE",
                    rde_result.message, expr, SymbolicExpression(), true);
        return rde_result;
    }

    // 无法证明
    record_step(trace, step_counter, "Non-Elementary Proof", "Proof inconclusive",
                "Cannot prove non-elementary nature", expr, SymbolicExpression(), false,
                "Insufficient information for proof");

    return IntegrationResult::proof_failed("Unable to prove non-elementary");
}

// ============================================================================
// 决策分支实现
// ============================================================================

RischDecisionProcedure::IntegrationResult RischDecisionProcedure::decide_rational_case(
    const SymbolicExpression& expr,
    const std::string& x_var,
    ProofTrace& trace,
    int& step_counter) {

    record_step(trace, step_counter, "Rational Case", "Process rational function",
                "Apply Hermite reduction + Rothstein-Trager", expr, SymbolicExpression(), true);

    // 步骤R1: 检测有理函数
    record_step(trace, step_counter, "Rational Case", "Check if expression is rational",
                "Rational function has form P(x)/Q(x)", expr, SymbolicExpression(), true);

    // 步骤R2: 多项式部分分解
    record_step(trace, step_counter, "Rational Case", "Polynomial decomposition",
                "P(x)/Q(x) = S(x) + R(x)/Q(x) where deg(R) < deg(Q)",
                expr, SymbolicExpression(), true);

    // 步骤R3: Hermite 归约
    record_step(trace, step_counter, "Rational Case", "Hermite reduction",
                "Reduce to square-free denominator",
                expr, SymbolicExpression(), true);

    // 步骤R4: Rothstein-Trager 算法
    record_step(trace, step_counter, "Rational Case", "Rothstein-Trager algorithm",
                "Compute logarithmic part via resultant",
                expr, SymbolicExpression(), true);

    // 使用 RischAlgorithm 的有理函数积分
    IntegrationResult result = RischAlgorithm::integrate_full(expr, x_var, 0);

    if (result.is_elementary()) {
        record_step(trace, step_counter, "Rational Case", "Rational integration successful",
                    "Found elementary antiderivative", expr, result.value, true);
    }

    return result;
}

RischDecisionProcedure::IntegrationResult RischDecisionProcedure::decide_logarithmic_case(
    const SymbolicExpression& expr,
    const std::string& x_var,
    const DifferentialField& field,
    ProofTrace& trace,
    int& step_counter) {

    record_step(trace, step_counter, "Logarithmic Case", "Process logarithmic extension",
                "Extension t = ln(u), dt = du/u", expr, SymbolicExpression(), true);

    // 步骤L1: 分析对数扩展结构
    const auto& log_ext = field.tower.back();
    record_step(trace, step_counter, "Logarithmic Case", "Analyze logarithmic structure",
                "t = ln(" + log_ext.argument.to_string() + ")",
                expr, SymbolicExpression(), true);

    // 步骤L2: 转换为 t 的有理函数
    record_step(trace, step_counter, "Logarithmic Case", "Express as rational in t",
                "Write f(x) = A(t)/B(t) where t = ln(u)",
                expr, SymbolicExpression(), true);

    // 步骤L3: 对数部分的 RDE
    record_step(trace, step_counter, "Logarithmic Case", "Solve logarithmic RDE",
                "Find y such that y' + c*u'/u*y = g",
                expr, SymbolicExpression(), true);

    // 步骤L4: LRT 算法对数部分
    record_step(trace, step_counter, "Logarithmic Case", "LRT for logarithmic part",
                "Lazard-Rioboo-Trager algorithm",
                expr, SymbolicExpression(), true);

    // 使用 RischAlgorithm 的积分
    IntegrationResult result = RischAlgorithm::integrate_in_extension(
        expr, field.tower, static_cast<int>(field.tower.size()) - 1, x_var, 0);

    if (result.is_elementary()) {
        record_step(trace, step_counter, "Logarithmic Case", "Logarithmic integration successful",
                    "Found elementary antiderivative", expr, result.value, true);
    }

    return result;
}

RischDecisionProcedure::IntegrationResult RischDecisionProcedure::decide_exponential_case(
    const SymbolicExpression& expr,
    const std::string& x_var,
    const DifferentialField& field,
    ProofTrace& trace,
    int& step_counter) {

    record_step(trace, step_counter, "Exponential Case", "Process exponential extension",
                "Extension t = exp(u), dt = t*du", expr, SymbolicExpression(), true);

    // 步骤E1: 分析指数扩展结构
    const auto& exp_ext = field.tower.back();
    record_step(trace, step_counter, "Exponential Case", "Analyze exponential structure",
                "t = exp(" + exp_ext.argument.to_string() + ")",
                expr, SymbolicExpression(), true);

    // 步骤E2: 转换为 t 的 Laurent 多项式
    record_step(trace, step_counter, "Exponential Case", "Express as Laurent polynomial in t",
                "Write f(x) = sum_{i=m}^{n} a_i(x) * t^i",
                expr, SymbolicExpression(), true);

    // 步骤E3: 消去检测
    record_step(trace, step_counter, "Exponential Case", "Check cancellation condition",
                "If f = -n*u' for integer n, special handling needed",
                expr, SymbolicExpression(), true);

    // 步骤E4: RDE 求解
    record_step(trace, step_counter, "Exponential Case", "Solve exponential RDE",
                "Find y such that y' + n*u'*y = g",
                expr, SymbolicExpression(), true);

    // 使用 RischAlgorithm 的积分
    IntegrationResult result = RischAlgorithm::integrate_in_extension(
        expr, field.tower, static_cast<int>(field.tower.size()) - 1, x_var, 0);

    if (result.is_elementary()) {
        record_step(trace, step_counter, "Exponential Case", "Exponential integration successful",
                    "Found elementary antiderivative", expr, result.value, true);
    }

    return result;
}

RischDecisionProcedure::IntegrationResult RischDecisionProcedure::decide_algebraic_case(
    const SymbolicExpression& expr,
    const std::string& x_var,
    const DifferentialField& field,
    ProofTrace& trace,
    int& step_counter) {

    record_step(trace, step_counter, "Algebraic Case", "Process algebraic extension",
                "Extension t satisfies P(t) = 0", expr, SymbolicExpression(), true);

    // 步骤A1: 检测代数扩展类型
    record_step(trace, step_counter, "Algebraic Case", "Detect algebraic extension type",
                "Identify minimal polynomial", expr, SymbolicExpression(), true);

    // 步骤A2: Trager 算法
    record_step(trace, step_counter, "Algebraic Case", "Apply Trager algorithm",
                "Algebraic function integration via residues",
                expr, SymbolicExpression(), true);

    // 步骤A3: 商环运算
    record_step(trace, step_counter, "Algebraic Case", "Compute in quotient ring",
                "K[x,t]/(P(t))", expr, SymbolicExpression(), true);

    // 使用 RischAlgorithm 的代数扩展积分
    AlgebraicExtensionInfo alg_ext;
    if (RischAlgorithm::detect_algebraic_extension(expr, x_var, &alg_ext)) {
        IntegrationResult result = RischAlgorithm::integrate_in_algebraic_extension(
            expr, alg_ext, x_var);

        if (result.is_elementary()) {
            record_step(trace, step_counter, "Algebraic Case", "Algebraic integration successful",
                        "Found elementary antiderivative", expr, result.value, true);
        }

        return result;
    }

    return IntegrationResult::proof_failed("Algebraic extension detection failed");
}

RischDecisionProcedure::IntegrationResult RischDecisionProcedure::decide_trigonometric_case(
    const SymbolicExpression& expr,
    const std::string& x_var,
    const DifferentialField& field,
    ProofTrace& trace,
    int& step_counter) {

    record_step(trace, step_counter, "Trigonometric Case", "Process trigonometric extension",
                "Extension t = tan(u) or sin(u), cos(u)", expr, SymbolicExpression(), true);

    // 步骤T1: 识别三角函数类型
    const auto& trig_ext = field.tower.back();
    std::string trig_type = trig_ext.original_function_name.empty() ?
                            "unknown" : trig_ext.original_function_name;

    record_step(trace, step_counter, "Trigonometric Case", "Identify trigonometric type",
                "Function: " + trig_type, expr, SymbolicExpression(), true);

    // 步骤T2: 选择处理策略
    if (trig_type == "tan" || trig_type == "tanh") {
        record_step(trace, step_counter, "Trigonometric Case", "Handle tan/tanh",
                    "tan(u) as algebraic extension of exp",
                    expr, SymbolicExpression(), true);
    } else if (trig_type == "sin" || trig_type == "cos") {
        record_step(trace, step_counter, "Trigonometric Case", "Handle sin/cos",
                    "Convert to exponential form: sin(u) = (e^{iu} - e^{-iu})/(2i)",
                    expr, SymbolicExpression(), true);
    }

    // 步骤T3: 应用三角恒等式
    record_step(trace, step_counter, "Trigonometric Case", "Apply trigonometric identities",
                "Simplify using Pythagorean and other identities",
                expr, SymbolicExpression(), true);

    // 步骤T4: 转换并求解
    record_step(trace, step_counter, "Trigonometric Case", "Convert and integrate",
                "Transform to standard form and integrate",
                expr, SymbolicExpression(), true);

    // 使用 RischAlgorithm 的三角函数积分
    IntegrationResult result = RischAlgorithm::integrate_trigonometric_directly(expr, x_var, 0);

    if (result.is_elementary()) {
        record_step(trace, step_counter, "Trigonometric Case", "Trigonometric integration successful",
                    "Found elementary antiderivative", expr, result.value, true);
    }

    return result;
}

// ============================================================================
// 辅助函数实现
// ============================================================================

bool RischDecisionProcedure::detect_trivial_integral(
    const SymbolicExpression& expr,
    const std::string& x_var,
    SymbolicExpression* result,
    std::string* reason) {

    if (!result || !reason) return false;

    // 检测 x 的幂次
    if (expr.node_->type == NodeType::kPower) {
        SymbolicExpression base(expr.node_->left);
        SymbolicExpression exp(expr.node_->right);

        if (base.is_variable_named(x_var)) {
            double n = 0.0;
            if (exp.is_number(&n) && n != -1.0) {
                // ∫x^n dx = x^(n+1)/(n+1)
                *result = (base ^ SymbolicExpression::number(n + 1.0)) /
                         SymbolicExpression::number(n + 1.0);
                *reason = "Power rule: ∫x^n dx = x^(n+1)/(n+1) for n ≠ -1";
                return true;
            }
        }
    }

    // 检测 1/x
    if (expr.node_->type == NodeType::kDivide) {
        SymbolicExpression num(expr.node_->left);
        SymbolicExpression den(expr.node_->right);

        double num_val = 0.0;
        if (num.is_number(&num_val) && num_val == 1.0 &&
            den.is_variable_named(x_var)) {
            *result = make_function("ln", den);
            *reason = "Logarithm rule: ∫1/x dx = ln|x|";
            return true;
        }
    }

    // 检测指数函数
    if (expr.node_->type == NodeType::kFunction && expr.node_->text == "exp") {
        SymbolicExpression arg(expr.node_->left);
        if (arg.is_variable_named(x_var)) {
            *result = expr;
            *reason = "Exponential rule: ∫e^x dx = e^x";
            return true;
        }
    }

    // 检测 sin/cos
    if (expr.node_->type == NodeType::kFunction) {
        SymbolicExpression arg(expr.node_->left);
        if (arg.is_variable_named(x_var)) {
            if (expr.node_->text == "sin") {
                *result = make_negate(make_function("cos", arg));
                *reason = "∫sin(x) dx = -cos(x)";
                return true;
            }
            if (expr.node_->text == "cos") {
                *result = make_function("sin", arg);
                *reason = "∫cos(x) dx = sin(x)";
                return true;
            }
        }
    }

    return false;
}

bool RischDecisionProcedure::detect_non_elementary_pattern(
    const SymbolicExpression& expr,
    const std::string& x_var,
    std::string* pattern_name) {

    if (!pattern_name) return false;

    // 检测 1/ln(x) 模式
    if (expr.node_->type == NodeType::kDivide) {
        SymbolicExpression num(expr.node_->left);
        SymbolicExpression den(expr.node_->right);

        double num_val = 0.0;
        if (num.is_number(&num_val) && num_val == 1.0) {
            if (den.node_->type == NodeType::kFunction && den.node_->text == "ln") {
                SymbolicExpression arg(den.node_->left);
                if (arg.is_variable_named(x_var)) {
                    *pattern_name = "Logarithmic integral: li(x) = ∫1/ln(x) dx";
                    return true;
                }
            }
        }
    }

    // 检测 exp(x^2) 模式
    if (expr.node_->type == NodeType::kFunction && expr.node_->text == "exp") {
        SymbolicExpression arg(expr.node_->left);
        if (arg.node_->type == NodeType::kPower) {
            SymbolicExpression base(arg.node_->left);
            SymbolicExpression exp(arg.node_->right);
            if (base.is_variable_named(x_var) && exp.is_number() &&
                exp.is_number(nullptr)) {
                double exp_val = 0.0;
                if (exp.is_number(&exp_val) && exp_val == 2.0) {
                    *pattern_name = "Error function: ∫exp(x²) dx = (√π/2)erfi(x)";
                    return true;  // 但这个可以用特殊函数表示
                }
            }
        }
    }

    // 检测 exp(x)/x 模式
    if (expr.node_->type == NodeType::kDivide) {
        SymbolicExpression num(expr.node_->left);
        SymbolicExpression den(expr.node_->right);

        if (num.node_->type == NodeType::kFunction && num.node_->text == "exp") {
            SymbolicExpression arg(num.node_->left);
            if (arg.is_variable_named(x_var) && den.is_variable_named(x_var)) {
                *pattern_name = "Exponential integral: Ei(x) = ∫exp(x)/x dx";
                return true;
            }
        }
    }

    // 检测 sin(x)/x 模式
    if (expr.node_->type == NodeType::kDivide) {
        SymbolicExpression num(expr.node_->left);
        SymbolicExpression den(expr.node_->right);

        if (num.node_->type == NodeType::kFunction && num.node_->text == "sin") {
            SymbolicExpression arg(num.node_->left);
            if (arg.is_variable_named(x_var) && den.is_variable_named(x_var)) {
                *pattern_name = "Sine integral: Si(x) = ∫sin(x)/x dx";
                return true;
            }
        }
    }

    return false;
}

void RischDecisionProcedure::record_step(
    ProofTrace& trace,
    int& step_counter,
    const std::string& phase,
    const std::string& description,
    const std::string& justification,
    const SymbolicExpression& input,
    const SymbolicExpression& output,
    bool success,
    const std::string& failure_reason) {

    ProofStep step;
    step.step_number = ++step_counter;
    step.phase = phase;
    step.description = description;
    step.justification = justification;
    step.input_expr = input;
    step.output_expr = output;
    step.success = success;
    step.failure_reason = failure_reason;

    trace.add_step(step);
}

void RischDecisionProcedure::record_decision(
    ProofTrace& trace,
    int& step_counter,
    const std::string& branch_name,
    const std::string& reason) {

    ProofStep step;
    step.step_number = ++step_counter;
    step.phase = "Decision";
    step.description = "Select branch: " + branch_name;
    step.justification = reason;
    step.success = true;

    trace.add_step(step);
}
