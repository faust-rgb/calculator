# 高精度计算集成计划 (High-Precision Integration Roadmap)

## 1. 目标
将计算器的核心能力（矩阵运算、数值分析、超越函数）从单一的 `double` 精度扩展为支持 `PreciseDecimal` 的多精度系统。在保证 50 位及以上精度的同时，通过算法优化最大程度减少软件模拟带来的性能损耗。

## 2. 现状分析
*   **优势**：已实现基于 $10^9$ 进制的大基数 `PreciseDecimal` 存储，基础算术性能优异。
*   **劣势**：核心模块（Matrix, Analysis）硬编码为 `double`；超越函数在高精度模式下仍退化为 `double` 计算；缺乏运算符重载，无法进行泛型编程。

## 3. 分阶段实施方案

### 第一阶段：增强 `PreciseDecimal` 基础设施 (当前阶段)
*   **目标**：使其在语法上与 `double` 兼容，打通泛型编程的基础。
*   **任务**：
    1.  实现完整的运算符重载 (`+`, `-`, `*`, `/`, `+=`, `*=`, `==`, `<`, `>` 等)。
    2.  支持与 `int`, `long long`, `double` 的混合运算（避免不必要的对象构造）。
    3.  实现高精度的 `abs`, `sqrt`（Newton 迭代），`pow`。
    4.  引入 `PrecisionContext` 管理全局目标精度（默认 40 位）。

### 第二阶段：矩阵模块模板化 (Matrix Templatization)
*   **目标**：使线性代数运算能够直接运行在 `PreciseDecimal` 之上。
*   **任务**：
    1.  将 `Matrix` 重构为 `template <typename T> struct Matrix`。
    2.  泛型化基础运算（加减乘、转置、行列式）。
    3.  泛型化核心分解（LU 分解、QR 分解）。
    4.  优化：针对 `PreciseDecimal` 块操作的特殊内存对齐。

### 第三阶段：数值分析泛型化 (Generic Analysis)
*   **目标**：提升微积分和方程求根的精度。
*   **任务**：
    1.  重构 `Integrator` 和 `Solver` 接口，接受 `std::function<T(T)>`。
    2.  实现高精度自适应步长策略。
    3.  在高精度模式下启用 Richardson 外推法以加速收敛。

### 第四阶段：超越函数与全局集成
*   **目标**：提供全链条的高精度体验。
*   **任务**：
    1.  为 `PreciseDecimal` 实现原生的高精度 `exp`, `ln`, `sin`, `cos`（基于 AGM 或二进制拆分级数）。
    2.  在 `UnifiedParser` 中引入 `precision` 指令。
    3.  实现 `StoredValue` 的全自动精度升级机制。

## 4. 性能优化策略
1.  **延迟规范化**：在连续运算（如向量点积）中，仅在最后一步执行 `normalize()`。
2.  **小整数优化**：对于能放入 `uint64_t` 的 `PreciseDecimal`，采用快速路径计算。
3.  **零拷贝机制**：引入移动语义和 `std::string_view` 解析。

## 5. 验收标准
1.  能够解出 $20 \times 20$ 的 Hilbert 矩阵而不产生显著精度偏差。
2.  `sin(pi)` 在高精度模式下返回的结果至少前 40 位为 0。
3.  全量单元测试通过。
