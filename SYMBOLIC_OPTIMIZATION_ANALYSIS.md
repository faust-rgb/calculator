# 符号计算优化分析

## 执行摘要

当前符号计算模块已实现了基础的微积分、简化、变换等功能。通过对代码和文档的分析，识别出 8 个核心优化方向，按优先级分为 3 个阶段。本分析基于已实现功能的边界和现有限制。

## 本次落地结果（2026-04-24）

本轮已优先完成文档中“在现有架构下可直接实现”的部分，并补齐回归测试，重点落在第一阶段的规范化与化简稳定性，以及第二阶段的符号微积分和多变量符号操作：

- 已完成：乘法因子按幂次归并，支持 `x * x * x -> x ^ 3`、`x ^ 3 / x ^ 2 -> x`
- 已完成：公因子提取，支持 `2 * x + 2 * y -> 2 * (x + y)`、`a * b + a * c -> a * (b + c)`
- 已完成：多项式商化简，支持 `(x ^ 2 - 1) / (x - 1) -> x + 1`
- 已完成：单变量多项式规范输出，改善 `x + x ^ 2 + 1` 这类表达式的稳定顺序
- 已完成：安全的 `exp(ln(...))` 简化，只在已知正值表达式上折叠，避免把 `exp(ln(x))` 误化简为 `x`
- 已完成：域感知根式化简，`sqrt(u ^ 2)` 规范为 `abs(u)`，避免在未知符号符号域中错误丢失符号条件
- 已完成：输出规范收敛，优先保留 `pi / 2`、`x ^ 3 / 3` 这类“纯数分母”形式，并修正因式顺序
- 已完成：反三角积分规则，支持 `integral(1 / (1 + x ^ 2)) -> atan(x) + C` 和 `integral(1 / sqrt(1 - x ^ 2)) -> asin(x) + C`
- 已完成：平方根特殊积分，支持 `integral(sqrt(1 - x ^ 2))`
- 已完成：有理函数积分基础层，支持多项式长除法后对线性/二次分母余项积分，例如 `integral(x / (1 + x ^ 2))`、`integral((x ^ 2 + 1) / (x + 1))`、`integral(1 / (x ^ 2 - 1))`
- 已完成：多项式 GCD 约分，支持非整除有理式中的公共多项式因子消去，例如 `(x ^ 3 - x) / (x ^ 2 - 2 * x + 1) -> (x ^ 2 + x) / (x - 1)`
- 已完成：基础部分分式积分，支持互异实线性因子、单一重复线性因子和常见重复不可约二次因子，例如 `integral(1 / (x ^ 3 - x))`、`integral(1 / (x - 1) ^ 2)`、`integral(1 / (x ^ 2 + 1) ^ 2)`
- 已完成：混合重复实线性因子的部分分式积分子集，例如 `integral(1 / ((x - 1) ^ 2 * (x + 1)))`
- 已完成：混合实线性因子 + 一个重复不可约二次因子的部分分式积分，例如 `integral(1 / ((x - 1) * (x ^ 2 + 1) ^ 2))`
- 已完成：三角恒等式积分基础扩展，支持 `sin(x)^2`、`cos(x)^2`、`tan(x)^2`、`sin(x)^3`、`cos(x)^3` 和 `sin(x)*cos(x)`
- 已完成：链式替换积分规则，支持 `g'(x) * F(g(x))` 形式，例如 `integral(2 * x * cos(x ^ 2)) -> sin(x ^ 2) + C`
- 已完成：显式多变量偏导入口，支持 `diff(expr, x)` 与 `diff(expr, x, y)` 这类连续偏导（命令层通过链式调用 `derivative()` 实现）
- 已完成：多变量符号 API 与命令入口，支持 `gradient(expr, vars...)`、`jacobian([exprs], vars...)`、`hessian(expr, vars...)`
- 已完成：连续多变量符号积分入口，支持 `integral(expr, x, y, ...)` 链式积分，并允许被积表达式中出现相对当前积分变量的符号常量
- 已完成：`critical(expr, vars...)` 命令入口，仿射梯度使用线性系统精确求解，非线性梯度使用有界 Newton 搜索，并基于 Hessian 对 1-3 维孤立驻点输出 `local min`、`local max`、`saddle` 或 `degenerate`
- 已完成：替换安全检查，禁止替换 `pi/e/i` 等保留符号或非法变量名；当前表达式树无绑定变量，因此不存在捕获型替换风险
- 已完成：变换框架第一步现代化，抽出 Fourier/Laplace/z 变换共用的线性、减法、取负、常数倍规则入口
- 已完成：`simplify()` 与 `derivative()` 的线程局部结构键缓存，缓存规模有上限以避免交互会话无限增长；Taylor/Pade 共享同一轮导数链缓存，避免重复构造中间导数
- 已完成：负有理数格式化修正，避免 `-1/8 * x ^ 2` 被错误显示为正项

## 性能优化落地结果（2026-04-25）

本轮按实际代码边界优先完成低风险、高收益优化：

- 已完成：表达式节点驻留（interning），复用结构相同的不可变节点，降低符号化简和求导中的重复分配
- 已完成：表达式节点结构键缓存，避免 `node_structural_key()` 对同一子树重复生成字符串
- 已完成：`simplify()` 与 `derivative()` 缓存从“满表清空”改为线程局部 LRU 淘汰，减少缓存抖动
- 已完成：多项式系数提取增加单次调用内的结构键记忆表，复用重复子树结果，缓解嵌套表达式中的重复递归
- 已完成：矩阵乘法改为小块平铺的 `i-k-j` 累加顺序，改善行优先存储下的缓存局部性
- 已完成：`dft/fft` 和 `idft/ifft` 对 2 的幂长度输入使用迭代 Cooley-Tukey FFT，非 2 的幂长度继续使用原 O(n^2) DFT
- 已完成：Jacobi SVD 底层对称特征分解增加尺度感知的 off-diagonal 收敛阈值，避免在已收敛矩阵上继续迭代
- 已完成：数值函数分析复用求值器并使用实例级 LRU 求值缓存，减少导数、极限、积分流程中相同采样点的重复解析和求值
- 已完成：节点 interning 改为 LRU 淘汰，避免容量满时整表清空导致复杂符号任务缓存抖动
- 已完成：Hessian 先缓存一阶梯度并只计算上三角，再按对称性填充下三角

未在本轮直接实现的优化：

- 完整跨工作流 CSE：当前通过节点驻留、结构键缓存和多项式提取记忆化获得共享收益，但没有引入可输出临时变量的全局 CSE pass
- 完整 AST 预解析函数分析：`Calculator` 的表达式解析器目前仍封装在命令求值路径中，本轮先用求值缓存降低重复解析成本
- 多个不同不可约高次因子的专用部分分式快速算法：当前仍复用通用线性系统路径

仍保留为后续工作：

- 更完整的有理函数体系：多个不同不可约二次因子、一般不可约高次因子、符号参数化因子
- 更广的积分规则：更完整的三角恒等变换、一般非线性分部积分策略、完整 Risch 风格判定
- 更强的非线性多变量 critical point 全局完备求解；当前已有 Hessian 分类但搜索仍是有界启发式
- 更完整的变换规则表和模式匹配编译器
- 跨命令/跨工作流的显式 CSE pass、更细粒度 Taylor/Pade 中间导数缓存，以及面向特殊部分分式系统的专用求解器

---

## 优化方向分类

### 第一阶段：核心基础设施优化（高优先级）

#### 1. **规范简化（Canonical Simplification）**

**当前状态：已完成第一批核心能力**

**现状：**
- 当前简化是启发式和代数式的
- 同等数学表达式可能打印为不同形式（如 `1 - x^2` vs `-(x^2) + 1`）
- 缺乏统一的标准形式
- 乘法因子收集不完整（不能跟踪多个相同因子的指数）

**优化空间：**
```
目标优化示例：
- x^3 / x^2 → x               (当前不能可靠地简化)
- x^5 * x^2 → x^7             (需要完整的指数追踪)
- x * x * x → x^3             (需要统一的幂次标准化)
- (x + 1)^2 / (x + 1) → x + 1 (复杂因子取消)
```

**实现策略：**
1. 在乘法/除法分解中扩展 `(base, exponent)` 对的追踪
2. 规范化重复相同因子为幂
3. 在分数取消时减去指数
4. 使用 `polynomial_coefficients()` 检测多项式并转换为规范输出

**预期收益：**
- 改进的符号表达式可读性
- 减少输出噪声
- 更好的重复简化稳定性

**代码位置：**
- `src/symbolic/algebra_helpers.cpp` - `try_combine_like_terms()`, `collect_division_factors()`
- `src/symbolic/symbolic_expression_internal.h` - 数据结构定义

**本次已落地：**
- `make_sorted_product()` 中按底数合并指数
- `try_canonical_factor_quotient()` 中统一处理幂次抵消
- `maybe_canonicalize_polynomial()` 中将单变量多项式重建为稳定顺序
- `to_string_impl()` 中修正负项显示，减少 `+ -...` 形式

**实现说明：**
- 幂次合并使用 `std::map<std::string, PowerGroup>` 按底数结构键分组，相同底数的指数相加
- 分数约分通过将分母指数取负后统一合并，根据最终指数符号决定放入分子或分母
- 单变量多项式通过 `polynomial_coefficients_from_simplified()` 提取系数后，按降幂顺序重建

---

#### 2. **公因子提取（Common Factor Extraction）**

**当前状态：已完成基础版**

**现状：**
- 仅支持"相同的完整剩余表达式"匹配
- 不能分解加法项以检测共享乘法因子
- 缺乏高级代数因式分解

**优化空间：**
```
目标优化示例：
- x*(y+1) + 2*x*(y+1) → 3*x*(y+1)      (识别公因子)
- 2*x + 2*y → 2*(x + y)                (常数提取)
- a*b + a*c → a*(b + c)                (分布律逆向)
- -x - x → -2*x                        (符号处理)
```

**实现策略：**
1. 在 `simplify()` 中添加加法项分解逻辑
2. 对每一项提取公共乘法因子
3. 收集相同基数的系数
4. 渐进式处理嵌套表达式

**预期收益：**
- 显著改进多项式简化
- 更自然的表达式形式
- 在微分和积分后更好的简化结果

**代码位置：**
- `src/symbolic/simplify.cpp` - `simplify_once()` / `simplify_impl()`
- `src/symbolic/algebra_helpers.cpp` - 公因子提取辅助函数

**本次已落地：**
- `try_factor_common_terms()` 已支持公共数值因子和公共符号因子提取
- 加法/减法分支已接入该规则
- 因式输出顺序已调整为更自然的 `a * (b + c)` 而不是 `(b + c) * a`

**实现说明：**
- 使用 `collect_multiplicative_terms()` 提取两边的数值系数和符号因子
- 使用 `expressions_match()`（基于 structural key）匹配公共符号因子
- 使用 `common_numeric_factor()` 计算公共数值因子（整数用欧几里得GCD）
- 重建表达式为 `outer * (left_inner +/- right_inner)` 形式

---

#### 3. **有理表达式化简（Rational Expression Reduction）**

**当前状态：已完成可直接落地的简化层，并补齐单变量多项式 GCD 约分**

**现状：**
- 基础数值因子取消工作
- 已支持多项式级别的整除化简和 GCD 约分
- 已支持多项式长除法
- 可处理如 `(x^2 - 1)/(x - 1) → x + 1` 和部分非整除公共因子约分

**优化空间：**
```
目标优化示例：
- (2*x)/(4*y) → x/(2*y)          (数值因子简化 ✓ 已有)
- (x*y)/x → y                    (变量取消 - 部分有)
- (x^2*y)/x → x*y                (幂次简化 - 缺失)
- (x+1)/(x+1) → 1                (完全表达式取消)
- (x^2 - 1)/(x - 1) → x + 1      (多项式分解 - 缺失)
```

**实现策略：**
1. 完成因子列表取消的幂次追踪（见第1项）
2. 扩展为多项式检测和因式分解
3. ✅ 实现 Euclidean 多项式GCD以检测可约因子
4. ✅ 应用多项式长除法化简可约分数

**预期收益：**
- 更智能的有理表达式处理
- 符号积分和变换中的更佳结果
- 提升用户表达式简洁性

**代码位置：**
- `src/algebra/polynomial.cpp` - 利用现有多项式工具
- `src/symbolic/polynomial_helpers.cpp` - 多项式表达式提取与重建
- `src/symbolic/simplify.cpp` - 集成到 `simplify()`

**本次已落地：**
- 幂次级别的变量约消
- 纯数因子约分
- 可整除多项式分式的长除法化简
- 非整除有理式的单变量多项式 GCD 约分
- 对”表达式 / 纯数字”保留分式形式，提升结果可读性

**实现说明：**
- 分子分母分别收集乘积项后统一处理
- 使用 `PowerGroup` 结构追踪 `(base, exponent)` 对
- 分母指数取负后与分子指数合并
- 最终根据指数正负分别重建分子分母

---

### 第二阶段：符号微积分扩展（中优先级）

#### 4. **符号积分规则扩展**

**当前状态：已完成反三角、特殊根式、有理函数基础规则、混合线性/重复不可约二次部分分式子集、三角恒等式基础扩展和链式替换积分**

**现状：**
- 基于规则而非完整 CAS
- 支持：基本多项式、三角、指数、分段函数、部分反三角形式、平方根特殊形式
- 不支持：完整 Risch 风格积分、多个不同不可约二次因子、不可约高次因子、一般非线性分部积分策略
- 局限于 `a*x + b` 形式的线性替换

**优化空间：**
```
支持的形式示例 ✓：
- ∫ x^n dx
- ∫ sin(ax+b) dx
- ∫ exp(ax+b) dx
- ∫ (polynomial) * exp(...) dx
- ∫ (polynomial) * sin(...) dx

新增支持的形式 ✓：
- ∫ 1/(1+x^2) dx                       (反三角)
- ∫ 1/sqrt(1-x^2) dx                   (反三角)
- ∫ sqrt(1-x^2) dx                     (平方根特殊形式)

新增支持的形式 ✓：
- ∫ (polynomial)/(linear/quadratic) dx             (有理函数基础层)
- ∫ 1/(x^3-x) dx                                   (互异实线性因子部分分式)
- ∫ 1/(x-1)^2 dx, ∫ 1/(x^2+1)^2 dx                 (重复线性/常见重复不可约二次)
- ∫ 1/((x-1)^2*(x+1)) dx                           (混合重复实线性因子)
- ∫ 1/((x-1)*(x^2+1)^2) dx                         (实线性 + 重复不可约二次)
- ∫ sin(x)^2 dx, ∫ cos(x)^2 dx, ∫ tan(x)^2 dx      (基础三角恒等式)
- ∫ sin(x)^3 dx, ∫ cos(x)^3 dx, ∫ sin(x)*cos(x) dx (新增三角幂/乘积恒等式)
- ∫ g'(x)*F(g(x)) dx                               (链式替换，支持 exp/sin/cos/tan)
- ∫ x*sin(x) dx                                    (递推分部积分形式，已有 `polynomial * sin/cos/exp`)

仍需支持的形式 ✗：
- ∫ rational functions with multiple distinct irreducible quadratic/high-order factors
- ∫ broader identity-driven trig forms and full algorithmic integration
```

**实现策略：**
1. 扩展 `integrate_polynomial_times_function()` 支持更多基础情况
2. 实现递推积分法（recurrence relations）
3. ✅ 添加基础有理函数部分分式分解
4. 增加三角和反三角特殊恒等式

**预期收益：**
- 覆盖更多教科书形式
- 改进用户体验
- 为多变量工作奠定基础

**代码位置：**
- `src/symbolic/symbolic_expression_calculus.cpp` - `integral()` 函数
- 需要新增：递推积分辅助函数

**本次相关修正：**
- 维持已有 `polynomial * exp/sin/cos` 递推积分
- 配合规范输出，积分结果更稳定地显示为 `x ^ 3 / 3 + C` 这类形式
- 新增 `1 / (1 + x ^ 2)`、`1 / sqrt(1 - x ^ 2)`、`sqrt(1 - x ^ 2)` 的符号原函数规则
- 新增多项式商积分：先长除法，再对线性/二次分母余项生成 `ln(abs(...))` 或 `atan(...)`
- 新增部分分式积分：覆盖互异实线性因子、单一重复线性因子、混合重复实线性因子和 `1/(x^2+1)^2`
- 新增基础三角恒等式积分：覆盖 `sin(x)^2`、`cos(x)^2`、`tan(x)^2`
- 新增链式替换积分：识别 `g'(x) * F(g(x))` 及常数比例因子

**实现说明：**
- 反三角积分通过 `is_one_plus_variable_squared()` 和 `is_sqrt_one_minus_variable_squared()` 辅助函数检测特定模式
- `sqrt(1-x^2)` 积分使用标准公式 `(x*sqrt(1-x^2) + asin(x))/2`
- 有理积分通过 `try_integrate_polynomial_quotient()` 统一入口处理长除法、一次/二次分母和基础部分分式
- 链式替换通过 `try_integrate_substitution_product()` 比对内层导数与乘法因子
- 积分结果自动调用 `simplify()` 进行化简

---

#### 5. **多变量符号操作（Multi-variable Symbolic Support）**

**当前状态：已完成基础偏导、矩阵型符号 API、连续符号积分和 nonlinear critical point 搜索**

**现状：**
- `diff(expr, var)` 支持显式按变量求偏导
- `diff(expr, var1, var2, ...)` 支持连续混合偏导（命令层通过循环调用 `derivative()` 实现链式求导）
- `gradient(expr, vars...)`、`jacobian([exprs], vars...)`、`hessian(expr, vars...)` 已提供命令入口
- `integral(expr, var)` 和 `integral(expr, var1, var2, ...)` 支持显式积分变量与链式多变量符号积分
- `critical(expr, vars...)` 支持梯度为仿射线性系统时求孤立驻点，并对非线性梯度执行有界 Newton 搜索
- `taylor()`, `limit()`, `extrema()` 仍以单变量工作流为主
- 自由变量分析已有 `identifier_variables()` 基础能力，但多变量替换安全检查仍有限

**优化空间：**
```
已支持的操作：

1. 按变量的显式偏导数：
   - derivative(expr, var1) vs derivative(expr, var2)
   - 多次混合偏导数：∂²f/∂x∂y

2. 多变量工具集：
   - gradient(expr, [x, y, z]) → [∂f/∂x, ∂f/∂y, ∂f/∂z]
   - jacobian(expr_list, var_list) → 矩阵
   - hessian(expr, var_list) → 对称矩阵

仍需完善的操作：

3. 安全的多变量替换：
   - 自由变量分析（哪些变量未被绑定）
   - 防止符号冲突
   - 链式替换安全

4. 多变量微积分：
   - ∫∫ f(x,y) dx dy 的解析支持 ✓（通过链式 `integral(expr, x, y, ...)`）
   - 二次型/仿射梯度的 critical point 分析 ✓
   - 非线性 critical point 的有限起点 Newton 搜索 ✓
   - 全局完备求解和分类仍需后续扩展
```

**实现策略：**
1. ✅ 重构 `derivative()` 支持可选变量列表（通过命令层循环调用实现多变量偏导）
2. ✅ 添加 `gradient()`, `jacobian()`, `hessian()` 包装函数
3. ✅ 实现自由变量分析辅助函数
4. 为多变量替换添加检查

**预期收益：**
- 支持更复杂的数学工作流
- 科学计算更强的适用性
- 为多变量优化打下基础

**代码位置：**
- `src/symbolic/symbolic_expression_calculus.cpp` - 核心微积分
- `src/symbolic/node_parser.cpp` - 变量分析与符号表达式入口
- 核心 API - `symbolic_expression.h` 中的新函数

**本次已落地：**
- `SymbolicExpression::gradient()`、`SymbolicExpression::hessian()`、`SymbolicExpression::jacobian()`
- 命令层 `gradient(...)`、`jacobian(...)`、`hessian(...)`
- `diff(expr, var...)` 支持显式变量和连续偏导（命令层通过链式调用 `derivative()` 实现）
- `integral(expr, var...)` 支持显式变量和连续符号积分
- `critical(expr, var...)` 支持仿射梯度系统和常见非线性梯度的驻点求解

---

### 第三阶段：高级性能与可维护性（长期优先级）

#### 6. **记忆化与性能优化（Memoization）**

**当前状态：已完成节点驻留、结构键缓存、LRU、局部多项式提取记忆化和 Taylor/Pade 导数链缓存，完整显式 CSE pass 尚未实现**

**现状：**
- `simplify()` 和 `derivative()` 已有线程局部 LRU 结构键缓存
- 表达式节点带结构键缓存，构造路径会驻留相同结构节点
- 多项式系数提取在单次调用内复用相同子树结果
- `derivative()` 仍可能产生复杂中间表达式
- Taylor 级数/Pade 逼近会在线程局部缓存中复用相同表达式、变量、展开中心下的中间导数和值
- 暂无生成临时变量/重写表达式的显式共同子表达式消除（CSE）

**优化空间：**
```
性能关键路径：
- simplify() 被频繁调用（解析期间和微积分后）
- 复杂表达式的导数增长迅速
- Taylor 展开需要重复计算相同的子式
```

**实现策略：**
1. ✅ 为 `simplify()` 结果添加记忆表（map<structural_key, simplified>）
2. ✅ 复用结构键以启用缓存
3. 部分完成：表达式节点驻留和局部记忆化共享重复子树；显式 CSE pass 尚未实现
4. ✅ 在 Taylor/Pade 工作流中缓存中间导数

**预期收益：**
- 显著加快复杂表达式的处理
- 改进用户交互响应性
- 支持更深的符号计算

**代码位置：**
- `src/symbolic/node_parser.cpp` - 节点构造、结构键、interning 与 simplify 入口
- `src/symbolic/simplify.cpp` - 主要化简规则与 fixpoint pass
- `src/symbolic/symbolic_expression_calculus.cpp` - 求导缓存、梯度、Jacobian、Hessian 与积分规则

---

#### 7. **符号常数与精确代数（Exact Constants）**

**现状：**
- 符号常数（pi, e）存储为 `double`
- 精确模式仅限于有理数
- 不支持根式形式（如 √2, ∛3）
- 无条件追踪和域感知简化

**优化空间：**
```
理想状态：
- 存储精确表示 pi, e, √2 等
- 支持根式代数运算
- 条件追踪（x > 0 时的简化差异）
- 改进 taylor() 的精确系数

示例：
- tan(pi/4) 应返回 1（精确）而非 0.99999...
- sqrt(2) * sqrt(3) = sqrt(6)（精确）
- ln(exp(x)) = x（有条件，当x实数时）
```

**实现策略：**
1. 定义内部"符号常数"表示（pi, e, 根式）
2. 构建有限的精确代数运算表
3. 为简化规则添加条件检查
4. 改进 `taylor()` 以保留精确系数

**预期收益：**
- 更好的符号计算准确性
- 改进的教学输出
- 为更高级功能奠定基础

**代码位置：**
- `src/symbolic/symbolic_expression_internal.h` - 新的常数表示
- `src/symbolic/node_parser.cpp` / `src/symbolic/simplify.cpp` - 精确符号节点与化简
- `src/core/exact_and_symbolic_render.cpp` - 与已有精确/符号显示模式集成

---

#### 8. **变换框架现代化（Transform Framework Refactor）**

**当前状态：已完成第一步公共线性规则抽取，完整规则表仍保留为后续工作**

**现状：**
- Fourier/Laplace/z 变换仍有各自的具体模式实现
- 线性、减法、取负、常数倍这些公共外壳规则已抽到 `apply_linear_transform_rules()`
- 难以添加新变换类型
- 尚未形成完整可声明的模式匹配表

**优化空间：**
```
改进的架构设想：

unified_transform(expr, rules_table, variable_mappings)
  ├─ pattern_match(expr, rules) → matched_patterns
  ├─ apply_rule(expr, pattern, template)
  └─ return combined_result

规则表结构：
[
  { input_pattern: "sin(w*t)", transform: "fourier" → "pi*delta(w-w0)" },
  { input_pattern: "exp(-a*t)", transform: "laplace" → "1/(s+a)" },
  ...
]
```

**实现策略：**
1. ✅ 抽象出公共线性和常数倍规则入口
2. 将现有变换规则转换为统一表格
3. 实现模式匹配编译器以提升效率
4. 支持用户定义规则（未来功能）

**预期收益：**
- 代码更易维护
- 更容易添加新变换
- 为可扩展性奠定基础
- 减少重复代码

**代码位置：**
- `src/symbolic/transforms.cpp` - 当前变换规则和 `apply_linear_transform_rules()`
- `src/symbolic/symbolic_expression_transforms.cpp` - 公共 API 包装

**本次已落地：**
- Fourier/Laplace/z 及其逆变换共享线性组合、取负和常数倍处理
- 保持现有变换输出不变，降低后续迁移到规则表的风险

---

## 优化优先级矩阵

| 优化项 | 影响范围 | 实现复杂度 | 推荐顺序 |
|--------|--------|---------|--------|
| 1. 规范简化 | 高 | 中 | **1** |
| 2. 公因子提取 | 高 | 中-高 | **2** |
| 3. 有理表达式化简 | 中-高 | 中-高 | **3** |
| 4. 积分规则扩展 | 中 | 中-高 | **4** |
| 5. 多变量支持 | 中 | 高 | **5** |
| 6. 记忆化 | 中 | 中 | **并行于3** |
| 7. 精确常数 | 中 | 中-高 | **6** |
| 8. 变换框架 | 低 | 高 | **7** (长期) |

---

## 当前测试覆盖与基准

**测试文件：** `test/tests.cpp`

**关键测试区域（已验证）：**
- 基础算术简化
- 相似项合并（`2*x + 3*x -> 5*x`）
- 变量因子取消（`(x*x)/x -> x`）
- 数值因子取消（`(2*x)/4 -> x/2`）
- 多项式样式清理
- 三角简化（部分）
- 指数/对数简化（部分）
- 反三角积分规则（`integral(1/(1+x^2))`, `integral(1/sqrt(1-x^2))`）
- 平方根特殊积分（`integral(sqrt(1-x^2))`）
- 有理积分规则（长除法、线性/二次分母、互异实线性因子部分分式）
- 重复线性因子、混合重复实线性因子、常见重复不可约二次因子和链式替换积分
- 基础三角恒等式积分（`sin^2/cos^2/tan^2`）
- 多变量偏导（`diff(x^2*y + y^3, y)`）
- gradient/jacobian/hessian 输出验证
- chained symbolic integral、affine-gradient 和 nonlinear `critical(...)` 输出验证

**建议的测试增强：**
```
优化 1 - 规范简化：
  - 幂次组合 (x^3 / x^2)
  - 重复因子规范化 (x*x*x)
  - 已在 `make_sorted_product()` 和 `try_canonical_factor_quotient()` 中实现

优化 2 - 公因子提取：
  - 多项式提取 (2x + 2y)
  - 因子识别 (x(y+1) + 2x(y+1))
  - 已在 `try_factor_common_terms()` 中实现

优化 3 - 有理表达式：
  - 多项式长除法验证
  - 多项式 GCD 约分验证
  - 已在 `try_canonical_factor_quotient()` 中实现

优化 4 - 符号积分：
  - 反三角积分规则 ✓ 已有测试
  - 平方根特殊形式 ✓ 已有测试
  - 有理函数长除法/部分分式 ✓ 已有测试
  - 重复因子和链式替换 ✓ 已有测试

优化 5 - 多变量符号操作：
  - 显式偏导和混合偏导 ✓ 已有测试
  - gradient/jacobian/hessian 输出 ✓ 已有测试
  - 链式多变量积分、仿射梯度和非线性 critical point ✓ 已有测试

仍待测试：
  - 混合重复不可约二次因子/一般不可约高次因子的完整部分分式
  - 更完整的三角恒等式驱动积分和非线性分部积分
  - nonlinear critical point 的全局完备性和分类
```

**重新进入计划：**
1. 运行 `make test` 获取当前基准
2. 检查 `src/symbolic/simplify.cpp`、`src/symbolic/algebra_helpers.cpp` 和 `src/symbolic/symbolic_expression_calculus.cpp` 中的关键函数
3. 如需新增优化规则，参考已实现的 `make_sorted_product()`、`try_canonical_factor_quotient()`、`try_factor_common_terms()` 模式
4. 避免回归：每次提交前检查完整测试套件

---

## 快速上手清单

对于下一个工作会话：

- [ ] 重新阅读本分析文档
- [ ] 读取 `src/symbolic/simplify.cpp`、`src/symbolic/algebra_helpers.cpp` 和 `src/symbolic/symbolic_expression_calculus.cpp`
- [ ] 搜索关键函数：`simplify()`, `try_combine_like_terms()`, `make_sorted_product()`, `try_canonical_factor_quotient()`, `try_factor_common_terms()`
- [ ] 运行 `make test` 确认基准
- [ ] 为新的优化项创建回归测试
- [ ] 实现单个优化规则（避免一次做太多）

---

## 相关文档引用

- `KNOWN_LIMITATIONS.md` - 系统范围限制
- `ROADMAP.md` - 长期路线图和指导原则
- `docs/archive/SIMPLIFY_IMPROVEMENTS.md` - 历史简化边界记录
- `src/symbolic/symbolic_expression.h` - 公共 API
- `src/symbolic/symbolic_expression_internal.h` - 内部结构

---

**文档创建时间：** 2026-04-24  
**版本：** 1.4  
**状态：** 阶段 1 已落地，阶段 2 已补齐有理积分基础层、混合重复实线性部分分式、基础三角恒等式积分、链式替换、链式多变量积分和 nonlinear critical point；阶段 3 已启动变换框架公共线性规则抽取与 simplify/derivative 局部缓存
