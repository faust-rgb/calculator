# 符号计算优化分析

## 执行摘要

当前符号计算模块已实现了基础的微积分、简化、变换等功能。通过对代码和文档的分析，识别出 5 个核心优化方向，按优先级分为3个阶段。本分析基于已实现功能的边界和现有限制。

## 本次落地结果（2026-04-24）

本轮已优先完成文档中“在现有架构下可直接实现”的部分，并补齐回归测试，重点落在第一阶段的规范化与化简稳定性：

- 已完成：乘法因子按幂次归并，支持 `x * x * x -> x ^ 3`、`x ^ 3 / x ^ 2 -> x`
- 已完成：公因子提取，支持 `2 * x + 2 * y -> 2 * (x + y)`、`a * b + a * c -> a * (b + c)`
- 已完成：多项式商化简，支持 `(x ^ 2 - 1) / (x - 1) -> x + 1`
- 已完成：单变量多项式规范输出，改善 `x + x ^ 2 + 1` 这类表达式的稳定顺序
- 已完成：安全的 `exp(ln(...))` 简化，只在已知正值表达式上折叠，避免把 `exp(ln(x))` 误化简为 `x`
- 已完成：输出规范收敛，优先保留 `pi / 2`、`x ^ 3 / 3` 这类“纯数分母”形式，并修正因式顺序

仍保留为后续工作：

- 有理函数的 GCD/部分分式体系
- 更广的积分规则（反三角、分部积分、特殊根式）
- 多变量符号 API（gradient/jacobian/hessian）
- simplify/derivative 的记忆化与结构缓存

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
- `src/symbolic/symbolic_expression_core.cpp` - `try_combine_like_terms()`, `collect_division_factors()`
- `src/symbolic/symbolic_expression_internal.h` - 数据结构定义

**本次已落地：**
- `make_sorted_product()` 中按底数合并指数
- `try_canonical_factor_quotient()` 中统一处理幂次抵消
- `maybe_canonicalize_polynomial()` 中将单变量多项式重建为稳定顺序
- `to_string_impl()` 中修正负项显示，减少 `+ -...` 形式

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
- `src/symbolic/symbolic_expression_core.cpp` - `simplify()` 函数
- 需要新增：`extract_common_factors()` 辅助函数

**本次已落地：**
- `try_factor_common_terms()` 已支持公共数值因子和公共符号因子提取
- 加法/减法分支已接入该规则
- 因式输出顺序已调整为更自然的 `a * (b + c)` 而不是 `(b + c) * a`

---

#### 3. **有理表达式化简（Rational Expression Reduction）**

**当前状态：已完成可直接落地的简化层，未做 GCD/部分分式**

**现状：**
- 基础数值因子取消工作
- 缺乏多项式级别的分解和化简
- 不支持多项式长除法
- 无法处理如 `(x^2 - 1)/(x - 1) → x + 1` 的形式

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
3. 实现 Euclidean 多项式GCD以检测可约因子
4. 应用多项式长除法化简可约分数

**预期收益：**
- 更智能的有理表达式处理
- 符号积分和变换中的更佳结果
- 提升用户表达式简洁性

**代码位置：**
- `src/algebra/polynomial.cpp` - 利用现有多项式工具
- `src/symbolic/symbolic_expression_core.cpp` - 集成到 `simplify()`

**本次已落地：**
- 幂次级别的变量约消
- 纯数因子约分
- 可整除多项式分式的长除法化简
- 对“表达式 / 纯数字”保留分式形式，提升结果可读性

---

### 第二阶段：符号微积分扩展（中优先级）

#### 4. **符号积分规则扩展**

**当前状态：部分完成，未扩展到分析文档列出的高级规则**

**现状：**
- 基于规则而非完整
- 支持：基本多项式、三角、指数、分段函数
- 不支持：分部积分、替换、高级恒等式
- 局限于 `a*x + b` 形式的线性替换

**优化空间：**
```
支持的形式示例 ✓：
- ∫ x^n dx
- ∫ sin(ax+b) dx
- ∫ exp(ax+b) dx
- ∫ (polynomial) * exp(...) dx
- ∫ (polynomial) * sin(...) dx

需要支持的形式 ✗：
- ∫ (polynomial)/(polynomial) dx      (有理函数)
- ∫ x*sin(x) dx                        (分部积分)
- ∫ 1/(1+x^2) dx                       (反三角)
- ∫ sqrt(1-x^2) dx                     (平方根特殊形式)
```

**实现策略：**
1. 扩展 `integrate_polynomial_times_function()` 支持更多基础情况
2. 实现递推积分法（recurrence relations）
3. 添加有理函数部分分式分解
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

---

#### 5. **多变量符号操作（Multi-variable Symbolic Support）**

**现状：**
- `diff()`, `integral()`, `taylor()`, `limit()`, `extrema()` 都仅限单变量
- 部分导数功能不存在
- 无梯度、Jacobian、Hessian 辅助函数
- 无多变量替换和自由变量分析

**优化空间：**
```
需要支持的操作：

1. 按变量的显式偏导数：
   - derivative(expr, var1) vs derivative(expr, var2)
   - 多次混合偏导数：∂²f/∂x∂y

2. 多变量工具集：
   - gradient(expr, [x, y, z]) → [∂f/∂x, ∂f/∂y, ∂f/∂z]
   - jacobian(expr_list, var_list) → 矩阵
   - hessian(expr, var_list) → 对称矩阵

3. 安全的多变量替换：
   - 自由变量分析（哪些变量未被绑定）
   - 防止符号冲突
   - 链式替换安全

4. 多变量微积分：
   - ∫∫ f(x,y) dx dy 的解析支持
   - 极值的 critical point 分析
```

**实现策略：**
1. 重构 `derivative()` 支持可选变量列表
2. 添加 `gradient()`, `jacobian()`, `hessian()` 包装函数
3. 实现自由变量分析辅助函数
4. 为多变量替换添加检查

**预期收益：**
- 支持更复杂的数学工作流
- 科学计算更强的适用性
- 为多变量优化打下基础

**代码位置：**
- `src/symbolic/symbolic_expression_calculus.cpp` - 核心微积分
- `src/symbolic/symbolic_expression_core.cpp` - 变量分析
- 核心 API - `symbolic_expression.h` 中的新函数

---

### 第三阶段：高级性能与可维护性（长期优先级）

#### 6. **记忆化与性能优化（Memoization）**

**现状：**
- 每次 `simplify()` 重新计算
- `derivative()` 可能产生复杂中间表达式
- Taylor 级数/Pade 逼近无共享子表达式重用
- 无结构哈希缓存

**优化空间：**
```
性能关键路径：
- simplify() 被频繁调用（解析期间和微积分后）
- 复杂表达式的导数增长迅速
- Taylor 展开需要重复计算相同的子式
```

**实现策略：**
1. 为 `simplify()` 结果添加记忆表（map<structural_key, simplified>）
2. 计算结构哈希以启用缓存键
3. 在表达式树上实现共同子表达式消除（CSE）
4. 在 Taylor/Pade 工作流中缓存中间导数

**预期收益：**
- 显著加快复杂表达式的处理
- 改进用户交互响应性
- 支持更深的符号计算

**代码位置：**
- `src/symbolic/symbolic_expression_core.cpp` - `simplify()`
- `src/symbolic/symbolic_expression_internal.h` - 新的缓存结构

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
- `src/symbolic/symbolic_expression_core.cpp` - 精确运算
- `src/core/calculator.cpp` - 与已有精确模式集成

---

#### 8. **变换框架现代化（Transform Framework Refactor）**

**现状：**
- Fourier/Laplace/z 变换有独立实现
- 代码重复较多（pattern matching 逻辑分散）
- 难以添加新变换类型
- 缺乏可重用的模式匹配表

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
1. 抽象出通用的模式匹配和规则应用引擎
2. 将现有变换规则转换为统一表格
3. 实现模式匹配编译器以提升效率
4. 支持用户定义规则（未来功能）

**预期收益：**
- 代码更易维护
- 更容易添加新变换
- 为可扩展性奠定基础
- 减少重复代码

**代码位置：**
- `src/symbolic/symbolic_expression_transforms.cpp` - 现有实现
- 需要新增：`transform_framework.h/.cpp`

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

**现有基准状态：**
```
Passed: 373 (最新简化工作后)
关键测试区域：
  ✓ 基础算术简化
  ✓ 相似项合并
  ✓ 变量因子取消
  ✓ 数值因子取消
  ✓ 多项式样式清理
  ✓ 三角简化（部分）
  ✓ 指数/对数简化（部分）
```

**建议的测试增强：**
```
优化 1 - 规范简化：
  - 幂次组合 (x^3 / x^2)
  - 重复因子规范化 (x*x*x)

优化 2 - 公因子提取：
  - 多项式提取 (2x + 2y)
  - 因子识别 (x(y+1) + 2x(y+1))

优化 3 - 有理表达式：
  - 多项式GCD测试
  - 长除法验证

等等...
```

**重新进入计划：**
1. 运行 `make test` 获取当前基准
2. 从优化 1 开始（规范简化）
3. 逐个添加简化规则，运行测试
4. 避免回归：每次提交前检查完整测试套件
5. 定期参考 `SIMPLIFY_IMPROVEMENTS.md`

---

## 快速上手清单

对于下一个工作会话：

- [ ] 重新阅读本分析文档
- [ ] 读取 `src/symbolic/symbolic_expression_core.cpp`
- [ ] 搜索关键函数：`simplify()`, `try_combine_like_terms()`, `collect_division_factors()`
- [ ] 运行 `make test` 确认基准
- [ ] 为选中的优化项创建回归测试
- [ ] 实现单个简化规则（避免一次做太多）
- [ ] 定期检查 `SIMPLIFY_IMPROVEMENTS.md` 中的风险区域

---

## 相关文档引用

- `SIMPLIFY_IMPROVEMENTS.md` - 当前简化边界和已验证示例
- `KNOWN_LIMITATIONS.md` - 系统范围限制
- `ROADMAP.md` - 长期路线图和指导原则
- `src/symbolic/symbolic_expression.h` - 公共 API
- `src/symbolic/symbolic_expression_internal.h` - 内部结构

---

**文档创建时间：** 2026-04-24  
**版本：** 1.0  
**状态：** 概念验证 (未实现任何建议)
