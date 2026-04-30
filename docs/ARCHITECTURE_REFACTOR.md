# 计算器架构优化方案

## 一、目录结构重组

### 当前结构
```
src/
├── core/           # 过于臃肿 (7645 行)
├── analysis/       # 职责过多 (积分、ODE、优化、级数)
├── symbolic/       # 符号计算
├── polynomial/     # 多项式 (与 symbolic 重叠)
├── matrix/         # 矩阵运算
├── math/           # 基础数学函数
├── script/         # 脚本引擎
├── statistics/     # 统计
├── dsp/            # 信号处理
├── plot/           # 绑图
└── app/            # 主程序
```

### 优化后结构
```
src/
├── api/                    # 公共 API 层
│   ├── calculator.h        # 主接口
│   ├── types.h             # 公共类型定义
│   └── exceptions.h        # 异常定义
│
├── core/                   # 核心运行时 (精简)
│   ├── context.h/cpp       # 计算上下文 (变量、函数表)
│   ├── registry.h/cpp      # 命令注册表
│   ├── dispatcher.h/cpp    # 命令分发器
│   └── state.h/cpp         # 状态持久化
│
├── parser/                 # 解析器层 (新建)
│   ├── base_parser.h       # 基类
│   ├── decimal_parser.h/cpp
│   ├── exact_parser.h/cpp
│   ├── matrix_parser.h/cpp
│   ├── symbolic_parser.h/cpp
│   └── parser_factory.h    # 解析器工厂
│
├── runtime/                # 运行时层 (新建)
│   ├── evaluator.h/cpp     # 表达式求值
│   ├── script_engine.h/cpp # 脚本执行
│   └── scope.h/cpp         # 作用域管理
│
├── types/                  # 类型层 (新建)
│   ├── rational.h/cpp      # 有理数
│   ├── precise_decimal.h/cpp
│   ├── stored_value.h/cpp  # 存储值
│   └── function.h/cpp      # 用户函数
│
├── math/                   # 基础数学 (保持)
│   ├── mymath.h/cpp
│   └── special_functions.h/cpp
│
├── numeric/                # 数值计算 (从 analysis 拆分)
│   ├── integration/
│   ├── optimization/
│   ├── ode/
│   └── linear_algebra/
│
├── symbolic/               # 符号计算 (保持，合并 polynomial)
│   ├── expression.h/cpp
│   ├── calculus.h/cpp
│   ├── simplify.h/cpp
│   └── polynomial.h/cpp    # 从 polynomial/ 移入
│
├── matrix/                 # 矩阵运算 (保持)
├── statistics/             # 统计 (保持)
├── signal/                 # 信号处理 (从 dsp 重命名)
├── plot/                   # 绑图 (保持)
└── app/                    # 主程序 (保持)
```

## 二、核心模块职责划分

### api/ - 公共接口层
```cpp
// api/calculator.h - 唯一对外头文件
class Calculator {
public:
    std::string evaluate(const std::string& expr);
    std::string execute_script(const std::string& source);
    void set_variable(const std::string& name, double value);
    // ...
};

// api/types.h - 公共类型
struct Matrix;  // 前向声明
struct Result;  // 统一返回类型

// api/exceptions.h
class CalculatorError;
class SyntaxError;
class DomainError;
```

### core/ - 核心运行时 (精简到 ~2000 行)
```cpp
// core/context.h - 计算上下文
struct ComputeContext {
    VariableTable* variables;
    FunctionTable* functions;
    Config* config;
};

// core/registry.h - 命令注册
struct CommandEntry {
    std::string name;
    CommandHandler handler;
    std::string help_text;
};
class CommandRegistry {
    void register_module(Module* module);
    bool dispatch(const std::string& cmd, Context& ctx);
};

// core/dispatcher.h - 命令分发
class Dispatcher {
    bool try_command(const std::string& input, Context& ctx);
};
```

### parser/ - 解析器层
```cpp
// parser/base_parser.h - 解析器基类
class BaseParser {
protected:
    void skip_spaces();
    bool match(char c);
    std::string parse_identifier();
    double parse_number();
    // ...
};

// parser/parser_factory.h - 统一入口
class ParserFactory {
public:
    enum class Mode { Decimal, Exact, Matrix, Symbolic };
    
    static std::unique_ptr<BaseParser> create(Mode mode, std::string source);
    
    template<typename T>
    T parse(const std::string& source, Mode mode);
};
```

### types/ - 类型层
```cpp
// types/stored_value.h - 存储值 (拆分自 calculator_internal_types.h)
struct StoredValue {
    double decimal;
    Matrix matrix;
    std::string string;
    ValueType type;
};

// types/rational.h - 有理数
struct Rational { ... };

// types/precise_decimal.h - 精确小数
struct PreciseDecimal { ... };
```

## 三、依赖关系优化

### 当前依赖问题
```
calculator_internal_types.h (被 27 个文件依赖)
    ├── 包含 Rational
    ├── 包含 PreciseDecimal
    ├── 包含 StoredValue
    ├── 包含 CustomFunction
    ├── 包含 Calculator::Impl
    └── 包含各种回调类型
```

### 优化后依赖
```
api/types.h          <- 只被对外接口依赖
types/stored_value.h <- 被需要存储的模块依赖
types/rational.h     <- 只被解析器依赖
types/function.h     <- 只被运行时依赖

依赖方向:
api -> core -> parser -> types -> math
                ↓
            runtime -> types
```

## 四、解析器统一方案

### 当前问题
- 6 个解析器类分散各处
- 词法分析代码重复
- 接口不统一

### 优化方案
```cpp
// parser/base_parser.h
class BaseParser {
public:
    virtual ~BaseParser() = default;
    
    // 统一接口
    virtual Value parse() = 0;
    
protected:
    // 共享词法分析
    std::string source_;
    size_t pos_;
    
    void skip_spaces();
    bool match(char c);
    bool match_string(const std::string& s);
    std::string parse_identifier();
    double parse_number();
    std::string parse_string_literal();
    // ...
};

// parser/decimal_parser.h
class DecimalParser : public BaseParser {
public:
    double parse() override;
    
private:
    double parse_expression();
    double parse_term();
    double parse_factor();
    // ...
};

// parser/parser_factory.h - 统一入口
class ParserFactory {
public:
    static double parse_decimal(const std::string& source, Context& ctx);
    static Rational parse_exact(const std::string& source, Context& ctx);
    static Matrix parse_matrix(const std::string& source, Context& ctx);
    static SymbolicExpr parse_symbolic(const std::string& source, Context& ctx);
};
```

## 五、模块注册机制优化

### 当前机制
```cpp
// calculator_commands.cpp 中硬编码
const std::vector<FunctionCommandRegistration> command_registry = {
    {polynomial_ops::is_polynomial_command, handle_polynomial_command},
    {series_ops::is_series_command, handle_series_command},
    // ... 12 个模块
};
```

### 优化方案 - 插件式注册
```cpp
// core/module.h
class Module {
public:
    virtual ~Module() = default;
    virtual std::string name() const = 0;
    virtual std::vector<CommandDef> commands() const = 0;
    virtual void initialize(Context& ctx) {}
    virtual void shutdown(Context& ctx) {}
};

// core/registry.h
class ModuleRegistry {
public:
    void register_module(std::unique_ptr<Module> module) {
        modules_.push_back(std::move(module));
        for (const auto& cmd : modules_.back()->commands()) {
            commands_[cmd.name] = cmd;
        }
    }
    
private:
    std::vector<std::unique_ptr<Module>> modules_;
    std::map<std::string, CommandDef> commands_;
};

// 各模块自注册
// numeric/integration_module.cpp
class IntegrationModule : public Module {
public:
    std::string name() const override { return "integration"; }
    std::vector<CommandDef> commands() const override {
        return {
            {"double_integral", handle_double_integral, "二重积分"},
            {"triple_integral", handle_triple_integral, "三重积分"},
            // ...
        };
    }
};

REGISTER_MODULE(IntegrationModule);  // 自动注册
```

## 六、实施步骤

### 第一阶段：拆分 core (优先级高)
1. 创建 `parser/` 目录，迁移所有解析器
2. 创建 `types/` 目录，拆分 `calculator_internal_types.h`
3. 创建 `runtime/` 目录，迁移 `script_runtime.cpp`
4. 精简 `core/` 到核心分发逻辑

### 第二阶段：合并模块
1. 合并 `polynomial/` 到 `symbolic/`
2. 重命名 `dsp/` 为 `signal/`
3. 拆分 `analysis/` 为 `numeric/` 下的子模块

### 第三阶段：接口统一
1. 创建 `api/` 目录，定义公共接口
2. 实现模块自注册机制
3. 统一解析器工厂

### 第四阶段：文档和测试
1. 更新架构文档
2. 添加模块开发指南
3. 补充单元测试

## 七、预期收益

| 指标 | 当前 | 优化后 |
|------|------|--------|
| core/ 行数 | 7645 | ~2000 |
| 单文件最大行数 | 3618 | <1500 |
| 解析器数量 | 6 个分散 | 6 个统一 |
| 依赖 calculator_internal_types.h | 27 文件 | ~10 文件 |
| 新增模块步骤 | 修改 3+ 文件 | 实现 Module 接口 |
