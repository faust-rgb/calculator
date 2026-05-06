# Architecture

## Overview

This project is a C++ command-line calculator with three overlapping roles:

- a normal scientific calculator
- a rational/exact-mode calculator
- a small programmer/CAS-style terminal tool

The implementation is intentionally lightweight and self-contained. Core math is
implemented in project code instead of relying on the standard math library for
the main functions.

## Module Layering

```
┌─────────────────────────────────────────────────────────────┐
│                        app/ (应用层)                         │
│                     main.cpp, CLI 入口                       │
├─────────────────────────────────────────────────────────────┤
│                       module/ (模块层)                       │
│              CalculatorModule, 外部功能扩展点                 │
├─────────────────────────────────────────────────────────────┤
│                      analysis/ (分析层)                      │
│            微积分、极限、ODE、向量场等高级数学功能              │
├─────────────────────────────────────────────────────────────┤
│    ┌──────────┐  ┌──────────┐  ┌──────────┐  ┌──────────┐   │
│    │  plot/   │  │  stats/  │  │  dsp/    │  │ symbolic/│   │
│    │  绘图    │  │  统计    │  │  信号处理 │  │  符号计算 │   │
│    └──────────┘  └──────────┘  └──────────┘  └──────────┘   │
├─────────────────────────────────────────────────────────────┤
│                       core/ (核心层)                         │
│     Calculator, Scope, Exceptions, Utils, Services          │
├─────────────────────────────────────────────────────────────┤
│                     execution/ (执行层)                      │
│    CommandRegistry, ScriptRuntime, VariableResolver         │
├─────────────────────────────────────────────────────────────┤
│                       parser/ (解析层)                       │
│   ExpressionParser, CommandParser, ScriptParser, AST        │
├─────────────────────────────────────────────────────────────┤
│                        io/ (IO层)                            │
│                  文件读写、导入导出                            │
├─────────────────────────────────────────────────────────────┤
│                       types/ (类型层)                        │
│          StoredValue, Function, Matrix, Rational            │
├─────────────────────────────────────────────────────────────┤
│                        math/ (数学层)                        │
│              基础数学函数、数值计算、辅助工具                    │
└─────────────────────────────────────────────────────────────┘
```

### Dependency Rules

1. **Upward dependency prohibited**: Lower layers cannot depend on upper layers
2. **Same-layer dependency**: Same-layer modules can depend on each other, but minimize
3. **Cross-layer dependency**: Can only depend downward, never upward

### Module Responsibilities

#### types/ (Type Layer)
- `stored_value.h` - Stored value types (scalar, matrix, string, etc.)
- `function.h` - Function type definitions

#### math/ (Math Layer)
- `mymath.h` - Basic math functions
- `helpers/` - Math helper tools (base conversion, integer operations, etc.)
- `matrix/` - Matrix operations
- `precise/` - Precise computation (rational, high-precision decimal)

#### parser/ (Parser Layer)
- `expression_ast.h` - Expression AST definition and compilation
- `command_parser.h` - Command parser
- `script_parser.h` - Script parser
- `unified_expression_parser.h` - Unified expression parsing
- `parser_utils.h` - Parsing utility functions
- `token_types.h` - Token type definitions
- `lazy_token_stream.h` - Lazy token stream

#### execution/ (Execution Layer)
- `command_registry.h` - Command registry
- `script_runtime.h` - Script runtime
- `variable_resolver.h` - Variable resolver
- `builtin_constants.h` - Built-in constants
- `inline_expander.h` - Inline function expansion
- `script_signal.h` - Script control flow signals

#### core/ (Core Layer)
- `calculator.h` - Main calculator class (Pimpl pattern)
- `calculator_impl.h` - Calculator implementation
- `scope.h` - Scope management
- `calculator_exceptions.h` - Exception types
- `service_interfaces.h` - Service interfaces
- `string_utils.h` - String utilities
- `format_utils.h` - Formatting utilities
- `expression_utils.h` - Expression utilities

#### module/ (Module Layer)
- `calculator_module.h` - Module interface definition

#### analysis/ (Analysis Layer)
- Calculus, limits, Taylor expansion, ODE, etc.

#### plot/ (Plot Layer)
- Function plotting, 3D plotting

#### statistics/ (Statistics Layer)
- Probability distributions, statistical functions

#### dsp/ (Signal Processing Layer)
- FFT, filter design, spectrum analysis

#### symbolic/ (Symbolic Layer)
- Symbolic expressions, symbolic integration/differentiation

#### io/ (IO Layer)
- File read/write, state persistence

## Main Files

- `src/app/main.cpp`
  Terminal interaction, history, autocomplete, and command dispatch
- `src/core/decimal_parser.cpp`
  Standard double-based expression parsing
- `src/core/precise_decimal_parser.cpp`
  High precision decimal parsing and arithmetic
- `src/core/exact_and_symbolic_render.cpp`
  Exact-mode parsing plus symbolic-constants rendering
- `src/core/calculator_commands.cpp`
  Command-style function processing and advanced calculator entry points
- `src/core/state_persistence.cpp`
  Save/load format handling
- `src/core/calculator_lifecycle.cpp`
  Calculator construction, mode toggles, completion lists, and lightweight runtime helpers
- `src/core/variable_resolver.cpp`
  Variable lookup with scoped resolution and FlatScopeStack support
- `src/core/expression_compiler.h`
  Expression cache for AST reuse in loops
- `src/core/expression_ast.cpp`
  Compiled expression AST for high-performance evaluation
- `src/core/calculator.h`
  Public calculator API
- `src/script/script_parser.cpp`
  Script lexer and parser (Python-style syntax)
- `src/script/script_runtime.cpp`
  Script execution engine with FlatScopeStack optimization
- `src/script/script_ast.h`
  Script AST node definitions
- `src/symbolic/node_parser.cpp`
  Symbolic expression parsing and node construction
- `src/symbolic/simplify.cpp`
  Symbolic simplification rules
- `src/symbolic/algebra_helpers.cpp`
  Symbolic substitution and algebra support helpers
- `src/symbolic/polynomial_helpers.cpp`
  Polynomial-oriented symbolic helpers
- `src/symbolic/symbolic_expression_calculus.cpp`
  Symbolic differentiation and integration rules
- `src/symbolic/symbolic_expression_transforms.cpp`
  Fourier/Laplace/z transform entry points
- `src/math/mymath.cpp`
  Core numerical algorithms and domain handling
- `src/math/mymath_special_functions.cpp`
  Trigonometric, inverse-trigonometric, gamma, and Bessel-related implementations
- `src/matrix/matrix.cpp`
  Matrix storage, core operations, statistics, signal helpers, and polynomial helpers
- `src/matrix/matrix_expression.cpp`
  Matrix expression parsing, matrix literals, and matrix function dispatch
- `src/matrix/matrix_linear_algebra.cpp`
  Inversion, decompositions, eigensolvers, RREF, and related linear-algebra routines
- `src/math/mymath.h`
  Math declarations and shared constants
- `test/main.cpp`
  Test entry point and suite orchestrator
- `test/suites/*.cpp`
  Modular test suites (Core, Analysis, Symbolic, etc.)
- `test/test_helpers.h`
  Shared test utilities and structures

## Directory Layout

- `src/app`
  CLI entry point
- `src/core`
  Calculator runtime and API
- `src/math`
  Custom numeric functions
- `src/matrix`
  Matrix representation and algorithms
- `src/analysis`
  Function analysis helpers
- `src/algebra`
  Polynomial operations
- `src/symbolic`
  Symbolic expression support
- `src/script`
  Script AST and parser
- `src/statistics`
  Statistics, probability distributions, and hypothesis testing
- `src/time`
  Time functions (timestamp, formatting, timers, sleep)
- `src/io`
  File I/O operations
- `test`
  Regression tests and runnable example scripts
- `bin`
  Ignored build outputs such as `calculator` and `calculator_tests`
- `build`
  Ignored object files and generated dependency files for incremental builds

Large implementation areas are split into private implementation `.cpp` files
with internal headers for shared declarations. Current internal split headers
include:

- `src/core/calculator_internal_types.h` (includes FlatScopeStack, VariableSlot, ExpressionCache)
- `src/core/variable_resolver.h`
- `src/core/expression_ast.h`
- `src/math/mymath_internal.h`
- `src/matrix/matrix_internal.h`
- `src/symbolic/symbolic_expression_internal.h`

## Execution Flow

User input goes through this rough path:

1. `src/app/main.cpp` reads a line from the terminal
2. command-style inputs such as `:help`, `:vars`, `:save`, `:load` are handled first
3. special display-only expressions such as `factor(...)` are handled explicitly
4. other expressions are passed to `Calculator::process_line(...)`
5. `Calculator` either:
   - treats the line as an assignment, or
   - evaluates it for display

## Parsing Model

Scalar parser implementations are split across focused implementation files
under `src/core`.

### `DecimalParser`

Used for standard evaluation with `double`.

Implementation file: `src/core/decimal_parser.cpp`.

This path supports:

- all ordinary math functions
- all programmer-style integer functions
- prefixed integer literals such as `0b1010`, `0o77`, `0xFF`
- variable lookup

### `ExactParser`

Used when exact fraction mode is enabled.

Implementation file: `src/core/exact_and_symbolic_render.cpp`.

This path tries to preserve expressions as `Rational` values.

It supports only operations/functions that still make sense in rational form.
If a feature cannot stay exact, the code throws `ExactModeUnsupported`, and the
caller falls back to decimal display.

### `PreciseDecimalParser`

Used for decimal text that needs more precision than `double` display can safely
preserve.

Implementation file: `src/core/precise_decimal_parser.cpp`.

## Stored Values

Calculator state lives in `Calculator::Impl`.

Variables are stored in:

- `std::map<std::string, StoredValue> variables`

`StoredValue` can hold either:

- an exact rational value
- or a decimal fallback value
- and optionally a scalar symbolic display string used by symbolic constants mode

This split allows exact mode to keep fractions when possible, while still
supporting assignments like `y = sin(pi / 2)`.

## Symbolic Constants Flow

Symbolic constants mode is a scalar display feature layered on top of the normal
numeric engine.

When enabled:

1. ordinary evaluation still computes the numeric result
2. the symbolic render path tries to rebuild a scalar expression using `pi` and
   `e`
3. if it succeeds, the symbolic form is displayed and can also be stored with
   scalar variables for later reuse
4. if it cannot represent the expression safely, the code falls back to normal
   numeric display

The current implementation lives mainly in:

- `src/core/calculator_lifecycle.cpp` and `src/core/exact_and_symbolic_render.cpp`
  mode flag, storage, and display dispatch
- `src/symbolic/node_parser.cpp`, `src/symbolic/simplify.cpp`, and
  `src/symbolic/algebra_helpers.cpp`
  symbolic parsing, simplification, and substitution rules

## Display-only Features

Some features are not plain numeric expressions and therefore are handled
outside normal parser return values.

These currently include:

- `factor(...)`
- `bin(...)`
- `oct(...)`
- `hex(...)`
- `base(...)`

These return formatted strings rather than just numbers.

## Symbolic Algebra Core

The symbolic engine (in `src/symbolic`) uses an **interning** pattern to ensure structural uniqueness of nodes. This enables:

- $O(1)$ structural comparison (via pointer equality or structural key checks).
- Efficient caching for simplification and derivatives.
- Reduced memory footprint for complex expressions.

**Eviction Strategy:** The interning pool uses an incremental LRU eviction strategy. When the pool (default 8192 nodes) is full, it performs a limited scan to prune expired weak references and evicts the oldest entries if necessary, ensuring stable $O(1)$ amortized insertion performance.

### Exact Constants and Normalization

Version 1.5 introduced dedicated node types for `pi` and `e` (`NodeType::kPi` and `NodeType::kE`). This prevents symbolic constants from being prematurely collapsed into floating-point numbers during simplification.

**Normalization Rules:**
- `e ^ x` is automatically normalized to `exp(x)` during the `simplify()` pass. This allows the engine to reuse established exponential rules for calculus and transforms while keeping the external display consistent.

## Terminal UX

`src/app/main.cpp` implements a custom REPL with raw terminal support:

- **Navigation:** Left/Right arrow keys for inline cursor movement, Home (`Ctrl+A`), and End (`Ctrl+E`).
- **History:** Up/Down arrow keys for command recall.
- **Editing:** `Backspace`, `Ctrl+D` (delete character at cursor), and `Ctrl+K` (clear to end of line).
- **Autocomplete:** Single `Tab` for completion, double `Tab` for listing candidates.
- **Fallback:** Transparently falls back to standard `getline` when input is not a TTY.

## Persistence

Variable state can be saved and restored with:

- `:save file`
- `:load file`

The persistence format is simple tab-separated text written by
`Calculator::save_state(...)` and read by `Calculator::load_state(...)`.
Current format version supports scalar variables, strings, matrices,
expression-style custom functions, and script functions.

## Numeric Design Notes

Important numeric behavior:

- near-zero results are normalized to `0`
- decimal results within `1e-10` of the nearest integer are displayed as that integer
- exact mode does not try symbolic algebra beyond rational arithmetic
- scientific notation input such as `1e-3`, `1e20`, and `1e-300` is supported
- negative-base fractional powers are supported only when they correspond to a
  real odd-denominator rational exponent
- base conversion output supports bases `2..16`
- bitwise functions only accept integers

## Script Engine Performance

The script engine (in `src/script`) uses several optimizations for efficient execution:

### FlatScopeStack

Traditional scope implementations use `std::vector<std::map<std::string, StoredValue>>`,
which has O(log n) lookup per scope level. The `FlatScopeStack` replaces this with
a flat array of `VariableSlot` structures:

- All variables stored in contiguous memory (`std::vector<VariableSlot>`)
- Scope boundaries tracked via `scope_starts` indices
- Linear search from innermost scope outward (fast for small local variable sets)
- Cache-friendly memory layout

### Expression Cache

Loop bodies avoid re-parsing expressions on each iteration:

- First iteration compiles the expression to AST
- AST stored in `ExpressionCache` keyed by source text
- Subsequent iterations reuse compiled AST
- Supports arithmetic, comparisons, function calls, and variable references

### Compile-time Variable Binding

When possible, AST variable nodes are bound to slot indices at compile time:

- Built-in constants (`pi`, `e`) folded directly into AST nodes
- Local variables resolved to slot indices for O(1) access
- Dynamic lookup fallback for variables that cannot be statically bound

### Stack Frame Allocation

Function calls use pre-allocated frame mechanism:

- Parameters stored directly in flat scope slots
- No per-call `std::map` allocation
- Recursive calls efficiently manage scope push/pop

## Best Re-entry Points

For future work, the fastest way to rebuild context is:

1. `README.md`
2. `ARCHITECTURE.md`
3. `test/main.cpp` and `test/suites/*.cpp`
4. `src/core/calculator_commands.cpp`
5. `src/script/script_runtime.cpp`
6. `src/symbolic/node_parser.cpp`
7. `src/math/mymath.cpp`
