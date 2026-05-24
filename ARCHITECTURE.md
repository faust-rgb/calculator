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
│                        app/ (Application)                   │
│                     main.cpp, CLI Entry                     │
├─────────────────────────────────────────────────────────────┤
│                       module/ (Module Layer)                │
│              CalculatorModule, Extension Points             │
├─────────────────────────────────────────────────────────────┤
│                      analysis/ (Analysis Layer)             │
│            Calculus, Limits, ODE, Optimization, etc.        │
├─────────────────────────────────────────────────────────────┤
│    ┌──────────┐  ┌──────────┐  ┌──────────┐  ┌──────────┐   │
│    │  plot/   │  │  stats/  │  │  dsp/    │  │ symbolic/│   │
│    │  Plotting│  │  Stats   │  │  Signal  │  │   CAS    │   │
│    └──────────┘  └──────────┘  └──────────┘  └──────────┘   │
├─────────────────────────────────────────────────────────────┤
│                       core/ (Core Layer)                    │
│     Calculator, Scope, Exceptions, Utils, Services          │
├─────────────────────────────────────────────────────────────┤
│                     execution/ (Execution Layer)            │
│    CommandRegistry, ScriptRuntime, VariableResolver         │
├─────────────────────────────────────────────────────────────┤
│                       parser/ (Parser Layer)                │
│   UnifiedParser, CommandParser, ScriptParser, AST           │
├─────────────────────────────────────────────────────────────┤
│                        io/ (IO Layer)                       │
│              Persistence and File Operations                │
├─────────────────────────────────────────────────────────────┤
│                       types/ (Type Layer)                   │
│          StoredValue, Function, Matrix, Scalar              │
├─────────────────────────────────────────────────────────────┤
│                        math/ (Math Layer)                   │
│              Base Functions, Numeric, Helpers               │
└─────────────────────────────────────────────────────────────┘
```

### Dependency Rules

1. **Upward dependency prohibited**: Lower layers cannot depend on upper layers
2. **Same-layer dependency**: Same-layer modules can depend on each other, but minimize
3. **Cross-layer dependency**: Can only depend downward, never upward

### Module Responsibilities

#### types/ (Type Layer)
- `scalar_type.h` - Core scalar type definitions
- `stored_value.h` - Stored value types (scalar, matrix, string, etc.)
- `function.h` - Function type definitions

#### math/ (Math Layer)
- `mymath.h` - Basic math functions and constants
- `helpers/` - Math helper tools (integer operations, base conversion)
- `types/float128.h` - High-precision math support (float128_t)

#### parser/ (Parser Layer)
- `unified_expression_parser.h` - Main entry point for all expression types
- `command_parser.h` - Calculator command parser
- `script_parser.h` - Script language parser
- `expression_ast.h` - Expression AST definition
- `token_types.h` - Lexical token definitions

#### execution/ (Execution Layer)
- `command_registry.h` - Command registration and lookup
- `script_runtime.h` - Script execution engine
- `variable_resolver.h` - Variable lookup with scoped resolution
- `builtin_constants.h` - Physical and mathematical constants

#### core/ (Core Layer)
- `calculator.h` - Public calculator API
- `calculator_impl.h` - Internal implementation state
- `scope.h` - Variable and function scope management
- `calculator_service_factory.h` - Dependency injection for services

#### analysis/ (Analysis Layer)
- `calculus/`, `integration/`, `optimization/`, `differential_equations/`, `rootfinding/`, `series/`

#### symbolic/ (Symbolic Layer)
- `core/`, `algebra/`, `calculus/`, `transformation/`, `solver/`

#### dsp/ (Signal Processing Layer)
- `fft.cpp`, `filter_design.cpp`, `signal_processing.h`

## Main Files

- `src/app/main.cpp`
  Terminal interaction, history, and command dispatch
- `src/core/api/calculator_core.cpp`
  Main calculator logic and service orchestration
- `src/parser/grammars/unified_expression_parser.cpp`
  Centralized parsing logic for decimal, exact, and symbolic expressions
- `src/parser/ast/expression_ast.cpp`
  Compiled expression AST for performance-sensitive evaluation
- `src/execution/engine/script_runtime.cpp`
  Script execution engine with scope optimization
- `src/execution/resolver/variable_resolver.cpp`
  Variable lookup and resolution strategy
- `src/symbolic/core/symbolic_expression.h`
  Symbolic expression representation
- `src/math/core/` and `src/math/transcendental/`
  Core numerical algorithms
- `src/matrix/matrix.cpp`
  Matrix storage and base operations
- `src/matrix/matrix_linear_algebra.cpp`
  Advanced linear algebra (SVD, QR, Eigen, etc.)

## Directory Layout

- `src/app`
  CLI entry point
- `src/core`
  Calculator runtime and service definitions
- `src/parser`
  All parsing logic (Expression, Command, Script)
- `src/execution`
  Runtime execution (Commands, Script, Variables)
- `src/math`
  Numeric functions and specialized math
- `src/matrix`
  Matrix representation and linear algebra
- `src/analysis`
  Advanced numerical analysis (split by domain)
- `src/symbolic`
  Symbolic expression engine (split by domain)
- `src/polynomial`
  Polynomial arithmetic and root-finding
- `src/statistics`
  Statistics and probability distributions
- `src/time`
  Time and benchmarking functions
- `src/io`
  Persistence and file operations
- `test`
  Regression suite and script examples

Large implementation areas are split into focused translation units with internal
headers for shared declarations. Current internal split headers include:

- `src/core/api/calculator_internal_types.h`
- `src/matrix/matrix_internal.h`
- `src/symbolic/core/symbolic_expression_internal.h`

## Execution Flow

User input goes through this rough path:

1. `src/app/main.cpp` reads a line from the terminal
2. expressions are passed to `Calculator::process_line(...)`
3. command-style inputs (including native commands like `plot(...)` and meta-commands like `:help`) are dispatched via `CommandRegistry` to the corresponding `CalculatorModule`.
4. The module executes the command using `execute_args_view`, resolving dependencies via `ServiceLocator`.
5. `Calculator` either:
   - treats the line as an assignment, or
   - evaluates it for display via the `UnifiedParser`.

## Parsing Model

Parsing is handled by the `UnifiedParser` in `src/parser/grammars/unified_expression_parser.cpp`.

### Decimal Path

Used for standard evaluation with `double` or `Scalar`.
Supports all ordinary math functions and variable lookup.

### Exact Path

Enabled via `:exact on`.
Handled within the unified parser.
Preserves expressions as `Rational` values; falls back to decimal if exactness
cannot be maintained.

### Symbolic Path

Used for CAS-style operations and enabled for scalar display via `:symbolic on`.
Rebuilds expressions using symbolic nodes (`src/symbolic/core/`).

## Stored Values

Calculator state lives in `Calculator::Impl`.
Variables are managed by `Scope` and resolved by `VariableResolver`.
`StoredValue` (`src/types/stored_value.h`) holds decimals, matrices, strings,
or symbolic expressions.

## Symbolic Constants Flow

Symbolic constants mode layers a symbolic render path on top of the numeric engine.
If successful, the symbolic form is displayed and can be stored for later reuse.

Implementation:
- `src/parser/grammars/` - Parser grammars for expression handling
- `src/symbolic/transformation/simplify/simplify.cpp`
- `src/symbolic/algebra/algebra_helpers.cpp`

## Terminal UX

`src/app/main.cpp` implements a custom REPL with raw terminal support:

- **Navigation:** Arrow keys and Emacs-style shortcuts.
- **History:** Up/Down arrow recall.
- **Autocomplete:** Tab-based completion for commands and functions.

## Persistence

Variable state is persisted via `:save` and `:load` commands, handled by
the `IoModule` and `Calculator::save_state/load_state`.

## Numeric Design Notes

- Near-zero normalization and integer-snapping for display.
- Exact mode for rational arithmetic.
- Bitwise functions require integer arguments.

## Script Engine Performance

The script engine in `src/execution/engine/script_runtime.cpp` uses several optimizations:

### FlatScopeStack

Uses contiguous memory for variable slots to improve cache locality and speed up
lookups compared to traditional map-based scopes.

### Expression Cache

In `src/parser/ast/`, loop bodies are compiled once to AST
and cached for subsequent iterations.

### Compile-time Variable Binding

Variables in AST nodes are bound to slot indices during compilation where possible
to achieve O(1) access during execution.

## Best Re-entry Points

1. `README.md`
2. `ARCHITECTURE.md`
3. `test/suites/`
4. `src/core/api/calculator_core.cpp`
5. `src/execution/engine/script_runtime.cpp`
6. `src/parser/grammars/unified_expression_parser.cpp`
