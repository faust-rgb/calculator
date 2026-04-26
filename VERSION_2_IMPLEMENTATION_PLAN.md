# Version 2.0 Implementation Plan

## Background

The current repository represents the 1.0 architecture. It has grown into a
capable command-line calculator, but the earliest core design centered on
`double`, `long long` rationals, and several later side paths for exact,
precise-decimal, symbolic, matrix, and complex behavior.

That design is now the main blocker for 2.0 goals:

- arbitrary-precision calculations across the whole system
- first-class complex-number support
- a real CAS core instead of a display-oriented symbolic layer
- complete symbolic differentiation
- broad and extensible symbolic integration
- the same user-facing function coverage as 1.0

Version 2.0 should therefore be implemented as a replacement of the numeric,
expression, and CAS cores while preserving the 1.0 command surface wherever
possible.

## Implementation Constraints

2.0 must remain self-contained. Do not use external mathematics libraries or
standard-library math implementations as the calculator's computational backend.

Forbidden as core implementations:

- Boost.Multiprecision or other Boost math/numeric components
- GMP, MPFR, MPC, FLINT, Arb, Eigen, SymEngine, GiNaC, or similar external math
  libraries
- `<cmath>` or `math.h` implementations for calculator math functions
- `std::complex` as the calculator complex-number implementation
- delegating CAS behavior to external systems such as SymPy or Mathematica

Allowed uses:

- the C++ standard library for containers, strings, algorithms, I/O, memory
  management, and general utilities
- compiler-provided primitive integer and floating-point types as internal
  building blocks, test references, or fast approximate fallback paths
- local test tooling that compares committed golden data against external tools,
  as long as the released calculator does not depend on those tools

All major mathematical functionality must be implemented in project code,
including:

- arbitrary-precision integer arithmetic
- arbitrary-precision rational arithmetic
- arbitrary-precision decimal/real arithmetic
- complex arithmetic and complex elementary functions
- elementary functions such as `sqrt`, `exp`, `ln`, `sin`, `cos`, `tan`, inverse
  trigonometric functions, and hyperbolic functions
- common special functions such as `gamma`, `beta`, `zeta`, `erf`, `erfc`, and
  `bessel`
- matrix algorithms
- symbolic simplification, differentiation, and integration

## Goals

2.0 must support:

- all functions, commands, matrix helpers, script features, and analysis
  workflows currently supported by 1.0
- arbitrary-precision integer, rational, decimal, and real calculations
- complex input, output, arithmetic, elementary functions, and matrix values
- a unified value system for numbers, expressions, matrices, strings, functions,
  and unevaluated symbolic forms
- a unified expression parser shared by numeric, symbolic, matrix, and script
  evaluation
- a real CAS expression tree with exact numeric nodes
- canonical simplification
- complete symbolic differentiation for supported expression forms
- an extensible symbolic integration framework
- graceful unevaluated symbolic output when an integral cannot be expressed in
  the currently supported closed forms

## Non-Goals

2.0 should not try to make every mathematical expression produce a closed-form
integral. A practical CAS must be able to preserve unresolved forms such as:

```text
integral(exp(x ^ 2), x)
```

The correct target is an extensible integration engine that covers common
elementary, rational, trigonometric, exponential, logarithmic, and special
function patterns, while returning a valid unevaluated `Integral` expression
when no rule applies.

## Current 1.0 Architecture Issues

- `DecimalParser` returns `double`, so the main evaluation path is fixed to
  binary floating-point.
- `Rational` uses `long long`, so exact arithmetic overflows quickly.
- `PreciseDecimal` is a local string-decimal path rather than a unified numeric
  type.
- `matrix::Matrix` stores `std::vector<double>`, so matrices cannot naturally
  contain high-precision, complex, or symbolic values.
- `StoredValue` uses several boolean flags to combine unrelated value kinds.
  This will become brittle as 2.0 adds more value categories.
- `SymbolicExpression::number(double)` and polynomial coefficients based on
  `std::vector<double>` keep the symbolic layer tied to approximate values.
- exact mode, symbolic constants mode, and precise decimal mode are side paths
  rather than parts of a single evaluation model.
- numeric functions in `mymath` expose `double` APIs and cannot be reused for
  arbitrary precision or complex values without significant redesign.

## Target Module Layout

Recommended new or reorganized modules:

```text
src/numeric/
  bigint.h/.cpp
  rational.h/.cpp
  decimal.h/.cpp
  complex.h/.cpp
  number.h/.cpp
  precision_context.h/.cpp

src/expression/
  expr.h/.cpp
  expr_node.h
  parser.h/.cpp
  printer.h/.cpp
  evaluator.h/.cpp
  pattern.h/.cpp

src/cas/
  simplify.h/.cpp
  differentiate.h/.cpp
  integrate.h/.cpp
  polynomial.h/.cpp
  rational_function.h/.cpp
  assumptions.h/.cpp
  transforms.h/.cpp

src/runtime/
  value.h/.cpp
  environment.h/.cpp
  function_registry.h/.cpp
  command_registry.h/.cpp

src/matrix/
  matrix_value.h/.cpp
  matrix_algorithms.h/.cpp

src/compat/
  legacy_commands.h/.cpp
  migration.h/.cpp
```

The current CLI, command handling, documentation, and tests can be migrated
gradually onto this new core.

## Core Value Model

Introduce first-class runtime values:

```cpp
class BigInt;
class Rational;
class BigDecimal;
class Complex;
class Number;
class Expr;
class MatrixValue;
class Value;
```

Recommended runtime value representation:

```cpp
using Value = std::variant<
    Number,
    Expr,
    MatrixValue,
    std::string,
    FunctionValue,
    UnevaluatedValue
>;
```

Recommended numeric kind model:

```cpp
enum class NumberKind {
    Integer,
    Rational,
    Decimal,
    RealApprox,
    Complex
};
```

Complex numbers should be represented with high-precision-capable components:

```cpp
struct Complex {
    Number real;
    Number imag;
};
```

Do not use `std::complex<double>` as the calculator complex type. Complex
arithmetic, formatting, promotion, branch cuts, and elementary functions must be
implemented in project code.

## Arbitrary-Precision Numeric System

### Requirements

- arbitrary-size integers
- arbitrary-size rationals
- arbitrary-precision decimal/real approximations
- configurable precision and rounding
- exact-first promotion rules
- deterministic display formatting

### Implementation Tasks

1. Implement `BigInt`.
   - This must be implemented in project code.
   - Use a limb representation such as base `1e9`, base `2^32`, or base
     `2^64`.
   - Start with schoolbook arithmetic, then add Karatsuba/Toom-Cook only if
     benchmarks justify the added complexity.
   - Support arithmetic, comparison, gcd, power, modular operations, and bitwise
     operations.

2. Implement `Rational`.
   - Use `BigInt` numerator and denominator.
   - Normalize after construction and arithmetic.
   - Support arithmetic, comparison, integer checks, floor/ceil/round, power,
     and conversion to decimal.

3. Implement `BigDecimal` or arbitrary-precision real.
   - Add a `PrecisionContext`.
   - Support configurable digit precision and rounding.
   - Support elementary functions through custom iterative algorithms, argument
     reduction, series expansions, Newton iterations, AGM-style methods where
     appropriate, and carefully bounded convergence checks.

4. Implement promotion rules.
   - `Integer + Rational -> Rational`
   - `Rational + Decimal -> Decimal`
   - `Real + Complex -> Complex`
   - exact values should remain exact unless an approximate function requires
     approximation.

5. Replace public `double`-based evaluation paths with `Number` or `Value`.

### Precision Context

```cpp
struct PrecisionContext {
    int digits;
    RoundingMode rounding;
    int max_iterations;
    bool exact_preferred;
};
```

User commands:

```text
:precision 100
:rounding nearest
:rounding floor
:rounding ceil
:rounding zero
```

### Numeric Acceptance Examples

```text
1 / 3
=> 1/3

:precision 80
sqrt(2)
=> 1.41421356237309504880168872420969807856967187537694807317667973799...

factorial(100)
=> exact integer

nCr(200, 100)
=> exact integer
```

## Complex Number System

Complex values should be part of the normal numeric tower.

### Requirements

- built-in constant `i`
- parser support for imaginary literals such as `2i`, `3 + 4i`, `-i`
- exact or high-precision real and imaginary components
- complex arithmetic
- complex elementary functions
- complex-aware matrix operations
- branch-cut behavior documented for `ln`, `sqrt`, `pow`, and inverse
  trigonometric functions

### Functions To Support

Initial complex function coverage:

- `sqrt`
- `exp`
- `ln`
- `pow`
- `sin`, `cos`, `tan`
- `sinh`, `cosh`, `tanh`
- `asin`, `acos`, `atan`
- `real`, `imag`, `abs`, `arg`, `conj`, `polar`

Special-function coverage can be staged:

- phase 1: `gamma`, `beta`, `zeta`
- phase 2: `erf`, `erfc`, `bessel`

### Complex Acceptance Examples

```text
sqrt(-1)  => i
i ^ 2     => -1
exp(i*pi) => -1
ln(-1)    => i*pi
sin(i)    => 1.175201...i
```

## Expression And CAS Core

`SymbolicExpression` should become part of the primary expression model instead
of a side subsystem.

### Expression Node Kinds

Recommended canonical node kinds:

```cpp
enum class ExprKind {
    Number,
    Symbol,
    Add,
    Mul,
    Pow,
    Function,
    Matrix,
    Integral,
    Derivative,
    Equation,
    Piecewise,
    List
};
```

Prefer canonical forms:

```text
a - b  -> a + (-1) * b
a / b  -> a * b ^ -1
```

This simplifies pattern matching, simplification, differentiation, and
integration.

### CAS Requirements

- exact numeric nodes
- complex expression nodes
- structural equality and hashing
- canonical ordering
- expression interning or hash-consing
- pattern matching
- assumption context
- substitution
- free-variable analysis
- expression printing independent from parsing

### Assumptions

Add a first-class assumption system:

```text
assume(x > 0)
assume(n integer)
assume(a real)
```

Domain-aware simplifications must depend on assumptions. For example:

```text
sqrt(x ^ 2) => x
```

is only valid when `x >= 0` is known.

## Simplification Plan

Simplification should run as a staged canonical pipeline.

### Stage 1: Local Algebraic Rules

```text
x + 0 => x
x * 1 => x
x * 0 => 0
x ^ 1 => x
x ^ 0 => 1
```

### Stage 2: Exact Constant Folding

```text
2 + 3 => 5
2/3 + 5/7 => 29/21
```

### Stage 3: Add/Mul Canonicalization

```text
y + x => x + y
2*x + 3*x => 5*x
x*x => x^2
```

### Stage 4: Power Rules

```text
x^a * x^b => x^(a+b)
(x^a)^b => x^(a*b)
```

Power rules must be guarded by domain assumptions where needed.

### Stage 5: Polynomial And Rational Simplification

- expand
- factor
- collect
- cancel common factors
- normalize rational functions

### Stage 6: Function Identities

Examples:

```text
sin(x)^2 + cos(x)^2 => 1
exp(ln(x)) => x
ln(exp(x)) => x
```

Again, some identities require assumptions.

## Symbolic Differentiation

Symbolic differentiation should be complete for all supported expression forms.

### Required Coverage

- constants
- variables
- addition and multiplication
- powers
- chain rule
- product rule
- quotient behavior through canonical multiplication by inverse powers
- elementary functions
- reciprocal trigonometric functions
- inverse trigonometric functions
- hyperbolic functions
- logarithmic and exponential functions
- `abs`, `sign`, and piecewise-aware derivatives where possible
- known special-function derivative rules
- unknown functions represented as `Derivative(...)`
- higher-order derivatives
- partial derivatives
- gradient, Jacobian, and Hessian

### Differentiation Acceptance Examples

```text
diff(sin(x^2), x)
=> 2*x*cos(x^2)

diff(x^x, x)
=> x^x*(ln(x)+1)

diff(ln(sin(x)), x)
=> cot(x)
```

## Symbolic Integration

Integration should be implemented as a layered solver.

### Integration Pipeline

```text
integrate(expr, x)
  -> normalize
  -> constant and linearity rules
  -> table rules
  -> polynomial integrator
  -> rational function integrator
  -> substitution detector
  -> integration by parts detector
  -> trigonometric integrator
  -> exponential/logarithmic integrator
  -> special-function integrator
  -> heurisch or Risch-lite layer
  -> return Integral(expr, x)
```

### Phase 1 Coverage

- constants
- powers
- `1/x`
- polynomials
- basic trigonometric functions
- basic reciprocal trigonometric forms
- exponentials of linear arguments
- logarithms
- simple chain-rule substitutions
- common integration-by-parts forms

Examples:

```text
integrate(x^3 + 2*x, x)
=> x^4/4 + x^2

integrate(2*x*cos(x^2), x)
=> sin(x^2)

integrate(x*exp(x), x)
=> exp(x)*(x - 1)
```

### Phase 2 Coverage

- rational-function integration
- polynomial division
- square-free decomposition
- partial fractions
- real and complex factorization
- trigonometric rational forms
- Weierstrass substitution with `t = tan(x/2)`
- polynomial times `exp`, `sin`, or `cos`

### Phase 3 Coverage

- heurisch-style integration
- Risch-lite extension framework
- special functions:
  - `erf`
  - `erfi`
  - `Si`
  - `Ci`
  - `Ei`

Examples:

```text
integrate(1/(x^2 + 1), x)
=> atan(x)

integrate(exp(-x^2), x)
=> sqrt(pi)/2 * erf(x)

integrate(exp(x^2), x)
=> sqrt(pi)/2 * erfi(x)
```

When unsupported:

```text
integrate(f(x), x)
=> integral(f(x), x)
```

## Parser Rewrite

2.0 should replace the current split parser model:

- `DecimalParser`
- exact parser
- precise decimal parser
- matrix expression parser
- symbolic parser

with one expression parser:

```cpp
Expr parse_expression(std::string_view input);
Value evaluate(const Expr& expr, Environment& env, EvalOptions options);
```

### Parser Requirements

- numeric literals
- scientific notation
- arbitrary-precision integers and decimals
- rational forms
- imaginary literals
- matrices
- function calls
- assignments
- function definitions
- equations
- lists
- unevaluated symbolic forms

The parser should build expressions. It should not directly perform `double`
evaluation.

## Function Registry

Function behavior should be centralized.

```cpp
struct FunctionSpec {
    std::string name;
    Arity arity;
    FunctionCategory category;
    EvalHandler numeric_eval;
    EvalHandler symbolic_eval;
    DiffRule diff_rule;
    IntegrationHints integration_hints;
};
```

Each function should be registered once with:

- name
- aliases
- arity
- numeric evaluator
- complex evaluator
- symbolic evaluator
- derivative rule
- integration hints
- help text
- autocomplete metadata

This registry should drive help text, autocomplete, documentation generation,
and evaluation dispatch.

## Runtime And Environment

Replace `StoredValue` boolean-flag composition with typed bindings.

```cpp
struct Binding {
    Value value;
    std::optional<Expr> original_expr;
    SourceLocation defined_at;
};

class Environment {
    std::map<std::string, Binding> variables;
    std::map<std::string, FunctionDef> functions;
    PrecisionContext precision;
    AssumptionContext assumptions;
};
```

Supported examples:

```text
x = 1/3
z = 3 + 4i
f(x) = sin(x)^2 + cos(x)^2
A = mat(2, 2, 1, 2, 3, 4)
assume(x > 0)
```

## Matrix Upgrade

Current matrices store only `double`. 2.0 matrices should support exact,
high-precision, complex, and symbolic elements.

Recommended representation:

```cpp
class MatrixValue {
    std::size_t rows;
    std::size_t cols;
    std::vector<Value> data;
};
```

Alternatively, use a templated matrix for numeric fast paths while exposing
`MatrixValue` at runtime.

### Matrix Requirements

- basic matrix arithmetic over `Value`
- exact rational determinant
- exact rational inverse
- exact rref
- complex matrix arithmetic
- symbolic matrix addition, multiplication, determinant, and small inverse
- approximate high-precision QR/SVD/eigen workflows
- clear errors or explicit approximation prompts when symbolic input is not
  supported by a numeric-only algorithm

### Matrix Acceptance Examples

```text
det(mat(2, 2, 1/3, 1/2, 2/3, 3/4))
=> exact rational

inverse(mat(2, 2, 1, 2, 3, 4))
=> exact rational matrix

mat(2, 2, 1, i, -i, 1) * mat(2, 1, 1, i)
=> complex matrix result
```

## 1.0 Compatibility

2.0 should preserve the 1.0 user-facing function surface.

Compatibility areas:

- basic arithmetic
- scientific functions
- integer and programmer utilities
- rational helpers
- aggregate and statistics helpers
- probability helpers
- unit conversion helpers
- base conversion helpers
- bitwise operations
- matrix creation, editing, arithmetic, and linear algebra
- complex helper functions
- symbolic constants mode
- exact mode
- variables
- custom functions
- history, help, save, and load
- script execution
- ODE solving
- root finding
- one-variable and multi-variable integration
- polynomial helpers
- LP, ILP, MILP, and binary planning helpers

Legacy commands should map onto the new runtime model:

```text
:exact on
```

should become equivalent to setting exact-preferred evaluation.

```text
:symbolic on
```

should become equivalent to symbolic-preferred output or evaluation mode.

## Numeric Function Migration

`src/math/mymath.*` should be split into backend-independent algorithms where
possible.

Recommended layers:

```text
numeric elementary:
  add/sub/mul/div/pow/sqrt/exp/ln/sin

arbitrary-precision real backend:
  BigDecimal algorithms

complex backend:
  complex functions based on real backend

legacy double backend:
  retained only for tests, comparison, and carefully isolated fast approximate
  paths; it must not be the authoritative implementation for 2.0 high-precision
  or complex behavior
```

Migration order:

1. `sqrt`, `exp`, `ln`, `pow`
2. `sin`, `cos`, `tan`
3. inverse trigonometric and hyperbolic functions
4. `gamma`, `beta`, `zeta`
5. `erf`, `erfc`, `bessel`
6. probability and statistics helpers

## Testing Plan

### 1.0 Regression

- preserve and run current `make test`
- preserve and run script validation
- turn README examples into golden tests
- create a compatibility matrix for all documented functions

### High-Precision Tests

Examples:

```text
sqrt(2) at 100 digits
exp(1) at 100 digits
pi at 100 digits
factorial(500)
nCr(500, 250)
```

### Exact Arithmetic Tests

```text
1/3 + 1/6 => 1/2
(2/3)^20 => exact rational
factorial(100) => exact integer
```

### Complex Tests

```text
i^2 => -1
sqrt(-4) => 2i
exp(i*pi) => -1
ln(-1) => i*pi
```

### CAS Tests

- simplification golden tests
- derivative golden tests
- integral golden tests
- `diff(integrate(f, x), x)` round-trip tests where valid
- tests for unsupported integrals returning unevaluated `Integral`

### Matrix Tests

- exact determinant
- exact inverse
- rational rref
- complex matrix multiplication
- symbolic determinant for small matrices
- high-precision numeric decompositions

### Property-Based Tests

Useful identities:

```text
(a + b) - b == a
a * inverse(a) == 1 when a != 0
diff(f + g, x) == diff(f, x) + diff(g, x)
transpose(transpose(A)) == A
```

### External Golden Generation

Use Python `mpmath` and `sympy` to generate optional golden data. The committed
test fixtures should not require Python at test runtime unless explicitly
documented.

## Migration Phases

### Phase 0: Freeze 1.0 Baseline

Tasks:

- run all current tests
- export the complete 1.0 function list
- create golden compatibility tests
- document current public behavior
- avoid adding major new features to the 1.0 core

Deliverables:

```text
docs/2.0_COMPATIBILITY_MATRIX.md
test/golden/v1/*.txt
```

### Phase 1: New Numeric Core

Tasks:

- implement or integrate `BigInt`
- implement arbitrary-precision `Rational`
- implement `BigDecimal`
- implement high-precision-capable `Complex`
- implement `Number`
- implement `PrecisionContext`
- add focused unit tests

Acceptance:

- exact large factorials
- exact rational arithmetic
- high-precision decimal arithmetic
- complex arithmetic

### Phase 2: Unified Value And Expression Core

Tasks:

- add `Value`
- add `Expr`
- add parser
- add evaluator
- add printer
- support arithmetic, variables, and constants

Acceptance:

```text
1 + 2*3
1/3 + 1/6
x = sqrt(2)
3 + 4i
```

### Phase 3: Function Registry And Basic Function Migration

Tasks:

- build `FunctionRegistry`
- migrate constants
- migrate arithmetic helpers
- migrate scientific functions
- migrate integer utilities
- migrate aggregate/statistics helpers
- migrate unit conversions

Acceptance:

- documented 1.0 scalar functions pass through the new evaluator
- high-precision versions of core functions pass
- complex versions of supported functions pass

### Phase 4: CAS Simplification And Differentiation

Tasks:

- implement canonical simplification
- implement polynomial collection
- implement symbolic differentiation
- migrate gradient, Jacobian, and Hessian
- remove `double` dependency from symbolic numbers

Acceptance:

```text
simplify(sin(x)^2 + cos(x)^2) => 1
diff(x^x, x) => x^x*(ln(x)+1)
hessian(x^2 + x*y + y^2, [x, y])
```

### Phase 5: CAS Integration

Tasks:

- implement table integrals
- implement linearity
- implement polynomial integration
- implement rational-function integration
- implement trigonometric integration
- implement substitution detection
- implement integration by parts
- implement special-function outputs
- implement unevaluated `Integral`

Acceptance:

```text
integrate(2*x*cos(x^2), x) => sin(x^2)
integrate(1/(x^2+1), x) => atan(x)
integrate(x*exp(x), x) => exp(x)*(x-1)
```

### Phase 6: Matrix System Upgrade

Tasks:

- add `MatrixValue`
- migrate matrix literals
- migrate matrix arithmetic
- migrate determinant, inverse, and rref
- support complex matrix values
- support selected symbolic matrix operations
- migrate numeric decompositions with high-precision support where feasible

Acceptance:

- exact rational determinant
- rational inverse
- complex matrix multiplication
- selected symbolic matrix workflows

### Phase 7: Commands, Scripts, And Persistence

Tasks:

- connect `Calculator::process_line` to the new runtime
- migrate `:vars`, `:clear`, `:save`, `:load`
- map `:exact` and `:symbolic` onto new evaluation settings
- add `:precision`
- migrate script runtime
- migrate custom functions
- update help and autocomplete generation

Acceptance:

- current script validation passes
- 1.0 state files can be loaded
- 2.0 state files use a versioned format

### Phase 8: Analysis Module Migration

Tasks:

- migrate root solving to high-precision values
- migrate numerical integration to precision-aware algorithms
- upgrade ODE solving from fixed-step RK4 to adaptive methods
- support complex polynomial roots
- keep LP, ILP, and MILP exact where appropriate
- make tolerances precision-aware

### Phase 9: Performance And Stability

Tasks:

- expression hash-consing
- simplification memoization
- BigInt small-integer optimization
- rational gcd optimization
- lazy decimal precision where possible
- dense numeric matrix fast paths
- benchmarks

### Phase 10: Documentation And Release

Tasks:

- update `README.md`
- update `FUNCTIONS_REFERENCE.md`
- update `COMMANDS_REFERENCE.md`
- add `CAS_GUIDE.md`
- add `PRECISION_GUIDE.md`
- add `COMPLEX_GUIDE.md`
- add `MIGRATION_1_TO_2.md`
- update `CHANGELOG.md`

## Suggested Milestones

```text
2.0-alpha.1
  BigInt/Rational/BigDecimal/Complex/Number complete

2.0-alpha.2
  new parser + Value + Expr + basic evaluation complete

2.0-alpha.3
  1.0 scalar function migration complete

2.0-alpha.4
  symbolic simplify + differentiation complete

2.0-alpha.5
  symbolic integration phase 1 complete

2.0-beta.1
  matrix Value migration complete

2.0-beta.2
  scripts, commands, persistence, and compatibility complete

2.0-rc.1
  full 1.0 regression and 2.0 feature regression pass

2.0
  documentation, performance checks, and release notes complete
```

## Main Risks

- Self-implemented arbitrary precision can consume a large amount of time.
  This is required by the project constraint, so keep the first implementation
  simple and heavily tested before adding faster multiplication/division
  algorithms.
- Avoid accidentally reintroducing external math behavior through `<cmath>`,
  Boost, GMP, MPFR, `std::complex`, or similar dependencies.
- Symbolic integration can grow without bound. Always support unevaluated
  `Integral` as a valid result.
- Matrix genericization will touch many files. Keep a numeric fast path while
  adding a generic runtime matrix type.
- Complex branch cuts must be documented and tested.
- Domain-aware simplification requires assumptions. Avoid unsafe
  transformations without them.
- 2.0 correctness improvements may change output compared with 1.0. Keep a
  legacy compatibility mode where practical.

## Recommended Development Strategy

Do not rewrite the existing 1.0 parser and matrix implementation in place at
the start.

Recommended approach:

1. Build the 2.0 numeric, expression, and CAS cores alongside the 1.0 core.
2. Add independent tests for the new core.
3. Gradually route `Calculator::process_line` through the new evaluator.
4. Keep legacy fallback paths while compatibility is incomplete.
5. Remove the old `double`-first paths only after full regression coverage is
   passing.
