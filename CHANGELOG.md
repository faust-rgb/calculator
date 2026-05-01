# Changelog

## Current State

The project has evolved from a minimal C++ hello-world style setup into a
feature-rich command-line calculator with exact rational mode, programmable
helpers, interactive terminal UX, and project documentation.

## Latest Script Engine Performance Optimizations

- **FlatScopeStack Implementation**:
  - Replaced `std::vector<std::map<std::string, StoredValue>>` with flat array storage
  - Variables stored in contiguous memory with O(n) linear search (fast for small scopes)
  - Cache-friendly memory layout improves loop iteration performance
- **Expression Cache**:
  - Loop bodies compile expressions to AST on first iteration
  - Subsequent iterations reuse cached AST, avoiding repeated parsing
  - Supports arithmetic, comparisons, function calls, and variable references
- **Compile-time Variable Binding**:
  - Built-in constants (`pi`, `e`) folded directly into AST nodes
  - Local variables can be bound to slot indices for O(1) access
  - Dynamic lookup fallback for variables that cannot be statically bound
- **Stack Frame Optimization**:
  - Function calls use pre-allocated frame mechanism
  - Parameters stored directly in flat scope slots
  - No per-call `std::map` allocation for local variables
- **New Files**:
  - `src/core/expression_ast.cpp` / `src/core/expression_ast.h`: Compiled expression AST
  - `src/core/variable_resolver.cpp` / `src/core/variable_resolver.h`: Scoped variable resolution
- **Documentation**:
  - Added `test/script/SYNTAX_GUIDE_CN.md`: Chinese script syntax guide
  - Updated `ARCHITECTURE.md` with script engine performance section

## Latest Vector Calculus Improvements

- **Added Advanced Differential Form and Vector Analysis Commands**:
  - `implicit_diff(F, y, x)`: Implicit differentiation using the formula dy/dx = -F_x / F_y
  - `param_deriv(x_t, y_t, t)`: Parametric curve derivatives (supports higher-order)
  - `directional(expr, vars..., direction...)`: Directional derivative with automatic normalization
  - `line_integral_scalar(f, curve, t, a, b)`: Scalar line integral ∫_C f(r) |r'| dt
  - `line_integral_vector(F, curve, t, a, b)`: Vector line integral (work) ∫_C F·r' dt
  - `surface_integral_scalar(f, surface, u, a, b, v, c, d)`: Scalar surface integral ∬_S f |r_u × r_v| du dv
  - `surface_integral_flux(F, surface, u, a, b, v, c, d)`: Flux integral ∬_S F·(r_u × r_v) du dv
- **Integration Engine Enhancements**:
  - Added `IntegrationEngine` class with strategy pattern for multi-method integration
  - Implemented LIATE rule for integration by parts
  - Added generic substitution with candidate collection
  - Added cyclic integration detection
- **Symbolic Polynomial Support**:
  - New `SymbolicPolynomial` class with symbolic coefficients
  - Square-free decomposition (Yun algorithm)
  - Coefficient identity method for partial fraction decomposition
  - Symbolic quadratic integration formulas

## Latest Series and Limit Improvements

- **Introduced Power Series Arithmetic (PSA) Engine**:
  - Implemented a dedicated engine for efficient Taylor series expansion of
    composite functions (e.g., `sin(exp(x))`, `exp(sin(x))`) using recursive
    coefficient formulas.
  - Significantly improved performance for high-degree expansions by avoiding
    symbolic derivative explosion.
- **Enhanced Limit Calculation**:
  - Integrated symbolic PSA for limit detection, allowing exact resolution of
    removable singularities like `sin(x)/x` or `(1-cos(x))/x^2`.
  - Added support for infinity limits via symbolic substitution ($x=1/t$) and
    Laurent-style series analysis.
  - Upgraded numerical fallback to a high-order Richardson extrapolation
    (up to 14th order), improving accuracy from $10^{-5}$ to $10^{-11}$ or better.
- **Advanced Series Summation**:
  - Added support for infinite series involving Riemann Zeta function values
    $\zeta(2k)$ (e.g., $\sum 1/n^2 = \pi^2/6$).
  - Implemented robust symbolic geometric series detection via term-ratio
    simplification.
- **Modular Test Architecture**:
  - Refactored the massive `test/tests.cpp` into a modular suite under
    `test/suites/` (`core`, `analysis`, `symbolic`).
  - Updated build system to automatically detect and compile new test suites.

## Latest Numerical Work

- Replaced one-variable definite integration with an adaptive Gauss-Kronrod
  G7-K15 path, including endpoint-singularity transforms for cases such as
  `integral(1 / sqrt(x), 0, 1)`
- Upgraded one-variable numerical differentiation to use adaptive central
  differences with four Richardson extrapolation layers
- Added symbolic derivative and integral coverage for reciprocal trigonometric
  functions including `sec`, `csc`, `cot`, `sec(x)^2`, `csc(x)^2`, and
  `sec(x)^2 * tan(x)`
- Added `round`, `trunc`, `clamp`, `log(x, base)`, `log2`, and `exp2`
- Added programmer helpers:
  - `rol`, `ror`
  - `popcount`, `bitlen`, `ctz`, `clz`, `parity`, `reverse_bits`
- Added `percentile(...)` and `quartile(...)` in both scalar-aggregate and vector forms
- Extracted calculator help generation into `src/core/calculator_help.cpp` so the
  main calculator implementation no longer carries the full help-text payload
- Added scientific-notation parsing for decimal, exact, and high-precision
  decimal paths
- Stabilized large-angle trigonometric reduction to avoid overflow-prone integer
  casts
- Switched `gamma`/`beta` critical paths to log-space evaluation and tightened
  special-function behavior for large arguments
- Reworked `bessel` evaluation to avoid the previous overflow/`nan` failures on
  large inputs
- Improved one-variable `limit(...)` sampling so removable singularities such as
  `sin(x)/x` and `(1-cos(x))/x^2` converge much more accurately
- Replaced the old rank-sensitive QR path with a Householder-based
  implementation and rebuilt reduced SVD/pseudo-inverse/least-squares support
  on top of it
- `pinv(...)`, `cond(...)`, `least_squares(...)`, `null(...)`, `svd_*`, and
  `eigvals(...)` now behave correctly on a wider set of singular and wide
  matrices
- Split oversized numeric implementation files by extracting
  `src/math/mymath_special_functions.cpp` and
  `src/matrix/matrix_linear_algebra.cpp` from the previous monolithic sources
- Expanded regression coverage for the reported numerical edge cases

## Major Additions Implemented

- Interactive command-line calculator REPL
- VS Code build/debug configuration cleanup
- Custom math implementation without relying on the standard math library for
  core functionality
- Expression parser with:
  - operator precedence
  - parentheses
  - unary operators
  - right-associative exponentiation
- Mathematical functions:
  - `sin`, `cos`, `tan`
  - `asin`, `acos`, `atan`
  - `sinh`, `cosh`, `tanh`
  - `exp`, `ln`, `log10`
  - `gamma`
  - `sqrt`, `cbrt`, `root`
  - `pow`
  - first-order ODE solving with `ode` and `ode_table`
  - multi-variable integration in Cartesian, cylindrical, and spherical coordinates
- Combinatorics helpers:
  - `factorial`
  - `nCr`
  - `nPr`
- Aggregate helpers:
  - `sum`
  - `avg`
  - `median`
  - `mean`
  - `mode`
  - `var`
  - `std`
  - `percentile`
  - `quartile`
- Unit conversion helpers:
  - `deg2rad`, `rad2deg`
  - `celsius`, `fahrenheit`, `kelvin`
- Result normalization for near-zero floating-point noise
- Terminal UX improvements:
  - up/down history navigation
  - `Tab` autocomplete
  - double-Tab completion candidate listing
  - context-aware completion for commands, help topics, variables, and custom functions
- Exact rational mode:
  - `:exact on`
  - `:exact off`
  - exact fraction display when possible
- Symbolic constants mode:
  - `:symbolic on`
  - `:symbolic off`
  - scalar symbolic display for expressions involving `pi` and `e`
  - symbolic propagation through scalar variable assignment and custom unary functions
  - common simplifications such as `sin(pi / 2)`, `sin(pi / 6)`, `sin(pi / 3)`,
    `cos(pi)`, `cos(pi / 3)`, `cos(pi / 6)`, `tan(pi / 4)`, `tan(pi / 6)`,
    `tan(pi / 3)`, `ln(e)`, and `exp(ln(pi))`
- Symbolic differentiation and integration improvements:
  - raw one-variable expressions such as `diff(x ^ 2)` and `integral(1 / x)` now infer the variable automatically
  - added symbolic derivative rules for `asin`, `acos`, `atan`, `sinh`, `cosh`, `tanh`, `abs`, and `cbrt`
  - added symbolic integral rules for `sqrt`, `cbrt`, `tan`, `1 / (a*x + b)`, linear powers `(a*x + b)^n`, and polynomial-times-`exp/sin/cos`
  - `integral(1 / x)` now returns `ln(abs(x)) + C`
  - `simplify(ln(exp(x)))` now reduces to `x`, while `exp(ln(x))` is only collapsed when the argument is known positive
  - multi-variable `simplify(...)` expressions continue to work without being mistaken for one-variable analysis input
  - replaced parse-based polynomial canonicalization with direct expression rebuilding to avoid recursive simplify stack overflows
- Display normalization for near-integer floating-point results
- Integer / utility functions:
  - `gcd`, `lcm`, `mod`
  - `abs`, `sign`, `floor`, `ceil`, `round`, `trunc`
  - `min`, `max`, `clamp`
  - `log`, `log2`, `exp2`
- Variable system:
  - assignments
  - variable reuse
  - variable listing/clearing
- Help system:
  - `help`
  - `:help`
  - topic help for commands/functions/examples
  - dedicated help for exact mode, variables, persistence, and programmer tools
- Persistence:
  - `:save file`
  - `:load file`
- Prime factorization:
  - `factor(n)`
- Base conversion:
  - `bin`, `oct`, `hex`, `base`
  - configurable `0x` prefix and uppercase/lowercase output
- Prefixed integer literal parsing:
  - `0b`
  - `0o`
  - `0x`
- Bitwise functions:
  - `and`, `or`, `xor`, `not`, `shl`, `shr`
  - `rol`, `ror`, `popcount`, `bitlen`, `ctz`, `clz`, `parity`, `reverse_bits`
- Expanded regression tests
- Handoff and architecture documentation

## Version 1.6 (2026-04-28)

- **ODE and Symbolic Calculus Enhancements:**
  - High-order ODE automatic reduction: `ode(y'' + y, x0, [y0, y'0], x1)` syntax now supported, automatically converting to first-order systems.
  - Symbolic ODE solver `dsolve(rhs, x_var, y_var)` for first-order linear ODEs and `y' = f(x)` forms.
  - Constrained optimization via Lagrange multipliers: `lagrange(f, [g1, g2...], vars...)`.
  - Field theory operators: `divergence`/`div`, `curl`, and `laplacian` with full CLI autocomplete support.
  - Identifier system extended to support derivative notation like `y'`, `y''` in ODE expressions.
- **Extended Trigonometric Integration:**
  - Added integral rules for `sec(x)^2`, `csc(x)^2`, `cot(x)^2`.
  - Added integral rules for products: `sec(x)*tan(x)`, `csc(x)*cot(x)`, `sec(x)^2*tan(x)`, `csc(x)^2*cot(x)`.
  - Implemented Weierstrass substitution for rational functions of `sin`/`cos`.
  - Extended chain-rule substitution to support `sec`, `csc`, `cot`, `ln`, `asin`, `acos`, `atan`, `sinh`, `cosh`, `tanh`.

## Version 1.5 (2026-04-27)

- **Performance Optimization:**
  - Improved symbolic node interning with an incremental LRU eviction strategy in `src/symbolic/node_parser.cpp`, eliminating $O(N)$ scan overhead during high-frequency node creation.
- **Terminal UX Enhancements:**
  - Upgraded REPL in `src/app/main.cpp` with a full-featured line editor:
    - Inline cursor movement (Left/Right arrow keys).
    - Quick navigation: Home (`Ctrl+A`), End (`Ctrl+E`).
    - Editing shortcuts: Delete character at cursor (`Ctrl+D`), Clear to end of line (`Ctrl+K`).
    - Precise terminal redrawing and cursor synchronization.
- **Symbolic Algebra Core Evolution:**
  - Introduced dedicated `NodeType` for exact constants `pi` and `e`, enabling more robust pattern matching and preventing premature numerical collapse.
  - Added `SymbolicExpression::common_subexpressions()` (CSE base) to identify and count repeated sub-trees within an expression.
  - Unified internal representation by normalizing `e ^ x` to `exp(x)` during simplification.
  - Updated derivative, integral, and transform engines to natively handle the new exact constant nodes.

## Current Quality Status

- `make test` currently passes all tests
- Current test summary:
  - `Passed: 867, Failed: 0`
  - `Planning passed: 6, Planning failed: 0`
- Added comprehensive symbolic-computing coverage for calculus, limits,
  matrices, equation solving, and planning helpers
