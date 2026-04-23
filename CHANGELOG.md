# Changelog

## Current State

The project has evolved from a minimal C++ hello-world style setup into a
feature-rich command-line calculator with exact rational mode, programmable
helpers, interactive terminal UX, and project documentation.

## Latest Numerical Work

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
  `src/math/mymath_special_functions.inc` and
  `src/matrix/matrix_linear_algebra.inc` from the previous monolithic sources
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

## Current Quality Status

- `make test` currently passes
- Expected test summary:
  - `Passed: 671, Failed: 0`
