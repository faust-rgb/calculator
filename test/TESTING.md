# Testing Guide

## Run Tests

```bash
make test
```

The test binaries are `bin/calculator_tests` and `bin/planning_tests`.

The automated regression source is `test/tests.cpp`.

Use `make script-test` for the CLI script validation, or `make check` to run
both the C++ regression suite and the script validation.

## What Is Covered

- Operator precedence and parentheses
- Power associativity
- Unary operators
- Exponential, logarithmic, square root, cube root, general root, trigonometric, and inverse trigonometric functions
- Exact fraction mode and simplest-form rational output
- `gcd(a, b)`, `lcm(a, b)`, and `mod(a, b)` evaluation
- `abs(x)`, `sign(x)`, `floor(x)`, and `ceil(x)` evaluation
- `min(a, b)` and `max(a, b)` evaluation
- `pow(a, b)` function-style exponentiation and `factor(n)` factorization
- Base conversion functions `bin`, `oct`, `hex`, and `base`
- Prefixed integer literal parsing for `0b`, `0o`, and `0x`
- Bitwise functions `and`, `or`, `xor`, `not`, `shl`, and `shr`
- Variable assignment and variable reuse in both decimal and exact modes
- Symbolic constants mode for scalar expressions using `pi` and `e`
- Variable listing and clearing commands
- Session history command
- Save/load of persisted calculator state, including matrix variables
- Matrix creation with `vec`, `mat`, `zeros`, and `eye`
- Matrix reshaping and editing with `resize`, `append_row`, `append_col`, `transpose`, `get`, and `set`
- Matrix arithmetic including matrix/matrix, matrix/scalar, explicit inverse, positive powers, negative powers for invertible matrices, dot products, and outer products
- Matrix analysis including `norm`, `trace`, `det`, `rank`, `rref`, `eigvals`, `eigvecs`, `null`, `least_squares`, `qr_q`, `qr_r`, `svd_u`, `svd_s`, `svd_vt`, and `solve`
- Matrix-specific error paths such as shape mismatch, non-square analysis, singular inverse attempts, invalid indices, and unsupported complex eigenvalue cases
- Symbolic calculus including product and chain-rule differentiation, mixed partial derivatives, rule-based integration, gradients, Jacobians, Hessians, and critical-point classification
- Limits including removable singularities, one-sided limits, trigonometric/exponential cancellation, and non-existent two-sided limit errors
- Equation solving including Newton solve, bisection, secant, fixed-point iteration, polynomial roots, and matrix linear-system solving
- Constants `pi` and `e`
- Error paths such as division by zero, invalid functions, and invalid domains

## Expected Result

`make test` should report zero failed cases and exit with status code `0`.

Current expected summary:

- `Passed: 792`
- `Failed: 0`
- `Planning passed: 6`
- `Planning failed: 0`

## Example Scripts

Runnable example inputs are stored in `test/script/`:

- `test/script/comprehensive_validation.calc`
  Broad script and calculator feature validation
- `test/script/SYNTAX_GUIDE.md`
  The dedicated script syntax guide

You can run them with:

```bash
bin/calculator test/script/comprehensive_validation.calc
printf ':run test/script/comprehensive_validation.calc\n' | bin/calculator
```
