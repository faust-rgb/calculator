# Testing Guide

## Run Tests

```bash
make test
```

The test binary is `calculator_tests`.

The automated regression source is `test/tests.cpp`.

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
- Save/load of persisted calculator state
- Matrix creation with `vec`, `mat`, `zeros`, and `eye`
- Matrix reshaping and editing with `resize`, `append_row`, `append_col`, `transpose`, `get`, and `set`
- Matrix arithmetic including matrix/matrix, matrix/scalar, explicit inverse, positive powers, negative powers for invertible matrices, dot products, and outer products
- Matrix analysis including `norm`, `trace`, `det`, `rank`, `rref`, `eigvals`, `eigvecs`, `null`, `least_squares`, `qr_q`, `qr_r`, `svd_u`, `svd_s`, `svd_vt`, and `solve`
- Matrix-specific error paths such as shape mismatch, non-square analysis, singular inverse attempts, invalid indices, and unsupported complex eigenvalue cases
- Constants `pi` and `e`
- Error paths such as division by zero, invalid functions, and invalid domains

## Expected Result

`make test` should report zero failed cases and exit with status code `0`.

Current expected summary:

- `Failed: 0`

## Example Scripts

Runnable example inputs are stored in `test/script/`:

- `test/script/test_variables.calc`
  Multi-line non-interactive REPL-style input
- `test/script/test_function.calc`
  Script functions and return values
- `test/script/test_control_flow.calc`
  Loops and conditional control flow
- `test/script/test_matrix_basic.calc`
  Matrix-oriented redirected-input workflow examples
- `test/script/SYNTAX_GUIDE.md`
  Script syntax notes and redirected-stdin behavior

You can run them with:

```bash
./calculator < test/script/test_variables.calc
./calculator < test/script/test_function.calc
./calculator < test/script/test_matrix_basic.calc
```
