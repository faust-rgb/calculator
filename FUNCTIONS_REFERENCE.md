# Functions Reference

## Arithmetic

- `a + b`
- `a - b`
- `a * b`
- `a / b`
- `a ^ b`
- `pow(a, b)`

## Exponential And Logarithmic

- `exp(x)`
- `ln(x)`
- `log10(x)`

## Roots

- `sqrt(x)`
- `cbrt(x)`
- `root(value, degree)`

Notes:

- `root(value, degree)` requires an integer degree
- negative values are allowed only for odd roots

## Trigonometric

- `sin(x)`
- `cos(x)`
- `tan(x)`
- `asin(x)`
- `acos(x)`
- `atan(x)`

Notes:

- all trigonometric functions use radians

## Numeric Helpers

- `abs(x)`
- `sign(x)`
- `floor(x)`
- `ceil(x)`
- `min(a, b)`
- `max(a, b)`

## Matrix Creation

- `vec(a, b, c, ...)`
- `mat(rows, cols, values...)`
- `zeros(rows, cols)`
- `eye(n)`
- `identity(n)`

## Matrix Shape And Editing

- `resize(m, rows, cols)`
- `append_row(m, values...)`
- `append_col(m, values...)`
- `transpose(m)`
- `inverse(m)`
- `get(m, row, col)`
- `set(m, row, col, value)`
- `get(v, index)`
- `set(v, index, value)`

Notes:

- matrix and vector indices are zero-based
- `get(v, index)` and `set(v, index, value)` only work on vectors

## Matrix Operations

- `A + B`
- `A - B`
- `A * B`
- `A / scalar`
- `A ^ n`
- `dot(a, b)`
- `outer(a, b)`

Notes:

- matrix addition and subtraction require the same shape
- matrix multiplication requires `lhs.cols == rhs.rows`
- matrix powers require a square matrix and an integer exponent
- negative integer powers are supported for invertible matrices
- division by a matrix is not supported

## Matrix Analysis

- `norm(m)`
- `trace(m)`
- `det(m)`
- `rank(m)`
- `rref(m)`
- `eigvals(m)`
- `eigvecs(m)`
- `null(m)`
- `least_squares(A, b)`
- `qr_q(A)`
- `qr_r(A)`
- `svd_u(A)`
- `svd_s(A)`
- `svd_vt(A)`
- `solve(A, b)`

Notes:

- `norm(m)` returns the Euclidean norm for vectors and the Frobenius norm for matrices
- `trace`, `det`, `eigvals`, and `eigvecs` require square matrices
- `eigvals` and `eigvecs` currently support real-valued results
- `dot(a, b)` and `outer(a, b)` require vector arguments
- `qr_q(A)` and `qr_r(A)` currently require square matrices
- `null(m)` returns a basis matrix whose columns span the nullspace
- `svd_u(A)`, `svd_s(A)`, and `svd_vt(A)` return the reduced SVD factors
- `solve(A, b)` solves `Ax = b` for square `A` and vector `b`

## Integer / Number Theory

- `gcd(a, b)`
- `lcm(a, b)`
- `mod(a, b)`
- `factor(n)`

Notes:

- `gcd`, `lcm`, and `mod` require integer arguments
- `factor(n)` returns a factorization string, not a plain numeric value

## Base Conversion

- `bin(n)`
- `oct(n)`
- `hex(n)`
- `base(n, b)`

Notes:

- these are display-oriented conversion functions
- arguments must be integers
- supported bases for `base(n, b)` are `2..16`

Examples:

- `bin(10)` -> `1010`
- `oct(83)` -> `123`
- `hex(255)` -> `FF`
- `base(-31, 16)` -> `-1F`

## Bitwise

- `and(a, b)`
- `or(a, b)`
- `xor(a, b)`
- `not(a)`
- `shl(a, n)`
- `shr(a, n)`

Notes:

- all bitwise functions require integer arguments
- shift counts must be non-negative

## Constants

- `pi`
- `e`

## Symbolic Constants Mode

When `:symbolic on` is enabled, scalar expressions involving `pi` or `e` are
displayed with those symbols instead of immediately formatting their decimal
approximations.

Examples:

- `pi / 2 + e` -> `pi / 2 + e`
- `x = pi / 2` followed by `x + 1` -> `pi / 2 + 1`
- `f(x) = x + pi`, then `f(e)` -> `e + pi`
- `sin(pi / 2)` -> `1`
- `cos(pi)` -> `-1`
- `tan(pi / 4)` -> `1`
- `ln(e)` -> `1`
- `exp(ln(pi))` -> `pi`

Notes:

- this is a display-oriented mode for scalar expressions
- common exact-looking identities for `pi` and `e` are simplified automatically
- unsupported symbolic forms fall back to the normal numeric result

## Decimal Display Notes

Decimal results are formatted from `double` values.

For readability:

- values with `abs(value) <= 1e-10` are displayed as `0`
- values with `abs(value - round(value)) <= 1e-10` are displayed as the nearest
  integer

Examples:

- `2.00000000005` -> `2`
- `2.0000000002` -> `2.0000000002`
- `0.00000000001` -> `0`

## Exact Mode Notes

Exact mode is enabled with `:exact on`.

In exact mode:

- rational expressions are preserved in simplest fractional form
- unsupported exact operations fall back to decimal display
- some functions still return exact integers/rationals when possible

Examples:

- `1/3 + 1/4` -> `7/12`
- `abs(-7/3)` -> `7/3`
- `floor(7/3)` -> `2`
- `min(7/3, 5/2)` -> `7/3`
