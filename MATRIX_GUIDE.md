# Matrix Guide

## Creation

- `[a, b; c, d]`
  Create a matrix literal, using `,` to separate columns and `;` to separate rows
- `vec(a, b, c)`
  Create a row vector
- `mat(rows, cols, values...)`
  Create a matrix in row-major order
- `zeros(rows, cols)`
  Create a zero matrix
- `eye(n)`
- `identity(n)`
  Create an identity matrix

Examples:

- `[1, 2; 3, 4]` -> `[[1, 2], [3, 4]]`
- `[1, 2; 3]` -> `[[1, 2], [3, 0]]`
- `vec(1, 2, 3)` -> `[1, 2, 3]`
- `mat(2, 2, 1, 2, 3, 4)` -> `[[1, 2], [3, 4]]`

## Resizing And Editing

- `resize(m, rows, cols)`
  Resize and zero-fill any new cells
- `append_row(m, values...)`
  Append one row
- `append_col(m, values...)`
  Append one column
- `transpose(m)`
  Return the transposed matrix
- `inverse(m)`
  Return the inverse matrix
- `get(m, row, col)`
  Read a matrix element using zero-based indices
- `set(m, row, col, value)`
  Return a modified matrix

Vector shortcuts:

- `get(v, index)`
- `set(v, index, value)`

Examples:

- `get(mat(2, 2, 1, 2, 3, 4), 1, 0)` -> `3`
- `set(vec(5, 6, 7), 1, 42)` -> `[5, 42, 7]`
- `set([1, 2; 3], 1, 2, 9)` -> `[[1, 2, 0], [3, 0, 9]]`

Notes:

- matrix literals pad missing elements with `0`
- `append_row(...)` and `append_col(...)` pad missing values with `0`
- if an appended row or column is longer than the current shape, the matrix expands and fills new cells with `0`
- `set(...)` expands the matrix with `0` fill if the target index is outside the current shape

## Operations

- `A + B`
- `A - B`
- `A * B`
- `A / scalar`
- `A ^ n`
- matrix-scalar mixing is supported for `+ - * /`
- `dot(a, b)`
- `outer(a, b)`

Notes:

- matrix addition and subtraction require the same shape
- matrix multiplication requires `lhs.cols == rhs.rows`
- matrix powers require a square matrix and an integer exponent
- negative integer powers are supported for invertible matrices, for example `A ^ -1`
- division by a matrix is not supported

Examples:

- `mat(2, 2, 1, 2, 3, 4) + eye(2)`
- `2 * mat(2, 2, 1, 2, 3, 4)`
- `mat(2, 3, 1, 2, 3, 4, 5, 6) * mat(3, 1, 7, 8, 9)`
- `mat(2, 2, 1, 2, 3, 4) ^ -1`
- `dot(vec(1, 2, 3), vec(4, 5, 6))`
- `outer(vec(1, 2), vec(3, 4, 5))`

## Analysis

- `norm(m)`
  Vector 2-norm or matrix Frobenius norm
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
- `lu_l(A)`
- `lu_u(A)`
- `svd_u(A)`
- `svd_s(A)`
- `svd_vt(A)`
- `solve(A, b)`

Notes:

- `trace`, `det`, `eigvals`, and `eigvecs` require square matrices
- `eigvals` and `eigvecs` are implemented for real-valued results; matrices with only complex eigenvalues are rejected
- `qr_q`, `qr_r`, `lu_l`, and `lu_u` currently require square matrices
- `lu_l` and `lu_u` use LU decomposition without pivoting, so leading pivots must stay non-zero
- `svd_u`, `svd_s`, and `svd_vt` return the reduced SVD factors
- `solve(A, b)` requires a square coefficient matrix and a vector right-hand side, and it returns the solution vector
- `least_squares(A, b)` currently uses the normal equations and expects a vector right-hand side

Examples:

- `norm(vec(3, 4))` -> `5`
- `det(mat(2, 2, 1, 2, 3, 4))` -> `-2`
- `det(mat(3, 3, 1, 2, 3, 0, 1, 4, 5, 6, 0))` -> `1`
- `rref(mat(2, 3, 1, 2, 3, 2, 4, 6))` -> `[[1, 2, 3], [0, 0, 0]]`
- `rref(mat(3, 4, 1, 2, -1, -4, 2, 3, -1, -11, -2, 0, -3, 22))` -> `[[1, 0, 0, -8], [0, 1, 0, 1], [0, 0, 1, -2]]`
- `eigvals(mat(2, 2, 2, 0, 0, 3))` -> `[3, 2]`
- `null(mat(2, 3, 1, 2, 3, 2, 4, 6))` -> `[[-2, -3], [1, 0], [0, 1]]`
- `least_squares(mat(2, 1, 1, 1), vec(2, 4))` -> `[3]`
- `qr_q(mat(2, 2, 2, 0, 0, 3))` -> `[[1, 0], [0, 1]]`
- `lu_l(mat(2, 2, 4, 3, 6, 3))` -> `[[1, 0], [1.5, 1]]`
- `lu_u(mat(2, 2, 4, 3, 6, 3))` -> `[[4, 3], [0, -1.5]]`
- `svd_s(mat(3, 2, 3, 0, 0, 2, 0, 0))` -> `[[3, 0], [0, 2]]`
- `solve(mat(2, 2, 2, 1, 5, 3), vec(1, 2))` -> `[[1], [-1]]`
- `solve(mat(3, 3, 3, 2, -1, 2, -2, 4, -1, 0.5, -1), vec(1, -2, 0))` -> `[[1], [-2], [-2]]`

## Variables

Matrix expressions can be assigned to variables:

- `m = mat(2, 2, 1, 2, 3, 4)`
- `m = set(m, 1, 0, 8)`
- `n = m + eye(2)`

Use `:vars` to inspect stored values.

## Persistence Note

Matrix variables are included in `:save` and `:load` state files alongside
scalar, string, custom-function, and script-function values.
