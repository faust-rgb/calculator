# Complex Validation Notes

This note records a manual validation pass against more complex expressions,
custom functions, matrix routines, and script features.

Validation input file:

- `test/script/complex_validation.calc`

Script-mode recursion check used separately:

- `print(fib(8));` with a recursive `fib(n)` script function

## Verified Cases

- exact rational expression
  `((1/3 + 2/5)^2) / (7/6 - 1/2)` -> expected `121/150`, actual `121/150`
- function identity + composition
  `f(x) = sin(x)^2 + cos(x)^2 + exp(ln(x + 1))`, `f(2)` -> expected `4`, actual `4`
- derivative of a composed function
  `g(x) = exp(-x) * sin(2 * x)`, `diff(g, 1.5)` -> expected about `-0.473282498921`, actual `-0.47328249854`
- definite integral
  `integral(g, 0, 2)` -> expected about `0.455868833837`, actual `0.455868833839`
- cubic roots
  `roots(x^3 - 6*x^2 + 11*x - 6)` -> expected `1, 2, 3`, actual `1, 2, 3`
- quartic roots
  `roots(x^4 - 5*x^2 + 4)` -> expected `-2, -1, 1, 2`, actual `-2, -1, 1, 2`
- eigenvalues
  `eigvals([[3, 1], [1, 3]])` -> expected `[4, 2]`, actual `[4, 2]`
- least squares
  `least_squares([[1, 1], [1, 2], [1, 3]], [1, 2, 2])` -> expected `[[0.666666666667], [0.5]]`, actual `[[0.666666666667], [0.5]]`
- inverse consistency
  `(A ^ -1) * A` for `A = [[1, 2], [3, 4]]` -> expected identity, actual `[[1, 0], [0, 1]]`
- recursive script function
  `fib(8)` -> expected `21`, actual `21`

## Observation

- `limit(u, 1)` for `u(x) = (x^2 - 1) / (x - 1)` now displays as `2`
  after near-integer display normalization
- current decimal display simplifies values within `1e-10` of `0` or the
  nearest integer
