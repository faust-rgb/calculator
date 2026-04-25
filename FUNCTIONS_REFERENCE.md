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
- `exp2(x)`
- `ln(x)`
- `log(x, base)`
- `log2(x)`
- `log10(x)`
- `gamma(x)`
- `beta(a, b)`
- `zeta(s)`
- `erf(x)`
- `erfc(x)`
- `bessel(n, x)`

## Hyperbolic

- `sinh(x)`
- `cosh(x)`
- `tanh(x)`
- `asinh(x)`
- `acosh(x)`
- `atanh(x)`

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
- `sec(x)`
- `csc(x)`
- `cot(x)`
- `asin(x)`
- `acos(x)`
- `atan(x)`
- `asec(x)`
- `acsc(x)`
- `acot(x)`
- `deg(x)`
- `rad(x)`
- `sin_deg(x)`
- `cos_deg(x)`

Notes:

- all trigonometric functions use radians

## Numeric Helpers

- `abs(x)`
- `sign(x)`
- `step(x)`
- `delta(x)`
- `heaviside(x)`
- `impulse(x)`
- `floor(x)`
- `ceil(x)`
- `round(x)`
- `trunc(x)`
- `min(a, b)`
- `max(a, b)`
- `clamp(x, min, max)`
- `sum(a, b, c, ...)`
- `mean(a, b, c, ...)`
- `avg(a, b, c, ...)`
- `median(a, b, c, ...)`
- `mode(a, b, c, ...)`
- `percentile(p, a, b, c, ...)`
- `quartile(q, a, b, c, ...)`
- `var(a, b, c, ...)`
- `std(a, b, c, ...)`
- `factorial(n)`
- `nCr(n, r)`
- `binom(n, r)`
- `nPr(n, r)`
- `fib(n)`
- `is_prime(n)`
- `next_prime(n)`
- `rand()`
- `randn()`
- `randint(a, b)`
- `pdf_normal(x, mu, sigma)`
- `cdf_normal(x, mu, sigma)`
- `rat(x)`
- `rat(x, max_denominator)`

Notes:

- `step(x)` / `heaviside(x)` return `1` for `x >= 0`, otherwise `0`
- `delta(x)` / `impulse(x)` are engineering shorthands; in numeric mode they return `1` at `0` and `0` elsewhere

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
- `pinv(m)`
- `get(m, row, col)`
- `set(m, row, col, value)`
- `get(v, index)`
- `set(v, index, value)`
- `reshape(m, rows, cols)`
- `vec(m)`
- `diag(m)`
- `diag(v)`

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
- `kron(a, b)`
- `hadamard(a, b)`

Notes:

- matrix addition and subtraction require the same shape
- matrix multiplication requires `lhs.cols == rhs.rows`
- matrix powers require a square matrix and an integer exponent
- negative integer powers are supported for invertible matrices
- division by a matrix is not supported

## Matrix Analysis

- `norm(m)`
- `cond(m)`
- `trace(m)`
- `det(m)`
- `rank(m)`
- `rref(m)`
- `eigvals(m)`
- `eigvecs(m)`
- `eig(m)`
- `null(m)`
- `least_squares(A, b)`
- `qr_q(A)`
- `qr_r(A)`
- `svd_u(A)`
- `svd_s(A)`
- `svd_vt(A)`
- `svd(A)`
- `solve(A, b)`
- `cholesky(A)`
- `schur(A)`
- `hessenberg(A)`
- `lu_l(A)`
- `lu_u(A)`
- `mean(v)`
- `median(v)`
- `mode(v)`
- `percentile(v, p)`
- `quartile(v, q)`
- `var(v)`
- `std(v)`
- `cov(a, b)`
- `corr(a, b)`
- `linear_regression(x, y)`
- `lagrange(x, y, xi)`
- `spline(x, y, xi)`
- `dft(x)`
- `fft(x)`
- `idft(X)`
- `ifft(X)`
- `conv(a, b)`
- `convolve(a, b)`
- `poly_fit(x, y, degree)`
- `polynomial_fit(x, y, degree)`
- `poly_eval(p, x)`
- `poly_deriv(p)`
- `poly_integ(p)`
- `poly_compose(p, q)`
- `poly_gcd(p, q)`
- `poly_add(p, q)`
- `poly_sub(p, q)`
- `poly_mul(p, q)`
- `poly_div(p, q)`
- `roots(p)`
- `complex(re, im)`
- `real(z)`
- `imag(z)`
- `abs(z)`
- `arg(z)`
- `conj(z)`
- `polar(r, theta)`

Notes:

- `norm(m)` returns the Euclidean norm for vectors and the Frobenius norm for matrices
- `trace`, `det`, `eigvals`, and `eigvecs` require square matrices
- `eigvals` and `eigvecs` currently support real-valued results
- `dot(a, b)` and `outer(a, b)` require vector arguments
- `qr_q(A)` and `qr_r(A)` currently require square matrices
- `lu_l(A)` and `lu_u(A)` currently require square matrices and use LU decomposition without pivoting
- `null(m)` returns a basis matrix whose columns span the nullspace
- `svd_u(A)`, `svd_s(A)`, and `svd_vt(A)` return the reduced SVD factors
- `solve(A, b)` solves `Ax = b` for square `A` and vector `b`
- `dft` / `fft` return an `N x 2` matrix whose rows are `[real, imag]`
- `idft` / `ifft` accept a real vector or an `N x 2` complex matrix
- `conv` / `convolve` perform linear convolution on real vectors or `N x 2` complex sequences
- `abs(z)` returns the magnitude of a complex number

## Function Analysis And ODE

- `diff(f)`
- `diff(f, x0)`
- `diff(expr, variable)`
- `diff(expr, variable1, variable2, ...)`
- `integral(f)`
- `integral(expr, variable)`
- `integral(f, x0)`
- `integral(f, a, b)`
- `gradient(expr, variable1, variable2, ...)`
- `jacobian([expr1; expr2; ...], variable1, variable2, ...)`
- `hessian(expr, variable1, variable2, ...)`
- `critical(expr, variable1, variable2, ...)`
- `taylor(expr, a, n)`
- `pade(expr, m, n)`
- `pade(expr, a, m, n)`
- `puiseux(expr, degree, denominator)`
- `puiseux(expr, a, degree, denominator)`
- `series_sum(expr, n, lower, upper)`
- `summation(expr, n, lower, upper)`
- `fourier(expr)`
- `fourier(expr, t, w)`
- `ifourier(expr)`
- `ifourier(expr, w, t)`
- `laplace(expr)`
- `laplace(expr, t, s)`
- `ilaplace(expr)`
- `ilaplace(expr, s, t)`
- `ztrans(expr)`
- `ztrans(expr, n, z)`
- `iztrans(expr)`
- `iztrans(expr, z, n)`
- `simplify(expr)`
- `limit(f, x0)`
- `limit(f, x0, direction)`
- `extrema(f, a, b)`
- `double_integral(expr, x0, x1, y0, y1)`
- `double_integral(expr, x0, x1, y0, y1, nx, ny)`
- `double_integral_cyl(expr, r0, r1, theta0, theta1)`
- `double_integral_cyl(expr, r0, r1, theta0, theta1, nr, ntheta)`
- `double_integral_polar(expr, r0, r1, theta0, theta1)`
- `double_integral_polar(expr, r0, r1, theta0, theta1, nr, ntheta)`
- `triple_integral(expr, x0, x1, y0, y1, z0, z1)`
- `triple_integral(expr, x0, x1, y0, y1, z0, z1, nx, ny, nz)`
- `triple_integral_cyl(expr, r0, r1, theta0, theta1, z0, z1)`
- `triple_integral_cyl(expr, r0, r1, theta0, theta1, z0, z1, nr, ntheta, nz)`
- `triple_integral_sph(expr, rho0, rho1, theta0, theta1, phi0, phi1)`
- `triple_integral_sph(expr, rho0, rho1, theta0, theta1, phi0, phi1, nrho, ntheta, nphi)`
- `ode(rhs, x0, y0, x1)`
- `ode(rhs, x0, y0, x1, steps)`
- `ode(rhs, x0, y0, x1, steps, params_vec)`
- `ode(rhs, x0, y0, x1, steps, event_expr)`
- `ode(rhs, x0, y0, x1, steps, event_expr, params_vec)`
- `ode_table(rhs, x0, y0, x1)`
- `ode_table(rhs, x0, y0, x1, steps)`
- `ode_table(rhs, x0, y0, x1, steps, params_vec)`
- `ode_table(rhs, x0, y0, x1, steps, event_expr)`
- `ode_table(rhs, x0, y0, x1, steps, event_expr, params_vec)`
- `ode_system(rhs_vec, x0, y0_vec, x1)`
- `ode_system(rhs_vec, x0, y0_vec, x1, steps)`
- `ode_system(rhs_vec, x0, y0_vec, x1, steps, params_vec)`
- `ode_system(rhs_vec, x0, y0_vec, x1, steps, event_expr)`
- `ode_system(rhs_vec, x0, y0_vec, x1, steps, event_expr, params_vec)`
- `ode_system_table(rhs_vec, x0, y0_vec, x1)`
- `ode_system_table(rhs_vec, x0, y0_vec, x1, steps)`
- `ode_system_table(rhs_vec, x0, y0_vec, x1, steps, params_vec)`
- `ode_system_table(rhs_vec, x0, y0_vec, x1, steps, event_expr)`
- `ode_system_table(rhs_vec, x0, y0_vec, x1, steps, event_expr, params_vec)`
- `lp_max(c, A, b, lower, upper)`
- `lp_max(c, A, b, Aeq, beq, lower, upper)`
- `lp_min(c, A, b, lower, upper)`
- `lp_min(c, A, b, Aeq, beq, lower, upper)`
- `ilp_max(c, A, b, lower, upper)`
- `ilp_max(c, A, b, Aeq, beq, lower, upper)`
- `ilp_min(c, A, b, lower, upper)`
- `ilp_min(c, A, b, Aeq, beq, lower, upper)`
- `milp_max(c, A, b, lower, upper, integrality)`
- `milp_max(c, A, b, Aeq, beq, lower, upper, integrality)`
- `milp_min(c, A, b, lower, upper, integrality)`
- `milp_min(c, A, b, Aeq, beq, lower, upper, integrality)`
- `bip_max(c, A, b)`
- `bip_max(c, A, b, Aeq, beq)`
- `bip_min(c, A, b)`
- `bip_min(c, A, b, Aeq, beq)`
- `solve(expr, x0)`
- `bisect(expr, a, b)`
- `secant(expr, x0, x1)`
- `fixed_point(expr, x0)`

Notes:

- `diff`, `integral`, `taylor`, `limit`, and `extrema` operate on one-variable functions
- symbolic `diff(expr)` and `integral(expr)` also accept raw one-variable expressions and infer the variable name automatically
- `taylor` also accepts raw one-variable symbolic expressions, not only named custom functions
- `pade` builds a rational approximant from the local Taylor coefficients of `expr`
- `puiseux` uses a denominator grid; for example `denominator = 2` allows half-integer powers
- `series_sum` / `summation` currently cover polynomial summands up to degree 3 and common geometric series
- use `inf`, `oo`, or `infinity` as the upper bound for supported infinite geometric sums
- symbolic transform commands default to `(t, w)` for Fourier, `(t, s)` for Laplace, and `(n, z)` for z transforms when variables are omitted
- current symbolic Fourier/Laplace/z support is rule-based and focuses on common signal-analysis forms such as constants, exponentials, `sin/cos`, `step`, and `delta`
- symbolic `simplify(expr)` may be used on multi-variable expressions, but symbolic `diff/integral/taylor/limit/extrema` still require a one-variable input
- `double_integral` and `triple_integral` use Cartesian coordinates
- `double_integral_cyl` / `double_integral_polar` use the planar polar Jacobian `r`
- `triple_integral_cyl` uses cylindrical coordinates `(r, theta, z)` with Jacobian `r`
- `triple_integral_sph` uses spherical coordinates `(rho, theta, phi)` with Jacobian `rho^2 sin(phi)`
- in cylindrical and spherical modes, Cartesian variables such as `x`, `y`, and `z` are also available inside the integrand
- `ode` solves the first-order initial value problem `y' = rhs(x, y)` with the classical RK4 method
- `ode` returns the approximated value at `x1`
- `ode_table` returns a two-column matrix whose rows are sampled `(x, y)` pairs
- `ode_system` solves nonlinear ODE systems whose right-hand side evaluates to a vector
- inside `ode_system`, the state is available as `y1`, `y2`, ... and also as the vector variable `y`
- `ode_system_table` returns a matrix with columns `[x, y1, y2, ...]`
- optional ODE parameters are exposed as `p`, and vector parameters also expose `p1`, `p2`, ...
- optional ODE event expressions are scalar expressions; integration stops when the event crosses zero
- the ODE step count must be a positive integer
- integration subdivision counts must be positive integers; odd counts are rounded up internally for Simpson integration
- `ode` defaults to `100` steps and `ode_table` defaults to `10` steps
- `ode_system` defaults to `100` steps and `ode_system_table` defaults to `10` steps
- `lp_max/lp_min` solve box-constrained linear programs with `A * x <= b`, optional `Aeq * x = beq`, and explicit lower/upper bounds
- `ilp_max/ilp_min` solve bounded integer programs and require integer bounds
- `milp_max/milp_min` use an `integrality` vector where nonzero entries mark integer variables
- `bip_max/bip_min` are binary-program shortcuts with implicit bounds `0 <= x <= 1`
- the current planning solvers are intended for small problems and favor reliability over scale
- current symbolic integration rules cover common powers, `1/(ax+b)`, `sin/cos/tan`, `exp/ln`, `sqrt/cbrt`, `asin/acos/atan`, `abs`, and some polynomial-times-`exp/sin/cos` products

## Integer / Number Theory

- `gcd(a, b)`
- `lcm(a, b)`
- `mod(a, b)`
- `factor(n)`
- `factorial(n)`
- `nCr(n, r)`
- `binom(n, r)`
- `nPr(n, r)`
- `fib(n)`
- `is_prime(n)`
- `next_prime(n)`

Notes:

- `gcd`, `lcm`, and `mod` require integer arguments
- `factor(n)` returns a factorization string, not a plain numeric value
- `factorial(n)` requires an integer `n >= 0`
- `nCr(n, r)` and `nPr(n, r)` require integer inputs with `0 <= r <= n`

## Unit Conversion

- `deg2rad(x)`
- `rad2deg(x)`
- `deg(x)`
- `rad(x)`
- `sin_deg(x)`
- `cos_deg(x)`
- `celsius(f)`
- `fahrenheit(c)`
- `kelvin(c)`
- `c2f(c)`
- `f2c(f)`

Notes:

- `deg2rad(x)` converts degrees to radians
- `rad2deg(x)` converts radians to degrees
- `celsius(f)` converts Fahrenheit to Celsius
- `fahrenheit(c)` converts Celsius to Fahrenheit
- `kelvin(c)` converts Celsius to Kelvin

## Base Conversion

- `bin(n)`
- `oct(n)`
- `hex(n)`
- `base(n, b)`

Notes:

- these are display-oriented conversion functions
- arguments must be integers
- supported bases for `base(n, b)` are `2..16`
- `:hexprefix` and `:hexcase` control how hexadecimal output is displayed

Examples:

- `bin(10)` -> `1010`
- `oct(83)` -> `123`
- `hex(255)` -> `FF`
- with `:hexprefix on`, `hex(255)` -> `0xFF`
- with `:hexprefix on` and `:hexcase lower`, `hex(255)` -> `0xff`
- `base(-31, 16)` -> `-1F`

## Bitwise

- `and(a, b)`
- `or(a, b)`
- `xor(a, b)`
- `not(a)`
- `shl(a, n)`
- `shr(a, n)`
- `rol(a, n)`
- `ror(a, n)`
- `popcount(n)`
- `bitlen(n)`
- `ctz(n)`
- `clz(n)`
- `parity(n)`
- `reverse_bits(n)`

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
