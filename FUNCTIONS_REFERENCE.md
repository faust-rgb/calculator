# Functions Reference

## Commands Reference

### General

- `help`
- `:help`
- `:help commands`
- `:help functions`
- `:help matrix`
- `:help examples`
- `:help exact`
- `:help variables`
- `:help persistence`
- `:help programmer`
- `exit`
- `quit`

### Exact Mode

- `:exact on`
  Enable exact fraction mode
- `:exact off`
  Disable exact fraction mode
- `:exact`
  Show current exact mode status

### Symbolic Constants Mode

- `:symbolic on`
  Preserve `pi` and `e` in scalar display results
- `:symbolic off`
  Return to normal numeric display
- `:symbolic`
  Show current symbolic constants mode status

### Programmer Formatting

- `:hexprefix on`
  Show hex results with a `0x` prefix
- `:hexprefix off`
  Hide the `0x` prefix
- `:hexprefix`
  Show current hex prefix mode status
- `:hexcase upper`
  Use uppercase hex digits such as `0xFF`
- `:hexcase lower`
  Use lowercase hex digits such as `0xff`
- `:hexcase`
  Show current hex letter-case mode

### Variables

- `name = expression`
  Assign a variable
- `:vars`
  List all stored variables
- `:clear name`
  Remove a single variable
- `:clear`
  Remove all variables

### History

- up arrow
  Recall previous command in interactive mode
- down arrow
  Move toward newer history entries
- `:history`
  Print the current session history

### Persistence

- `:save file`
  Save variables to a file
- `:load file`
  Load variables from a file
- `:run file.calc`
  Execute a `.calc` script file

Notes:

- persistence supports scalar variables, strings, matrices, custom functions,
  and script functions

### Scripting

Use `bin/calculator file.calc` or `:run file.calc` to execute scripts. The
dedicated syntax guide is `test/script/SYNTAX_GUIDE.md`.

### Autocomplete

- `Tab`
  Autocomplete commands, functions, variables, and custom functions
- double `Tab`
  Show the current candidate list when multiple completions match

Examples:

- type `:he` then press `Tab` -> `:help`
- type `sq` then press `Tab` -> `sqrt(`
- type `:help ma` then press `Tab` -> `:help matrix`
- type `:help ` then double `Tab` -> show help topics
- type `g` inside `diff(g` then press `Tab` -> complete a matching custom function

### Prefixed Integer Literals

These are expression features rather than commands, but they are commonly used
like shell-style numeric shortcuts:

- `0b1010`
- `0o77`
- `0xFF`

They can be mixed directly into expressions:

- `0b1010 + 0xF`
- `and(0xF, 0b1010)`
- `rol(1, 3)`
- `popcount(0xF0)`

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

### Plotting

- `plot(f(x), x, start, end)`
- `plot(f(x), x, start, end, points)`
- `plot(f(x), start, end)` (defaults to variable `x`)

Examples:
- `plot(sin(x), -pi, pi)`
- `plot(x^2, x, -5, 5, 50)`

---

### Matrix Operations


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

## Signal Processing

### FFT Commands

- `fft(signal)`
  Fast Fourier Transform
- `ifft(spectrum)`
  Inverse FFT
- `rfft(signal)` (planned)
  Real FFT (optimized for real signals)

### Convolution & Correlation

- `conv(s1, s2)`
  Linear convolution
- `cconv(s1, s2 [, n])` (planned)
  Circular convolution
- `xcorr(s1, s2)` (planned)
  Cross-correlation
- `autocorr(signal)` (planned)
  Auto-correlation

### Window Functions

- `hann(n)` / `hanning(n)`
  Hanning window
- `hamming(n)`
  Hamming window
- `blackman(n)`
  Blackman window

### Filter Design (planned)

- `fir_design(order, cutoff, type [, window])`
  Design FIR filter (types: `lowpass`, `highpass`, `bandpass`, `bandstop`)
- `iir_design(order, cutoff, type)`
  Design Butterworth IIR filter
- `filter(b, a, signal)`
  Apply digital filter
- `freqz(b, a [, n])`
  Frequency response

### Time-Frequency Analysis (planned)

- `psd(signal [, nfft])`
  Power spectral density (Welch method)
- `stft(signal [, nfft])`
  Short-time Fourier transform
- `spectrogram(signal [, nfft])`
  Spectrogram

Notes:

- `fft` supports arbitrary-length signals using mixed-radix and Bluestein algorithms
- `rfft` (planned) will return only positive frequencies for real signals
- `hann`, `hamming`, `blackman` return normalized window vectors
- `fir_design` (planned) returns numerator `b` and denominator `a` coefficients
- `psd` (planned) uses Welch's method with 50% overlap and Hanning window by default

Examples:

- `fft([1, 2, 3, 4])` -> `[10+0j, -2+2j, -2+0j, -2-2j]`
- `conv([1, 2, 3], [1, 1])` -> `[1, 3, 5, 3]`
- `hann(5)` -> `[0, 0.5, 1, 0.5, 0]`
- `hamming(8)` -> `[0.08, 0.253, 0.642, 0.954, 0.954, 0.642, 0.253, 0.08]`
- `fir_design(16, 0.2, lowpass)` (planned) -> returns `b` and `a` coefficients
- `psd([1, 2, 3, 4, 5, 6, 7, 8], 8)` (planned) -> power spectral density

Notes:

- `norm(m)` returns the Euclidean norm for vectors and the Frobenius norm for matrices
- `trace`, `det`, `eigvals`, and `eigvecs` require square matrices
- `eigvals` returns real eigenvalues as a vector and complex eigenvalues as
  an `N x 2` matrix with rows `[real, imag]`; `eigvecs` currently supports
  real-valued eigenvectors only
- `dot(a, b)` and `outer(a, b)` require vector arguments
- `qr_q(A)` and `qr_r(A)` currently require square matrices
- `lu_l(A)` and `lu_u(A)` currently require square matrices and use LU decomposition without pivoting
- `null(m)` returns a basis matrix whose columns span the nullspace
- `svd_u(A)`, `svd_s(A)`, and `svd_vt(A)` return the reduced SVD factors
- `solve(A, b)` solves `Ax = b` for square `A` and vector `b`
- 3x3 inverse, determinant, RREF, and linear solve workflows are covered by the
  regression suite in addition to the smaller examples below
- `dft` / `fft` return an `N x 2` matrix whose rows are `[real, imag]`
- `idft` / `ifft` accept a real vector or an `N x 2` complex matrix
- `conv` / `convolve` perform linear convolution on real vectors or `N x 2` complex sequences
- `abs(z)` returns the magnitude of a complex number (for real inputs, returns absolute value)
- `complex(re, im)` creates a complex number from real and imaginary parts
- `real(z)` returns the real part of a complex number
- `imag(z)` returns the imaginary part of a complex number
- `arg(z)` returns the argument (phase angle) of a complex number
- `conj(z)` returns the complex conjugate
- `polar(r, theta)` creates a complex number from polar coordinates

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
- `lagrange(f, constraint, variable1, variable2, ...)`
- `lagrange(f, [constraint1, constraint2, ...], variable1, variable2, ...)`
- `dsolve(rhs, x_var, y_var)`
- `divergence([expr1, expr2, ...], variable1, variable2, ...)`
- `div([expr1, expr2, ...], variable1, variable2, ...)`
- `curl([expr1, expr2, expr3], variable1, variable2, variable3)`
- `curl_2d([expr1, expr2], variable1, variable2)`
- `laplacian(expr, variable1, variable2, ...)`
- `implicit_diff(F, y, x)`
- `param_deriv(x_t, y_t, t)`
- `param_deriv(x_t, y_t, t, order)`
- `directional(expr, var1, var2, ..., v1, v2, ...)`
- `line_integral(f_or_F, [x_t, y_t, ...], t, a, b)`
- `surface_integral(f_or_F, [x_u_v, y_u_v, z_u_v], u, a, b, v, c, d)`
- `greens_theorem([F_x, F_y], [x_t, y_t], t, a, b)`
- `stokes_theorem([F_x, F_y, F_z], [x_u_v, y_u_v, z_u_v], u, a, b, v, c, d)`
- `divergence_theorem([F_x, F_y, F_z], [x_u_v, y_u_v, z_u_v], u, a, b, v, c, d)`
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
- `double_integral(expr, r0, r1, theta0, theta1, "polar")`
- `triple_integral(expr, x0, x1, y0, y1, z0, z1)`
- `triple_integral(expr, r0, r1, theta0, theta1, z0, z1, "cyl")`
- `triple_integral(expr, rho0, rho1, theta0, theta1, phi0, phi1, "sph")`
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
- symbolic `diff(expr, variable1, variable2, ...)` applies chained partial
  differentiation in the provided variable order
- symbolic `gradient`, `jacobian`, and `hessian` accept multi-variable
  expressions and return vector/matrix displays
- `taylor` also accepts raw one-variable symbolic expressions, not only named custom functions
- `pade` builds a rational approximant from the local Taylor coefficients of `expr`
- `critical` appends a Hessian-based classification (`local min`, `local max`, `saddle`, or `degenerate`) for isolated 1-3 variable points
- `lagrange` solves constrained optimization problems using Lagrange multipliers; it constructs the Lagrangian and delegates to `critical`
- `dsolve` attempts symbolic solution of first-order ODEs; currently supports `y' = f(x)` and linear ODEs `y' + P(x)y = Q(x)`
- `divergence`/`div` computes the divergence of a vector field
- `curl` computes the curl of a 3D vector field
- `curl_2d([F_x, F_y], x, y)` computes the scalar 2D curl ∂Fy/∂x - ∂Fx/∂y
- `laplacian` computes the Laplacian of a scalar field
- `implicit_diff(F, y, x)` computes the implicit derivative dy/dx for F(x,y) = 0 using the formula dy/dx = -F_x / F_y
- `param_deriv(x_t, y_t, t)` computes the derivative dy/dx for a parametric curve (x(t), y(t)); higher-order derivatives are supported via the optional `order` parameter
- `directional(expr, vars..., direction...)` computes the directional derivative of `expr` in the given direction (automatically normalized)
- `line_integral(f_or_F, curve, t, a, b)` computes the scalar or vector line integral along a parametric curve
- `surface_integral(f_or_F, surface, u, a, b, v, c, d)` computes the scalar or flux surface integral over a parametric surface
- both automatically detect whether the integrand is a scalar field or a vector field
- `greens_theorem([F_x, F_y], [x_t, y_t], t, a, b)` computes the line integral ∮_C F·dr using Green's theorem (2D)
- `stokes_theorem([F_x, F_y, F_z], [x_u_v, y_u_v, z_u_v], u, a, b, v, c, d)` computes the surface integral ∬_S curl(F)·n dS using Stokes' theorem
- `divergence_theorem([F_x, F_y, F_z], [x_u_v, y_u_v, z_u_v], u, a, b, v, c, d)` computes the flux ∬_S F·n dS using the divergence theorem
- symbolic `integral` includes rule-based support for mixed real-linear plus repeated irreducible quadratic rational factors and common trig power/product identities
- `puiseux` uses a denominator grid; for example `denominator = 2` allows half-integer powers
- `series_sum` / `summation` covers polynomial summands (via Faulhaber formulas), geometric series, and infinite series involving Riemann Zeta values $\zeta(2k)$ up to $k=5$
- use `inf`, `oo`, or `infinity` as the upper bound for supported infinite sums
- `limit` uses a hybrid approach: symbolic Power Series Arithmetic (PSA) for exact results and removable singularities, plus high-order Richardson extrapolation (14th order) for numerical fallback
- infinity limits (`inf` or `-inf`) are handled via symbolic substitution $x=1/t$ followed by series analysis
- symbolic transform commands default to `(t, w)` for Fourier, `(t, s)` for Laplace, and `(n, z)` for z transforms when variables are omitted
- current symbolic Fourier/Laplace/z support is rule-based and focuses on common signal-analysis forms such as constants, exponentials, `sin/cos`, `step`, and `delta`
- symbolic `simplify(expr)` may be used on multi-variable expressions, but symbolic `diff/integral/taylor/limit/extrema` still require a one-variable input
- `double_integral` and `triple_integral` use Cartesian coordinates by default
- `double_integral` with `"polar"` uses the planar polar Jacobian `r`
- `triple_integral` with `"cyl"` uses cylindrical coordinates `(r, theta, z)` with Jacobian `r`
- `triple_integral` with `"sph"` uses spherical coordinates `(rho, theta, phi)` with Jacobian `rho^2 sin(phi)`
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
- `divisors(n)`
- `extended_gcd(a, b)`
- `xgcd(a, b)`

Notes:

- `gcd`, `lcm`, and `mod` require integer arguments
- `factor(n)` returns a factorization string, not a plain numeric value
- `factorial(n)` requires an integer `n >= 0`
- `nCr(n, r)` and `nPr(n, r)` require integer inputs with `0 <= r <= n`
- `divisors(n)` returns a vector of all positive divisors of `n` in ascending order
- `extended_gcd(a, b)` / `xgcd(a, b)` returns `[g, x, y]` where `g = gcd(a, b)` and `a*x + b*y = g`

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
- `c` (Speed of Light)
- `G` (Gravitational Constant)
- `h` (Planck Constant)
- `k` (Boltzmann Constant)
- `NA` (Avogadro Number)

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
- `eigvals` returns real eigenvalues as a vector and complex eigenvalues as an `N x 2` matrix of `[real, imag]` rows
- `eigvecs` currently supports real-valued eigenvectors only
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
