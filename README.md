# Command Line Calculator

This project is a C++ command line calculator that evaluates expressions without
using the standard math library implementations from `<cmath>` or `math.h`.

## Project Layout

- `src/app`
  REPL entry point and command dispatch
- `src/core`
  Calculator API and main expression/script execution logic
- `src/math`
  Custom numeric helpers and math primitives, split into standard `.cpp` translation units
- `src/matrix`
  Matrix types and operations, including a dedicated linear-algebra implementation unit
- `src/analysis`
  One-variable function analysis
- `src/algebra`
  Polynomial helpers
- `src/symbolic`
  Symbolic expression support
- `src/script`
  Script AST and parser
- `test`
  Regression tests and runnable example scripts

## Additional Docs

- `ARCHITECTURE.md`
- `FUNCTIONS_REFERENCE.md`
- `MATRIX_GUIDE.md`
- `COMMANDS_REFERENCE.md`
- `HANDOFF.md`
- `TESTING.md`
- `CHANGELOG.md`
- `ROADMAP.md`
- `KNOWN_LIMITATIONS.md`
- `SIMPLIFY_IMPROVEMENTS.md`

## Source Organization

Large subsystems are now split across normal `.cpp` files instead of implementation
fragments:

- `src/math/mymath.cpp` and `src/math/mymath_special_functions.cpp`
- `src/matrix/matrix.cpp` and `src/matrix/matrix_linear_algebra.cpp`
- `src/core/calculator.cpp` and `src/core/calculator_lifecycle.cpp`
- `src/symbolic/symbolic_expression_core.cpp`,
  `src/symbolic/symbolic_expression_calculus.cpp`, and
  `src/symbolic/symbolic_expression_transforms.cpp`

Shared internal declarations for these splits live in private headers such as
`src/math/mymath_internal.h`, `src/matrix/matrix_internal.h`,
`src/core/calculator_internal_types.h`, and
`src/symbolic/symbolic_expression_internal.h`.

## Features

- Basic arithmetic: `+`, `-`, `*`, `/`
- Power operator: `^`
- Scientific notation input such as `1e-3`, `1e20`, and `1e-300`
- Functions: `sin(x)`, `cos(x)`, `tan(x)`, `asin(x)`, `acos(x)`, `atan(x)`,
  `sec(x)`, `csc(x)`, `cot(x)`, `asec(x)`, `acsc(x)`, `acot(x)`,
  `sinh(x)`, `cosh(x)`, `tanh(x)`, `asinh(x)`, `acosh(x)`, `atanh(x)`,
  `exp(x)`, `exp2(x)`, `ln(x)`, `log(x, base)`, `log2(x)`, `log10(x)`,
  `gamma(x)`, `beta(a, b)`, `zeta(s)`, `erf(x)`, `erfc(x)`, `bessel(n, x)`,
  `sqrt(x)`, `cbrt(x)`, `root(a, n)`, `abs(x)`, `sign(x)`, `floor(x)`,
  `ceil(x)`, `round(x)`, `trunc(x)`, `pow(a, b)`
- Comparison and integer utilities: `min(a, b)`, `max(a, b)`, `clamp(x, min, max)`, `gcd(a, b)`, `lcm(a, b)`, `mod(a, b)`, `factorial(n)`, `nCr(n, r)`, `binom(n, r)`, `nPr(n, r)`, `fib(n)`, `is_prime(n)`, `next_prime(n)`
- Rational approximation helper: `rat(x)` and `rat(x, max_denominator)`
- Aggregate helpers: `sum(...)`, `mean(...)`, `avg(...)`, `median(...)`, `mode(...)`, `percentile(...)`, `quartile(...)`, `var(...)`, `std(...)`
- Probability helpers: `rand()`, `randn()`, `randint(a, b)`, `pdf_normal(x, mu, sigma)`, `cdf_normal(x, mu, sigma)`
- Unit conversion helpers: `deg(x)`, `rad(x)`, `deg2rad(x)`, `rad2deg(x)`, `sin_deg(x)`, `cos_deg(x)`, `celsius(f)`, `fahrenheit(c)`, `kelvin(c)`, `c2f(c)`, `f2c(f)`
- Prime factorization with `factor(n)`
- Base conversion with `bin(n)`, `oct(n)`, `hex(n)`, `base(n, b)`
- Hex formatting controls with `:hexprefix` and `:hexcase`
- Bitwise operations with `and(a, b)`, `or(a, b)`, `xor(a, b)`, `not(a)`, `shl(a, n)`, `shr(a, n)`, `rol(a, n)`, `ror(a, n)`, `popcount(n)`, `bitlen(n)`, `ctz(n)`, `clz(n)`, `parity(n)`, `reverse_bits(n)`
- Matrix creation with `vec(...)`, `mat(...)`, `zeros(...)`, `eye(...)`
- Matrix editing with `resize(...)`, `append_row(...)`, `append_col(...)`, `get(...)`, `set(...)`
- Matrix transpose with `transpose(...)`
- Matrix explicit inverse with `inverse(...)`
- Matrix pseudo-inverse, condition number, and structured products with `pinv(...)`, `cond(...)`, `kron(...)`, `hadamard(...)`
- Matrix operations with matrix/matrix and matrix/scalar `+`, `-`, `*`, `/`, `^`
- Invertible matrices support negative integer powers such as `A ^ -1`
- Matrix analysis with `norm(...)`, `trace(...)`, `det(...)`, `rank(...)`, `rref(...)`, `eigvals(...)`, `eigvecs(...)`, `eig(...)`, `svd(...)`
- Linear system solving with `solve(A, b)` for `Ax = b`
- Vector/matrix helpers with `dot(...)`, `outer(...)`, `null(...)`, `least_squares(...)`, `qr_q(...)`, `qr_r(...)`, `diag(...)`, `reshape(...)`, `vec(m)`
- Reduced singular value decomposition with `svd_u(...)`, `svd_s(...)`, `svd_vt(...)`
- Wide-matrix QR/SVD support, SVD-based pseudo-inverse for singular matrices,
  and minimum-norm `least_squares(...)` solutions for underdetermined systems
- Matrix decompositions with `cholesky(...)`, `schur(...)`, `hessenberg(...)`
- Vector statistics and fitting with `cov(...)`, `corr(...)`, `lagrange(...)`, `spline(...)`, `linear_regression(...)`, `poly_fit(...)`
- Lightweight complex helpers with `complex(re, im)`, `real(z)`, `imag(z)`, `abs(z)`, `arg(z)`, `conj(z)`, `polar(r, theta)`
- Constants: `pi`, `e`
- Symbolic constants display mode with `:symbolic on` and `:symbolic off`
- Parentheses and unary plus/minus
- Variable assignment, for example `x = 1/3 + 1/4`
- Toggleable exact fraction mode with `:exact on` and `:exact off`
- Variable management commands: `:vars`, `:clear name`, `:clear`
- Custom function commands: `:funcs`, `:clearfunc name`, `:clearfuncs`
- Session history command: `:history`
- Interactive help topics for exact mode, variables, persistence, and programmer tools
- State persistence commands: `:save file`, `:load file`
- Script execution with `:run file` or redirected stdin
- Script language support for `fn`, `if/else`, `while`, `for`, `return`, `break`, `continue`, strings, and `print(...)`
- Separate one-variable custom function analysis module with evaluation,
  derivative, definite integral, indefinite integral value, interval
  extrema solving, polynomial arithmetic, Taylor expansion, polynomial roots,
  first-order ODE initial-value solving with `ode(...)` / `ode_table(...)`,
  nonlinear ODE system solving with `ode_system(...)` / `ode_system_table(...)`,
  optional ODE event stopping and parameter passing,
  root-finding with `solve(...)`, `bisect(...)`, `secant(...)`, `fixed_point(...)`,
  and multi-variable integration in Cartesian, cylindrical, and spherical coordinates
- Box-constrained linear and integer planning with `lp_max(...)`, `lp_min(...)`,
  `ilp_max(...)`, mixed-integer planning with `milp_max(...)` / `milp_min(...)`,
  optional equality constraints `Aeq * x = beq`, and binary shortcuts `bip_max(...)` / `bip_min(...)`
- Improved numerical stability for large-angle trigonometric reduction, `beta`,
  `gamma`, `bessel`, and removable-singularity style `limit(...)` evaluations

## Build

```bash
make
```

This generates the executable `calculator`.

The regression suite lives in `test/tests.cpp`.

Current validation status:

- `make test`
- expected result: `Failed: 0`

## Run

```bash
./calculator
```

Example session:

```text
> 1 + 2 * 3
7
> :exact on
Exact fraction mode: ON
> 1/3 + 1/4
7/12
> gcd(48, 18)
6
> mod(17, 5)
2
> pow(3, 4)
81
> cbrt(-8)
-2
> root(27, 3)
3
> min(7/3, 5/2)
7/3
> sign(-7/3)
-1
> factor(360)
2^3 * 3^2 * 5
> rat(pi)
355/113
> hex(255)
FF
> m = mat(2, 2, 1, 2, 3, 4)
m = [[1, 2], [3, 4]]
> m + eye(2)
[[2, 2], [3, 5]]
> get(m, 1, 0)
3
> m = set(m, 1, 0, 8)
m = [[1, 2], [8, 4]]
> det(m)
-12
> mat(2, 2, 1, 2, 3, 4) ^ -1
[[-2, 1], [1.5, -0.5]]
> inverse(mat(2, 2, 1, 2, 3, 4))
[[-2, 1], [1.5, -0.5]]
> dot(vec(1, 2, 3), vec(4, 5, 6))
32
> svd_s(mat(3, 2, 3, 0, 0, 2, 0, 0))
[[3, 0], [0, 2]]
> solve(mat(2, 2, 2, 1, 5, 3), vec(1, 2))
[[1], [-1]]
> rref(mat(2, 3, 1, 2, 3, 2, 4, 6))
[[1, 2, 3], [0, 0, 0]]
> and(0xF, 0b1010)
10
> x = 1/3 + 1/4
x = 7/12
> x + 1/6
3/4
> :vars
x = 7/12
> lcm(12, 18)
36
> :exact off
Exact fraction mode: OFF
> :symbolic on
Symbolic constants mode: ON
> pi / 2 + e
pi / 2 + e
> sin(pi / 2)
1
> log10(1000)
3
> atan(1)
0.785398163397
> f(x) = sin(x) + x ^ 2
f(x) = sin(x) + x ^ 2
> f(2)
4.90929742683
> p(x) = x ^ 2 + 2 * x + 1
p(x) = x ^ 2 + 2 * x + 1
> q(x) = x - 1
q(x) = x - 1
> poly_mul(p, q)
x ^ 3 + x ^ 2 - x - 1
> roots(p)
-1
> diff(f)
cos(x) + 2 * x
> diff(abs(x))
sign(x)
> diff(f, 0)
1
> integral(f)
-cos(x) + x ^ 3 / 3 + C
> integral(1 / x)
ln(abs(x)) + C
> integral(x * exp(x))
exp(x) * x - exp(x) + C
> taylor(f, 0, 3)
x + x ^ 2 - 0.166666666667 * x ^ 3
> h(x) = sin(x) / x
h(x) = sin(x) / x
> limit(h, 0)
1
> integral(f, 0, 1)
0.793031027466
> double_integral(x + y, 0, 1, 0, 2)
3
> triple_integral_sph(1, 0, 1, 0, 2 * pi, 0, pi)
4.18879020479
> ode(y - x ^ 2 + 1, 0, 0.5, 2, 20)
5.30547195055
> ode_table(y, 0, 1, 1, 4)
[0, 1; 0.25, 1.28402541668; 0.5, 1.64872127069; 0.75, 2.11700001658; 1, 2.7182818284]
> ode_system(vec(y2, -y1), 0, vec(0, 1), pi / 2, 40)
[1, 0]
> lp_max(vec(3, 2), mat(3, 2, 1, 1, 1, 0, 0, 1), vec(4, 2, 3), vec(0, 0), vec(10, 10))
x = [2, 2]
objective = 10
> extrema(f, -1, 1)
min: x = -0.450183689594, f(x) = -0.232465575151
> exit
```

## Supported Expression Examples

```text
2 ^ 8
sqrt(25) + exp(1)
ln(e)
sin(pi / 6) * cos(pi / 3)
asin(0.5) + acos(0.5)
gcd(48, 18)
lcm(12, 18)
mod(17, 5)
abs(-7/3)
sign(-7/3)
floor(7/3)
ceil(7/3)
min(7/3, 5/2)
max(7/3, 5/2)
pow(3, 4)
sinh(1)
gamma(5)
factorial(5)
nCr(5, 2)
nPr(5, 2)
sum(1, 2, 3, 4)
avg(1, 2, 3, 4)
median(9, 1, 5, 2)
cbrt(-8)
root(27, 3)
deg2rad(180)
rad2deg(pi / 2)
celsius(212)
fahrenheit(100)
kelvin(0)
factor(360)
bin(10)
oct(83)
hex(255)
base(31, 2)
and(6, 3)
or(6, 3)
xor(6, 3)
not(0)
shl(3, 2)
shr(16, 2)
vec(1, 2, 3)
mat(2, 2, 1, 2, 3, 4)
zeros(2, 3)
eye(3)
resize(mat(2, 2, 1, 2, 3, 4), 3, 3)
append_row(mat(1, 2, 1, 2), 3, 4)
append_col(mat(2, 1, 1, 2), 3, 4)
transpose(mat(2, 3, 1, 2, 3, 4, 5, 6))
inverse(mat(2, 2, 1, 2, 3, 4))
mat(2, 2, 1, 2, 3, 4) ^ -1
dot(vec(1, 2, 3), vec(4, 5, 6))
outer(vec(1, 2), vec(3, 4, 5))
null(mat(2, 3, 1, 2, 3, 2, 4, 6))
least_squares(mat(2, 1, 1, 1), vec(2, 4))
qr_q(mat(2, 2, 2, 0, 0, 3))
qr_r(mat(2, 2, 2, 0, 0, 3))
svd_u(mat(3, 2, 3, 0, 0, 2, 0, 0))
svd_s(mat(3, 2, 3, 0, 0, 2, 0, 0))
svd_vt(mat(3, 2, 3, 0, 0, 2, 0, 0))
get(mat(2, 2, 1, 2, 3, 4), 1, 0)
set(vec(5, 6, 7), 1, 42)
norm(vec(3, 4))
trace(mat(2, 2, 1, 2, 3, 4))
det(mat(2, 2, 1, 2, 3, 4))
rank(mat(2, 2, 1, 2, 2, 4))
rref(mat(2, 3, 1, 2, 3, 2, 4, 6))
eigvals(mat(2, 2, 2, 0, 0, 3))
eigvecs(mat(2, 2, 2, 0, 0, 3))
solve(mat(2, 2, 2, 1, 5, 3), vec(1, 2))
f(x) = sin(x) + x ^ 2
f(2)
p(x) = x ^ 2 + 2 * x + 1
q(x) = x - 1
poly_add(p, q)
poly_sub(p, q)
poly_mul(p, q)
poly_div(p, q)
roots(p)
diff(f)
diff(f, 2)
integral(f)
taylor(f, 0, 3)
limit(f, 0)
integral(f, 0, 3)
integral(f, 3)
double_integral(x + y, 0, 1, 0, 2)
double_integral_cyl(1, 0, 1, 0, 2*pi)
triple_integral(x * y * z, 0, 1, 0, 1, 0, 1)
triple_integral_sph(1, 0, 1, 0, 2*pi, 0, pi)
ode(y - x, 0, 1, 2)
ode_table(y, 0, 1, 1, 4)
ode_system(vec(y2, -y1), 0, vec(0, 1), pi / 2, 40)
lp_max(vec(3, 2), mat(3, 2, 1, 1, 1, 0, 0, 1), vec(4, 2, 3), vec(0, 0), vec(10, 10))
ilp_max(vec(3, 2), mat(3, 2, 1, 1, 1, 0, 0, 1), vec(4, 2, 3), vec(0, 0), vec(10, 10))
extrema(f, -2, 2)
```

You can also enter prefixed integer literals directly:

```text
0b1010
0o77
0xFF
```

## Exact Fraction Mode

When exact mode is enabled, expressions that stay in the rational-number domain
are reduced to simplest form.

```text
> :exact on
> 1/3 + 1/4
7/12
> (2/3) * (9/4)
3/2
```

For expressions that are not naturally rational, such as `sin(pi)` or
`sqrt(2)`, the calculator falls back to decimal output.

## Symbolic Constants Mode

Use symbolic constants mode when you want scalar results that involve `pi` or
`e` to stay in symbolic form instead of being displayed as decimal
approximations.

```text
> :symbolic on
Symbolic constants mode: ON
> x = pi / 2
x = pi / 2
> x + 1
pi / 2 + 1
> f(x) = x + pi
f(x) = x + pi
> f(e)
e + pi
> :symbolic off
Symbolic constants mode: OFF
```

Notes:

- this mode affects scalar display and variable reuse
- it is intended for symbolic-looking output, not exact algebraic simplification
- common constant identities such as `sin(pi / 2)`, `cos(pi)`, `tan(pi / 4)`,
  `ln(e)`, and `exp(ln(pi))` are simplified automatically
- unsupported expressions fall back to the normal numeric result

## Variables

You can store results in variables and reuse them later:

```text
> :exact on
> x = 1/3 + 1/4
x = 7/12
> x + 1/6
3/4
> y = sin(pi / 2)
y = 1
> y + 1/2
1.5
```

Variable names must start with a letter and may contain letters, digits, and
underscores.

## Custom Functions

You can define a unary custom function with `name(x) = expression`, where the
expression is composed from the existing operators and supported `mymath`
functions.

```text
> f(x) = sin(x) + x ^ 2
f(x) = sin(x) + x ^ 2
> f(2)
4.90929742683
> p(x) = x ^ 2 + 2 * x + 1
p(x) = x ^ 2 + 2 * x + 1
> q(x) = x - 1
q(x) = x - 1
> poly_add(p, q)
x ^ 2 + 3 * x
> poly_div(p, q)
quotient: x + 3, remainder: 4
> roots(p)
-1
> diff(f)
cos(x) + 2 * x
> diff(f, 2)
3.58385316348
> integral(f)
-cos(x) + x ^ 3 / 3 + C
> taylor(f, 0, 3)
x + x ^ 2 - 0.166666666667 * x ^ 3
> h(x) = sin(x) / x
h(x) = sin(x) / x
> limit(h, 0)
1
> limit(h, 0, 1)
1
> integral(f, 0, 3)
10.9899924966
> integral(f, 3)
10.9899924966
> extrema(f, -2, 2)
min: x = -0.450183689594, f(x) = -0.232465575151
```

Use `:funcs` to list custom functions, `:clearfunc name` to remove one, and
`:clearfuncs` to clear all of them.

`poly_add/poly_sub/poly_mul/poly_div/roots` operate on custom functions that
can be recognized as one-variable polynomials. `roots(p)` returns real roots
only. `diff(f)` returns a symbolic derivative expression, while `diff(f, value)`
returns the numeric derivative at a point. `integral(f)` returns a symbolic
indefinite integral, while `integral(...)` with more arguments keeps the
existing numeric integration behavior. Symbolic `diff(...)` and `integral(...)`
also accept raw one-variable expressions such as `diff(x ^ 2)` or
`integral(1 / x)`, and infer the variable automatically. Current symbolic rules
cover common algebraic powers, `sin/cos/tan`, `exp/ln`, `sqrt/cbrt`, `abs`, the
inverse trigonometric basics `asin/acos/atan`, and some polynomial-times-
`exp/sin/cos` cases such as `integral(x * exp(x))`. `taylor(f, a, n)` returns
the Taylor polynomial around `a` up to degree `n`. `ode(rhs, x0, y0, x1[, steps])`
solves the first-order initial value problem `y' = rhs(x, y)`, while
`ode_system(rhs_vec, x0, y0_vec, x1[, steps])` solves nonlinear ODE systems
using state names such as `y1`, `y2`, and `y3`. `ode_table(...)` and
`ode_system_table(...)` return sampled trajectories. Planning helpers
`lp_max/lp_min` and `ilp_max/ilp_min` solve box-constrained linear and integer
programs of the form `A * x <= b` with explicit lower and upper bounds.
`double_integral(...)` / `triple_integral(...)` compute Cartesian
multi-integrals, `double_integral_cyl(...)` / `triple_integral_cyl(...)` use
cylindrical coordinates, and `triple_integral_sph(...)` uses spherical
coordinates.

## Scripting

The calculator can also execute small scripts. Use `:run file.calc` in the REPL
or redirect a file into stdin with `./calculator < file.calc`.

Runnable script-related inputs are provided in `test/script/`:

- `test/script/test_variables.calc`
  Non-interactive variable and expression workflow
- `test/script/test_function.calc`
  Script functions and return values
- `test/script/test_control_flow.calc`
  `if`, `for`, and `while` control flow examples
- `test/script/test_matrix_basic.calc`
  Matrix creation, analysis, and solving via redirected input
- `test/script/SYNTAX_GUIDE.md`
  Detailed scripting syntax notes

Example commands:

```bash
./calculator < test/script/test_variables.calc
./calculator < test/script/test_function.calc
./calculator < test/script/test_matrix_basic.calc
```

For a dedicated syntax and behavior guide, see `test/script/SYNTAX_GUIDE.md`.

Supported script constructs:

- `fn name(args) { ... }`
- `if (cond) { ... } else { ... }`
- `while (cond) { ... }`
- `for (init; cond; step) { ... }`
- `return expr;`
- `break;`
- `continue;`
- string literals such as `"hello"`
- `print(a, b, c);`

Example:

```text
fn fact(n) {
  if (n <= 1) {
    return 1;
  } else {
    return n * fact(n - 1);
  }
}

x = 0;
sum = 0;
while (x < 5) {
  sum = sum + x;
  x = x + 1;
}

print(sum);
fact(5);
```

`save/load` now persists scalar variables, string variables, expression-style
custom functions, and script functions. Matrix variables are still excluded.

## Matrices

Matrix support is available directly in expressions.

```text
> [1, 2; 3]
[[1, 2], [3, 0]]
> m = mat(2, 2, 1, 2, 3, 4)
m = [[1, 2], [3, 4]]
> m + eye(2)
[[2, 2], [3, 5]]
> get(m, 1, 0)
3
> m = set(m, 1, 0, 8)
m = [[1, 2], [8, 4]]
> set([1, 2; 3], 1, 2, 9)
[[1, 2, 0], [3, 0, 9]]
> append_row([1, 2], 3)
[[1, 2], [3, 0]]
> append_col([1; 2], 3, 4, 5)
[[1, 3], [2, 4], [0, 5]]
> transpose(m)
[[1, 8], [2, 4]]
> inverse(m)
[[-0.5, 0.25], [1, -0.125]]
> m ^ -1
[[-0.5, 0.25], [1, -0.125]]
> null(mat(2, 3, 1, 2, 3, 2, 4, 6))
[[-2, -3], [1, 0], [0, 1]]
> svd_s(mat(3, 2, 3, 0, 0, 2, 0, 0))
[[3, 0], [0, 2]]
> solve(mat(2, 2, 2, 1, 5, 3), vec(1, 2))
[[1], [-1]]
> eigvals(mat(2, 2, 2, 0, 0, 3))
[3, 2]
```

Indices are zero-based. For a fuller matrix-specific reference, see
`MATRIX_GUIDE.md`.

`solve(A, b)` currently expects a square coefficient matrix and a vector-shaped
right-hand side.

Matrix variables are currently displayable and usable in expressions, but they
cannot yet be persisted through `:save` / `:load`.

## Variable Management

You can inspect or clear stored variables with command-style inputs:

```text
> :vars
x = 7/12
y = 1
> :clear x
Cleared variable: x
> :clear
Cleared all variables.
```

## History

The current session history can be inspected with:

```text
> :history
1. :exact on
2. x = 1/3 + 1/4
3. x + 1/6
```

## Save And Load

You can persist variables to a file and restore them in a later session:

```text
> :save state.txt
Saved variables to: state.txt
> :clear
Cleared all variables.
> :load state.txt
Loaded variables from: state.txt
```

## Notes

- Trigonometric functions use radians.
- The implementation uses `double` and does not aim for arbitrary precision.
- Exact fraction mode is intended for rational arithmetic and integer helpers.
- Domain errors such as `ln(0)` and division by zero are reported at runtime.
