# Keywords and Reserved Names

This document lists all keywords and reserved names that cannot be used as
variable or function names in calculator expressions and scripts.

## Script Keywords

These keywords are used by the script parser and cannot be used as variable
names in `.calc` scripts:

| Keyword | Purpose |
|---------|---------|
| `def` | Function definition |
| `fn` | Function definition (alternative) |
| `if` | Conditional statement |
| `elif` | Else-if branch |
| `else` | Else branch |
| `while` | While loop |
| `for` | For loop |
| `in` | For-in iteration |
| `return` | Function return |
| `break` | Break out of loop |
| `continue` | Continue to next iteration |
| `pass` | No-op statement |

## Reserved Function Names

These names are reserved for built-in functions and cannot be used for
user-defined functions:

### Mathematical Functions

| Name | Description |
|------|-------------|
| `abs` | Absolute value |
| `acos` | Arc cosine |
| `acosh` | Inverse hyperbolic cosine |
| `asin` | Arc sine |
| `asinh` | Inverse hyperbolic sine |
| `atan` | Arc tangent |
| `atanh` | Inverse hyperbolic tangent |
| `cos` | Cosine |
| `cosh` | Hyperbolic cosine |
| `exp` | Exponential |
| `ln` | Natural logarithm |
| `log` | Logarithm |
| `pow` | Power |
| `sin` | Sine |
| `sinh` | Hyperbolic sine |
| `sqrt` | Square root |
| `tan` | Tangent |
| `tanh` | Hyperbolic tangent |

### Rounding and Integer Functions

| Name | Description |
|------|-------------|
| `ceil` | Ceiling |
| `floor` | Floor |
| `round` | Round to nearest |
| `trunc` | Truncate |

### Matrix Functions

| Name | Description |
|------|-------------|
| `cholesky` | Cholesky decomposition |
| `cond` | Condition number |
| `det` | Determinant |
| `diag` | Diagonal matrix |
| `eig` | Eigenvalues and eigenvectors |
| `eigvals` | Eigenvalues |
| `eigvecs` | Eigenvectors |
| `eye` | Identity matrix |
| `hessenberg` | Hessenberg decomposition |
| `identity` | Identity matrix |
| `inverse` | Matrix inverse |
| `kron` | Kronecker product |
| `least_squares` | Least squares solution |
| `lu_l` | LU decomposition (L matrix) |
| `lu_u` | LU decomposition (U matrix) |
| `norm` | Matrix/vector norm |
| `null` | Null space |
| `outer` | Outer product |
| `pinv` | Pseudo-inverse |
| `qr_q` | QR decomposition (Q matrix) |
| `qr_r` | QR decomposition (R matrix) |
| `rank` | Matrix rank |
| `rref` | Reduced row echelon form |
| `schur` | Schur decomposition |
| `solve` | Solve linear system |
| `svd` | Singular value decomposition |
| `svd_s` | SVD singular values |
| `svd_u` | SVD U matrix |
| `svd_vt` | SVD V^T matrix |
| `trace` | Matrix trace |
| `transpose` | Matrix transpose |

### Calculus and Analysis

| Name | Description |
|------|-------------|
| `bisect` | Bisection method |
| `critical` | Critical points |
| `curl` | Curl |
| `curl_2d` | 2D curl |
| `diff` | Derivative |
| `directional` | Directional derivative |
| `div` | Divergence |
| `divergence` | Divergence |
| `dsolve` | Symbolic ODE solver |
| `extrema` | Extrema finding |
| `fixed_point` | Fixed point iteration |
| `gradient` | Gradient |
| `greens_theorem` | Green's theorem |
| `hessian` | Hessian matrix |
| `implicit_diff` | Implicit differentiation |
| `integral` | Definite integral |
| `jacobian` | Jacobian matrix |
| `lagrange` | Lagrange multipliers |
| `laplacian` | Laplacian |
| `limit` | Limit |
| `line_integral_scalar` | Scalar line integral |
| `line_integral_vector` | Vector line integral |
| `ode` | ODE solver |
| `ode_system` | ODE system solver |
| `ode_table` | ODE table output |
| `ode_system_table` | ODE system table output |
| `param_deriv` | Parametric derivative |
| `roots` | Polynomial roots |
| `secant` | Secant method |
| `solve` | Solve equation or system |
| `stokes_theorem` | Stokes' theorem |
| `surface_integral_flux` | Surface flux integral |
| `surface_integral_scalar` | Scalar surface integral |
| `divergence_theorem` | Divergence theorem |

### Transforms

| Name | Description |
|------|-------------|
| `dft` | Discrete Fourier transform |
| `fft` | Fast Fourier transform |
| `fourier` | Fourier transform |
| `idft` | Inverse DFT |
| `ifft` | Inverse FFT |
| `ifourier` | Inverse Fourier transform |
| `ilaplace` | Inverse Laplace transform |
| `iztrans` | Inverse Z transform |
| `laplace` | Laplace transform |
| `ztrans` | Z transform |

### Series and Approximation

| Name | Description |
|------|-------------|
| `pade` | Pade approximant |
| `puiseux` | Puiseux series |
| `series_sum` | Series summation |
| `summation` | Summation |
| `taylor` | Taylor series |

### Root Finding and Solving

| Name | Description |
|------|-------------|
| `bisect` | Bisection method |
| `fixed_point` | Fixed point iteration |
| `roots` | Polynomial roots |
| `secant` | Secant method |
| `solve` | Solve equation or system |

### Polynomial Functions

| Name | Description |
|------|-------------|
| `poly_add` | Polynomial addition |
| `poly_compose` | Polynomial composition |
| `poly_deriv` | Polynomial derivative |
| `poly_div` | Polynomial division |
| `poly_eval` | Polynomial evaluation |
| `poly_gcd` | Polynomial GCD |
| `poly_integ` | Polynomial integral |
| `poly_mul` | Polynomial multiplication |
| `poly_sub` | Polynomial subtraction |

### Statistics

| Name | Description |
|------|-------------|
| `avg` | Average |
| `corr` | Correlation |
| `cov` | Covariance |
| `lagrange` | Lagrange interpolation |
| `linear_regression` | Linear regression |
| `max` | Maximum |
| `mean` | Mean |
| `median` | Median |
| `min` | Minimum |
| `mode` | Mode |
| `percentile` | Percentile |
| `poly_fit` | Polynomial fit |
| `quartile` | Quartile |
| `spline` | Spline interpolation |
| `std` | Standard deviation |
| `var` | Variance |

### Combinatorics

| Name | Description |
|------|-------------|
| `combination` | Combination count |
| `factorial` | Factorial |
| `nCr` | Combinations |
| `nPr` | Permutations |

### Number Theory

| Name | Description |
|------|-------------|
| `divisors` | All divisors |
| `extended_gcd` | Extended Euclidean algorithm |
| `xgcd` | Extended Euclidean algorithm (alias) |

### Random Numbers

| Name | Description |
|------|-------------|
| `rand` | Random number |
| `randint` | Random integer |
| `randn` | Normal random |

### Complex Numbers

| Name | Description |
|------|-------------|
| `abs` | Absolute value / magnitude |
| `arg` | Argument (phase angle) |
| `complex` | Create complex number |
| `conj` | Complex conjugate |
| `imag` | Imaginary part |
| `polar` | Create from polar coordinates |
| `real` | Real part |

### Window Functions

| Name | Description |
|------|-------------|
| `blackman` | Blackman window |
| `hamming` | Hamming window |
| `hann` | Hanning window |
| `hanning` | Hanning window (alias) |

### Vector Operations

| Name | Description |
|------|-------------|
| `append_col` | Append column |
| `append_row` | Append row |
| `cross` | Cross product |
| `dot` | Dot product |
| `get` | Get element |
| `hadamard` | Hadamard product |
| `mat` | Create matrix |
| `resize` | Resize matrix |
| `set` | Set element |
| `vec` | Create vector |
| `zeros` | Zero matrix |

### Signal Processing

| Name | Description |
|------|-------------|
| `conv` | Convolution |
| `convolve` | Convolution |
| `dft` | Discrete Fourier transform |
| `fft` | Fast Fourier transform |
| `freqz` | Frequency response |
| `idft` | Inverse DFT |
| `ifft` | Inverse FFT |
| `residue` | Partial fraction residues |

### Plotting

| Name | Description |
|------|-------------|
| `plot` | Plot function |

### Linear/Integer Programming

| Name | Description |
|------|-------------|
| `bip_max` | Binary integer programming max |
| `bip_min` | Binary integer programming min |
| `ilp_max` | Integer linear programming max |
| `ilp_min` | Integer linear programming min |
| `lp_max` | Linear programming max |
| `lp_min` | Linear programming min |
| `milp_max` | Mixed integer linear programming max |
| `milp_min` | Mixed integer linear programming min |

### Multi-Variable Integration

| Name | Description |
|------|-------------|
| `double_integral` | Double integral |
| `double_integral_cyl` | Double integral (cylindrical) |
| `double_integral_polar` | Double integral (polar) |
| `triple_integral` | Triple integral |
| `triple_integral_cyl` | Triple integral (cylindrical) |
| `triple_integral_sph` | Triple integral (spherical) |

## Built-in Constants

These constants are predefined and should not be reassigned:

| Constant | Value |
|----------|-------|
| `pi` | π ≈ 3.14159... |
| `e` | e ≈ 2.71828... |
| `i` | Imaginary unit (in expressions like `3+4i`) |

## Variable Naming Rules

Valid variable names must satisfy:

1. First character must be a letter (`A-Z`, `a-z`) or underscore (`_`)
2. Subsequent characters can be letters, digits (`0-9`), or underscores
3. Cannot be a keyword or reserved function name

### Examples

**Valid:**
- `x`, `y`, `z`
- `count`, `total_sum`
- `_temp`, `value1`

**Invalid:**
- `1var` (starts with digit)
- `my-var` (contains hyphen)
- `if`, `while`, `for` (keywords)
- `sin`, `cos`, `det` (reserved functions)

## Scope of Restrictions

| Context | Keywords | Reserved Functions | Constants |
|---------|----------|-------------------|-----------|
| Script variable names | Blocked | Warning | Warning |
| Script function names | Blocked | Blocked | Blocked |
| Function parameters | Should block | Should block | Should block |
| REPL variable names | Allowed | Warning | Warning |
| Expression-style functions | Allowed | Blocked | Blocked |

Note: Some restrictions are enforced at parse time, others at runtime.
The behavior may vary slightly between contexts.
