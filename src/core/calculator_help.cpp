#include "calculator.h"

#include <stdexcept>
#include <string>

namespace {

std::string build_help_topic(const std::string& topic) {
    if (topic == "commands") {
        return
        "Commands:\n"
        "  help, :help         Show this help message\n"
        "  :exact on|off       Toggle exact fraction mode\n"
        "  :exact              Show current exact mode status\n"
        "  :symbolic on|off    Preserve pi/e in scalar display results\n"
        "  :symbolic           Show current symbolic constants mode\n"
        "  :hexprefix on|off   Toggle 0x/0X prefix for hex output\n"
        "  :hexprefix          Show current hex prefix mode\n"
        "  :hexcase upper|lower Set hex letter case\n"
        "  :hexcase            Show current hex letter case\n"
        "  :vars               List stored variables\n"
        "  :funcs              List custom functions\n"
        "  :clear name         Clear one variable\n"
        "  :clearfunc name     Clear one custom function\n"
        "  :clear              Clear all variables\n"
        "  :clearfuncs         Clear all custom functions\n"
        "  :history            Show session input history\n"
        "  :save file          Save variables to a file\n"
        "  :load file          Load variables from a file\n"
        "  :run file.calc      Execute a .calc script file\n"
        "  exit, quit          Exit the calculator";
    }

    if (topic == "examples") {
        return
        "Examples:\n"
        "  x = 1/3 + 1/4       Assign a variable\n"
        "  2 ^ 10              Power operator\n"
        "  pow(3, 4)           Function-style power\n"
        "  v = vec(1, 2, 3)    Create a vector\n"
        "  m = mat(2, 2, 1, 2, 3, 4)  Create a matrix\n"
        "  m + eye(2)          Add matrices directly\n"
        "  2 * m               Matrix-scalar multiplication\n"
        "  transpose(m)        Matrix transpose\n"
        "  inverse(m)          Matrix inverse\n"
        "  dot(a, b)           Vector dot product\n"
        "  outer(a, b)         Vector outer product\n"
        "  null(m)             Nullspace basis\n"
        "  least_squares(A, b) Least-squares solution\n"
        "  qr_q(A), qr_r(A)    QR decomposition parts\n"
        "  lu_l(A), lu_u(A)    LU decomposition parts\n"
        "  svd_u/s/vt(A)       Reduced SVD factors\n"
        "  solve(A, b)         Solve Ax = b\n"
        "  get(m, 1, 0)        Read one element\n"
        "  m = set(m, 1, 0, 8) Update one element\n"
        "  det(m)              Determinant\n"
        "  rref(m)             Reduced row echelon form\n"
        "  resize(m, 3, 3)     Resize with zero-fill\n"
        "  factor(360)         Prime factorization\n"
        "  factorial(5)        Integer factorial\n"
        "  nCr(5, 2)           Combination count\n"
        "  nPr(5, 2)           Permutation count\n"
        "  rat(pi, 1000)       Best rational approximation with bounded denominator\n"
        "  round(2.6)          Round to nearest integer\n"
        "  clamp(12, 0, 10)    Clamp into a closed interval\n"
        "  log(8, 2)           Logarithm with an arbitrary base\n"
        "  exp2(5)             Power of two\n"
        "  percentile(75, 1, 2, 3, 4)  Scalar percentile\n"
        "  quartile(vec(1, 2, 3, 4), 1)  Vector quartile\n"
        "  rol(1, 3)           Rotate bits to the left\n"
        "  popcount(15)        Count set bits\n"
        "  sum(1, 2, 3, 4)     Aggregate sum\n"
        "  mean(1, 2, 3, 4)    Aggregate mean\n"
        "  avg(1, 2, 3, 4)     Aggregate average\n"
        "  median(1, 5, 2, 9)  Aggregate median\n"
        "  mode(1, 2, 2, 3)    Aggregate mode\n"
        "  sinh(1)             Hyperbolic sine\n"
        "  asinh(1)            Inverse hyperbolic sine\n"
        "  sec(pi/3)           Reciprocal trig function\n"
        "  gamma(5)            Gamma function\n"
        "  erf(1)              Error function\n"
        "  beta(2, 3)          Beta function\n"
        "  zeta(2)             Riemann zeta function\n"
        "  fib(10)             Fibonacci number\n"
        "  is_prime(17)        Prime test\n"
        "  deg2rad(180)        Angle conversion\n"
        "  sin_deg(30)         Degree-based sine\n"
        "  fahrenheit(25)      Celsius to Fahrenheit\n"
        "  c2f(100)            Celsius to Fahrenheit alias\n"
        "  f(x) = sin(x)+x^2   Define a custom unary function\n"
        "  f(2)                Evaluate a custom function\n"
        "  :symbolic on        Preserve pi/e symbolically in scalar output\n"
        "  :hexprefix on       Show hex results like 0xFF\n"
        "  :hexcase lower      Show hex results like 0xff\n"
        "  pi / 2 + e          Symbolic constants mode example\n"
        "  :run demo.calc      Run a script file\n"
        "  def fact(n):        Define a script function in a script\n"
        "  print(a, b, c)      Print script values, including strings\n"
        "  poly_add(p, q)      Polynomial addition\n"
        "  poly_sub(p, q)      Polynomial subtraction\n"
        "  poly_mul(p, q)      Polynomial multiplication\n"
        "  poly_div(p, q)      Polynomial division\n"
        "  roots(p)            Real roots of a polynomial\n"
        "  simplify(x^2 + x^2) Simplify a symbolic expression\n"
        "  simplify(expr)      Simplify a symbolic expression\n"
        "  diff(f)             Symbolic derivative expression\n"
        "  diff(f, 2)          Derivative at x = 2\n"
        "  integral(f)         Symbolic indefinite integral expression\n"
        "  integral(f, x, y)   Chained symbolic integral over x then y\n"
        "  critical(f, x, y)   Symbolic critical point search\n"
        "  step(t - 1)         Unit step / Heaviside function\n"
        "  delta(t - 1)        Unit impulse / Dirac delta shorthand\n"
        "  laplace(exp(-2*t), t, s)  Symbolic Laplace transform\n"
        "  ilaplace(1 / (s + 2), s, t)  Symbolic inverse Laplace transform\n"
        "  fourier(delta(t - 1), t, w)  Symbolic Fourier transform\n"
        "  ifourier(delta(w - 3), w, t) Symbolic inverse Fourier transform\n"
        "  ztrans(step(n), n, z)  Symbolic z transform\n"
        "  iztrans(z / (z - 1), z, n)  Symbolic inverse z transform\n"
        "  dft([1, 0, 0, 0])   Discrete Fourier transform\n"
        "  idft([[1, 0], [1, 0], [1, 0], [1, 0]])  Inverse DFT\n"
        "  convolve([1, 2], [3, 4, 5])  Linear convolution\n"
        "  pade(exp(x), 0, 2, 2)  Pade approximant around a point\n"
        "  puiseux((1 + x) ^ (1 / 2), 0, 4, 2)  Puiseux-style local series\n"
        "  series_sum(n^2, n, 1, N)  Symbolic finite series sum\n"
        "  taylor(f, 0, 5)     Taylor expansion up to degree 5\n"
        "  limit(f, 0)         Two-sided limit as x -> 0\n"
        "  integral(f, 0, 3)   Definite integral on [0, 3]\n"
        "  integral(f, 3)      Indefinite integral value at x = 3\n"
        "  double_integral(x + y, 0, 1, 0, 2)  Cartesian double integral\n"
        "  double_integral_cyl(x^2 + y^2, 0, 1, 0, 2*pi)  Polar/cylindrical double integral\n"
        "  triple_integral(x*y*z, 0, 1, 0, 1, 0, 1)  Cartesian triple integral\n"
        "  triple_integral_sph(1, 0, 1, 0, 2 * pi, 0, pi)  Spherical triple integral\n"
        "  ode(y - x, 0, 1, 2) Solve y' = y - x with y(0) = 1\n"
        "  ode(p1 * y, 0, 1, 1, 50, mat(1, 1, 2))  ODE with parameter vector p\n"
        "  ode_table(y, 0, 1, 1, 4)  Return sampled ODE trajectory\n"
        "  ode_system(vec(y2, -y1), 0, vec(0, 1), 1)  Solve a nonlinear ODE system\n"
        "  ode_system_table(mat(1, 1, p1 * y1), 0, mat(1, 1, 1), 2, 40, y1 - 2, mat(1, 1, 1))  Stop on an event\n"
        "  lp_max(vec(3, 2), mat(3, 2, 1, 1, 1, 0, 0, 1), vec(4, 2, 3), vec(0, 0), vec(10, 10))  Box-constrained linear planning\n"
        "  ilp_max(vec(3, 2), mat(3, 2, 1, 1, 1, 0, 0, 1), vec(4, 2, 3), vec(0, 0), vec(10, 10))  Integer planning\n"
        "  milp_max(vec(3, 1), mat(1, 2, 2, 1), vec(5), vec(0, 0), vec(2, 10), vec(1, 0))  Mixed-integer planning\n"
        "  bip_max(vec(5, 4, 3), mat(1, 3, 2, 1, 1), vec(2))  Binary planning shortcut\n"
        "  solve(x^2 - 2, 1)   Newton root solve\n"
        "  bisect(x^2 - 2, 1, 2)  Bisection root solve\n"
        "  extrema(f, -2, 2)   Solve extrema on an interval\n"
        "  :run script.calc    Execute a script file\n"
        "  root(27, 3)         General root\n"
        "  cbrt(-8)            Cube root\n"
        "  hex(255)            Base conversion\n"
        "  and(6, 3)           Bitwise and\n"
        "  min(7/3, 5/2)       Smaller of two values\n"
        "  :save state.txt     Save variables";
    }

    if (topic == "matrix") {
        return
        "Matrix guide:\n"
        "  Create:  [a,b;c,d] vec mat zeros eye identity\n"
        "  Shape:   resize append_row append_col transpose\n"
        "  Elem:    get set\n"
        "  Extra:   inverse dot outer null least_squares qr_q qr_r lu_l lu_u svd_u svd_s svd_vt pinv kron hadamard\n"
        "  Ops:     + - * / ^ with scalars and matrices\n"
        "  Anal.:   norm trace det rank rref eigvals eigvecs solve cond diag reshape cholesky schur hessenberg percentile quartile\n"
        "  Signal:  dft fft idft ifft conv convolve\n"
        "  Notes:   indices are zero-based\n"
        "  Notes:   dft/idft accept a real vector or an N x 2 complex matrix\n"
        "  Notes:   matrix literals pad missing elements with 0\n"
        "  Notes:   append_row/append_col also pad or expand with 0\n"
        "  Example: m = mat(2, 2, 1, 2, 3, 4)\n"
        "  Example: [1, 2; 3]\n"
        "  Example: get(m, 1, 0)\n"
        "  Example: m = set(m, 1, 0, 8)\n"
        "  Example: append_row([1, 2], 3)\n"
        "  Example: transpose(m)\n"
        "  Example: inverse(m)\n"
        "  Example: dot(vec(1, 2), vec(3, 4))\n"
        "  Example: lu_l(mat(2, 2, 4, 3, 6, 3))\n"
        "  Example: svd_s(mat(3, 2, 3, 0, 0, 2, 0, 0))\n"
        "  Example: solve(mat(2, 2, 2, 1, 5, 3), vec(1, 2))\n"
        "  Example: det(m)\n"
        "  Example: rref(m)";
    }

    if (topic == "symbolic") {
        return
        "Symbolic tools:\n"
        "  Simplify:       simplify(expr)\n"
        "  Derivative:     diff(expr), diff(expr, x), diff(expr, x, y)\n"
        "  Integral:       integral(expr), integral(expr, x), integral(expr, x, y)\n"
        "  Vector calc:    gradient(expr, x, y), jacobian([f; g], x, y), hessian(expr, x, y)\n"
        "  Critical:       critical(expr, x, y) with Hessian classification when available\n"
        "  Series:         taylor(expr, a, n), pade(expr, m, n), pade(expr, a, m, n), puiseux(expr, a, degree, denominator)\n"
        "  Sums:           series_sum(expr, n, lower, upper), summation(expr, n, lower, upper)\n"
        "  Signals:        step(t), delta(t), heaviside(t), impulse(t)\n"
        "  Transforms:     laplace(expr, t, s), ilaplace(expr, s, t), fourier(expr, t, w), ifourier(expr, w, t)\n"
        "  Z-transform:    ztrans(expr, n, z), iztrans(expr, z, n)\n"
        "  Examples:       integral(1 / (s + 2), s), ilaplace(1 / (s + 2), s, t)\n"
        "  Examples:       gradient(x ^ 2 + y ^ 2, x, y), critical(x ^ 2 - y ^ 2, x, y)\n"
        "  Notes:          symbolic integration is rule-based, not a full Risch integrator";
    }

    if (topic == "analysis") {
        return
        "Analysis and solving:\n"
        "  Roots:          solve(expr, guess), bisect(expr, left, right), secant(expr, x0, x1), fixed_point(expr, x0)\n"
        "  Limits:         limit(expr, x0), limit(expr, x0, direction)\n"
        "  Extrema:        extrema(f, left, right[, scan_segments])\n"
        "  Integrals:      integral(f, x0), integral(f, a, b)\n"
        "  Multi-var:      double_integral, double_integral_polar, double_integral_cyl\n"
        "  Multi-var:      triple_integral, triple_integral_cyl, triple_integral_sph\n"
        "  ODE scalar:     ode(rhs, x0, y0, x1[, steps]), ode_table(rhs, x0, y0, x1[, steps])\n"
        "  ODE system:     ode_system(rhs_vec, x0, y0_vec, x1[, steps]), ode_system_table(...)\n"
        "  Events/params:  ODE commands accept optional parameter vectors and event expressions\n"
        "  Examples:       solve(x ^ 2 - 2, 1), limit(sin(x) / x, 0)\n"
        "  Examples:       ode(y - x, 0, 1, 2), ode_system(vec(y2, -y1), 0, vec(0, 1), 1)";
    }

    if (topic == "planning") {
        return
        "Planning and optimization:\n"
        "  Linear:         lp_max(c, A, b, lo, hi), lp_min(c, A, b, lo, hi)\n"
        "  Integer:        ilp_max(c, A, b, lo, hi), ilp_min(c, A, b, lo, hi)\n"
        "  Mixed integer:  milp_max(c, A, b, lo, hi, integrality), milp_min(...)\n"
        "  Binary:         bip_max(c, A, b), bip_min(c, A, b)\n"
        "  Binary aliases: binary_max(c, A, b), binary_min(c, A, b)\n"
        "  Optional args:  equality constraints can be supplied as Aeq and beq\n"
        "  Examples:       lp_max(vec(3, 2), mat(3, 2, 1, 1, 1, 0, 0, 1), vec(4, 2, 3), vec(0, 0), vec(10, 10))\n"
        "  Examples:       bip_max(vec(5, 4, 3), mat(1, 3, 2, 1, 1), vec(2))";
    }

    if (topic == "functions") {
        return
        "Common functions:\n"
        "  Trigonometric: sin cos tan sec csc cot asin acos atan asec acsc acot sinh cosh tanh asinh acosh atanh\n"
        "  Exponential:   exp exp2 ln log log2 log10 pow gamma beta zeta erf erfc bessel\n"
        "  Roots:         sqrt cbrt root\n"
        "  Numeric:       abs sign floor ceil round trunc min max clamp sum mean avg median mode percentile quartile var std factorial nCr binom nPr fib is_prime next_prime rand randn randint\n"
        "  Legacy nums:   factorial nCr nPr\n"
        "  Signals:       step delta heaviside impulse fourier ifourier laplace ilaplace ztrans iztrans\n"
        "  Discrete sig:  dft fft idft ifft conv convolve\n"
        "  Series:        taylor pade puiseux series_sum summation\n"
        "  Convert:       deg rad deg2rad rad2deg sin_deg cos_deg celsius fahrenheit kelvin c2f f2c\n"
        "  Legacy conv:   deg2rad rad2deg celsius fahrenheit kelvin\n"
        "  Matrix create: vec mat zeros eye identity\n"
        "  Matrix shape:  resize append_row append_col transpose\n"
        "  Matrix elem:   get set\n"
        "  Matrix extra:  inverse dot outer null least_squares qr_q qr_r lu_l lu_u svd_u svd_s svd_vt pinv kron hadamard\n"
        "  Matrix ops:    + - * / ^ with scalars and matrices\n"
        "  Matrix anal.:  norm trace det rank rref eigvals eigvecs solve cond diag reshape cholesky schur hessenberg\n"
        "  Integer:       gcd lcm mod factor\n"
        "  Base convert:  bin oct hex base\n"
        "  Aggregate:     sum mean avg median mode percentile quartile var std cov corr\n"
        "  Legacy aggr:   sum avg median\n"
        "  Multi-var:     double_integral double_integral_cyl double_integral_polar triple_integral triple_integral_cyl triple_integral_sph\n"
        "  Bitwise:       and or xor not shl shr rol ror popcount bitlen ctz clz parity reverse_bits\n"
        "  Script:        def if elif else while for range return break continue print strings\n"
        "  Custom:        f(x)=...  poly_add poly_sub poly_mul poly_div roots poly_eval poly_deriv poly_integ poly_fit poly_compose poly_gcd\n"
        "  Symbolic:      simplify diff integral gradient jacobian hessian critical taylor pade puiseux series_sum summation\n"
        "  Analysis:      limit extrema ode ode_table ode_system ode_system_table "
        "lp_max lp_min ilp_max ilp_min milp_max milp_min bip_max bip_min solve bisect secant fixed_point eig svd";
    }

    if (topic == "exact") {
        return
        "Exact mode:\n"
        "  :exact on           Prefer rational results like 7/12 over decimals\n"
        "  :exact off          Return to normal decimal-first display\n"
        "  :exact              Show the current exact mode status\n"
        "  Works best with:    + - * / pow integer-exponent min max clamp sum avg median mean factorial nCr binom nPr\n"
        "  Integer helpers:    gcd lcm mod and programmer bitwise helpers stay exact when possible\n"
        "  Falls back to decimal display for non-rational functions like sin, cos, exp, ln, sqrt, percentile";
    }

    if (topic == "variables") {
        return
        "Variables and functions:\n"
        "  x = 3/4             Assign a scalar variable\n"
        "  v = vec(1, 2, 3)    Assign a matrix/vector value\n"
        "  :vars               List all stored variables\n"
        "  :clear x            Clear one variable\n"
        "  :clear              Clear all variables\n"
        "  f(x) = x^2 + 1      Define a custom expression function\n"
        "  :funcs              List custom expression/script functions\n"
        "  :clearfunc f        Clear one custom function\n"
        "  :clearfuncs         Clear all custom functions\n"
        "  Custom functions are available in expressions, analysis commands, and scripts";
    }

    if (topic == "persistence") {
        return
        "Persistence:\n"
        "  :save state.txt     Save current scalar variables and custom functions\n"
        "  :load state.txt     Load a previously saved state file\n"
        "  Save/load keeps:    scalar variables, exact values, string values, custom functions, script functions\n"
        "  Current limit:      matrix variables are not saved yet\n"
        "  Tip:                use separate files for different sessions or experiments";
    }

    if (topic == "programmer") {
        return
        "Programmer tools:\n"
        "  Base convert:       bin oct hex base\n"
        "  Bitwise:            and or xor not shl shr rol ror\n"
        "  Bit metrics:        popcount bitlen ctz clz parity reverse_bits\n"
        "  Hex formatting:     :hexprefix on|off, :hexcase upper|lower\n"
        "  Examples:           hex(255), base(255, 16), and(6, 3), shl(5, 2), rol(1, 3), popcount(15)\n"
        "  Notes:              bit helpers use 64-bit two's-complement integer semantics\n"
        "  Notes:              base conversion only accepts integers and bases 2..16\n"
        "  Notes:              hex formatting applies to hex(...) and base(..., 16)";
    }

    throw std::runtime_error("unknown help topic: " + topic);
}

}  // namespace

std::string Calculator::help_text() const {
    return
        "Help topics:\n"
        "  :help commands      Show command reference\n"
        "  :help functions     Show supported functions\n"
        "  :help matrix        Show matrix usage guide\n"
        "  :help symbolic      Show symbolic algebra and transform help\n"
        "  :help analysis      Show calculus, root solving, and ODE help\n"
        "  :help planning      Show linear/integer planning help\n"
        "  :help examples      Show example inputs\n"
        "  :help exact         Show exact fraction mode help\n"
        "  :help variables     Show variable and function usage help\n"
        "  :help persistence   Show save/load help\n"
        "  :help programmer    Show bitwise/base-conversion help\n"
        "\n" +
        build_help_topic("commands") + "\n\n" +
        build_help_topic("matrix");
}

std::string Calculator::help_topic(const std::string& topic) const {
    return build_help_topic(topic);
}
