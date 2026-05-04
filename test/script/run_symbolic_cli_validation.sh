#!/usr/bin/env bash
set -u

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
CALC="$ROOT_DIR/bin/calculator"

passed=0
failed=0

run_calc() {
    printf ":symbolic on\n%s\n" "$1" | "$CALC" 2>&1 | tail -n 1
}

expect_exact() {
    local expr="$1"
    local expected="$2"
    local actual
    actual="$(run_calc "$expr")"
    if [[ "$actual" == "$expected" ]]; then
        passed=$((passed + 1))
    else
        failed=$((failed + 1))
        printf 'FAIL exact: %s\n  expected: %s\n  actual:   %s\n' "$expr" "$expected" "$actual"
    fi
}

expect_one_of() {
    local expr="$1"
    shift
    local actual
    actual="$(run_calc "$expr")"
    local expected
    for expected in "$@"; do
        if [[ "$actual" == "$expected" ]]; then
            passed=$((passed + 1))
            return
        fi
    done

    failed=$((failed + 1))
    printf 'FAIL one-of: %s\n  expected one of:\n' "$expr"
    for expected in "$@"; do
        printf '    %s\n' "$expected"
    done
    printf '  actual:   %s\n' "$actual"
}

expect_numeric() {
    local expr="$1"
    local expected="$2"
    local tolerance="$3"
    local actual
    actual="$(run_calc "$expr")"
    
    if [[ "$actual" == Error* ]]; then
        failed=$((failed + 1))
        printf 'FAIL numeric (error): %s\n  expected: %s\n  actual:   %s\n' "$expr" "$expected" "$actual"
        return
    fi
    
    local val
    val=$(echo "$actual" | tr -d '[]' | cut -d',' -f1 | xargs)
    
    if awk -v a="$val" -v e="$expected" -v t="$tolerance" 'BEGIN { d = a - e; if (d < 0) d = -d; exit !(d <= t) }'; then
        passed=$((passed + 1))
    else
        failed=$((failed + 1))
        printf 'FAIL numeric: %s\n  expected: %s +/- %s\n  actual:   %s (parsed: %s)\n' "$expr" "$expected" "$tolerance" "$actual" "$val"
    fi
}

expect_error_contains() {
    local expr="$1"
    local expected_msg="$2"
    local actual
    actual="$(run_calc "$expr")"
    if [[ "$actual" == *"$expected_msg"* ]]; then
        passed=$((passed + 1))
    else
        failed=$((failed + 1))
        printf 'FAIL error: %s\n  expected to contain: %s\n  actual:              %s\n' "$expr" "$expected_msg" "$actual"
    fi
}

# =============================================================================
# 1. Symbolic Differentiation and Calculus
# =============================================================================
expect_exact 'diff(sin(x) * exp(x), x)' 'exp(x) * (cos(x) + sin(x))'
expect_exact 'diff(x ^ 2 * y + sin(y), x, y)' '2 * x'
expect_one_of 'integral(cos(3 * x), x)' 'sin(3 * x) / 3 + C' '1/3 * sin(3 * x) + C'
expect_exact 'integral(1 / (x + 2), x)' 'ln(abs(x + 2)) + C'
expect_exact 'integral(tan(x) ^ 2)' 'tan(x) - x + C'
expect_exact 'simplify(sec(x) ^ 2 - tan(x) ^ 2)' '1'
expect_exact 'diff(sin(x) * exp(x), x)' 'exp(x) * (cos(x) + sin(x))'
expect_exact 'diff(x ^ 2 * y + sin(y), x, y)' '2 * x'
expect_one_of 'diff(sin(x) / x)' \
    '(cos(x) * x - sin(x)) / x ^ 2' \
    '(cos(x) * x - sin(x)) / (x * x)'
expect_exact 'expand((x + y) ^ 3)' 'x ^ 3 + 3 * y ^ 2 * x + y ^ 3 + 3 * x ^ 2 * y'
expect_one_of 'integral(cos(3 * x), x)' 'sin(3 * x) / 3 + C' '1/3 * sin(3 * x) + C'
expect_exact 'integral(cos(pi * x), x)' 'sin(pi * x) / pi + C'
expect_one_of 'integral(exp(e * x), x)' 'exp(e * x + -1) + C' 'exp(e * x) / e + C'
expect_exact 'integral(3 * x ^ 2 * exp(x ^ 3), x)' 'exp(x ^ 3) + C'
expect_exact 'integral(1 / (x + 2), x)' 'ln(abs(x + 2)) + C'
expect_exact 'integral(1 / (x ^ 2 - 1))' '1/2 * (-ln(abs(x + 1)) + ln(abs(x - 1))) + C'
expect_one_of 'integral(sin(x) ^ 2)' \
    'x / 2 - sin(2 * x) / 4 + C' \
    '1/2 * x - 1/4 * sin(2 * x) + C'
expect_one_of 'integral(cos(x) ^ 2)' \
    'sin(2 * x) / 4 + x / 2 + C' \
    '1/4 * sin(2 * x) + 1/2 * x + C'
expect_exact 'integral(tan(x) ^ 2)' 'tan(x) - x + C'
expect_exact 'integral(sec(x) * tan(x), x)' 'sec(x) + C'
expect_exact 'integral(x * cos(x), x)' 'cos(x) + sin(x) * x + C'
expect_exact 'integral(x * sin(x), x)' '-(cos(x) * x) + sin(x) + C'
expect_exact 'integral(x ^ 2 * sin(x), x)' '-(cos(x) * x ^ 2) + 2 * (cos(x) + sin(x) * x) + C'
expect_one_of 'integral(x * ln(x), x)' \
    'x ^ 2 * (ln(x) / 2 - 1/4) + C' \
    'x ^ 2 * (1/2 * ln(x) - 1/4) + C'
expect_exact 'integral(x * atan(x), x)' '1/2 * (atan(x) * (x ^ 2 + 1) - x) + C'
expect_one_of 'integral(x / sqrt(1 - x ^ 2), x)' \
    'sqrt(-(x ^ 2) + 1) / -1 + C' \
    '-sqrt(1 - x ^ 2) + C'
expect_exact 'gradient(x ^ 2 + x * y + y ^ 2, x, y)' '[2 * x + y, x + 2 * y]'
expect_exact 'hessian(x ^ 2 + x * y + y ^ 2, x, y)' '[[2, 1], [1, 2]]'
expect_exact 'jacobian([x ^ 2 + y; sin(x * y)], x, y)' '[[2 * x, 1], [cos(x * y) * y, cos(x * y) * x]]'
expect_exact 'critical(x ^ 2 + y ^ 2, x, y)' '[x = 0, y = 0] (local min)'
expect_exact 'critical(x ^ 2 - y ^ 2, x, y)' '[x = 0, y = 0] (saddle)'
expect_exact 'critical(x ^ 4 + y ^ 4, x, y)' '[x = 0, y = 0] (degenerate)'
expect_exact 'taylor(sin(x), 0, 5)' 'x - 1/6 * x ^ 3 + 1/120 * x ^ 5'
# =============================================================================
# 2. Multivariable Calculus
# =============================================================================
expect_exact 'gradient(x ^ 2 + y * z + exp(z), x, y, z)' '[2 * x, z, exp(z) + y]'
expect_exact 'hessian(x * y * z, x, y, z)' '[[0, z, y], [z, 0, x], [y, x, 0]]'
expect_exact 'critical(x ^ 2 + y ^ 2 - 4 * x + 6 * y, x, y)' '[x = 2, y = -3] (local min)'

# =============================================================================
# 3. Numeric Integration & Differentiation
# =============================================================================
expect_numeric 'integral(sin(x) ^ 2, 0, 3.1415926535)' '1.570796' '0.001'
expect_numeric 'diff(sin(x) * exp(x), 0)' '1' '0.0001'
expect_numeric 'integral(exp(-x), 0, 50)' '1' '0.0001'
expect_numeric 'integral(exp(-(x ^ 2)), -inf, inf)' '1.77245' '0.0001'

# =============================================================================
# 4. Limits and Root Finding
# =============================================================================
expect_numeric 'limit(sin(x) / x, 0)' '1' '0.000001'
expect_numeric 'solve(cos(x) - x, 0.5)' '0.739085' '0.0001'
expect_numeric 'limit(sin(x) / x, 0)' 1 0.00000001
expect_numeric 'limit((exp(x) - 1) / x, 0)' 1 0.00000001
expect_numeric 'limit((1 - cos(x)) / x ^ 2, 0)' 0.5 0.00000001
expect_exact 'limit((x ^ 3 - 1) / (x - 1), 1)' '3'
expect_exact 'limit(abs(x) / x, 0, 1)' '1'
expect_exact 'limit(abs(x) / x, 0, -1)' '-1'
expect_numeric 'limit((1 + 1 / x) ^ x, inf)' 2.71828182842 0.00000001
expect_numeric 'limit(x * sin(1 / x), inf)' 1 0.00000001
expect_error_contains 'limit(1 / x, 0)' 'limit did not converge'
expect_numeric 'solve(cos(x) - x, 0.5)' 0.739085133215 0.0000000001
expect_numeric 'bisect(x ^ 3 - x - 2, 1, 2)' 1.5213797068 0.0000000001
expect_numeric 'secant(x ^ 3 - x - 2, 1, 2)' 1.5213797068 0.0000000001
expect_numeric 'fixed_point(cos(x), 0.5)' 0.73908513325 0.0000000001
# =============================================================================
# 5. Matrix and Linear Algebra
# =============================================================================
expect_exact 'det([1, 2, 3; 0, 1, 4; 5, 6, 0])' '1'
expect_exact 'rank([1, 2, 3; 2, 4, 6; 1, 1, 1])' '2'
expect_exact 'pinv([1, 2; 2, 4])' '[[0.04, 0.08], [0.08, 0.16]]'
expect_exact 'eigvals([2, 0; 0, 3])' '[3, 2]'
expect_exact 'inverse([1, 2, 3; 0, 1, 4; 5, 6, 0])' '[[-24, 18, 5], [20, -15, -4], [-5, 4, 1]]'
expect_exact 'det([1, 2, 3; 0, 1, 4; 5, 6, 0])' '1'
expect_exact 'rref([1, 2, -1, -4; 2, 3, -1, -11; -2, 0, -3, 22])' '[[1, 0, 0, -8], [0, 1, 0, 1], [0, 0, 1, -2]]'
expect_exact 'solve([3, 2, -1; 2, -2, 4; -1, 0.5, -1], [1, -2, 0])' '[[1], [-2], [-2]]'
expect_exact 'rank([1, 2, 3; 2, 4, 6; 1, 1, 1])' '2'
expect_exact 'null([1, 2, 3; 2, 4, 6])' '[[-2, -3], [1, 0], [0, 1]]'
expect_exact 'pinv([1, 2; 2, 4])' '[[0.04, 0.08], [0.08, 0.16]]'
expect_exact 'eigvals([2, 0; 0, 3])' '[3, 2]'
# =============================================================================
# 6. Complex Number Arithmetic and Analysis
# =============================================================================
expect_numeric 'real(exp(complex(0, 1) * 3.1415926535))' '-1' '0.000001'
# Complex multiplication (1+2i)*(3+4i) = -5+10i
expect_exact 'complex(1, 2) * complex(3, 4)' 'complex(-5, 10)'
# Division 1/i = -i
expect_exact '1 / complex(0, 1)' 'complex(0, -1)'
# Residues
expect_exact 'residue(1/z, z, 0)' '[1, 0]'
expect_exact 'residue(1/(z^2+1), z, complex(0,1))' '[0, -0.5]'

# =============================================================================
# 7. Transform Rules
# =============================================================================
expect_one_of 'laplace(exp(-2 * t), t, s)' '1 / (s + 2)' '1 / (s - -2)'
expect_exact 'fourier(exp(-abs(t)), t, w)' '2 / (w ^ 2 + 1)'
expect_one_of 'ztrans(step(n - 2), n, z)' \
    '1 / (z * (z - 1))' \
    'z ^ -1 / (z - 1)'
expect_one_of 'ztrans(n ^ 2, n, z)' \
    'z * (z + 1) / (z - 1) ^ 3' \
    'z * (z + 1) / ((z - 1) * (z - 1) * (z - 1))'

printf '\nCLI comprehensive validation passed: %d\n' "$passed"
printf 'CLI comprehensive validation failed: %d\n' "$failed"

if [[ "$failed" -ne 0 ]]; then
    exit 1
fi
