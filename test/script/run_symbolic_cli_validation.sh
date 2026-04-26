#!/usr/bin/env bash
set -u

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
CALC="$ROOT_DIR/bin/calculator"

passed=0
failed=0

run_calc() {
    printf '%s\n' "$1" | "$CALC" 2>&1
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

expect_numeric() {
    local expr="$1"
    local expected="$2"
    local tolerance="$3"
    local actual
    actual="$(run_calc "$expr")"
    if awk -v a="$actual" -v e="$expected" -v t="$tolerance" 'BEGIN { d = a - e; if (d < 0) d = -d; exit !(d <= t) }'; then
        passed=$((passed + 1))
    else
        failed=$((failed + 1))
        printf 'FAIL numeric: %s\n  expected: %s +/- %s\n  actual:   %s\n' "$expr" "$expected" "$tolerance" "$actual"
    fi
}

expect_error_contains() {
    local expr="$1"
    local expected="$2"
    local actual
    actual="$(run_calc "$expr")"
    if [[ "$actual" == *"$expected"* ]]; then
        passed=$((passed + 1))
    else
        failed=$((failed + 1))
        printf 'FAIL error: %s\n  expected substring: %s\n  actual:             %s\n' "$expr" "$expected" "$actual"
    fi
}

# Symbolic differentiation and calculus.
expect_exact 'diff(sin(x) * exp(x), x)' 'exp(x) * (cos(x) + sin(x))'
expect_exact 'diff(x ^ 2 * y + sin(y), x, y)' '2 * x'
expect_exact 'diff(sin(x) / x)' '(cos(x) * x - sin(x)) / x ^ 2'
expect_exact 'integral(cos(3 * x), x)' 'sin(3 * x) / 3 + C'
expect_exact 'integral(3 * x ^ 2 * exp(x ^ 3), x)' 'exp(x ^ 3) + C'
expect_exact 'integral(1 / (x + 2), x)' 'ln(abs(x + 2)) + C'
expect_exact 'integral(1 / (x ^ 2 - 1))' '1/2 * ln(abs((2 * x - 2) / (2 * x + 2))) + C'
expect_exact 'integral(sin(x) ^ 2)' 'x / 2 - sin(2 * x) / 4 + C'
expect_exact 'integral(cos(x) ^ 2)' 'sin(2 * x) / 4 + x / 2 + C'
expect_exact 'integral(tan(x) ^ 2)' 'tan(x) - x + C'
expect_exact 'gradient(x ^ 2 + x * y + y ^ 2, x, y)' '[2 * x + y, 2 * y + x]'
expect_exact 'hessian(x ^ 2 + x * y + y ^ 2, x, y)' '[[2, 1], [1, 2]]'
expect_exact 'jacobian([x ^ 2 + y; sin(x * y)], x, y)' '[[2 * x, 1], [cos(x * y) * y, cos(x * y) * x]]'
expect_exact 'critical(x ^ 2 + y ^ 2, x, y)' '[x = 0, y = 0] (local min)'
expect_exact 'critical(x ^ 2 - y ^ 2, x, y)' '[x = 0, y = 0] (saddle)'
expect_exact 'taylor(sin(x), 0, 5)' 'x - 1/6 * x ^ 3 + 1/120 * x ^ 5'

# Limits and root finding.
expect_numeric 'limit(sin(x) / x, 0)' 1 0.00000001
expect_numeric 'limit((exp(x) - 1) / x, 0)' 1 0.00000001
expect_numeric 'limit((1 - cos(x)) / x ^ 2, 0)' 0.5 0.00000001
expect_exact 'limit((x ^ 3 - 1) / (x - 1), 1)' '3'
expect_exact 'limit(abs(x) / x, 0, 1)' '1'
expect_exact 'limit(abs(x) / x, 0, -1)' '-1'
expect_error_contains 'limit(1 / x, 0)' 'limit did not converge'
expect_numeric 'solve(cos(x) - x, 0.5)' 0.739085133215 0.0000000001
expect_numeric 'bisect(x ^ 3 - x - 2, 1, 2)' 1.5213797068 0.0000000001
expect_numeric 'secant(x ^ 3 - x - 2, 1, 2)' 1.5213797068 0.0000000001
expect_numeric 'fixed_point(cos(x), 0.5)' 0.73908513325 0.0000000001

# Matrix and linear algebra workflows.
expect_exact 'inverse(mat(3, 3, 1, 2, 3, 0, 1, 4, 5, 6, 0))' '[[-24, 18, 5], [20, -15, -4], [-5, 4, 1]]'
expect_exact 'det(mat(3, 3, 1, 2, 3, 0, 1, 4, 5, 6, 0))' '1'
expect_exact 'rref(mat(3, 4, 1, 2, -1, -4, 2, 3, -1, -11, -2, 0, -3, 22))' '[[1, 0, 0, -8], [0, 1, 0, 1], [0, 0, 1, -2]]'
expect_exact 'solve(mat(3, 3, 3, 2, -1, 2, -2, 4, -1, 0.5, -1), vec(1, -2, 0))' '[[1], [-2], [-2]]'
expect_exact 'rank(mat(3, 3, 1, 2, 3, 2, 4, 6, 1, 1, 1))' '2'
expect_exact 'null(mat(2, 3, 1, 2, 3, 2, 4, 6))' '[[-2, -3], [1, 0], [0, 1]]'
expect_exact 'pinv(mat(2, 2, 1, 2, 2, 4))' '[[0.04, 0.08], [0.08, 0.16]]'
expect_exact 'eigvals(mat(2, 2, 2, 0, 0, 3))' '[3, 2]'

# Transform rules used by symbolic workflows.
expect_exact 'laplace(exp(-2 * t), t, s)' '1 / (s + 2)'
expect_exact 'ilaplace(1 / (s + 2), s, t)' 'exp(-2 * t) * step(t)'
expect_exact 'fourier(exp(-2 * t) * step(t), t, w)' '1 / (i * w + 2)'
expect_exact 'ztrans(step(n - 2), n, z)' '1 / (z * (z - 1))'
expect_exact 'iztrans(z / (z - 1), z, n)' 'step(n)'

printf 'CLI symbolic validation passed: %d\n' "$passed"
printf 'CLI symbolic validation failed: %d\n' "$failed"

if [[ "$failed" -ne 0 ]]; then
    exit 1
fi
