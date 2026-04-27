# Simplify Status And Next Steps

## Purpose

This note captures the current symbolic simplification boundary and a practical
roadmap for future work. It is intended as the first stop for the next session
that continues simplification improvements.

## Current Command

- `simplify(expr)`
  - Simplifies a symbolic expression using the current symbolic expression
    parser and simplifier
  - Can also simplify nested command results such as `diff(...)`,
    `integral(...)`, and `poly_*` expressions when they can be converted into a
    symbolic expression first

## Current Support Matrix

| Category | Supported Now | Not Reliably Supported Yet |
| --- | --- | --- |
| Basic arithmetic simplification | `x * x / x -> x`, `2 * 4 * x / 4 -> 2 * x`, zero/one cleanup | Deep normalization across large mixed expressions |
| Like-term merging | `2*x + 3*x -> 5*x`, `x + x -> 2*x`, `5*sin(x) - 2*sin(x) -> 3*sin(x)` | Broad equivalence beyond currently normalized structures |
| Commuted grouped terms | `(x + 1) + 2 * (1 + x) -> 3 * (x + 1)` | More complex nested sums/products with multiple normalization layers |
| Variable-factor cancellation | `(x*x)/x -> x`, `(x^2)/x -> x` | Full power arithmetic like `x^5 / x^2 -> x^3` in canonical display form |
| Numeric factor cancellation | `(2*x)/4 -> x/2`, `(2*6)/(4*y)` evaluates correctly | Rich rational-expression normalization across larger symbolic fractions |
| Polynomial-style cleanup | Many expanded polynomial combinations simplify reasonably | Full polynomial canonicalization, term sorting, and factor extraction |
| Trig simplification | Repeated identical trig terms merge; special angles simplify | General trig identities like `sin(x)^2 + cos(x)^2 -> 1` |
| Exponential/log simplification | `ln(e) -> 1`, `exp(ln(x)) -> x` | General laws like `exp(a)*exp(b) -> exp(a+b)`, `ln(a)+ln(b) -> ln(a*b)` |
| Mixed command simplification | `simplify(diff(diff(g2))) -> 6 * x` | All nested symbolic command combinations in every display path |

## Examples That Work Today

| Input | Current Output |
| --- | --- |
| `simplify(x*x/x)` | `x` |
| `simplify((x + 1) + 2 * (1 + x))` | `3 * (x + 1)` |
| `simplify(diff(diff(g2)))` with `g2(x)=x^3` | `6 * x` |
| `simplify(5*sin(x) - 2*sin(x))` | `3 * sin(x)` |
| `simplify((2*x)/4)` | `x / 2` |

## Important Display Constraints

The simplifier must preserve human-friendly output. Recent work specifically
protected these forms from regressing into less readable decimal-factor forms:

- `pi / 2` should stay `pi / 2`, not `0.5 * pi`
- `sqrt(3) / 2` should stay `sqrt(3) / 2`, not `0.5 * sqrt(3)`
- `x ^ 3 / 3` should stay `x ^ 3 / 3`, not `0.333333333333 * x ^ 3`

When extending simplification, always check both:

1. algebraic correctness
2. display quality

## Good Next Improvements

### 1. Power combination and cancellation

High-value targets:

- `x^3 / x^2 -> x`
- `x^5 * x^2 -> x^7`
- `x * x * x -> x^3`
- `(x + 1)^2 / (x + 1) -> x + 1`

Likely implementation direction:

- Extend multiplicative factor collection to track `(base, exponent)` pairs
- Normalize repeated identical factors into powers
- Subtract exponents during fraction cancellation

### 2. Stronger common-factor extraction

High-value targets:

- `x*(y+1) + 2*x*(y+1) -> 3*x*(y+1)`
- `2*x + 2*y -> 2*(x + y)`
- `a*b + a*c -> a*(b + c)`
- `-x - x -> -2*x`

Likely implementation direction:

- Move beyond “same whole rest-expression” matching
- Introduce additive term decomposition with shared multiplicative-factor
  detection

### 3. Better canonical polynomial output

High-value targets:

- Stable term ordering by degree
- Constant, linear, quadratic term collection into standard polynomial form
- More predictable output for mixed expanded polynomial expressions

Likely implementation direction:

- Reuse `polynomial_coefficients(...)` when an expression is polynomial
- Render through a dedicated canonical polynomial output path

### 4. Stronger rational-expression reduction

High-value targets:

- `(2*x)/(4*y) -> x/(2*y)`
- `(x*y)/x -> y`
- `(x^2*y)/x -> x*y`
- `(x+1)/(x+1) -> 1`
- `(x^2 - 1)/(x - 1) -> x + 1`

Likely implementation direction:

- Continue factor-list cancellation
- Later add polynomial division / factorization for reducible polynomial ratios

### 5. Broader symbolic identities

Potential future targets:

- `sin(x)^2 + cos(x)^2 -> 1`
- `exp(a) * exp(b) -> exp(a + b)`
- `ln(a) + ln(b) -> ln(a*b)` with domain caveats

These should be treated cautiously because display quality and mathematical
validity can depend on domains and branch behavior.

## Risk Areas

- Over-simplifying into less readable decimal-factor output
- Infinite rewrite loops between equivalent forms
- Applying identities that are not universally valid over all domains
- Expanding and re-factoring the same expression repeatedly

## Suggested Re-entry Plan

1. Read this file
2. Read `src/symbolic/symbolic_expression_core.cpp`
3. Search for:
   - `simplify(expr)`
   - `canonical_expression_key`
   - `collect_division_factors`
   - `try_combine_like_terms`
4. Run `make test`
5. Add one targeted simplification rule at a time, with regression tests

## Current Baseline

- Test command: `make test`
- Expected status after the latest simplification work:
  - `Passed: 373`
  - `Failed: 0`
