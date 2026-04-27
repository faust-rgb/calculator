# CAS Guide

The calculator now has two symbolic layers:

- legacy 1.0 symbolic commands such as `diff(f)`, `integral(f)`, `limit(...)`,
  transforms, series, and custom unary functions
- 2.0 expression/CAS commands through `:v2`, including exact-first
  `simplify(expr)`, `diff(expr, x)`, `integrate(expr, x)`, `gradient`,
  `jacobian`, and `hessian`

## 2.0 Examples

```text
:v2 simplify(sin(x)^2 + cos(x)^2)
:v2 diff(x^x, x)
:v2 integrate(2*x*cos(x^2), x)
:v2 gradient(x^2 + x*y + y^2, [x, y])
:v2 hessian(x^2 + x*y + y^2, [x, y])
```

Unsupported 2.0 integrals return an unevaluated `integral(expr, x)` instead of
falling back to an approximate numeric answer.

## Current Coverage

- simplification: exact constant folding, identity removal, same-base power
  merging, repeated-factor powers, basic trig identity folding, and
  `exp(ln(x))` / `ln(exp(x))`
- differentiation: sums, products, powers including `x^x`, elementary
  trig/log/exp/hyperbolic rules, `sqrt`, and unevaluated `Derivative`
- integration: constants, linearity, powers, polynomials, selected rational
  forms, trig/exponential table rules, substitution detection, integration by
  parts, and selected special-function outputs

## Performance Notes

2.0 expression nodes are structurally interned, and CAS simplification keeps a
bounded memoization cache. This helps repeated script, matrix, and CAS workflows
reuse identical expression structure without changing public output.

## Limitations

The legacy symbolic tree still stores numeric nodes as `double`; use `:v2` when
you need exact rational or high-precision decimal expression arithmetic.
