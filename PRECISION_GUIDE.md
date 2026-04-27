# Precision Guide

The calculator has two numeric paths during the 2.0 transition.

## 2.0 High-Precision Path

Migrated 2.0 forms route to the exact-first runtime without requiring `:v2`.
The prefix is still accepted when you want to force or debug the 2.0 path.

```text
:precision 60
1000000000000000000000000000000 + 1
[[1, 2], [3, 4]]
simplify(sin(x)^2 + cos(x)^2)
```

Expected behavior:

- integer literals use `numeric::BigInt`
- exact divisions like `1 / 7` use `numeric::Rational`
- decimal divisions like `1.0 / 7` use `numeric::BigDecimal` and the current
  `PrecisionContext`
- complex literals such as `3 + 4i` use the 2.0 complex backend

The 2.0 arithmetic evaluator does not round ordinary arithmetic through
`double`. Unsupported transcendental numeric cases use the legacy compatibility
runtime unless explicitly forced with `:v2`; forced 2.0 evaluation keeps them
symbolic instead of silently replacing exact values with double approximations.

## Legacy 1.0 Path

Inputs that are not yet migrated still use the compatibility runtime for broad
1.0 support. That path intentionally retains `double` for scientific functions,
analysis commands, and the legacy matrix backend.

Use the legacy path for mature 1.0 features such as:

- numerical root solving and ODE commands
- matrix decompositions such as SVD/eigenvalue helpers
- probability helpers and random distributions
- broad scientific-function approximations

## Persistence

New saves use `STATE_V5` and preserve 2.0 precision plus 2.0 variables. Older
`STATE_V1` through `STATE_V4` files remain loadable.

## Validation

Regression coverage checks that 60-digit `:v2 1.0 / 7` and very large integer
addition keep their exact/high-precision text, which would fail if the path
rounded through `double`.
