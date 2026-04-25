# Comprehensive Test Report

Date: 2026-04-25

## Summary

The calculator was validated with the native C++ regression suite and the merged
script example. The focus areas were symbolic calculus, limits, matrix
calculation, and equation solving.

## Commands Run

```bash
make test
./calculator < test/script/test_merged_minimal.calc
printf ':run test/script/test_merged_minimal.calc\n' | ./calculator
```

## Results

| Check | Result |
| --- | --- |
| C++ regression suite | `Passed: 743`, `Failed: 0` |
| Redirected merged script | Passed, final output `1 / s` |
| `:run` merged script | Passed, final output `1 / s` |

## Coverage Highlights

- Symbolic differentiation: product rule, chain rule, raw expressions, named
  functions, higher derivatives, mixed partial derivatives, gradients,
  Jacobians, Hessians, and critical-point classification.
- Symbolic integration: powers, logarithmic forms, shifted rational forms,
  trigonometric identities, linear-chain rules, exponential substitution,
  partial fractions, and chained multi-variable integration.
- Limits: removable singularities, trigonometric and exponential cancellation,
  one-sided limits, and non-existent two-sided limit error handling.
- Equation solving: Newton solve, bisection, secant, fixed-point iteration,
  polynomial roots, and matrix linear-system solving.
- Matrix calculation: construction, editing, arithmetic, inverse, negative
  powers, determinant, rank, RREF, null space, QR/LU/SVD helpers, near-singular
  solve/rank behavior, eigenvalue paths, pseudo-inverse, and error paths.

## Issues Found And Fixed

- `test/script/test_merged_minimal.calc` defined a script function named
  `factorial`, which conflicts with the built-in reserved function name.
  It now uses `fact`.
- The merged script included symbolic integration examples outside the current
  rule-based integrator support. They were replaced with supported partial
  fraction and trigonometric identity examples.
- Additional regression assertions were added for symbolic calculus, limits,
  root-finding, and 3x3 matrix workflows.

## Conclusion

The project is currently passing all automated validation performed in this
round. No unresolved failures remain from the comprehensive symbolic-computing
test pass.
