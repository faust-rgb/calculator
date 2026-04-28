# Known Limitations

## Parsing

- The parser is expression-oriented and not a full symbolic algebra parser
- Function argument handling is intentionally simple and explicit

## Exact Mode

- Exact mode is limited to rational arithmetic and related exact-friendly helpers
- Expressions that cannot remain rational fall back to decimal display
- Exact mode is not a full symbolic CAS

## Base Conversion

- `bin`, `oct`, `hex`, and `base` only accept integer values
- `base(n, b)` currently supports bases `2..16`

## Bitwise Operations

- Bitwise functions only accept integers
- Shift counts must be non-negative
- Behavior follows normal C++ integer bitwise semantics on signed values

## Numeric Model

- Internal decimal calculations use `double`
- Arbitrary precision is not supported
- Some operations are approximate by nature and rely on tolerance handling
- Numerical stability for extremely high-order derivatives or limits of highly oscillatory functions may still be limited by double-precision floating point (though PSA/Richardson optimizations significantly improve this)
- Finite one-variable integrals use adaptive Gauss-Kronrod G7-K15 and include
  endpoint-singularity transforms, but infinite integration bounds are not yet
  exposed through command syntax
- ODE solving uses adaptive RKF45 with configurable internal tolerances, but it
  is not a dedicated stiff solver
- Polynomial fitting still uses the normal equations path and can be sensitive
  on ill-conditioned datasets

## Persistence

- Save/load format is simple text and not versioned
- Loading invalid files raises an error instead of attempting partial recovery

## Terminal UX

- Autocomplete uses a fixed word list
- Line editing provides basic navigation (arrows, Home/End) and deletion (Ctrl+D/K), but is not as comprehensive as `readline` or `zsh`.

## Scope

This project is now a capable calculator and mini CAS-style tool, but it is not:

- a full computer algebra system
- a big-number arbitrary-precision engine
- a full programming language REPL

## Symbolic Algebra

- Symbolic integration is rule-based rather than complete
  - full Risch-style integration is not implemented
  - multiple distinct irreducible quadratic factors, general irreducible high-degree factors, and symbolic-parameter factorization are still limited
  - reciprocal trigonometric rules cover common `sec`/`csc`/`cot` cases
  - Weierstrass substitution (universal tangent substitution) is implemented for
    rational functions of `sin`/`cos`, but general trigonometric reduction is not
    complete
- Nonlinear multi-variable `critical(...)` uses bounded numeric search
  - Hessian-based classification is available for isolated 1-3 variable solutions, but global completeness is not guaranteed
- Common subexpression elimination (CSE)
  - Nodes are interned and structural keys are cached.
  - Basic CSE extraction logic is available via `common_subexpressions()`, but the core simplifier does not yet automatically rewrite expressions into variable-sharing forms.
- Simplification is algebraic and heuristic
  - mathematically equivalent expressions may print in a different but still correct form, such as reordered products or `1 - x ^ 2` vs `-(x ^ 2) + 1`
- Domain-aware simplification is still partial
  - some expressions are intentionally not collapsed unless positivity is known; full condition-tracking is not implemented
- Symbolic Constants
  - `pi` and `e` are now first-class exact constant nodes. They are preserved during symbolic algebra passes unless explicit numerical evaluation is requested.
