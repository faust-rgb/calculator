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
- Interactive editing is intentionally minimal and does not provide full shell-like editing

## Scope

This project is now a capable calculator and mini CAS-style tool, but it is not:

- a full computer algebra system
- a big-number arbitrary-precision engine
- a full programming language REPL

## Symbolic Algebra

- Symbolic integration is rule-based rather than complete
  - full Risch-style integration is not implemented
  - multiple distinct irreducible quadratic factors, general irreducible high-degree factors, and symbolic-parameter factorization are still limited
  - reciprocal trigonometric rules cover common `sec`/`csc`/`cot` cases, but
    general trigonometric reduction and universal tangent substitution are not
    implemented
- Nonlinear multi-variable `critical(...)` uses bounded numeric search
  - Hessian-based classification is available for isolated 1-3 variable solutions, but global completeness is not guaranteed
- Common subexpression elimination is still partial
  - expression nodes are interned and structural keys are cached, but there is no explicit CSE pass that rewrites expressions with shared temporary subexpressions
- Simplification is algebraic and heuristic
  - mathematically equivalent expressions may print in a different but still correct form, such as reordered products or `1 - x ^ 2` vs `-(x ^ 2) + 1`
- Domain-aware simplification is still partial
  - some expressions are intentionally not collapsed unless positivity is known; full condition-tracking is not implemented
- Internal symbolic coefficients still rely on `double`
  - this keeps the implementation lightweight, but it can limit robustness for more advanced exact algebra
