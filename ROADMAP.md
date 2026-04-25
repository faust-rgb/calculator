# Roadmap

## Short Term

- More programmer-style functions:
  - bit-mask helpers
  - signed/unsigned display helpers
- More statistical helpers:
  - `skewness`
  - `kurtosis`
  - `histogram`
- More matrix reporting helpers:
  - decomposition summary commands
  - richer formatting for factorized results

## Medium Term

- Better autocomplete:
  - more context-sensitive completion
- Better output formatting options:
  - configurable digit grouping like `1,000,000`
  - scientific notation display toggles
  - aligned tabular matrix/vector output
- More aggregate helpers:
  - weighted statistics
  - streaming summaries
- Symbolic engine foundation upgrades:
  - stronger canonical simplify for additive and multiplicative expressions
  - exact rational constants instead of relying only on `double`
  - stable handling of symbolic constants such as `pi`, `e`, and radical forms
  - repeated simplify-to-fixed-point passes with regression tests

## Long Term

- Symbolic calculus expansion:
  - broader symbolic integral coverage for trigonometric, hyperbolic, and less common irreducible rational forms
  - stronger substitution pattern recognition beyond current chain-rule product matching
  - broader nonlinear integration-by-parts strategies beyond current polynomial times `exp`/`sin`/`cos` recurrences
- Multi-variable symbolic workflows:
  - stronger nonlinear critical point solving beyond bounded Newton search
  - safer multi-variable substitution and free-variable analysis
- Symbolic transformation quality:
  - stronger factor/expand/collect style helpers
  - better fraction cancellation and power merging
  - function identity simplifications such as `ln(exp(x))` and symmetry rules
- Symbolic transform framework:
  - continue moving Fourier/Laplace/z-transform support toward reusable pattern-matching tables
  - improve linearity decomposition and inverse-transform recognition
- Symbolic performance and maintainability:
  - structural hashing and memoization for `simplify`/`derivative`/`substitute`
  - common subexpression reuse for Taylor/Pade/Puiseux workflows
  - clearer separation between parser forms and canonical symbolic forms
- Better persistence format with versioning
- Multi-session history persistence

## Symbolic Priorities

Recommended implementation order for symbolic work:

1. Canonical simplification and exact constant representation
2. Symbolic integral rule expansion
3. Multi-variable symbolic differentiation interfaces
4. Memoization and transform-framework cleanup

## Guiding Principle

Keep the project small, understandable, and self-contained.

Prefer:

- predictable behavior
- strong regression coverage
- transparent implementation

over:

- maximal feature count
- hidden parser magic
- difficult-to-maintain abstractions
