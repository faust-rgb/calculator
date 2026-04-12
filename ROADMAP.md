# Roadmap

## Short Term

- Scientific notation parsing such as `1e-3`
- More programmer-style functions:
  - `rol`
  - `ror`
  - `popcount`
  - `bitlen`
- More numeric helpers:
  - `round`
  - `trunc`
  - `clamp`
- More logarithm helpers:
  - `log(x, base)`
  - `log2`
  - `exp2`

## Medium Term

- Better autocomplete:
  - show candidates on double-Tab
  - more context-sensitive completion
- More help topics:
  - `:help exact`
  - `:help variables`
  - `:help persistence`
  - `:help programmer`
- Better output formatting options:
  - prefixed hex output like `0xFF`
  - lowercase/uppercase hex mode
- More aggregate helpers:
  - `sum`
  - `avg`
  - `median`

## Long Term

- Lightweight symbolic simplification
- Symbolic substitution
- Small expression transformation helpers
- Better persistence format with versioning
- Multi-session history persistence

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
