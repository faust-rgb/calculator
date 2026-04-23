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
