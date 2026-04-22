# Known Limitations

## Parsing

- Scientific notation is not supported yet
  - example: `1e-3` is not valid input
- The parser is expression-oriented and not a full symbolic algebra parser
- Function argument handling is intentionally simple and explicit

## Exact Mode

- Exact mode is limited to rational arithmetic and related exact-friendly helpers
- Expressions that cannot remain rational fall back to decimal display
- Exact mode is not a full symbolic CAS

## Base Conversion

- `bin`, `oct`, `hex`, and `base` only accept integer values
- `base(n, b)` currently supports bases `2..16`
- No prefixed output mode toggle yet

## Bitwise Operations

- Bitwise functions only accept integers
- Shift counts must be non-negative
- Behavior follows normal C++ integer bitwise semantics on signed values

## Numeric Model

- Internal decimal calculations use `double`
- Arbitrary precision is not supported
- Some operations are approximate by nature and rely on tolerance handling
- Near-zero values are normalized to `0` for readability
- Decimal display also simplifies values that are within `1e-10` of the nearest
  integer
- In practice, this means values with `abs(value) <= 1e-10` display as `0`, and
  values with `abs(value - round(value)) <= 1e-10` display as that integer

## Persistence

- Save/load format is simple text and not versioned
- Loading invalid files raises an error instead of attempting partial recovery

## Terminal UX

- Autocomplete uses a fixed word list
- Double-Tab candidate listing is not implemented
- Interactive editing is intentionally minimal and does not provide full shell-like editing

## Scope

This project is now a capable calculator and mini CAS-style tool, but it is not:

- a full computer algebra system
- a big-number arbitrary-precision engine
- a full programming language REPL

## Symbolic Algebra

- Symbolic `diff`, `integral`, `taylor`, `limit`, and `extrema` are still limited to one-variable inputs
- Symbolic integration is rule-based rather than complete
  - many textbook forms now work, but integration by parts, substitutions, and identities are not exhaustive
- Simplification is algebraic and heuristic
  - mathematically equivalent expressions may print in a different but still correct form, such as reordered products or `1 - x ^ 2` vs `-(x ^ 2) + 1`
- Domain-aware simplification is still partial
  - some expressions are intentionally not collapsed unless positivity is known, and full condition-tracking is not implemented
- Internal symbolic coefficients still rely on `double`
  - this keeps the implementation lightweight, but it can limit robustness for more advanced exact algebra
