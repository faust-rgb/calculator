# Codex Handoff

## Project

This is a C++ command-line calculator / mini CAS in `/home/roselia/ai-code/calculator`.

Main capabilities already implemented:

- Arithmetic, parentheses, unary signs, power operator `^`
- `pow(a, b)`, `sqrt(x)`, `cbrt(x)`, `root(a, n)`
- `exp`, `ln`, `log10`
- `sin`, `cos`, `tan`, `asin`, `acos`, `atan`
- Exact rational mode via `:exact on` / `:exact off`
- Symbolic constants mode via `:symbolic on` / `:symbolic off`
- Rational simplification in exact mode
- Variables: assignment and reuse
- Variable commands: `:vars`, `:clear`, `:clear name`
- Session history command: `:history`
- Help system: `help`, `:help`, `:help commands|functions|examples`
- Save/load variables: `:save file`, `:load file`
- Integer helpers: `gcd`, `lcm`, `mod`, `factor`
- Numeric helpers: `abs`, `sign`, `floor`, `ceil`, `min`, `max`
- Base conversion output: `bin`, `oct`, `hex`, `base`
- Prefixed integer parsing: `0b...`, `0o...`, `0x...`
- Bitwise functions: `and`, `or`, `xor`, `not`, `shl`, `shr`
- Input UX: arrow-up/down history and `Tab` autocomplete

Current test status:

- `make test`
- Expected: `Passed: 743, Failed: 0`
- Comprehensive script validation: `./calculator test/script/comprehensive_validation.calc`

## Important Files

- `/home/roselia/ai-code/calculator/src/app/main.cpp`
  REPL / terminal UX / commands / history / autocomplete
- `/home/roselia/ai-code/calculator/src/core/calculator.cpp`
  Core parser + exact mode + symbolic constants mode + variables + help + factor/base conversion + persistence
- `/home/roselia/ai-code/calculator/src/core/calculator.h`
  Public calculator interface
- `/home/roselia/ai-code/calculator/src/symbolic/symbolic_expression_core.cpp`
  Symbolic expression parsing, rendering, and symbolic-constants simplification rules
- `/home/roselia/ai-code/calculator/src/symbolic/symbolic_expression_calculus.cpp`
  Symbolic differentiation and integration rules
- `/home/roselia/ai-code/calculator/src/symbolic/symbolic_expression_transforms.cpp`
  Fourier/Laplace/z transform entry points
- `/home/roselia/ai-code/calculator/src/math/mymath.cpp`
  Custom math implementations, no standard math core dependency for main functions
- `/home/roselia/ai-code/calculator/src/math/mymath.h`
  Math declarations
- `/home/roselia/ai-code/calculator/test/tests.cpp`
  Main regression suite and best overview of supported behavior
- `/home/roselia/ai-code/calculator/SIMPLIFY_IMPROVEMENTS.md`
  Current simplification boundary, support matrix, and next-step roadmap
- `/home/roselia/ai-code/calculator/README.md`
  User-facing usage doc
- `/home/roselia/ai-code/calculator/test/TESTING.md`
  Testing summary
- `/home/roselia/ai-code/calculator/TEST_REPORT.md`
  Latest comprehensive validation report

## Directory Layout

- `src/app`
  CLI entry point
- `src/core`
  Calculator runtime and public API
- `src/math`
  Numeric helpers
- `src/matrix`
  Matrix implementation
- `src/analysis`
  Function analysis helpers
- `src/algebra`
  Polynomial helpers
- `src/symbolic`
  Symbolic expression support
- `src/script`
  Script AST and parser
- `test`
  Regression tests and example scripts

## Architecture Notes

### Parsers

`src/core/calculator.cpp` has two parser paths:

- `DecimalParser`
  General floating-point evaluation path
- `ExactParser`
  Rational-only path for exact mode

Exact mode falls back to decimal display if expression cannot remain rational.

### Stored Variables

Variables live in `Calculator::Impl::variables`.

Stored type is `StoredValue`:

- `exact = true` with `Rational`
- or decimal fallback with `double`
- optional symbolic display text when symbolic constants mode captures a scalar expression involving `pi` or `e`

### Display-only Features

These are not plain numeric evaluation features; they are handled specially:

- `factor(...)` via `factor_expression(...)`
- `bin/oct/hex/base(...)` via `base_conversion_expression(...)`

### Input / UX

`src/app/main.cpp` implements:

- raw terminal mode
- up/down history recall
- `Tab` completion from a fixed word list
- REPL commands for exact mode and symbolic constants mode

Non-interactive input still falls back to `getline`.
Script files can be run directly with `./calculator file.calc`.

### Symbolic Constants Mode

There is now a display-oriented scalar mode controlled by:

- `:symbolic on`
- `:symbolic off`
- `:symbolic`

This mode preserves `pi` and `e` in scalar output when possible.

Currently verified examples include:

- `pi / 2 + e` -> `pi / 2 + e`
- `x = pi / 2`, then `x + 1` -> `pi / 2 + 1`
- `f(x) = x + pi`, then `f(e)` -> `e + pi`
- `sin(pi / 2)` -> `1`
- `sin(pi / 6)` -> `1 / 2`
- `sin(pi / 3)` -> `sqrt(3) / 2`
- `cos(pi)` -> `-1`
- `cos(pi / 3)` -> `1 / 2`
- `cos(pi / 6)` -> `sqrt(3) / 2`
- `tan(pi / 4)` -> `1`
- `tan(pi / 6)` -> `1 / sqrt(3)`
- `tan(pi / 3)` -> `sqrt(3)`
- `ln(e)` -> `1`
- `exp(ln(pi))` -> `pi`

Important boundary:

- this is primarily a display mode for scalar expressions
- unsupported symbolic forms fall back to the normal numeric result
- save/load now persists scalar symbolic display text as well as exact/decimal/string/function state

## Important Behavioral Notes

- Scientific notation like `1e-3` is supported by the parser.
- Prefixed literals support integers only: `0b`, `0o`, `0x`.
- Bitwise functions only accept integers.
- Base conversion currently supports bases `2..16`.
- Near-zero floating-point results are normalized to `0`.
- Decimal display simplifies values within `1e-10` of `0` or the nearest integer.
- `-2^2` is parsed as `-(2^2)` and returns `-4`.
- Negative-base fractional exponent support exists for odd-denominator rational exponents in numeric power handling.

## Useful Commands

```bash
make test
make
./calculator
```

## Recommended Re-entry Flow

If continuing work on another machine/session:

1. Read `/home/roselia/ai-code/calculator/HANDOFF.md`
2. Read `/home/roselia/ai-code/calculator/README.md`
3. Read `/home/roselia/ai-code/calculator/test/tests.cpp`
4. Read `/home/roselia/ai-code/calculator/src/core/calculator.cpp`
5. Read `/home/roselia/ai-code/calculator/src/symbolic/symbolic_expression_core.cpp`
6. Run `make test`

## Good Next Extensions

High-value next steps:

- Scientific notation parsing
- Extend symbolic constants mode to more unit-circle angles
- Decide whether symbolic scalar display text should be persisted by save/load
- More programmer functions: `rol`, `ror`, `popcount`, `bitlen`
- More math functions: `log(x, base)`, `log2`, `exp2`, `round`, `trunc`, `clamp`
- Better help topics: `:help exact`, `:help variables`, `:help persistence`
- Better autocomplete: show candidates on double-Tab
- Continue symbolic simplification work from `SIMPLIFY_IMPROVEMENTS.md`
