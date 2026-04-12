# Architecture

## Overview

This project is a C++ command-line calculator with three overlapping roles:

- a normal scientific calculator
- a rational/exact-mode calculator
- a small programmer/CAS-style terminal tool

The implementation is intentionally lightweight and self-contained. Core math is
implemented in project code instead of relying on the standard math library for
the main functions.

## Main Files

- `/home/roselia/code/src/app/main.cpp`
  Terminal interaction, history, autocomplete, and command dispatch
- `/home/roselia/code/src/core/calculator.cpp`
  Expression parsing, exact mode, symbolic constants mode, variables, display-only features, persistence
- `/home/roselia/code/src/core/calculator.h`
  Public calculator API
- `/home/roselia/code/src/symbolic/symbolic_expression.cpp`
  Symbolic expression parsing, simplification, and rendering
- `/home/roselia/code/src/math/mymath.cpp`
  Numerical algorithms and domain handling
- `/home/roselia/code/src/math/mymath.h`
  Math declarations and shared constants
- `/home/roselia/code/test/tests.cpp`
  Regression suite covering supported behavior

## Directory Layout

- `src/app`
  CLI entry point
- `src/core`
  Calculator runtime and API
- `src/math`
  Custom numeric functions
- `src/matrix`
  Matrix representation and algorithms
- `src/analysis`
  Function analysis helpers
- `src/algebra`
  Polynomial operations
- `src/symbolic`
  Symbolic expression support
- `src/script`
  Script AST and parser
- `test`
  Regression tests and runnable example scripts

## Execution Flow

User input goes through this rough path:

1. `src/app/main.cpp` reads a line from the terminal
2. command-style inputs such as `:help`, `:vars`, `:save`, `:load` are handled first
3. special display-only expressions such as `factor(...)` are handled explicitly
4. other expressions are passed to `Calculator::process_line(...)`
5. `Calculator` either:
   - treats the line as an assignment, or
   - evaluates it for display

## Parsing Model

There are two parser implementations in `src/core/calculator.cpp`.

### `DecimalParser`

Used for standard evaluation with `double`.

This path supports:

- all ordinary math functions
- all programmer-style integer functions
- prefixed integer literals such as `0b1010`, `0o77`, `0xFF`
- variable lookup

### `ExactParser`

Used when exact fraction mode is enabled.

This path tries to preserve expressions as `Rational` values.

It supports only operations/functions that still make sense in rational form.
If a feature cannot stay exact, the code throws `ExactModeUnsupported`, and the
caller falls back to decimal display.

## Stored Values

Calculator state lives in `Calculator::Impl`.

Variables are stored in:

- `std::map<std::string, StoredValue> variables`

`StoredValue` can hold either:

- an exact rational value
- or a decimal fallback value
- and optionally a scalar symbolic display string used by symbolic constants mode

This split allows exact mode to keep fractions when possible, while still
supporting assignments like `y = sin(pi / 2)`.

## Symbolic Constants Flow

Symbolic constants mode is a scalar display feature layered on top of the normal
numeric engine.

When enabled:

1. ordinary evaluation still computes the numeric result
2. the symbolic render path tries to rebuild a scalar expression using `pi` and
   `e`
3. if it succeeds, the symbolic form is displayed and can also be stored with
   scalar variables for later reuse
4. if it cannot represent the expression safely, the code falls back to normal
   numeric display

The current implementation lives mainly in:

- `src/core/calculator.cpp`
  mode flag, storage, and display dispatch
- `src/symbolic/symbolic_expression.cpp`
  symbolic parsing and simplification rules

## Display-only Features

Some features are not plain numeric expressions and therefore are handled
outside normal parser return values.

These currently include:

- `factor(...)`
- `bin(...)`
- `oct(...)`
- `hex(...)`
- `base(...)`

These return formatted strings rather than just numbers.

## Terminal UX

`src/app/main.cpp` implements a raw terminal reader with:

- up/down arrow history recall
- `Tab` autocomplete
- fallback to `getline` when input is not a TTY

Autocomplete currently uses a fixed dictionary of common commands and functions.

## Persistence

Variable state can be saved and restored with:

- `:save file`
- `:load file`

The persistence format is simple tab-separated text written by
`Calculator::save_state(...)` and read by `Calculator::load_state(...)`.

## Numeric Design Notes

Important numeric behavior:

- near-zero results are normalized to `0`
- decimal results within `1e-10` of the nearest integer are displayed as that integer
- exact mode does not try symbolic algebra beyond rational arithmetic
- scientific notation input such as `1e-3` is not currently supported
- negative-base fractional powers are supported only when they correspond to a
  real odd-denominator rational exponent
- base conversion output supports bases `2..16`
- bitwise functions only accept integers

## Best Re-entry Points

For future work, the fastest way to rebuild context is:

1. `/home/roselia/code/HANDOFF.md`
2. `/home/roselia/code/README.md`
3. `/home/roselia/code/test/tests.cpp`
4. `/home/roselia/code/src/core/calculator.cpp`
5. `/home/roselia/code/src/symbolic/symbolic_expression.cpp`
6. `/home/roselia/code/src/math/mymath.cpp`
