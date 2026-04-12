# Changelog

## Current State

The project has evolved from a minimal C++ hello-world style setup into a
feature-rich command-line calculator with exact rational mode, programmable
helpers, interactive terminal UX, and project documentation.

## Major Additions Implemented

- Interactive command-line calculator REPL
- VS Code build/debug configuration cleanup
- Custom math implementation without relying on the standard math library for
  core functionality
- Expression parser with:
  - operator precedence
  - parentheses
  - unary operators
  - right-associative exponentiation
- Mathematical functions:
  - `sin`, `cos`, `tan`
  - `asin`, `acos`, `atan`
  - `exp`, `ln`, `log10`
  - `sqrt`, `cbrt`, `root`
  - `pow`
- Result normalization for near-zero floating-point noise
- Terminal UX improvements:
  - up/down history navigation
  - `Tab` autocomplete
- Exact rational mode:
  - `:exact on`
  - `:exact off`
  - exact fraction display when possible
- Symbolic constants mode:
  - `:symbolic on`
  - `:symbolic off`
  - scalar symbolic display for expressions involving `pi` and `e`
  - symbolic propagation through scalar variable assignment and custom unary functions
  - common simplifications such as `sin(pi / 2)`, `sin(pi / 6)`, `sin(pi / 3)`,
    `cos(pi)`, `cos(pi / 3)`, `cos(pi / 6)`, `tan(pi / 4)`, `tan(pi / 6)`,
    `tan(pi / 3)`, `ln(e)`, and `exp(ln(pi))`
- Display normalization for near-integer floating-point results
- Integer / utility functions:
  - `gcd`, `lcm`, `mod`
  - `abs`, `sign`, `floor`, `ceil`
  - `min`, `max`
- Variable system:
  - assignments
  - variable reuse
  - variable listing/clearing
- Help system:
  - `help`
  - `:help`
  - topic help for commands/functions/examples
- Persistence:
  - `:save file`
  - `:load file`
- Prime factorization:
  - `factor(n)`
- Base conversion:
  - `bin`, `oct`, `hex`, `base`
- Prefixed integer literal parsing:
  - `0b`
  - `0o`
  - `0x`
- Bitwise functions:
  - `and`, `or`, `xor`, `not`, `shl`, `shr`
- Expanded regression tests
- Handoff and architecture documentation

## Current Quality Status

- `make test` currently passes
- Expected test summary:
  - `Passed: 339, Failed: 0`
