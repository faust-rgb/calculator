# Migration 1.0 To 2.0

2.0 is currently shipped as an additive runtime beside the 1.0 calculator. This
keeps compatibility while the exact/high-precision path expands.

## What To Use

- Use regular expressions for mature 1.0 behavior and broad function coverage.
- Migrated 2.0 forms now route to the 2.0 runtime by default: large integer
  arithmetic, complex literals, nested matrix literals, and 2.0 CAS commands.
- Use `:v2 expr` when you want to force 2.0 behavior for a partially migrated
  or symbolic expression.
- Use `:precision digits` to control 2.0 decimal division precision.

## Command Changes

- `:vars` lists both 1.0 and 2.0 variables.
- `:clear name` and `:clear` clear both variable stores.
- `:save` writes `STATE_V5`, including 2.0 precision and variables.
- `:load` still accepts old `STATE_V1` through `STATE_V4` files.

## Output Differences

2.0 prefers exact results:

```text
:v2 1/3 + 1/6
1/2
```

Unsupported exact transcendental cases stay symbolic:

```text
:v2 sin(1)
sin(1)
```

This is intentional: the 2.0 path avoids quietly replacing exact values with
legacy approximate double values.

## Remaining Compatibility Layer

The following areas still rely on the legacy runtime for full feature coverage:

- most transcendental numeric approximation
- probability helpers
- matrix decompositions
- broad analysis commands outside the migrated tolerance and complex-root work
