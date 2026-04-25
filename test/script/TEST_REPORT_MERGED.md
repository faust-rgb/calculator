# Merged Script Test Report

Date: 2026-04-25

## Scope

`test_merged_minimal.calc` consolidates representative examples from the older
standalone script inputs. It covers variable assignment, control flow, script
functions, selected math helpers, matrix arithmetic, symbolic calculus, critical
points, transforms, and simplification.

## Validation

```bash
./calculator < test/script/test_merged_minimal.calc
printf ':run test/script/test_merged_minimal.calc\n' | ./calculator
```

Both execution modes passed. The script's final expression is
`laplace(step(t), t, s)`, so the expected final output is:

```text
1 / s
```

## Notes

- The script uses `fact(n)` for recursion because `factorial` is a built-in
  reserved function name.
- Symbolic integration examples are intentionally limited to expressions
  supported by the current rule-based integrator.
