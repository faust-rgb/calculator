# Comprehensive Script Test Report

Date: 2026-04-25

## Scope

`comprehensive_validation.calc` covers representative script syntax and
calculator features in one executable regression script.

## Validation

```bash
bin/calculator test/script/comprehensive_validation.calc
printf ':run test/script/comprehensive_validation.calc\n' | bin/calculator
```

Both execution modes passed. The script uses an internal score counter; the
expected final output is:

```text
11
```

## Notes

- Detailed syntax rules are maintained only in `test/script/SYNTAX_GUIDE.md`.
- The script keeps command-only symbolic expressions as standalone statements.
