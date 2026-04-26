# Test Scripts

This directory contains runnable calculator script examples and script syntax
notes.

## Running The Comprehensive Script

From the project root:

```bash
bin/calculator test/script/comprehensive_validation.calc
printf ':run test/script/comprehensive_validation.calc\n' | bin/calculator
```

The expected final output is:

```text
11
```

## Files

- `comprehensive_validation.calc`
  Broad script and calculator feature validation.
- `SYNTAX_GUIDE.md`
  The dedicated script syntax guide.
- `TEST_REPORT_MERGED.md`
  Historical validation notes for the older merged script.
