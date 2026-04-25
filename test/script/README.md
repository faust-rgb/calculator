# Test Scripts

This directory contains runnable calculator script examples and script syntax
notes.

## Running The Merged Script

From the project root:

```bash
./calculator < test/script/test_merged_minimal.calc
printf ':run test/script/test_merged_minimal.calc\n' | ./calculator
```

The expected final output is:

```text
1 / s
```

## Files

- `test_merged_minimal.calc`
  Consolidated representative workflow covering variables, `if`, `for`,
  `while`, recursive script functions, selected math helpers, matrix arithmetic,
  supported symbolic calculus examples, critical points, transforms, and
  simplification.
- `SYNTAX_GUIDE.md`
  Script syntax and redirected-stdin behavior notes.
- `TEST_REPORT_MERGED.md`
  Current validation notes for the merged script.

## Script Syntax Notes

1. Use `#` for line comments.
2. End simple statements with `;`.
3. Use `^` for exponentiation in calculator expressions.
4. Do not define script functions with built-in reserved names such as
   `factorial`, `sin`, or `solve`.
5. Do not put `;` after the closing brace of `if`, `while`, `for`, or `fn`
   blocks.
