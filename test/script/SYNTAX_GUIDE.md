# Calculator Script Syntax Guide

This is the single source of truth for `.calc` script syntax. Scripts now use a
Python-like layout: `#` comments, colon-started blocks, indentation for scope,
and no statement-ending semicolons.

## Running Scripts

From the project root:

```bash
./calculator test/script/comprehensive_validation.calc
```

From the interactive prompt:

```text
> :run test/script/comprehensive_validation.calc
```

Script files should use the `.calc` extension.

## File Structure

- Use `#` for line comments.
- Write one statement per line.
- Do not end statements with `;`.
- Start `if`, `elif`, `else`, `while`, `for`, and `def` blocks with `:`.
- Use consistent indentation to define block bodies. Spaces are recommended.
- Matrix literals still use `;` inside brackets to separate rows, for example
  `[1, 2; 3, 4]`.

```calc
# comment
x = 10
print("x = ", x)
```

## Values And Variables

Scripts can use numbers, strings, scalar variables, expression-style custom
functions, and matrices.

```calc
count = 3
label = "ready"
m = [1, 2; 3, 4]
n = mat(2, 2, 5, 6, 7, 8)
```

Variables are case-sensitive. Constants such as `pi` and `e`, built-in function
names, and active function names should not be reused as ordinary variables.

When a variable must be read after a block, initialize it before the block and
assign to it inside the block.

## Expressions

Calculator expressions use the same operators and built-ins as the REPL.

```calc
total = 1 + 2 * 3
power = 2 ^ 10
distance = sqrt(3 ^ 2 + 4 ^ 2)
ok = abs(sin(pi / 2) - 1) < 0.000001
```

Use nested `if` statements instead of `&&` or `||` in script conditions.

## Control Flow

```calc
if score >= 90:
    grade = "A"
elif score >= 80:
    grade = "B"
else:
    grade = "other"

i = 1
sum = 0
while i <= 10:
    sum = sum + i
    i = i + 1

odd_sum = 0
for j in range(0, 10):
    if j == 7:
        break
    if mod(j, 2) == 0:
        continue
    odd_sum = odd_sum + j
```

`range(stop)`, `range(start, stop)`, and `range(start, stop, step)` are
supported. The loop stops before `stop`, matching Python's convention.

## Script Functions

Use `def` for script functions. `return` is required when the function is meant
to produce a value.

```calc
def fact_script(n):
    if n <= 1:
        return 1
    return n * fact_script(n - 1)

fact_script(6)
```

Avoid defining script functions with names already used by built-ins.

## Expression-Style Functions

Expression-style functions can be defined in scripts and then used by analysis,
polynomial, and calculus commands.

```calc
quad_fn(x) = x ^ 2 + 2 * x + 1
diff(quad_fn)
integral(quad_fn, 0, 1)
taylor(quad_fn, 0, 3)
```

## Output

`print(...)` accepts script values, strings, and formatted matrices. Symbolic
commands such as `diff(x ^ 2)` or `laplace(step(t), t, s)` should usually be
written as standalone statements instead of being nested inside `print(...)`.

```calc
print("det = ", det([1, 2; 3, 4]))
laplace(step(t), t, s)
```

## Compatibility Note

The parser still accepts the older brace-and-semicolon form internally so saved
state files containing script functions remain loadable. New scripts and
examples should use the Python-like form described above.

## Feature Coverage Example

Use `test/script/comprehensive_validation.calc` as the broad regression script.
It covers:

- comments, assignments, strings, arithmetic, and comparisons
- `if/elif/else`, nested conditionals, `while`, `for range`, `break`, and
  `continue`
- recursive `def` functions and expression-style functions
- numeric math helpers and constants
- matrices, determinants, trace, rank, inverse, solve, transpose, and rref
- polynomial vector helpers
- derivative, integral, limit, Taylor, root solving, simplification, critical
  points, Pade, Laplace, Fourier, and Z transforms
- programmer-style bit helpers

The expected final output is:

```text
11
```
