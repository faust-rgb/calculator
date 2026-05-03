# Calculator Script Syntax Guide

This is the single source of truth for `.calc` script syntax. Scripts now use a
Python-like layout: `#` comments, colon-started blocks, indentation for scope,
and no statement-ending semicolons.

## Running Scripts

From the project root:

```bash
bin/calculator test/script/comprehensive_validation.calc
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
- Start `if`, `elif`, `else`, `while`, `for`, `def`, and `match` blocks with `:`.
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

### Keywords and Reserved Names

Script keywords cannot be used as variable names:

```
def, fn, if, elif, else, while, for, in, match, case, return, break, continue, pass
```

Reserved function names cannot be used for user-defined functions. See
`KEYWORDS_REFERENCE.md` for the complete list. Common examples include:

```
sin, cos, tan, exp, ln, sqrt, abs, det, trace, inv, transpose, diff, integral
```

Variable naming rules:

- First character must be a letter or underscore
- Subsequent characters can be letters, digits, or underscores
- Cannot be a keyword or reserved function name

## Expressions

Calculator expressions use the same operators and built-ins as the REPL.

```calc
total = 1 + 2 * 3
power = 2 ^ 10
distance = sqrt(3 ^ 2 + 4 ^ 2)
ok = abs(sin(pi / 2) - 1) < 0.000001
```

Logical operators `&&` (and) and `||` (or) are supported with short-circuit
evaluation:

```calc
if x > 0 && x < 10:
    print("x is in range (0, 10)")
if a == 0 || b == 0:
    print("at least one is zero")
```

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

### For-In Loops

The `for-in` loop supports iterating over lists, matrices (row-by-row), and
strings (character-by-character):

```calc
# Iterate over a list
for item in [1, 2, 3, 4, 5]:
    print(item)

# Iterate over matrix rows
m = [1, 2; 3, 4; 5, 6]
for row in m:
    print("row: ", row)

# Iterate over string characters
for ch in "hello":
    print(ch)
```

### Match-Case Pattern Matching

Pattern matching similar to Python 3.10+ is supported:

```calc
match value:
    case 0:
        print("zero")
    case 1:
        print("one")
    case _:
        print("other")
```

Cases can have guard conditions with `if`:

```calc
match x:
    case 0 if y > 0:
        print("zero with positive y")
    case 0:
        print("zero")
    case _:
        print("other")
```

Match works with scalars, strings, and matrices:

```calc
match status:
    case "ok":
        print("success")
    case "error":
        print("failed")
    case _:
        print("unknown")
```

## Indexing and Assignment

### List and Dictionary Indexing

Lists and dictionaries support Python-style indexing:

```calc
lst = [1, 2, 3, 4, 5]
print(lst[0])      # first element
print(lst[-1])     # last element
print(lst[1:3])    # slice [2, 3]
lst[2] = 10        # assignment

d = {"a": 1, "b": 2}
print(d["a"])
d["c"] = 3
```

### Matrix Indexing

Matrices support bracket indexing for both reading and assignment:

```calc
m = [1, 2; 3, 4]
print(m[0, 0])     # first element (1)
print(m[1])        # linear index (2)
m[0, 0] = 5        # assignment
print(m)           # [5, 2; 3, 4]
```

Negative indices are supported:

```calc
print(m[-1, -1])   # last element (4)
```

## File I/O

Scripts can read and write files using the I/O module:

```calc
# Write to a file
fd = open("output.txt", "w")
write(fd, "Hello, World!")
close(fd)

# Read from a file
fd = open("input.txt", "r")
content = read(fd)
close(fd)
print(content)

# Read line by line
fd = open("data.txt", "r")
lines = read_lines(fd)
close(fd)
for line in lines:
    print(line)

# File positioning
fd = open("data.txt", "r")
seek(fd, 10)        # Move to position 10
pos = tell(fd)      # Get current position
line = readline(fd) # Read single line
close(fd)

# File management
if exists("data.txt"):
    delete("data.txt")
```

### File Functions

| Function | Description |
|----------|-------------|
| `open(path, mode)` | Open a file, returns file descriptor |
| `close(fd)` | Close a file |
| `read(fd)` | Read entire file content as string |
| `write(fd, text)` | Write text to file |
| `read_lines(fd)` | Read all lines into a list |
| `readline(fd)` | Read single line |
| `seek(fd, pos)` | Set file position |
| `tell(fd)` | Get current position |
| `exists(path)` | Check if file exists (returns 1 or 0) |
| `delete(path)` | Delete a file |

### Open Modes

| Mode | Description |
|------|-------------|
| `"r"` | Read only (default) |
| `"w"` | Write (truncates existing file) |
| `"a"` | Append |
| `"rw"` or `"r+"` | Read and write |

### CSV and JSON Support

```calc
# CSV file operations
m = [1, 2, 3; 4, 5, 6]
write_csv("matrix.csv", m)

loaded = read_csv("matrix.csv")
print("Loaded matrix: ", loaded)

# JSON file operations
data = {"name": "test", "values": [1, 2, 3]}
write_json("data.json", data)

loaded_data = read_json("data.json")
print("Loaded data: ", loaded_data)
```

| Function | Description |
|----------|-------------|
| `read_csv(path)` | Read matrix from CSV file |
| `write_csv(path, matrix)` | Write matrix to CSV file |
| `read_json(path)` | Read data from JSON file |
| `write_json(path, data)` | Write data to JSON file |

### State Persistence

Use `:save` and `:load` commands to persist calculator state:

```text
> :save state.txt        # Save all variables and functions
> :load state.txt        # Restore saved state
```

### Exporting Data

Export matrices to CSV files:

```text
> m = [1, 2; 3, 4]
> :export "matrix.csv" m
Exported m to matrix.csv
```

The matrix will be saved in CSV format with comma-separated values.

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
