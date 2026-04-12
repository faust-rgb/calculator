# Calculator Script Syntax Guide

This document describes the syntax for writing calculator scripts (`.calc` files).

## Basic Structure

- Scripts are plain text files with `.calc` extension
- Each statement ends with a semicolon `;`
- Whitespace is generally ignored (except within string literals)
- Comments start with `#` and continue to end of line

## Comments

```calc
# This is a single-line comment
x = 10;  # Comments can appear after statements
```

**Note**: Use `#` for comments, not `//`

## Variables

### Naming Rules
- Must start with a letter or underscore
- Can contain letters, digits, and underscores
- Case-sensitive

### Reserved Names (Cannot be used as variables)
- `e` - Euler's number (2.71828...)
- `pi` - Pi (3.14159...)
- Mathematical function names: `sin`, `cos`, `exp`, `ln`, etc.

### Assignment
```calc
x = 10;
y = 3.14;
result = x + y;
```

## Data Types

### Numbers
- Integers: `42`, `-10`
- Floating point: `3.14`, `-0.5`
- Scientific notation: `1.5e3` (1500), `2.5e-2` (0.025)

### Matrices
```calc
# Row-major format: [row1; row2; ...]
m = [1, 2, 3; 4, 5, 6];  # 2x3 matrix
v = vec(1, 2, 3);         # Column vector
z = zeros(2, 3);          # 2x3 zero matrix
i = eye(3);               # 3x3 identity matrix
```

## Operators

### Arithmetic
| Operator | Description | Example |
|----------|-------------|---------|
| `+` | Addition | `x + y` |
| `-` | Subtraction | `x - y` |
| `*` | Multiplication | `x * y` |
| `/` | Division | `x / y` |
| `^` | XOR (bitwise) | `x ^ y` |
| `%` | Not supported | Use `mod(a, b)` |

**Important**: For power/exponentiation, use `pow(base, exponent)` not `^`
```calc
# Correct
result = pow(2, 10);  # 1024

# Wrong - this is XOR!
result = 2 ^ 10;      # 8 (binary XOR)
```

### Comparison (for conditions)
| Operator | Description |
|----------|-------------|
| `==` | Equal to |
| `!=` | Not equal to |
| `<` | Less than |
| `>` | Greater than |
| `<=` | Less than or equal |
| `>=` | Greater than or equal |

## Control Flow

### If-Else Statement
```calc
if (condition) {
    # then branch
} else {
    # else branch
}
```

**Important**: Do NOT put a semicolon after the closing brace!
```calc
# Correct
if (x > 0) {
    result = 1;
} else {
    result = 0;
}

# Wrong - semicolon after if-else
if (x > 0) {
    result = 1;
} else {
    result = 0;
};  # ERROR!
```

### While Loop
```calc
while (condition) {
    # loop body
}
```

Example:
```calc
sum = 0;
i = 1;
while (i <= 10) {
    sum = sum + i;
    i = i + 1;
}
```

### For Loop
```calc
for (initializer; condition; step) {
    # loop body
}
```

Example:
```calc
sum = 0;
for (i = 1; i <= 10; i = i + 1) {
    sum = sum + i;
}
```

### Break and Continue
```calc
while (condition) {
    if (some_condition) {
        break;      # Exit loop
    }
    if (other_condition) {
        continue;   # Skip to next iteration
    }
}
```

## Functions

### Defining Functions
```calc
fn function_name(param1, param2) {
    # function body
    return expression;
}
```

Example:
```calc
fn square(x) {
    return x * x;
}

fn hypotenuse(a, b) {
    return sqrt(a*a + b*b);
}

fn factorial(n) {
    if (n <= 1) {
        return 1;
    }
    return n * factorial(n - 1);
}
```

**Note**: Function definitions do not need a trailing semicolon.

### Calling Functions
```calc
result = square(5);
h = hypotenuse(3, 4);
f = factorial(10);
```

## Mathematical Functions

### Trigonometric
- `sin(x)`, `cos(x)`, `tan(x)` - Arguments in radians
- `asin(x)`, `acos(x)`, `atan(x)` - Inverse trig functions

### Exponential and Logarithmic
- `exp(x)` - e^x
- `ln(x)` - Natural logarithm
- `log10(x)` - Base-10 logarithm
- `sqrt(x)` - Square root
- `cbrt(x)` - Cube root
- `root(x, n)` - n-th root
- `pow(x, y)` - x^y

### Other
- `abs(x)` - Absolute value
- `floor(x)` - Floor function
- `ceil(x)` - Ceiling function
- `sign(x)` - Sign function (-1, 0, or 1)
- `gcd(a, b)` - Greatest common divisor
- `lcm(a, b)` - Least common multiple
- `mod(a, b)` - Modulo operation
- `min(a, b)`, `max(a, b)` - Minimum/maximum

## Matrix Operations

### Creation
```calc
m = [1, 2; 3, 4];           # Literal
v = vec(1, 2, 3);           # Column vector
z = zeros(2, 3);            # Zero matrix
e = eye(3);                 # Identity matrix
```

### Operations
```calc
# Arithmetic
m1 + m2;        # Addition
m1 - m2;        # Subtraction
m1 * m2;        # Matrix multiplication
m * 2;          # Scalar multiplication

# Properties
det(m);         # Determinant
trace(m);       # Trace
rank(m);        # Rank
norm(m);        # Frobenius norm

# Transformations
transpose(m);   # Transpose
inverse(m);     # Inverse
rref(m);        # Row-reduced echelon form

# Decompositions
qr_q(m);        # QR decomposition: Q matrix
qr_r(m);        # QR decomposition: R matrix
lu_l(m);        # LU decomposition: L matrix
lu_u(m);        # LU decomposition: U matrix
svd_u(m);       # SVD: U matrix
svd_s(m);       # SVD: S (singular values) matrix
svd_vt(m);      # SVD: V^T matrix

# Solving
solve(A, b);            # Solve Ax = b
least_squares(A, b);    # Least squares solution
```

## Variable Scope

Variables defined inside blocks (`{ }`) are not visible outside:

```calc
# Correct - declare outside, assign inside
result = 0;
if (x > 5) {
    result = 1;
}
result;

# Wrong - result is not visible outside
if (x > 5) {
    result = 1;  # result only exists here
}
result;  # ERROR: unknown variable
```

## Common Mistakes

### 1. Using `//` for comments
```calc
// Wrong - use # instead
# Correct
```

### 2. Using `^` for power
```calc
# Wrong - this is XOR
result = 2 ^ 10;  # = 8

# Correct
result = pow(2, 10);  # = 1024
```

### 3. Using `e` or `pi` as variables
```calc
# Wrong - e is a constant
e = 10;  # ERROR

# Correct
my_e = 10;
```

### 4. Semicolon after if-else
```calc
# Wrong
if (x > 0) { ... } else { ... };  # ERROR

# Correct
if (x > 0) { ... } else { ... }
```

### 5. Missing semicolons inside blocks
```calc
# Wrong
if (x > 0) {
    result = 1  # Missing semicolon!
}

# Correct
if (x > 0) {
    result = 1;
}
```

## Complete Example

```calc
# Calculate factorial using a loop
fn factorial(n) {
    result = 1;
    i = 1;
    while (i <= n) {
        result = result * i;
        i = i + 1;
    }
    return result;
}

# Calculate combination C(n, k)
fn combination(n, k) {
    return factorial(n) / (factorial(k) * factorial(n - k));
}

# Test
combination(10, 5);  # Returns 252
```

## Running Scripts

From the calculator prompt:
```
:run filename.calc
```

From command line:
```bash
./calculator <<< ":run filename.calc"
```
