# Test Scripts

This directory contains test scripts for the calculator. All scripts use the `.calc` extension and can be executed using the `:run` command.

## Running Tests

From the project root directory:

```bash
./calculator
:run test/script/test_variables.calc
```

Or run all tests:

```bash
for f in test/script/*.calc; do
    echo "Testing: $f"
    ./calculator <<< ":run $f"
done
```

## Test Scripts

### test_variables.calc
Tests variable assignment and basic expressions.
- Expected result: `500` (10² + 20²)

### test_math_functions.calc
Tests mathematical functions including trigonometry, logarithms, and power functions.
- Expected result: `2.71828182846` (value of e)

### test_matrix_basic.calc
Tests basic matrix operations (addition).
- Expected result: `[[6, 8], [10, 12]]`

### test_control_flow.calc
Tests if-else conditional statements.
- Expected result: `1` (true branch)

### test_while_loop.calc
Tests while loop with accumulation.
- Expected result: `55` (sum of 1 to 10)

### test_for_loop.calc
Tests for loop with accumulation.
- Expected result: `55` (sum of 1 to 10)

### test_function.calc
Tests custom function definition and recursion.
- Defines `square(x)` and `factorial(n)` functions
- Expected result: `145` (square(5) + factorial(5) = 25 + 120)

## Script Syntax Notes

1. **Comments**: Use `#` for line comments (not `//`)
2. **Statements**: End each statement with `;`
3. **Power operator**: Use `pow(base, exp)` instead of `^` (which is XOR)
4. **Reserved names**: Cannot use `e` or `pi` as variable names (they are constants)
5. **If-else**: Do not put `;` after the closing brace of if-else statements
6. **Variable scope**: Variables defined inside blocks are not visible outside
