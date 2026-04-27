# Migration Plan: 1.0 to 2.0 High-Precision Path

## Context

The user wants FULL migration of all 1.0 double-precision code to use 2.0's `numeric::Number` arbitrary precision path. After migration, legacy 1.0 code can be deleted.

Current state:
- 810 tests passing
- 2.0 runtime uses `numeric::Number` with arbitrary precision
- 1.0 code uses `double`/`long double` throughout
- Node structure already has hybrid storage (double + Number)

## Implementation Plan

### Phase 1: Foundation Layer

#### Step 1.1: Add Missing Functions to numeric::functions
**Files:** `src/numeric/functions.h`, `src/numeric/functions.cpp`

Add functions that mymath provides but numeric lacks:
- `sec`, `csc`, `cot`, `asec`, `acsc`, `acot`
- `sign`, `floor`, `ceil`, `round` (already have via function_registry)
- Utility: `is_finite`, `infinity`

**Verification:** Build and run tests.

#### Step 1.2: Create Number Conversion Utilities
**File:** `src/numeric/conversion.h` (new)

```cpp
namespace numeric {
double to_double(const Number& value, const PrecisionContext& context = {});
bool can_convert_to_double(const Number& value);
}
```

**Verification:** Build and run tests.

### Phase 2: Symbolic System Migration

#### Step 2.1: Update Node Accessors
**File:** `src/symbolic/symbolic_expression_internal.h`

Add unified accessor:
```cpp
numeric::Number get_number_value() const;
```

**Verification:** Build and run tests.

#### Step 2.2: Migrate simplify.cpp
**File:** `src/symbolic/simplify.cpp`

Replace all `mymath::sin/cos/exp/ln/sqrt/abs` calls with `numeric::` equivalents.

Pattern change:
```cpp
// Before:
return SymbolicExpression::number(mymath::sin(numeric));

// After:
return SymbolicExpression::number(numeric::sin(numeric::Number(numeric), ctx));
```

**Verification:** Build and run tests.

#### Step 2.3: Migrate transforms.cpp
**File:** `src/symbolic/transforms.cpp`

- Replace `factorial_double` with `factorial_number`
- Replace `mymath::exp/sqrt/kPi` with `numeric::` equivalents
- Change double coefficients to Number

**Verification:** Build and run tests.

#### Step 2.4: Migrate symbolic_expression_calculus.cpp
**File:** `src/symbolic/symbolic_expression_calculus.cpp`

- Change `vector<double>` coefficients to `vector<numeric::Number>`
- Replace mymath calls with numeric equivalents
- Update `solve_dense_linear_system` to use Number

**Verification:** Build and run tests.

#### Step 2.5: Migrate algebra_helpers.cpp and polynomial_helpers.cpp
**Files:** `src/symbolic/algebra_helpers.cpp`, `src/symbolic/polynomial_helpers.cpp`

- Update `decompose_linear` to return Number
- Update polynomial coefficient functions
- Replace mymath calls

**Verification:** Build and run tests.

#### Step 2.6: Migrate node_parser.cpp
**File:** `src/symbolic/node_parser.cpp`

- Update number formatting to use Number
- Replace mymath calls in function evaluation

**Verification:** Build and run tests.

### Phase 3: Analysis Module Migration

#### Step 3.1: Migrate function_analysis.cpp
**File:** `src/analysis/function_analysis.cpp`

- Replace `mymath::abs/sqrt` with numeric equivalents
- Keep numerical algorithms (Gauss-Kronrod, Richardson) using long double internally
- Return Number at API boundary

**Verification:** Build and run tests.

#### Step 3.2: Migrate ode_solver.cpp
**File:** `src/analysis/ode_solver.cpp`

- Keep RKF45 algorithm using double internally (numerical method)
- Convert Number input/output at API boundary
- Replace mymath calls

**Verification:** Build and run tests.

### Phase 4: Matrix Module Migration

#### Step 4.1: Create NumberMatrix
**File:** `src/matrix/number_matrix.h` (new)

```cpp
struct NumberMatrix {
    std::size_t rows = 0;
    std::size_t cols = 0;
    std::vector<numeric::Number> data;
    numeric::PrecisionContext context;

    static NumberMatrix from_double_matrix(const Matrix& m);
    Matrix to_double_matrix() const;

    // Basic operations
    static NumberMatrix add(const NumberMatrix& lhs, const NumberMatrix& rhs);
    static NumberMatrix multiply(const NumberMatrix& lhs, const NumberMatrix& rhs);
    // ... etc
};
```

**Verification:** Build.

#### Step 4.2: Implement NumberMatrix Operations
**File:** `src/matrix/number_matrix.cpp` (new)

Implement basic operations using `numeric::add/multiply/divide`.

**Verification:** Build and add unit tests.

#### Step 4.3: Implement NumberMatrix Linear Algebra
**File:** `src/matrix/number_matrix.cpp`

- QR decomposition
- SVD
- Eigenvalues
- Matrix inverse

Use precision-controlled iteration.

**Verification:** Build and run tests.

#### Step 4.4: Update evaluator.cpp
**File:** `src/expression/evaluator.cpp`

Replace `RealMatrix` (long double) with `NumberMatrix` for exact computation paths.

**Verification:** Build and run tests.

### Phase 5: Cleanup

#### Step 5.1: Remove mymath Dependencies
Remove all `#include "mymath.h"` from:
- `src/symbolic/*.cpp`
- `src/analysis/*.cpp`
- `src/matrix/*.cpp`
- `src/expression/*.cpp`

**Verification:** Build and run tests.

#### Step 5.2: Delete Legacy Files
Delete:
- `src/math/mymath.h`
- `src/math/mymath.cpp`
- `src/math/mymath_special_functions.cpp`
- `src/core/decimal_parser.cpp`
- `src/core/exact_and_symbolic_render.cpp`
- `src/core/precise_decimal_parser.cpp`

**Verification:** Build and run tests.

#### Step 5.3: Remove Double Fields from Node
**File:** `src/symbolic/symbolic_expression_internal.h`

Remove `double number_value` field, keep only `numeric::Number exact_number_value`.

**Verification:** Build and run tests.

## Critical Files

| File | Purpose |
|------|---------|
| `src/symbolic/symbolic_expression_internal.h` | Node storage, hybrid already present |
| `src/symbolic/polynomial_helpers.cpp` | Partial Number polynomial support exists |
| `src/numeric/functions.cpp` | Mathematical functions, needs completion |
| `src/matrix/matrix.h` | Matrix struct, needs Number variant |
| `src/symbolic/simplify.cpp` | Primary mymath user |

## Verification

After each step:
1. `make clean && make`
2. `make test` - verify 810 tests pass
3. Check for performance regressions

Final verification:
1. All 810 tests pass
2. No mymath dependencies remain
3. Legacy files deleted
4. Build succeeds without warnings
