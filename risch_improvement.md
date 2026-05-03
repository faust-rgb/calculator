# Implementation Plan - Complete Risch Decision Procedure

This plan aims to improve the current Risch algorithm implementation in the calculator project to make it more robust, recursive, and closer to a complete decision procedure for transcendental extensions.

## Objective
Fix identified issues in the Risch algorithm:
1.  Implement a truly recursive integration framework.
2.  Ensure algebraic independence in the differential tower.
3.  Improve Risch Differential Equation (RDE) and Parametric RDE solvers.
4.  Robustly handle logarithmic and exponential extensions.

## Key Files & Context
- `src/symbolic/risch_algorithm.h`: Interface for the Risch algorithm.
- `src/symbolic/risch_algorithm.cpp`: Implementation of the algorithm.
- `src/symbolic/symbolic_polynomial.h/cpp`: Polynomial operations (GCD, Resultant, etc.).
- `src/symbolic/integration_engine.cpp`: Integration engine that calls the Risch algorithm.

## Implementation Steps

### 1. Refactor `RischAlgorithm::integrate` to be Recursive
- Current implementation has some recursion but it's not consistent.
- The new `integrate` will:
    a. Build the differential tower $x = t_0, t_1, \dots, t_n$.
    b. Define a recursive helper `integrate_in_extension(expression, tower_index)`.
    c. If `tower_index == 0`, it's the rational case (base case).
    d. Otherwise, treat expression as $f \in K(t_n)$ where $K = \mathbb{Q}(t_0, \dots, t_{n-1})$.
    e. Apply Risch reduction for $t_n$ (Hermite reduction, Rothstein-Trager, and polynomial part integration).
    f. Polynomial part integration for $t_n$ will recursively call `integrate_in_extension` for coefficients in $K$ and solve RDEs.

### 2. Improve Differential Tower Construction
- Implement `check_algebraic_independence(extension, current_tower)`.
- For $t = \ln(u)$, check if $u'/u$ has an integral in $K$. If $\int u'/u = v \in K$, then $t = v + c$ is not a new transcendental extension.
- For $t = \exp(u)$, check if $u' = v'/v$ for some $v \in K$. If $u' = v'/v$, then $t = c \cdot v$ is not a new transcendental extension.

### 3. Enhance RDE Solver (`solve_rde`)
- Implement proper degree bounding for $y' + fy = g$.
- For $t = \exp(u)$, handle the case where $f$ has a term $\alpha u'$.
- Implement the "Special Polynomial" case for exponential extensions (poles at $t=0$).

### 4. Improve Parametric RDE Solver (`solve_parametric_rde`)
- Current implementation is only for polynomials.
- Extend it to handle rational functions if possible, or at least be more robust in the polynomial case.

### 5. Robust Transcendental Reductions
- **Logarithmic Extension ($t = \ln(u)$)**:
    - $\int \sum a_i t^i dx$: $\int a_n t^n = \frac{a_n}{n+1} t^{n+1} + \dots$ is wrong. 
    - Correct: $\int a t^n = y t^n + \int (a - y' - n y \frac{u'}{u}) t^n$? No.
    - Actually, the standard Risch says: $\int a t^n = (\int a) t^n - \int (\int a) n t^{n-1} \frac{u'}{u} dx$. This requires $(\int a) \in K$.
    - If $(\int a) \notin K$, then the integral might not be elementary or might require a higher power of $t$.
- **Exponential Extension ($t = \exp(u)$)**:
    - $\int a t^n dx = y t^n$ where $y' + n u' y = a$. Solve RDE for each $i$.

### 6. Verification & Testing
- Add more complex test cases in `test/suites/test_risch.cpp`.
- Test nested extensions: $\ln(\ln(x))$, $\exp(\exp(x))$.
- Test mixed extensions: $\exp(x) \ln(x)$, $\frac{\exp(x)}{\ln(x)}$ (non-elementary).
- Test algebraic independence: $\ln(x^2)$ should be reduced to $2\ln(x)$.

## Migration & Rollback
- The changes are internal to `RischAlgorithm`.
- Revert by restoring `src/symbolic/risch_algorithm.cpp/h` from git history.

## Success Criteria
- Passes all existing tests.
- Successfully integrates $\ln(\ln(x))$ and $\exp(\exp(x))$ (if they have elementary integrals).
- Correctly identifies non-elementary integrals like $\int \frac{e^x}{x} dx$ and returns `false`.
