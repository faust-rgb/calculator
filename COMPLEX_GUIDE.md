# Complex Guide

Complex support is split between the legacy compatibility layer and the 2.0
numeric runtime.

## 2.0 Complex Values

Use `i` directly:

```text
:v2 z = 3 + 4i
:v2 conj(z)
:v2 real(z)
:v2 imag(z)
:v2 mat(2, 2, 1, i, 2i, 3) * mat(2, 2, 1, 2, 3, 4)
```

Any input containing a standalone `i` is routed to the 2.0 runtime, so `z = 3 +
4i` also works from the normal prompt.

2.0 complex arithmetic stores real and imaginary parts as `BigDecimal` values
rather than `std::complex<double>`.

## Polynomial Complex Roots

Use `roots_complex(name)` or `complex_roots(name)` on a polynomial custom
function:

```text
c(x) = x ^ 2 + 1
roots_complex(c)
```

`roots(name)` keeps the legacy behavior and returns real roots only.

## Legacy Complex Helpers

The 1.0 path also has lightweight helpers such as `complex(re, im)`, `real`,
`imag`, `abs`, `arg`, `conj`, and `polar`. Those helpers are still backed by the
legacy matrix/scalar representation and may round through `double`.
