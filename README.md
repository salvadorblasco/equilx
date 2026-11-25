# EQUILX 
EQ FORTRAN LIB FOR EQUILIBRIA RESOLUTION
by Salvador Blasco <salvador.blasco@uv.es>

## Non-negative least squares

The Fortran `linalg` module provides `nnls` for both a single right-hand side
and a matrix of spectra. To estimate extinction coefficients, use:

```fortran
use linalg, only: nnls
! C(NPT,NSP): concentrations; A(NPT,NWL): absorbances
! E(NSP,NWL): fitted extinction coefficients
call nnls(C, A, E, info, failed_column)
```

This minimizes the sum of squared residuals in `A - matmul(C, E)` subject to
`E >= 0`, by solving an independent NNLS problem for each wavelength. Zero
coefficients are allowed. The original vector form `call nnls(A, b, x, info)`
remains available.

`info` is 0 on success; -1 for invalid dimensions; -2 for non-finite inputs;
1 for an iteration limit; 2 for numerical rank failure; and 3 for a LAPACK or
numerical failure. The matrix form sets all of `E` to zero on failure. The
optional `failed_column` is the 1-based wavelength index for a solver failure,
or zero on success or input validation failure. Do not use `E` when `info` is
nonzero.
