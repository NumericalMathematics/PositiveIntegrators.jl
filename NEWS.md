# Changelog

PositiveIntegrators.jl.jl follows the interpretation of
[semantic versioning (semver)](https://julialang.github.io/Pkg.jl/dev/compatibility/#Version-specifier-format-1)
used in the Julia ecosystem. Notable changes will be documented in this file
for human readability.


## Changes in the v0.2 lifecycle

### Changed

- The minimum required Julia version was updated to v1.10
- A failed linear solve now rejects the time step. With `adaptive = true`, the
  step is repeated with a smaller time step; with `adaptive = false`, the
  integration stops with `ReturnCode.ConvergenceFailure`. Previously, the
  output of the failed linear solve was used as the solution.
- The parameters of `MPRK22`, `MPRKO22`, `MPRK43I`, `MPRK43II`, and `SSPMPRK22`
  are converted to the floating point type of the problem before the
  coefficients of the method are computed. Previously, parameters of a type
  different from the problem's (e.g., `MPRK22(1.0)` for a `Float32` or
  `BigFloat` problem) caused a `MethodError`.
- `MPRK43I` and `SSPMPRK22` promote their two parameters to a common type, so
  that, e.g., `SSPMPRK22(0, 0.75)` no longer throws an `InexactError`.
- For problems with dual numbers of ForwardDiff.jl, the coefficients and small
  constants of all schemes have the type of the primal values (e.g., `Float64`),
  so that the primal values of the solution agree with the solution of the
  corresponding `Float64` problem. Previously, `MPRK22`, `MPRKO22`, `MPRK43I`,
  `MPRK43II`, `SSPMPRK22`, and `SSPMPRK43` threw a `MethodError`, and `MPDeC`
  used a different small constant.
- The in-place `SSPMPRK43` converts its coefficients to the floating point type
  of the problem. Previously, it threw a `MethodError` for problems that do not
  use `Float64` (e.g., `Float32`).
- The in-place `SSPMPRK43` no longer adds the diagonal of the production
  matrix of a `ConservativePDSProblem` to the right-hand side of its third
  stage, so that it agrees with the out-of-place version. The diagonal of the
  production matrix of a `ConservativePDSProblem` is not used.
- `MPDeC` no longer uses the diagonal of the production matrix of a
  `ConservativePDSProblem`, which violated conservation.
- All schemes solve their linear systems for the Patankar weights `u ./ σ`
  (with the Patankar weight denominators `σ`) instead of `u`. The system
  matrices are assembled without divisions by `σ`, so tiny denominators no
  longer produce infinite matrix entries. Before the solve, columns of the
  system matrices are scaled by powers of two: all columns of sparse matrices
  and the columns of dense, `Tridiagonal`, and static matrices with tiny or
  huge diagonal entries or with Patankar weights that could overflow. The
  results agree with the previous formulation up to rounding errors.
- A step attempt is also rejected if a linear solve returns non-finite values.
  With `adaptive = true`, the step is repeated with a smaller time step.
- If a component of the scaled error estimate of an adaptive scheme is not
  finite, the error estimate is set to `Inf`, so that the step size controller
  rejects the step attempt and reduces the time step.
- With `linsolve = nothing`, the schemes use `LUFactorization()` for
  `Tridiagonal` matrices: the default solver of LinearSolve.jl throws an
  exception for some singular `Tridiagonal` matrices instead of reporting
  the failure.

### Added

- MPDeC methods of Öffner and Torlo
- MPRKO22 methods of Ávila et al.
- HIRES test problem `prob_pds_hires`

## Changes when updating to v0.2 from v0.1.x

#### Removed

- The optional keyword argument `d_prototype` has been removed from `PDSProblem`


## Changes in the v0.1 lifecycle

### Added

- MPRK methods `MPRK43I` and `MPRK43II` of Kopecz and Meister
- SSP MPRK methods `SSPMPRK22` and `SSPMPRK43` of Huang and Shu


## Initial release v0.1.0

- Production-destruction problems `PDSProblem` and `ConservativePDSProblem`,
  including conversions to standard `ODEProblem`s from OrdinaryDiffEq.jl
- Some default problems such as `prob_pds_bertolazzi`, `prob_pds_brusselator`,
  `prob_pds_linmod`, `prob_pds_nonlinmod`, `prob_pds_npzd`, `prob_pds_robertson`,
  `prob_pds_sir`, `prob_pds_stratreac`
- Modified Patankar methods `MPE` and `MPRK22`
