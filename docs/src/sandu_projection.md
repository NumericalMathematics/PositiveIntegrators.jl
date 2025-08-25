# [Tutorial: Positive-projection method](@id tutorial-sandu)

This tutorial is about solving an ODE using the positive-projection method introduced by Adrian Sandu in [Positive Numerical Integration Methods for Chemical Kinetic Systems](https://doi.org/10.1006%2Fjcph.2001.6750). It guarantees positivity through the solution of an optimization problem while preserving all linear invariants.


## Definition of the positive-projection method

The idea is to project the numerical solutions ``y^n`` of an ODE system obtained by a linear-preserving-one step integration method if they are negative. To guarantee that the projected value ``z^{n+1}`` preserves all linear invariants, defined by ``A^T y^{n} = b`` for all ``n`` and is near to the true solution, we set up the following quadratic optimization problem
```math
\begin{equation*}
\mathrm{min} \frac12 \||z^{n+1}-y^{n+1}\||_\mathrm{G}^2 \quad \text{subject to} \ A^T z^{n+1} = b, z^{n+1} \geq 0
\end{equation*}
```
with the norm ```$\||y\||_G = \sqrt(y^T G y)$```, where ```G```is a positive definite matrix that depends on the relative and absolute error tolerances used for the underlying method as well as ``y^{n}``.

## Solution of the minimization problem

Now we are ready to define an optimization problem using [JuMP.jl](https://jump.dev/) and solve it with a JuMP.jl [supported solver](https://jump.dev/JuMP.jl/stable/installation/#Supported-solvers). 

First recall the model with its linear invariants. Here, we use the NPZD model [`prob_pds_npzd`](@ref) as an example. The corresponding linear invariant matrix is already provided for each problem in [Example Problems](https://NumericalMathematics.github.io/PositiveIntegrators.jl/dev/api_reference/#Example-problems). Then, we set up the necessary callback. 

```@example Sandu_NPZD
using PositiveIntegrators
using JuMP
using Clarabel

prob = prob_pds_npzd
cb = SanduProjection(Model(Clarabel.Optimizer), prob.f.linear_invariants, prob.f.linear_invariants * prob.u0, 0.0)
nothing
```
In some cases it might be senseful to guarantee the projected solution stays away from zero. This can be set in the last input parameter of the callback.
Finally solve the ODE with an algorithm out of [OrdinaryDiffEq.jl](https://docs.sciml.ai/OrdinaryDiffEq/stable/).
```@example Sandu_NPZD
using OrdinaryDiffEqLowOrderRK
sol = solve(prob, Heun(); abstol = 5e-2, reltol = 1e-1, dt = 0.1,
            save_everystep = false, callback = cb)
nothing
```

```@example Sandu_NPZD
using Plots

plot(sol; label = ["N" "P" "Z" "D"], xguide = "t")
nothing
```


### Performance comparison

Finally, we use [BenchmarkTools.jl](https://github.com/JuliaCI/BenchmarkTools.jl)
to show the benefit of using static arrays.

```@example Sandu_NPZD
using BenchmarkTools
@benchmark solve(prob, Heun(); abstol = 5e-2, reltol = 1e-1, dt = 0.1,
            save_everystep = false, callback = cb)
```

```@example Sandu_NPZD
using BenchmarkTools
@benchmark solve(prob, MPRK22(1.0); abstol = 5e-2, reltol = 1e-1, dt = 0.1)
```

```@example Sandu_NPZD
@benchmark solve(prob, Heun(); abstol = 1e-10, reltol = 1e-9)
```

## Package versions

These results were obtained using the following versions.
```@example NPZD
using InteractiveUtils
versioninfo()
println()

using Pkg
Pkg.status(["PositiveIntegrators", "JuMP", "Clarabel", "OrdinaryDiffEq", "Plots"],
           mode=PKGMODE_MANIFEST)
nothing # hide
```
