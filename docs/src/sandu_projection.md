# [Tutorial: Positive-projection method](@id tutorial-sandu)

This tutorial is about solving an ODE using the positive-projection method introduced by Adrian Sandu in [Positive Numerical Integration Methods for Chemical Kinetic Systems](https://doi.org/10.1006%2Fjcph.2001.6750). It guarantees positivity by solving an optimization problem while preserving all linear invariants.

The Sandu projection is a post-processing technique which can be used in combination with any ODE solver.
If the ODE solver computes a negative approximation at any time step, the projection method calculates a positive approximation, also taking into account the linear invariants.

## Solution of the ODE system

As an example we want to the solve the NPZD problem [`prob_pds_npzd`](@ref), which is an ODE system in which negative approximations quickly lead to unacceptable solutions. First, we solve the problem without Sandu projection and select `ROS2` form [`OrdinaryDiffEq.jl`](https://docs.sciml.ai/OrdinaryDiffEq/stable/) as ODE solver.

```@example Sandu_NPZD
using PositiveIntegrators
using OrdinaryDiffEqRosenbrock
using Plots

prob = prob_pds_npzd

ref_sol = solve(prob, ROS2(); abstol = 1e-12, reltol = 1e-10); # reference solution for plotting

sol = solve(prob, ROS2(); abstol = 5e-2, reltol = 1e-1, dt = 0.1) 

plot(ref_sol, linestyle = :dash, label = "", color = palette(:default)[1:4]')
plot!(sol, ylims = (-2.5, 12.5), denseplot = false,  markers = :circle, linewith = 2, color = palette(:default)[1:4]', label = ["N" "P" "Z" "D"], legend = :right)
nothing
```

The plot shows the solution obtained by `ROS2` compared to a reference solution (dashed lines).
We see that as soon as negative values of the ``N`` species occur, ``N`` continues to decrease and the solution becomes completely unacceptable.

Now, we want to avoid negative approximations by using [`SanduProjection`](@ref). For this, we need to choose one of the [supported optimization solvers](https://jump.dev/JuMP.jl/stable/installation/#Supported-solvers) of [JuMP.jl](https://jump.dev/JuMP.jl/stable/) and we select [Clarabel.jl](https://clarabel.org/stable/) for this tutorial.

In addition, we need to specify the linear invariants of the problem. 
The only linear invariant of the NPZD problem is ``N(t)+P(t)+Z(t)+D(t)=N(0)+P(0)+Z(0)D(0)=15`` for all times ``t≥0``.
This can be written in the form ``\\mathbf{A}^T\\begin{pmatrix}N(t)\\\\ P(t)\\\\ Z(t)\\\\ D(t)\\end{pmatrix} = \\mathbf{b}`` with ``\\mathbf{A}^T = [1.0  1.0  1.0  1.0]`` and ``\\mathbf{b} = [15]``.

The projection method [`SanduProjection`](@ref) is implemented as a callback and hence, must be passed as an argument to the keyword `callback`. In addition, we must also use `save_everystep = false`.

```@example Sandu_NPZD
using JuMP, Clarabel

AT = [1.0 1.0 1.0 1.0]
b = [15.0]
cb = SanduProjection(Model(Clarabel.Optimizer), AT, b)

sol_cb = solve(prob, ROS2(); abstol = 5e-2, reltol = 1e-1, dt = 0.1,
            save_everystep = false, callback = cb)

plot(ref_sol, linestyle = :dash, label = "", color = palette(:default)[1:4]')
plot!(sol_cb, ylims = (-2.5, 12.5), denseplot = false,  markers = :circle, linewith = 2, color = palette(:default)[1:4]', label = ["N" "P" "Z" "D"], legend = :right)            
nothing
```

We see that negative approximations no longer occur.

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
