# [Tutorial: Positive-projection method](@id tutorial-sandu)

This tutorial is about solving an ODE using the projection method introduced by Adrian Sandu in [Positive Numerical Integration Methods for Chemical Kinetic Systems](https://doi.org/10.1006%2Fjcph.2001.6750). It guarantees positivity by solving an optimization problem while preserving all linear invariants.

The Sandu projection is a post-processing technique that can be used in combination with any ODE solver.
If the ODE solver computes a negative approximation at any time step, the projection method calculates a positive approximation, also taking into account the linear invariants.

## Solution of the ODE system

As an example we want to the solve the NPZD problem [`prob_pds_npzd`](@ref), which is an ODE system in which negative approximations quickly lead to unacceptable solutions. First, we solve the problem without Sandu projection and select `ROS2` form [`OrdinaryDiffEq.jl`](https://docs.sciml.ai/OrdinaryDiffEq/stable/) as ODE solver.

```@example Sandu_NPZD
using PositiveIntegrators
using OrdinaryDiffEqRosenbrock
using Plots

prob = prob_pds_npzd

ref_sol = solve(prob, ROS2(); abstol = 1e-8, reltol = 1e-6); # reference solution for plotting

sol = solve(prob, ROS2(); abstol = 5e-2, reltol = 1e-1)

plot(ref_sol, linestyle = :dash, label = "", color = palette(:default)[1:4]')
plot!(sol, ylims = (-2.5, 12.5), denseplot = false,  markers = :circle, linewidth = 2, color = palette(:default)[1:4]', label = ["N" "P" "Z" "D"], legend = :right)
```

The plot shows the numerical solution obtained with `ROS2` compared to a reference solution (dashed lines).
We see that the `ROS2` method produces negative approximations, which can occur because Rosenbrock methods are not positivity-preserving. For the NPZD problem, however, this is fatal and leads to a completely unacceptable numerical solution. It is therefore particularly important to use techniques that guarantee positivity of the numerical approximations for this problem. We achieve this below with the [`SanduProjection`](@ref)

To apply the [`SanduProjection`](@ref) we need to choose an [optimization solver](https://jump.dev/JuMP.jl/stable/installation/#Supported-solvers) which is supported by [JuMP.jl](https://jump.dev/JuMP.jl/stable/) and can handle quadratic optimization problems (QP). In this tutorial we select [Clarabel.jl](https://clarabel.org/stable/) as optimization solver.

In addition, we need to specify the linear invariants of the problem. 
The only linear invariant of the NPZD problem is ``N(t)+P(t)+Z(t)+D(t)=N(0)+P(0)+Z(0)+D(0)=15`` for all times ``t≥0``.
This can be written in the form 
```math
\mathbf{A}^T \begin{pmatrix} N(t)\\ P(t)\\ Z(t)\\ D(t) \end{pmatrix} = \mathbf{b}
```
with ``\mathbf{A}^T = [1.0,\  1.0,\  1.0,\  1.0]`` and ``\mathbf{b} = [15]``.

The projection method [`SanduProjection`](@ref) is implemented as a callback and hence, must be passed as an argument to the keyword `callback`. In addition, we must also use `save_everystep = false`.

```@example Sandu_NPZD
using JuMP, Clarabel

AT = [1.0 1.0 1.0 1.0]
b = [15.0]
proj = SanduProjection(Model(Clarabel.Optimizer), AT, b)

sol_proj = solve(prob, ROS2(); abstol = 5e-2, reltol = 1e-1
                 save_everystep = false, callback = proj);

plot(ref_sol, linestyle = :dash, label = "", color = palette(:default)[1:4]')
plot!(sol_proj, ylims = (-2.5, 12.5), denseplot = false,  markers = :circle, linewidth = 2, color = palette(:default)[1:4]', label = ["N" "P" "Z" "D"], legend = :right)            
```

As intended, negative approximations no longer occur and we obtain an acceptable approximation.

The [`SanduProjection`](@ref) is implemented as a [`DiscreteCallback`](https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions/#SciMLBase.DiscreteCallback) and we can display the number of projection steps in the following way.

```@example Sandu_NPZD
@show get_numsteps_SanduProjection(proj) 
```

We can see that in this example, a single projection step was already sufficient.

## Package versions

These results were obtained using the following versions.
```@example NPZD
using InteractiveUtils
versioninfo()
println()

using Pkg
Pkg.status(["PositiveIntegrators", "JuMP", "Clarabel", "OrdinaryDiffEqRosenbrock", "Plots"],
           mode=PKGMODE_MANIFEST)
nothing # hide
```
