# [Tutorial: Solving a scalar production-destruction equation with MPRK schemes](@id tutorial-scalar-pds)

Originally, modified Patankar-Runge-Kutta (MPRK) schemes were designed to solve positive and conservative systems of ordinary differential equations.
The conservation property requires that the system consists of at least two scalar differential equations. 
Nevertheless, we can also apply the idea of the Patankar trick to a scalar production-destruction system (PDS)

```math
u'(t)=p(u(t))-d(u(t)),\quad u(0)=u_0>0
```

with nonnegative functions ``p`` and ``d``.
Since conservation is not an issue here, we can apply the Patankar trick to the destruction term ``d`` to ensure positivity and leave the production term ``p`` unweighted. 
A first-order scheme of this type, based on the forward Euler method, reads

```math
u^{n+1}= u^n + Δ t p(u^n) - Δ t d(u^n)\frac{u^{n+1}}{u^n}
```
and this idea can easily be generalized to higher-order explicit Runge-Kutta schemes. 

By closer inspection we realize that this is exactly the approach the MPRK schemes of [PositiveIntegrators.jl](https://github.com/SKopecz/PositiveIntegrators.jl) use to solve non-conservative PDS for which the production matrix is diagonal. 
Hence, we can use the existing schemes to solve a scalar PDS by regarding the production term as a ``1×1``-matrix and the destruction term as a ``1``-vector.

# [Example 1](@id scalar-example-1)

We want to solve

```math
u' =  u^2 - u,\quad u(0) = 0.95
```

for ``0≤ t≤ 10``.
Here,  we can choose ``p(u)=u^2`` as production term and ``d(u)=u`` as destruction term.
The exact solution of this problem is

```math
u(t) = \frac{19}{19+e^{t}}.
```

Next, we show how to solve this scalar PDS in the way discussed above.
Please note that we must use [`PDSProblem`](@ref) to create the problem.
Furthermore, we use static matrices and vectors from [StaticArrays.jl](https://juliaarrays.github.io/StaticArrays.jl/stable/) instead of standard arrays for efficiency.


```@example scalar_example_1
using PositiveIntegrators, StaticArrays, Plots

u0 = @SVector [0.95] # 1-vector
tspan = (0.0, 10.0)

# Attention: Input u is a 1-vector
prod(u, p, t) = @SMatrix [u[1]^2] # create static 1x1-matrix
dest(u, p, t) = @SVector [u[1]] # create static 1-vector
prob = PDSProblem(prod, dest, u0, tspan) 

sol = solve(prob, MPRK22(1.0))

# plot
tt = 0:0.1:10
f(t) = 19.0 / (19.0 + exp(t)) # exact solution
plot(tt, f.(tt), label="exact")
plot!(sol, label="u")
```

# [Example 2] (@id scalar-example-2)

Next, we want to compute positive solutions of a more challenging scalar PDS. 
In [Example 1](@ref scalar-example-1), we could have also used standard schemes from [OrdinaryDiffEq.jl](https://docs.sciml.ai/OrdinaryDiffEq/stable/) and use the solver option `isoutofdomain` to ensure positivity.
But this is not always the case as the following example will show.

We want to compute the nonnegative solution of 

```math
u'(t) = -\sqrt{\lvert u(t)\rvert },\quad u(0)=1. 
```

for ``t≥ 0``.
Please note that this initial value problem has infinitely many solutions

```math 
u(t) = \begin{cases} \frac{1}{4}(t-2)^2, & 0≤ t< 2,\\ 0, & 2≤ t < t^*,\\ -\frac{1}{4}(t-2)^2, & t^*≤  t, \end{cases}
```

where ``t^*≥ 2`` is arbitrary.
But among these, the only nonnegative solution is

```math 
u(t) = \begin{cases} \frac{1}{4}(t-2)^2, & 0≤ t< 2,\\ 0, & 2≤ t. \end{cases}
```

This is the solution we want to compute.

First, we try this using a standard solver from [OrdinaryDiffEq.jl](https://docs.sciml.ai/OrdinaryDiffEq/stable/).
We try to enforce positivity with the solver option `isoutofdomain` by specifying that negative solution components are not acceptable.

```@example
using OrdinaryDiffEqRosenbrock

tspan = (0.0, 3.0)
u0 = 1.0

f(u, p, t) = -sqrt(abs(u))
prob = ODEProblem(f, u0, tspan)

sol = solve(prob, Rosenbrock23(); isoutofdomain = (u, p, t) -> any(<(0), u))
```

We see that `isoutofdomain` cannot be used to ensure nonnegative solutions in this case, as the computation stops at about ``t≈ 2`` before the desired final time is reached. 
For at least first- and second-order explicit Runge-Kutta schemes, this can also be shown analytically. A brief computation reveals that to ensure nonnegative solutions, the time step size must tend to zero if the numerical solution tends to zero.

Next, we want to use an MPRK scheme. 
We can choose ``p(u)=0`` as the production term and ``d(u)=\sqrt{\lvert u\rvert }`` as the destruction term. 
Furthermore, we create the [`PDSProblem`](@ref) in the same way as in [Example 1](@ref scalar-example-1).

```@example
using PositiveIntegrators, StaticArrays, Plots

tspan = (0.0, 3.0)
u0 = @SVector [1.0]

prod(u, p, t) = @SMatrix zeros(1,1)
dest(u, p, t) = @SVector [sqrt(abs(first(u)))]
prob = PDSProblem(prod, dest, u0, tspan)

sol = solve(prob, MPRK22(1.0))

# plot
tt = 0:0.03:3
f(t) = 0.25 * (t - 2)^2 * (t <= 2) # exact solution
plot(tt, f.(tt), label="exact")
plot!(sol, label="u")
```

We can see that the MPRK scheme used is well suited to solve the problem. 


## Package versions

These results were obtained using the following versions.
```@example scalar_example_1
using InteractiveUtils
versioninfo()
println()

using Pkg
Pkg.status(["PositiveIntegrators", "StaticArrays", "OrdinaryDiffEqRosenbrock"],
           mode = PKGMODE_MANIFEST)
nothing # hide
```
