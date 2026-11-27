# Troubleshooting and frequently asked questions

## Sparse matrices

You can use sparse matrices for the linear systems arising in
[PositiveIntegrators.jl](https://github.com/NumericalMathematics/PositiveIntegrators.jl),
as described, e.g., in the [tutorial on linear advection](@ref tutorial-linear-advection).
However, you need to make sure that you do not change the sparsity pattern
of the production term matrix since we assume that the structural nonzeros
are kept fixed. This is a [known issue](https://github.com/JuliaSparse/SparseArrays.jl/issues/190).
For example, you should avoid something like

```@repl
using SparseArrays
p = spdiagm(0 => ones(4), 1 => zeros(3))
p .= 2 * p
```

Instead, you should be able to use a pattern like the following, where the function `nonzeros` is used to modify the values of a sparse matrix.

```@repl
using SparseArrays
p = spdiagm(0 => ones(4), 1 => zeros(3))
for j in axes(p, 2)
    for idx in nzrange(p, j)
        i = rowvals(p)[idx]
        nonzeros(p)[idx] = 10 * i + j # value p[i, j]
    end
end; p
```

## Why was a time step rejected or the integration stopped?

The schemes reject a step attempt if

- a linear solve does not report success (its return code is neither
  successful nor `ReturnCode.Default`), or
- a linear solve returns non-finite values.

Rejected attempts are counted in `sol.stats.nreject`, and the solution keeps
the last accepted value. With `adaptive = true`, the step is repeated with a
smaller time step (divided by the keyword argument `failfactor` of the solver
interface, default 2). The integration stops if the time step falls below the
keyword argument `dtmin` (`ReturnCode.DtLessThanMin`) or becomes too small to
advance the time (`ReturnCode.Unstable`). With the default `dtmin = 0`, the
latter takes more than 1000 rejected attempts at `t = 0`; set `dtmin` to stop
earlier. With `adaptive = false`, the integration stops at the first rejected
attempt with `ReturnCode.ConvergenceFailure`. The messages follow the setting
of the keyword argument `verbose` of `solve`, e.g., the integration stops
silently with `verbose = DEVerbosity(SciMLLogging.None())`, where `DEVerbosity`
and `SciMLLogging` are exported by OrdinaryDiffEq.jl.

In addition, the step size controller of an adaptive scheme rejects a step
attempt if the error estimate is larger than one. If a component of the scaled
error estimate is not finite, e.g., because a Patankar weight denominator of
the embedded solution overflows, the error estimate is set to `Inf`. Such an
attempt is repeated with a time step reduced by the maximal factor allowed by
the keyword argument `qmin`, and the integration stops under the same
conditions on the time step as above.

Check `sol.retcode` after solving. A failure often indicates a too large
initial time step, production or destruction terms that become non-finite for
some states, or a linear solver that is not suitable for the system matrices.

## How can I set up fair performance comparisons between PDS and standard SciML solvers?

When benchmarking PDS algorithms (such as [`MPRK22`](@ref)) against standard implicit SciML integrators (such as `ROS2`), it is crucial for a fair comparison to ensure that both algorithms receive equivalent structural information, since using fallback settings for standard solvers can lead to misleading runtime comparisons.

Below, we demonstrate step-by-step how setup choices affect runtime using a 1D heat equation discretized via finite differences ($N=500$).

### Setup (1D Heat Equation)

First, we define the spatial grid, boundary conditions, and right-hand side functions for $N = 500$ grid points. In particular, we define both the production-matrix function `heat_eq_P!` and the ODE right-hand side `heat_eq_f!`. Implementing both `heat_eq_P!` and `heat_eq_f!` in-place provides the essential prerequisite for minimizing memory allocations.

```@example heat_benchmark
using PositiveIntegrators
using OrdinaryDiffEqRosenbrock
using LinearAlgebra
using BenchmarkTools

# Standard ODE right-hand side
function heat_eq_f!(du, u, μ, t)
    fill!(du, 0)
    N = length(u)
    Δx = 1 / N
    μ_Δx2 = μ / Δx^2

    du[1] = (-u[1] + u[2]) * μ_Δx2
    for i in 2:(N - 1)
        du[i] = (u[i - 1] - 2 * u[i] + u[i + 1]) * μ_Δx2
    end
    du[N] = (u[N - 1] - u[N]) * μ_Δx2
    return nothing
end

# Production matrix function
function heat_eq_P!(P, u, μ, t)
    fill!(P, 0)
    N = length(u)
    Δx = 1 / N
    μ_Δx2 = μ / Δx^2

    P[1, 2] = u[2] * μ_Δx2
    for i in 2:(N - 1)
        P[i, i - 1] = u[i - 1] * μ_Δx2
        P[i, i + 1] = u[i + 1] * μ_Δx2
    end
    P[end, end - 1] = u[end - 1] * μ_Δx2
    return nothing
end

# Problem parameters
N = 500
x_boundaries = range(0, 1, length = N + 1)
x = x_boundaries[1:(end - 1)] .+ step(x_boundaries) / 2
u0 = @. cospi(x)^2
tspan = (0.0, 1.0)
μ = 1.0e-2

# Tridiagonal prototype matching physical coupling
p_prototype = Tridiagonal(ones(eltype(u0), length(u0) - 1),
                          ones(eltype(u0), length(u0)),
                          ones(eltype(u0), length(u0) - 1));

# Algorithms to compare
alg1 = MPRK22(1.0)
alg2 = ROS2()      
nothing #hide                    
```  
### Comparisons

In the initial setup, we only use the production matrix function `heat_eq_P!` to create the PDS. The package automatically generates the standard ODE right-hand side, necessary for standard solvers, under the hood by summing over the production terms.

```@example heat_benchmark
prob1 = ConservativePDSProblem(heat_eq_P!, u0, tspan, μ)
t11 = @belapsed solve($prob1, $alg1; save_everystep = false)
t12 = @belapsed solve($prob1, $alg2; save_everystep = false)
(t11, t12)
```
In this configuration, `ROS2()` performs significantly worse than `MPRK22(1.0)`.
One reason is that the ODE right-hand side, necessary for standard solvers, must be constructed implicitly by summing over the production matrix elements. Another reason is that this implicit right-hand side is used for automatic differentiation (`ForwardDiff`).

A first step to increase the performance of `ROS2()` is to provide the ODE right-hand side explicitly.

```@example heat_benchmark
prob2 = ConservativePDSProblem(heat_eq_P!, u0, tspan, μ; std_rhs = heat_eq_f!)
t21 = @belapsed solve($prob2, $alg1; save_everystep = false)
t22 = @belapsed solve($prob2, $alg2; save_everystep = false)
(t21, t22)
```
By explicitly supplying `heat_eq_f!`, the runtime of `ROS2()` drops drastically and `ROS2()` becomes even faster than `MPRK22(1.0)`.

However, `MPRK22(1.0)` as well as `ROS2()` are still operating on dense matrices when solving linear systems, since no sparsity information is provided.
First, we specify `p_prototype` to speed up linear solves involving the production matrix. 

```@example heat_benchmark
prob3 = ConservativePDSProblem(heat_eq_P!, u0, tspan, μ; p_prototype = p_prototype, std_rhs = heat_eq_f!)
t31 = @belapsed solve($prob3, $alg1; save_everystep = false)
t32 = @belapsed solve($prob3, $alg2; save_everystep = false)
(t31, t32)
```
With the `p_prototype` supplied, the execution time of `MPRK22(1.0)` drops several orders of magnitude. However, `ROS2()` remains unaffected, since the production matrices play no role in the solution process of `ROS2()`. Instead, `ROS2()` requires the solution of linear systems containing the Jacobian of `std_rhs` and it is the sparsity structure of this Jacobian we must provide.
To do so, we must provide `std_rhs` as an `ODEFunction` which in addition allows us to pass the sparsity information by specifying `jac_prototype`.

```@example heat_benchmark
prob4 = ConservativePDSProblem(heat_eq_P!, u0, tspan, μ; p_prototype = p_prototype, std_rhs = ODEFunction(heat_eq_f!; jac_prototype = p_prototype))
t41 = @belapsed solve($prob4, $alg1; save_everystep = false)
t42 = @belapsed solve($prob4, $alg2; save_everystep = false)
(t41, t42)
```
Providing the `jac_prototype` also drops the runtime of `ROS2()` by several orders of magnitude.

We now have reached a point where the PDS solver efficiently exploits the sparsity of the production matrix using `p_prototype`, while the standard solver operates with both an explicit `std_rhs` and the corresponding Jacobian sparsity pattern, specified by `jac_prototype`, which allows a fair comparison of both schemes.