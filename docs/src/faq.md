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
end; 
```

## How can I set up fair performance comparisons between PDS and standard SciML solvers?

When benchmarking PDS algorithms (such as `MPRK22`) against standard implicit SciML integrators (such as `ROS2`), it is crucial for a fair comparison to ensure that both algorithms receive equivalent structural information, since using fallback settings for standard solvers can lead to misleading runtime comparisons.

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

alg1 = MPRK22(1.0)
alg2 = ROS2()                          
```  
### Comparisons

In the initial setup, we only use the production matrix function `heat_eq_P!` to create the PDS. The package automatically generates the standard ODE right-hand side, necessary for standard solvers, under the hood by summing over the production terms.

```@example heat_benchmark
prob1 = ConservativePDSProblem(heat_eq_P!, u0, tspan, μ)
t11 = @belapsed solve($prob1, $alg1; save_everystep = false)
t12 = @belapsed solve($prob1, $alg2; save_everystep = false)
(t11, t12)
```
In this configuration, `ROS2()` performs significantly worse than MPRK22(1.0).
This severe slowdown is driven by two compounding factors: The standard right-hand side is automatically constructed by summing over the production matrix elements and the automatic differentiation (`ForwardDiff`) of this auto-generated fallback right-hand side.

A first step to increase the performance of `ROS2()` is to provide an ODE right-hand side explicitly.

```@example heat_benchmark
prob2 = ConservativePDSProblem(heat_eq_P!, u0, tspan, μ; std_rhs = heat_eq_f!)
t21 = @belapsed solve($prob2, $alg1; save_everystep = false)
t22 = @belapsed solve($prob2, $alg2; save_everystep = false)
(t21, t22)
```
By explicitly supplying `heat_eq_f!`, the runtime of `ROS2()` drops drastically, since providing a standard right-hand side removes the overhead of the auto-generated summation fallback. 

However, `MPRK22(1.0)` as well as `ROS2()` are still solving linear systems without any structural sparsity information.
There is room for optimization by providing sparsity information. First, we define `p_prototype`.

```@example heat_benchmark
prob3 = ConservativePDSProblem(heat_eq_P!, u0, tspan, μ; p_prototype = p_prototype, std_rhs = heat_eq_f!)
t31 = @belapsed solve($prob3, $alg1; save_everystep = false)
t32 = @belapsed solve($prob3, $alg2; save_everystep = false)
(t31, t32)
```
With the `p_prototype` supplied, the execution time of `MPRK22(1.0)` drops significantly.

However, `ROS2()` remains unaffected. In order to let `ROS2()` benefit from sparsity, we must provide `std_rhs` as an `ODEFunction` which provides the sparsity information by specifying `jac_prototype`.

```@example heat_benchmark
prob4 = ConservativePDSProblem(heat_eq_P!, u0, tspan, μ; p_prototype = p_prototype, std_rhs = ODEFunction(heat_eq_f!; jac_prototype = p_prototype))
t41 = @belapsed solve($prob4, $alg1; save_everystep = false)
t42 = @belapsed solve($prob4, $alg2; save_everystep = false)
(t41, t42)
```
Now both solvers receive equivalent structural information and the runtime of `ROS2()` also drops dramatically.