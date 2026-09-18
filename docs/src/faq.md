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

## How can I set up fair performance comparisons between PDS and standard SciML solvers?

When benchmarking PDS algorithms (such as `MPRK22`) against standard implicit SciML integrators (such as `ROS2`), providing structural Jacobian information—like a `Tridiagonal` matrix prototype—is crucial for a fair comparison.

Without a `jac_prototype`, implicit solvers fall back to dense $N \times N$ finite-difference Jacobians, scaling with $\mathcal{O}(N^3)$ computational complexity.

### Setup (1D Heat Equation)

First, we define the spatial grid, boundary conditions, and right-hand side functions for $N = 2000$ grid points:


```@example heat_benchmark
using PositiveIntegrators
using OrdinaryDiffEqRosenbrock
using LinearAlgebra
using SciMLBase
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
    for i in 2:(length(u) - 1)
        P[i, i - 1] = u[i - 1] * μ_Δx2
        P[i, i + 1] = u[i + 1] * μ_Δx2
    end
    P[end, end - 1] = u[end - 1] * μ_Δx2
    return nothing
end

# Problem parameters
x_boundaries = range(0, 1, length = 2001)
x = x_boundaries[1:(end - 1)] .+ step(x_boundaries) / 2
u0 = @. cospi(x)^2
tspan = (0.0, 1.0)
μ = 1.0e-2

# Tridiagonal prototype matching physical coupling
p_prototype = Tridiagonal(ones(eltype(u0), length(u0) - 1),
                          ones(eltype(u0), length(u0)),
                          ones(eltype(u0), length(u0) - 1));
```                          
```@example heat_benchmark
# 1. PDS Problem with explicit MPRK22
prob1 = ConservativePDSProblem(heat_eq_P!, u0, tspan, μ; 
                               p_prototype = p_prototype)
t1 = @belapsed solve($prob1, MPRK22(1.0))

# 2. PDS Problem solved with ROS2 (plain function std_rhs / Dense fallback)
t2 = @belapsed solve($prob1, ROS2(autodiff = AutoFiniteDiff()))

# 3. PDS Problem solved with ROS2 (plain function std_rhs / Dense fallback)
prob2 = ConservativePDSProblem(heat_eq_P!, u0, tspan, μ; 
                                        p_prototype = p_prototype, std_rhs = heat_eq_f!)
t3 = @belapsed solve($prob2, ROS2())

# 4. PDS Problem solved with ROS2 (ODEFunction with jac_prototype delegation)
prob3 = ConservativePDSProblem( heat_eq_P!, u0, tspan, μ; 
                                p_prototype = p_prototype, 
                                std_rhs = ODEFunction(heat_eq_f!; jac_prototype = p_prototype))
t4 = @belapsed solve($prob3, ROS2())

println("1. MPRK22 (PDS):                               $(round(t1, digits=4))s")
println("2. ROS2 (PDS without std_rhs):                 $(round(t2, digits=4))s")
println("3. ROS2 (PDS with std_rhs):                    $(round(t3, digits=4))s")
println("4. ROS2 (PDS with std_rhs and jac_prototype):  $(round(t4, digits=4))s")
```