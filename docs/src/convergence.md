# [Experimental convergence order of MPRK schemes](@id convergence_mprk)

In this tutorial, we check that all implemented MPRK schemes can achieve their expected order of convergence.
We also address the issue that some methods suffer from order reduction when the solution gets too close to zero.

## Conservative production-destruction systems

First, we consider conservative production-destruction systems (PDS). To investigate the convergence order, we define the non-autonomous test problem 

```math
\begin{aligned}
u_1' &= \cos(\pi t)^2 u_2 - \sin(2\pi t)^2 u_1, & u_1(0)&=0.9, \\
u_2' & = \sin(2\pi t)^2 u_1 - \cos(\pi t)^2 u_2, & u_2(0)&=0.1,
\end{aligned}
```
for ``0≤ t≤ 1``.
The PDS is conservative since the sum of the right-hand side terms equals zero. 
An implementation of this problem is given next.


```@example eoc
using PositiveIntegrators

# define problem
P(u, p, t) = [0.0 cos.(π * t) .^ 2 * u[2]; sin.(2 * π * t) .^ 2 * u[1] 0.0]
prob = ConservativePDSProblem(P, [0.9; 0.1], (0.0, 1.0))

nothing # hide
```

To use `analyticless_test_convergence` from [DiffEqDevTools.jl](https://github.com/SciML/DiffEqDevTools.jl), we need to pick a solver to compute the reference solution and specify tolerances.
Since the problem is not stiff, we use the high-order explicit solver `Vern9()` from [OrdinaryDiffEqVerner.jl](https://docs.sciml.ai/OrdinaryDiffEq/stable/).
```@example eoc
using OrdinaryDiffEqVerner
using DiffEqDevTools: analyticless_test_convergence

# solver and tolerances to compute reference solution
test_setup = Dict(:alg => Vern9(), :reltol => 1e-14, :abstol => 1e-14)
nothing # hide
```

To keep the code short, we also define an auxiliary function that outputs a convergence table, which lists the errors obtained with the respective time step size ``Δ t`` as well as the estimated order of convergence in parentheses.

```@example eoc
using Printf: @sprintf
using PrettyTables: pretty_table

# auxiliary function
function convergence_table(dts, prob, algs, labels, test_setup)
    # compute errors and estimated convergence orders
    err_eoc = []
    for i in eachindex(algs)
        sim = analyticless_test_convergence(dts, prob, algs[i], test_setup)

        err = sim.errors[:l∞]
        eoc = [NaN; -log2.(err[2:end] ./ err[1:(end - 1)])]

        push!(err_eoc, tuple.(err, eoc))
    end

    # gather data for table
    data = hcat(dts, reduce(hcat,err_eoc))

    # print table
    formatter = (v, i, j) ->  (j>1) ? (@sprintf "%5.2e (%4.2f) " v[1] v[2]) : (@sprintf "%5.2e " v)
    pretty_table(data; formatters = [formatter], column_labels = ["Δt"; labels]) 
end

nothing # hide
```

### Second-order and third-order MPRK schemes

First, we test several second-order and third-order MPRK schemes.

```@example eoc
# choose step sizes
dts = 0.5 .^ (5:10)

# select 2nd order schemes
algs2a = [MPRK22(0.5); MPRK22(2.0 / 3.0); MPRK22(1.0)]
labels2a = ["MPRK22(0.5)"; "MPRK22(2.0/3.0)"; "MPRK22(1.0)"]
algs2b = [SSPMPRK22(0.5, 1.0); MPDeC(2)]
labels2b = ["SSPMPRK22(0.5, 1.0)"; "MPDeC(2)"]

# select 3rd order schemes
algs3a = [MPRK43I(1.0, 0.5); MPRK43I(0.5, 0.75)]
labels3a = ["MPRK43I(1.0,0.5)"; "MPRK43I(0.5, 0.75)"]
algs3b = [MPRK43II(0.5); MPRK43II(2.0 / 3.0)]
labels3b = [ "MPRK43II(0.5)"; "MPRK43II(2.0/3.0)"]
algs3c = [SSPMPRK43(); MPDeC(3)]
labels3c = ["SSPMPRK43()"; "MPDeC(3)"]

convergence_table(dts, prob, algs2a, labels2a, test_setup)
convergence_table(dts, prob, algs2b, labels2b, test_setup)

convergence_table(dts, prob, algs3a, labels3a, test_setup)
convergence_table(dts, prob, algs3b, labels3b, test_setup)
convergence_table(dts, prob, algs3c, labels3c, test_setup)
```

The tables show that all schemes converge as expected.

### Higher-order MPRK schemes

To actually see the order of higher-order methods we need to use more accurate floating-point numbers. Here, we use [`DoubleFloats`](https://github.com/JuliaMath/DoubleFloats.jl).

```@example eoc
using DoubleFloats 

# define problem using Double64
P(u, p, t) = [0 cospi(t)^2 * u[2]; sinpi(2 * t)^2 * u[1] 0]
u0 = [Double64(9) / 10; Double64(1) / 10]
tspan = (Double64(0), Double64(1))
prob_d64 = ConservativePDSProblem(P, u0, tspan)

# choose step sizes
dts_d64 = Double64(1/2) .^ (5:9)

# select higher-order schemes
algs4a = [MPDeC(4); MPDeC(5); MPDeC(6)]
labels4a = ["MPDeC(4)"; "MPDeC(5)"; "MPDeC(6)"]
algs4b = [MPDeC(7); MPDeC(8)]
labels4b = ["MPDeC(7)"; "MPDeC(8)"]
algs4c = [MPDeC(9); MPDeC(10)]
labels4c = ["MPDeC(9)"; "MPDeC(10)"]

# solver and tolerances to compute reference solution
test_setup_d64 = Dict(:alg => Vern9(), :reltol => 1e-30, :abstol => 1e-30)

# compute errors and experimental order of convergence
convergence_table(dts_d64, prob_d64, algs4a, labels4a, test_setup_d64)
convergence_table(dts_d64, prob_d64, algs4b, labels4b, test_setup_d64)
convergence_table(dts_d64, prob_d64, algs4c, labels4c, test_setup_d64)
```

Again, all schemes show the expected converge order.

## Non-conservative PDS

Next, we consider the non-autonomous but also non-conservative test problem 

```math
\begin{aligned}
u_1' &= \cos(\pi t)^2 u_2 - \sin(2\pi t)^2 u_1 - \cos(2\pi t)^2 u_1, & u_1(0)&=0.9,\\
u_2' & = \sin(2\pi t)^2 u_1 - \cos(\pi t)^2 u_2 - \sin(\pi t)^2 u_2, & u_2(0)&=0.1,
\end{aligned}
```

for ``0≤ t≤ 1``.
Since the sum of the right-hand side terms does not vanish, the PDS is indeed non-conservative.
Hence, we need to use [`PDSProblem`](@ref) for its implementation.

```@example eoc
# PDS
P(u, p, t) = [0.0 cospi(t)^2 * u[2]; sinpi(2 * t)^2 * u[1] 0.0]
D(u, p, t) = [cospi(2 * t)^2 * u[1]; sinpi(t)^2 * u[2]]
prob = PDSProblem(P, D, [0.9; 0.1], (0.0, 1.0))

nothing # hide
```

The following tables demonstrate that the chosen MPRK schemes converge as expected also for this non-conservative PDS.

### Second-order and third-order MPRK schemes

```@example eoc
convergence_table(dts, prob, algs2a, labels2a, test_setup)    
convergence_table(dts, prob, algs2b, labels2b, test_setup)

convergence_table(dts, prob, algs3a, labels3a, test_setup)
convergence_table(dts, prob, algs3b, labels3b, test_setup)
convergence_table(dts, prob, algs3c, labels3c, test_setup)
```

### Higher-order MPRK schemes

```@example eoc
# problem implementation using DoubleFloats
P(u, p, t) = [0 cospi(t)^2 * u[2]; sinpi(2 * t)^2 * u[1] 0]
D(u, p, t) = [cospi(2 * t)^2 * u[1]; sinpi(t)^2 * u[2]]
prob_d64 = PDSProblem(P, D, [Double64(9)/10; Double64(1)/10], (Double64(0), Double64(1)))

convergence_table(dts_d64, prob_d64, algs4a, labels4a, test_setup_d64)
convergence_table(dts_d64, prob_d64, algs4b, labels4b, test_setup_d64)
convergence_table(dts_d64, prob_d64, algs4c, labels4c, test_setup_d64)
```

## Order reduction

It was shown in [Torlo, Öffner, Ranocha: Issues with positivity-preserving Patankar-type schemes with positivity-preserving Patankar-type schemes](https://doi.org/10.1016/j.apnum.2022.07.014) that some MPRK methods 
suffer from order reduction if the solution of the PDS is too close to zero.
We demonstrate this by solving a problem where one component of the initial condition is equal to zero. 

The problem is

```math
\begin{aligned}
u_1' &= -u_1, & u_1(0)&=1, \\
u_2' & = u_1, & u_2(0)&=0,
\end{aligned}
```

for ``0≤ t≤ 1`` and can be implemented as follows.


```@example eoc
# PDS
P(u, p, t) = [0 0; u[1] 0]
prob = ConservativePDSProblem(P, [1.0; 0.0], (0.0, 1.0))
nothing # hide
```

Next, we generate the corresponding convergence tables as in the sections above.

```@example eoc
test_setup = Dict(:alg => Vern9(), :reltol => 1e-14, :abstol => 1e-14)

dts = 0.5 .^ (6:12)

convergence_table(dts, prob, algs2a, labels2a, test_setup)
convergence_table(dts, prob, algs2b, labels2b, test_setup)

convergence_table(dts, prob, algs3a, labels3a, test_setup)
convergence_table(dts, prob, algs3b, labels3b, test_setup)
convergence_table(dts, prob, algs3c, labels3c, test_setup)

convergence_table(dts, prob, algs4a, labels4a, test_setup)
convergence_table(dts, prob, algs4b, labels4b, test_setup)
convergence_table(dts, prob, algs4c, labels4c, test_setup) 
nothing # hide
```

We find that all methods apart from MPDeC(``K``) methods with ``K ≥ 3`` converge as expected.
The MPDeC(``K``) methods with ``K ≥ 3`` suffer from order reduction and show convergence order 2 instead of ``K``.
