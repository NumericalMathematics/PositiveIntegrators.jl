# [Benchmark: Solution of the HIRES problem](@id benchmark-hires)

We use the stiff HIRES problem [`prob_pds_hires`](@ref) to assess the efficiency of different solvers from [OrdinaryDiffEq.jl](https://docs.sciml.ai/OrdinaryDiffEq/stable/) and [PositiveIntegrators.jl](https://github.com/NumericalMathematics/PositiveIntegrators.jl).

```@example HIRES
using OrdinaryDiffEqFIRK, OrdinaryDiffEqRosenbrock, OrdinaryDiffEqSDIRK
using PositiveIntegrators

# select HIRES problem
prob = prob_pds_hires
nothing # hide
```

To keep the following code as clear as possible, we define a helper function `hires_plot` that we use for plotting.

```@example HIRES
using Plots

function hires_plot(sols, labels = fill("", length(sols)), sol_ref = nothing)
    if !(sols isa Vector)
        sols = [sols]
    end
    if !(labels isa Vector)
        labels = [labels]
    end

    layout = (4, 2)
    size = (800, 800)
    linewidth = 2
    tickfontsize = 7
    xguide = ["" "" "" "" "" "" "t" "t"]
    yguide = ["u₁" "u₂" "u₃" "u₄" "u₅" "u₆" "u₇" "u₈"]
    legend = [:topright false false false false false false false]
    widen = true

    # Time points for plotting interpolated solutions. They resolve the fast
    # initial transient while keeping the figure files small.
    tspan = prob_pds_hires.tspan
    plotat = [range(first(tspan), 3.0, length = 101);
              range(3.0, last(tspan), length = 301)[2:end]]

    p = plot(; layout, size)
    if !isnothing(sol_ref)
        plot!(p, sol_ref; linestyle = :dash, label = "Ref.", linewidth,
              plotat, color = :black, legend)
    end
    for (sol, label) in zip(sols, labels)
        # Solutions with more time steps than `plotat` are only plotted at
        # `plotat`, which reduces the file size without visible changes.
        points = length(sol.t) > length(plotat) ? (; plotat) : (; denseplot = false)
        plot!(p, sol; xguide, tickfontsize, yguide, legend, widen, linewidth, label,
              points...)
    end
    return p
end
nothing # hide
```

First, we show approximations of `Rosenbrock23()` and `MPRK22(1.0)` using loose tolerances.

```@example HIRES
# compute reference solution for plotting
ref_sol = solve(prob, Rodas4P(); abstol = 1e-14, reltol = 1e-13);

# compute solutions with loose tolerances
abstol = 1e-2
reltol = 1e-1
sol_Ros23 = solve(prob, Rosenbrock23(); abstol, reltol);
sol_MPRK = solve(prob, MPRK22(1.0); abstol, reltol);

# plot solutions
hires_plot([sol_Ros23, sol_MPRK], ["Rosenbrock23", "MPRK22(1.0)"], ref_sol)
```

The `Rosenbrock23` solution contains negative values in the components ``u_5`` and ``u_6``.

```@example HIRES
isnonnegative(sol_Ros23)
```

The `MPRK22(1.0)` solution is nonnegative, but the components ``u_5`` and ``u_6`` decay too slowly.
Neither solution resolves the rapid change of ``u_7`` and ``u_8`` near the final time.

[OrdinaryDiffEq.jl](https://docs.sciml.ai/OrdinaryDiffEq/stable/) provides the solver option `isoutofdomain`, which can be used in combination with [`isnegative`](@ref) to guarantee nonnegative solutions.

```@example HIRES
# compute solution with isoutofdomain = isnegative
sol_Ros23 = solve(prob, Rosenbrock23(); abstol, reltol,
                  isoutofdomain = isnegative); #reject negative solutions

isnonnegative(sol_Ros23)
```

With stricter tolerances, MPRK schemes also capture the behavior near the final time.

```@example HIRES
sol1 = solve(prob, MPRK43I(1.0, 0.5); abstol = 1e-4, reltol = 1e-3);
sol2 = solve(prob, MPRK43I(1.0, 0.5); abstol = 1e-6, reltol = 1e-5);

# plot solutions
hires_plot([sol1, sol2], ["MPRK43I(1.0,0.5), abstol=1e-4", "MPRK43I(1.0,0.5), abstol=1e-6"],
           ref_sol)
```

## Work-Precision diagrams

In the following we show several work-precision diagrams, which compare different methods with respect to computing time and the respective error.
We focus solely on adaptive methods.

Since the HIRES problem is stiff, we need to use a suited implicit scheme to compute a reference solution, see the [solver guide](https://docs.sciml.ai/DiffEqDocs/dev/solvers/ode_solve/#Stiff-Problems).
Note that we cannot use the recommended method `radau()`, since [`prob_pds_hires`](@ref) uses [StaticArrays.jl](https://juliaarrays.github.io/StaticArrays.jl/stable/) instead of standard `Array`s.

```@example HIRES
# select solver to compute reference solution
alg_ref = Rodas4P()
nothing # hide
```

We use the functions [`work_precision_adaptive`](@ref) and [`work_precision_adaptive!`](@ref) to compute the data for the diagrams.
Furthermore, the following absolute and relative tolerances are used.

```@example HIRES
# set absolute and relative tolerances
abstols = 1.0 ./ 10.0 .^ (2:1:9)
reltols = abstols .* 10.0
nothing # hide
```

The error chosen to compare the performances of different solvers is the relative maximum error at the final time ``t = 321.8122``.
We do not use the relative maximum error over all time steps ([`rel_max_error_overall`](@ref)), since most components of the initial value vanish.
The reference values of some components are therefore tiny after the first time step, and the relative error over all time steps is close to one for all tolerances considered here.

```@example HIRES
# select relative maximum error at the end of the problem's time span.
compute_error = rel_max_error_tend
nothing # hide
```

We start with a comparison of different adaptive MPRK schemes.

```@example HIRES
# choose methods to compare
algs = [MPRK22(0.5); MPRK22(2.0 / 3.0); MPRK22(1.0); SSPMPRK22(0.5, 1.0);
        MPRK43I(1.0, 0.5); MPRK43I(0.5, 0.75); MPRK43II(0.5); MPRK43II(2.0 / 3.0)]
labels = ["MPRK22(0.5)"; "MPRK22(2/3)"; "MPRK22(1.0)"; "SSPMPRK22(0.5,1.0)";
          "MPRK43I(1.0,0.5)"; "MPRK43I(0.5,0.75)"; "MPRK43II(0.5)"; "MPRK43II(2.0/3.0)"]

# compute work-precision data
wp = work_precision_adaptive(prob, algs, labels, abstols, reltols, alg_ref;
                             compute_error)

# plot work-precision diagram
plot(wp, labels; title = "HIRES benchmark", legend = :outerright,
     color = permutedims([repeat([1], 3)..., 2, repeat([3], 2)..., repeat([4], 2)...]),
     xlims = (10^-4, 10^2), xticks = 10.0 .^ (-4:1:2),
     ylims = (10^-5, 10^0), yticks = 10.0 .^ (-5:1:0), minorticks = 10)
```

The third-order schemes are more efficient than the second-order schemes, and `MPRK43I(1.0, 0.5)` performs best.
Next, we compare the `MPDeC` methods to `MPRK43I(1.0, 0.5)`.

```@example HIRES
# choose methods to compare
algs = [MPRK43I(1.0, 0.5);
        MPDeC(2); MPDeC(3); MPDeC(4); MPDeC(5); MPDeC(6); MPDeC(7); MPDeC(8); MPDeC(9); MPDeC(10)]
labels = ["MPRK43I(1.0,0.5)";
          "MPDeC(2)"; "MPDeC(3)"; "MPDeC(4)"; "MPDeC(5)"; "MPDeC(6)"; "MPDeC(7)"; "MPDeC(8)"; "MPDeC(9)"; "MPDeC(10)"]

# compute work-precision data
wp = work_precision_adaptive(prob, algs, labels, abstols, reltols, alg_ref;
                             compute_error)

# plot work-precision diagram
plot(wp, labels; title = "HIRES benchmark", legend = :outerright,
     color = permutedims([3, repeat([5], 9)...]),
     xlims = (10^-4, 10^2), xticks = 10.0 .^ (-4:1:2),
     ylims = (10^-5, 10^0), yticks = 10.0 .^ (-5:1:0), minorticks = 10)
```

The higher-order `MPDeC` methods yield smaller errors for loose tolerances, but their errors stop decreasing at about ``10^{-3}`` for the stricter tolerances.
Overall, the `MPDeC` methods do not provide a clear benefit compared to `MPRK43I(1.0, 0.5)`.

For comparisons with schemes from [OrdinaryDiffEq.jl](https://docs.sciml.ai/OrdinaryDiffEq/stable/), we choose the second-order scheme `MPRK22(1.0)`, the third-order scheme `MPRK43I(1.0, 0.5)`, and the 10th order scheme `MPDeC(10)`.
We compare these schemes with a selection of second- and third-order stiff solvers from [OrdinaryDiffEq.jl](https://docs.sciml.ai/OrdinaryDiffEq/stable/).
To guarantee nonnegative solutions, we use the solver option `isoutofdomain = isnegative`.

```@example HIRES
# select reference MPRK methods
algs1 = [MPRK22(1.0); MPRK43I(1.0, 0.5); MPDeC(10)]
labels1 = ["MPRK22(1.0)"; "MPRK43I(1.0,0.5)"; "MPDeC(10)"]

# select methods from OrdinaryDiffEq
algs2 = [TRBDF2(); Kvaerno3(); KenCarp3(); Rodas3(); ROS2(); ROS3(); Rosenbrock23()]
labels2 = ["TRBDF2"; "Kvaerno3"; "KenCarp3"; "Rodas3"; "ROS2"; "ROS3"; "Rosenbrock23"]

# compute work-precision data
wp = work_precision_adaptive(prob, algs1, labels1, abstols, reltols, alg_ref;
                             compute_error)
# add work-precision data with isoutofdomain = isnegative
work_precision_adaptive!(wp, prob, algs2, labels2, abstols, reltols, alg_ref;
                         compute_error, isoutofdomain = isnegative)

# plot work-precision diagram
plot(wp, [labels1; labels2]; title = "HIRES benchmark", legend = :outerright,
     color = permutedims([1, 3, 5, repeat([4], 3)..., repeat([6], 4)...]),
     xlims = (10^-8, 10^2), xticks = 10.0 .^ (-8:1:2),
     ylims = (10^-5, 10^0), yticks = 10.0 .^ (-5:1:0), minorticks = 10)
```

For this problem, the schemes from [OrdinaryDiffEq.jl](https://docs.sciml.ai/OrdinaryDiffEq/stable/) are clearly more efficient than the MPRK schemes for all tolerances considered here.
Among them, only `ROS2()` is noticeably slower. For loose tolerances, it rejects many steps due to `isoutofdomain = isnegative`, since it produces negative values otherwise.

Finally, we compare the selected MPRK schemes to some [recommended solvers](https://docs.sciml.ai/DiffEqDocs/dev/solvers/ode_solve/) of higher order from [OrdinaryDiffEq.jl](https://docs.sciml.ai/OrdinaryDiffEq/stable/).
Again, to guarantee nonnegative solutions we select the solver option `isoutofdomain = isnegative`.

```@example HIRES
algs3 = [Rodas5P(); Rodas4P(); RadauIIA5()]
labels3 = ["Rodas5P"; "Rodas4P"; "RadauIIA5"]

# compute work-precision data
wp = work_precision_adaptive(prob, algs1, labels1, abstols, reltols, alg_ref;
                             compute_error)
# add work-precision data with isoutofdomain = isnegative
work_precision_adaptive!(wp, prob, algs3, labels3, abstols, reltols, alg_ref;
                         compute_error, isoutofdomain = isnegative)

# plot work-precision diagram
plot(wp, [labels1; labels3]; title = "HIRES benchmark", legend = :outerright,
     color = permutedims([1, 3, 5, repeat([6], 2)..., 7]),
     xlims = (10^-8, 10^2), xticks = 10.0 .^ (-8:1:2),
     ylims = (10^-5, 10^0), yticks = 10.0 .^ (-5:1:0), minorticks = 10)
```

Again, the schemes from [OrdinaryDiffEq.jl](https://docs.sciml.ai/OrdinaryDiffEq/stable/) are more efficient than the MPRK schemes.

## Package versions

These results were obtained using the following versions.
```@example HIRES
using InteractiveUtils
versioninfo()
println()

using Pkg
Pkg.status(["PositiveIntegrators", "StaticArrays", "LinearSolve",
            "OrdinaryDiffEqFIRK", "OrdinaryDiffEqRosenbrock",
            "OrdinaryDiffEqSDIRK"],
           mode = PKGMODE_MANIFEST)
nothing # hide
```
