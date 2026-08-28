# [Benchmark: Solution of the Diffusion problem](@id benchmark-diffusion)

We consider the test problem [`prob_pds_diffusion`](@ref) of the spatially heterogeneous diffusion equation to assess the efficiency of different solvers from [OrdinaryDiffEq.jl](https://docs.sciml.ai/OrdinaryDiffEq/stable/) and [PositiveIntegrators.jl](https://github.com/NumericalMathematics/PositiveIntegrators.jl), especially on larger-scale problems.

```@example DIFFU
using OrdinaryDiffEqFIRK, OrdinaryDiffEqRosenbrock, OrdinaryDiffEqSDIRK
using PositiveIntegrators

# select spatially heterogeneous diffusion problem
prob_pds = prob_pds_diffusion
prob_ode = prob_ode_diffusion
nothing # hide
```

To keep the following code as clear as possible, we define a helper function `diffusion_plot` that we use for plotting.

```@example DIFFU
using Plots
using Interpolations

function plot_diffusion_xt(sol, title_str = "Diffusion Solution in (x,t)-Plane")
    # Extract original x and t grids
    x_orig = range(0.0, 1.0; length = length(sol.u[1]))
    t_orig = sol.t
    U_orig = hcat(sol.u...)'  # Dimensions: (length(t), length(x))

    # Create a continuous interpolation object over the data grid
    itp = linear_interpolation((t_orig, x_orig), U_orig, extrapolation_bc = Flat())

    # Generate a dense evaluation grid for smooth rendering (e.g., 300x300 points)
    t_fine = range(minimum(t_orig), maximum(t_orig); length = 300)
    x_fine = range(0.0, 1.0; length = 300)

    # Evaluate interpolated matrix on the fine grid
    U_fine = [itp(t, x) for t in t_fine, x in x_fine]

    # Create 2D heatmap using the fine interpolated grid
    p = heatmap(x_fine, t_fine, U_fine,
                c = :jet,             # Rainbow colormap matching the original paper
                clims = (0, 4),       # Color limits corresponding to original paper scale
                xlabel = "x",
                ylabel = "t",
                xticks = 0.0:0.2:1.0, # Set x-ticks in steps of 0.2
                colorbar = false,     # Hide the colorbar
                title = title_str,
                expand_limits = true)

    return p
end

nothing # hide
```

## Work-Precision diagrams

In the following we show several work-precision diagrams, which compare different methods with respect to computing time and the respective error.

Since the diffusion problem is stiff, we need to use a suited implicit scheme to compute a reference solution, see the [solver guide](https://docs.sciml.ai/DiffEqDocs/dev/solvers/ode_solve/#Stiff-Problems).

```@example DIFFU
# select solver to compute reference solution
alg_ref = Rodas5P()
nothing # hide
```

### Adaptive schemes

We use the functions [`work_precision_adaptive`](@ref) and [`work_precision_adaptive!`](@ref) to compute the data for the diagrams.
Furthermore, the following absolute and relative tolerances are used.

```@example DIFFU
# set absolute and relative tolerances
abstols = 1.0 ./ 10.0 .^ (2:1:7)
reltols = abstols .* 10.0
nothing # hide
```

#### Relative maximum error at the final time

In this section the chosen error is the relative maximum error at the final time ``t = 60.0``.

```@example DIFFU
# select relative maximum error at the end of the problem's time span.
compute_error = rel_max_error_tend
nothing # hide
```

We start with a comparison of different adaptive MPRK schemes.
```@example DIFFU
# choose methods to compare
algs = [MPRK22(0.5); MPRK22(2.0 / 3.0); MPRK22(1.0); SSPMPRK22(0.5, 1.0);
        MPRK43I(1.0, 0.5); MPRK43I(0.5, 0.75); MPRK43II(0.5); MPRK43II(2.0 / 3.0);
        MPDeC(2); MPDeC(3); MPDeC(4); MPDeC(5); MPDeC(6); MPDeC(7); MPDeC(8); MPDeC(9); MPDeC(10)]
labels = ["MPRK22(0.5)"; "MPPRK22(2/3)"; "MPRK22(1.0)"; "SSPMPRK22(0.5,1.0)";
          "MPRK43I(1.0, 0.5)"; "MPRK43I(0.5, 0.75)"; "MPRK43II(0.5)"; "MPRK43II(2.0/3.0)";
          "MPDeC(2)"; "MPDeC(3)"; "MPDeC(4)"; "MPDeC(5)"; "MPDeC(6)"; "MPDeC(7)"; "MPDeC(8)"; "MPDeC(9)"; "MPDeC(10)"]

# compute work-precision data
wp = work_precision_adaptive(prob_pds, algs, labels, abstols, reltols, alg_ref;
                            adaptive_ref = true, compute_error)

# plot work-precision diagram
plot(wp, labels; title = "Diffusion benchmark", legend = :outerright,
     color = permutedims([
    repeat([1], 3)..., 2, repeat([3], 2)..., repeat([4], 2)..., repeat([5], 9)...  ]),
     xlims = (1*10^-9, 10^-1), xticks = 10.0 .^ (-9:1:-1),
     ylims = (2*10^-3, 1*10^0), yticks = 10.0 .^ (-3:1:0), minorticks = 10)
```

For comparisons with other schemes we choose `MPRK22(1.0)`, `MPRK43I(1.0, 0.5)` and `MPDeC(10)`.

```@example DIFFU
# compute reference solution for plotting
ref_sol = solve(prob_ode, alg_ref; abstol = 1e-14, reltol = 1e-13);

# compute solutions
sol_MPRK22 = solve(prob_pds, MPRK22(1.0))
sol_MPRK43 = solve(prob_pds, MPRK43I(1.0, 0.5))
sol_MPDeC10 = solve(prob_pds, MPDeC(10))

p1 = plot_diffusion_xt(ref_sol, "Reference solution");
p2 = plot_diffusion_xt(sol_MPRK22, "MPRK22(1.0)");
p3 = plot_diffusion_xt(sol_MPRK43, "MPRK43I(0.5, 0.75)");
p4 = plot_diffusion_xt(sol_MPDeC10, "MPDeC(10)");
plot(p1, p2, p3, p4, layout=(2,2))
```

Next, we compare these three schemes with a selection of second- and third-order stiff solvers from [OrdinaryDiffEq.jl](https://docs.sciml.ai/OrdinaryDiffEq/stable/). 

```@example DIFFU
# select reference MPRK methods
algs1 = [MPRK22(1.0); MPRK43I(1.0, 0.5); MPDeC(10)]
labels1 = ["MPRK22(1.0)"; "MPRK43I(1.0,0.5)"; "MPDeC(10)"]

# select methods from OrdinaryDiffEq
algs2 = [TRBDF2(); Kvaerno3(); KenCarp3(); Rodas3(); ROS2(); ROS3(); Rosenbrock23()]
labels2 = ["TRBDF2"; "Kvearno3"; "KenCarp3"; "Rodas3"; "ROS2"; "ROS3"; "Rosenbrock23"]

# compute work-precision data
wp = work_precision_adaptive(prob_pds, algs1, labels1, abstols, reltols, alg_ref;
                               adaptive_ref = true, compute_error)
# add work-precision data
work_precision_adaptive!(wp, prob_ode, algs2, labels2, abstols, reltols, alg_ref;
                               adaptive_ref = true, compute_error)

# plot work-precision diagram
plot(wp, [labels1; labels2]; title = "Diffusion benchmark", legend = :outerright,
     color = permutedims([1, 3, 5, repeat([6], 3)..., repeat([7], 4)...]),
     xlims = (10^-9, 10^-1), xticks = 10.0 .^ (-11:1:-1),
     ylims = (10^-4, 10^1), yticks = 10.0 .^ (-4:1:1), minorticks = 10)
```

In addition,  we compare the selected MPRK schemes to some [recommended solvers](https://docs.sciml.ai/DiffEqDocs/dev/solvers/ode_solve/) of higher order from [OrdinaryDiffEq.jl](https://docs.sciml.ai/OrdinaryDiffEq/stable/). 

```@example DIFFU
algs3 = [Rodas5P(); Rodas4P(); RadauIIA5()]
labels3 = ["Rodas5P"; "Rodas4P"; "RadauIIA5"]

# compute work-precision data
wp = work_precision_adaptive(prob_pds, algs1, labels1, abstols, reltols, alg_ref;
                               adaptive_ref = true, compute_error)
# add work-precision data with isoutofdomain = isnegative
work_precision_adaptive!(wp, prob_ode, algs3, labels3, abstols, reltols, alg_ref;
                               adaptive_ref = true, compute_error)

# plot work-precision diagram
plot(wp, [labels1; labels3]; title = "Diffusion benchmark", legend = :topright,
     color = permutedims([1, 3, 5, repeat([4], 2)..., 6]),
     xlims = (10^-9, 10^-2), xticks = 10.0 .^ (-9:1:-2),
     ylims = (10^-4, 10^0), yticks = 10.0 .^ (-4:1:0), minorticks = 10)
```

#### Relative maximum error over all time steps

In this section we do not compare the relative maximum errors at the final time ``t = 60.0}``, but the relative maximum errors over all time steps.

```@example DIFFU
# select relative maximum error at the end of the problem's time span.
compute_error = rel_max_error_overall
nothing # hide
```

First, we compare different MPRK schemes.

```@example DIFFU
# compute work-precision data
wp = work_precision_adaptive(prob_pds, algs, labels, abstols, reltols, alg_ref;
                            adaptive_ref = true, compute_error)

# plot work-precision diagram
plot(wp, labels; title = "Diffusion benchmark", legend = :outerright,
     color = permutedims([repeat([1], 3)..., 2, repeat([3], 2)..., repeat([4], 2)..., repeat([5], 9)...  ]),
     xlims = (10^-7, 10^0), xticks = 10.0 .^ (-7:1:2),
     ylims = (3*10^-3, 10^0), yticks = 10.0 .^ (-3:1:0), minorticks = 10)
```

We choose the schemes `MPRK22(1.0)`, `SSPMPKR22(0.5,1.0)` and `MPDeC(7)` for comparison with solvers from [OrdinaryDiffEq.jl](https://docs.sciml.ai/OrdinaryDiffEq/stable/).

```@example DIFFU
# select reference MPRK methods
algs1 = [MPRK22(1.0); SSPMPRK22(0.5, 1.0); MPDeC(2)]
labels1 = ["MPRK22(1.0)"; "SSPMPRK22(0.5, 1.0)"; "MPDeC(2)"]

# compute work-precision data
wp = work_precision_adaptive(prob_pds, algs1, labels1, abstols, reltols, alg_ref;
                               adaptive_ref = true, compute_error)
# add work-precision data with isoutofdomain = isnegative
work_precision_adaptive!(wp, prob_ode, algs2, labels2, abstols, reltols, alg_ref; adaptive_ref = true, compute_error)

# plot work-precision diagram
plot(wp, [labels1; labels2]; title = "Diffusion benchmark", legend = :topright,
     color = permutedims([1, 2, 3, repeat([5], 3)..., repeat([6], 4)...]),
     xlims = (10^-7, 10^0), xticks = 10.0 .^ (-7:1:3),
     ylims = (2*10^-4, 10^0), yticks = 10.0 .^ (-4:1:0), minorticks = 10)
```

Finally, we compare `MPRK43I(0.5, 0.75)` and `MPRK22(1.0)` to [recommended solvers](https://docs.sciml.ai/DiffEqDocs/dev/solvers/ode_solve/) of higher order from [OrdinaryDiffEq.jl](https://docs.sciml.ai/OrdinaryDiffEq/stable/).

```@example DIFFU
# compute work-precision data
wp = work_precision_adaptive(prob_pds, algs1, labels1, abstols, reltols, alg_ref;
                               adaptive_ref = true, compute_error)
# add work-precision data with isoutofdomain = isnegative
work_precision_adaptive!(wp, prob_ode, algs3, labels3, abstols, reltols, alg_ref;adaptive_ref = true, compute_error)

# plot work-precision diagram
plot(wp, [labels1; labels3]; title = "Diffusion benchmark", legend = :topright,
     color = permutedims([1, 3, repeat([4], 2)..., 5]),
     xlims = (10^-10, 10^0), xticks = 10.0 .^ (-10:1:0),
     ylims = (10^-3, 2*10^1), yticks = 10.0 .^ (-3:1:1), minorticks = 10)
```

### Fixed time step sizes

Here we use fixed time step sizes instead of adaptive time stepping.
We use the functions [`work_precision_fixed`](@ref) and [`work_precision_fixed!`](@ref) to compute the data for the diagrams.
Please note that these functions set error and computing time to `Inf`, whenever a solution contains negative elements.
Consequently, such cases are not visible in the work-precision diagrams.

Within the work-precision diagrams we use the following time step sizes.

```@example DIFFU
# set time step sizes
dts = 6.0 ./ 2.0 .^ (0:1:5)
nothing # hide
```

#### Relative maximum error at the end of the problem's time span

Again, we start with the relative maximum error at the final time ``t = 60.0``.

```@example DIFFU
# select relative maximum error at the end of the problem's time span.
compute_error = rel_max_error_tend
nothing  # hide
```

First, we compare different MPRK methods.
For fixed time step sizes we can also consider `MPE()`, `SSPMPRK43()` and the various MPLM schemes.

```@example DIFFU
# choose MPRK methods to compare
algs = [MPE(); algs; SSPMPRK43(); 
        MPLM22(); MPLM33(); MPLM43();
        MPLM54(); MPLM75(); MPLM106()]
labels = ["MPE"; labels; "SSPMPRK43";
          "MPLM22"; "MPLM33"; "MPLM43";
          "MPLM54"; "MPLM75"; "MPLM106"]

# compute work-precision data
wp = work_precision_fixed(prob_pds, algs, labels, dts, alg_ref;
                          compute_error)

# plot work-precision diagram
plot(wp, labels; title = "Diffusion benchmark", legend = :outerright,
     color = permutedims([5,repeat([1], 3)..., 2, repeat([3], 2)..., repeat([4], 2)..., repeat([5], 9)..., 2, repeat([13], 6)...]),
     xlims = (10^-10, 10^2), xticks = 10.0 .^ (-10:1:2),
     ylims = (10^-3, 10^1), yticks = 10.0 .^ (-6:1:1), minorticks = 10)
```

Besides `MPLM22()` all other MPLM schemes show instablities.
As an example, we plot the numerical solutions computed with `MPLM22()` and `MPLM43()`.

```@example DIFFU
# compute solutions
sol_MPLM22 = solve(prob_pds, MPLM22(); dt = 1.0)
sol_MPLM43 = solve(prob_pds, MPLM43(); dt = 1.0)

p1 = plot_diffusion_xt(ref_sol, "Reference solution");
p2 = plot_diffusion_xt(sol_MPLM22, "MPLM22");
p3 = plot_diffusion_xt(sol_MPLM43, "MPLM43");
plot(p1, p2, p3, layout=(1,3))
```

We choose `MPRK22(1.0)`, `MPRK43I(1.0, 0.5)` and `MPDeC(4)` for comparisons with other schemes from [OrdinaryDiffEq.jl](https://docs.sciml.ai/OrdinaryDiffEq/stable/).
First, we compare these methods with other second- and third-order schemes.

```@example DIFFU
algs1 = [MPRK22(1.0); MPRK43I(1.0, 0.5); MPDeC(4)]
labels1 = ["MPRK22(1.0)"; "MPRK43I(1.0,0.5)"; "MPDeC(4)"]

# compute work-precision data
wp = work_precision_fixed(prob_pds, algs1, labels1, dts, alg_ref;
                               compute_error)
work_precision_fixed!(wp, prob_ode, algs2, labels2, dts, alg_ref;
                     compute_error)
# plot work-precision diagram
plot(wp, [labels1; labels2]; title = "Diffusion benchmark", legend = :topright,
     color = permutedims([1, 3, 5, repeat([4], 3)..., repeat([13],4)...,repeat([6],4)...]),
     xlims = (10^-9, 10^-1), xticks = 10.0 .^ (-9:2:-1),
     ylims = (10^-4, 10^0), yticks = 10.0 .^ (-5:1:0), minorticks = 10)
```
Finally, we show a comparison with [recommended solvers](https://docs.sciml.ai/DiffEqDocs/dev/solvers/ode_solve/) from [OrdinaryDiffEq.jl](https://docs.sciml.ai/OrdinaryDiffEq/stable/).

```@example DIFFU
# compute work-precision data
wp = work_precision_fixed(prob_pds, algs1, labels1, dts, alg_ref; compute_error)
work_precision_fixed!(wp, prob_ode, algs3, labels3, dts, alg_ref; compute_error)

# plot work-precision diagram
plot(wp, [labels1; labels3]; title = "Diffusion benchmark", legend = :bottomright,
     color = permutedims([1, 3, 5, repeat([4], 3)..., repeat([13],4)...,repeat([6],4)...]),
     xlims = (10^-12, 10^-2), xticks = 10.0 .^ (-12:2:-2),
     ylims = (10^-4, 10^0), yticks = 10.0 .^ (-4:1:0), minorticks = 10)
```

#### Relative maximum error over all time steps

As for the adaptive schemes, we also show work-precisions diagrams where the error is the relative maximum error over all time steps.

```@example DIFFU
# select relative maximum error over all time steps
compute_error = rel_max_error_overall
nothing  # hide
```

```@example DIFFU
# compute work-precision
wp = work_precision_fixed(prob_pds, algs, labels, dts, alg_ref;
                               compute_error)

#plot work-precision diagram
plot(wp, labels; title = "Diffusion benchmark", legend = :outerright,
     color = permutedims([5,repeat([1], 3)..., 2, repeat([3], 2)..., repeat([4], 2)..., repeat([5], 9)..., 2, repeat([13], 6)...]),
     xlims = (10^-3, 10^2), xticks = 10.0 .^ (-3:1:2),
     ylims = (10^-3, 10^1), yticks = 10.0 .^ (-3:1:1), minorticks = 10)
```

We choose `MPRK22(1.0)`, `MPDeC(2)` and `MPLM22()` for comparisons with other schemes from [OrdinaryDiffEq.jl](https://docs.sciml.ai/OrdinaryDiffEq/stable/).

```@example DIFFU
algs1 = [MPRK22(1.0); MPDeC(2); MPLM22()]
labels1 = ["MPRK22(1.0)"; "MPDeC(4)"; "MPLM22"]

wp = work_precision_fixed(prob_pds, algs1, labels1, dts, alg_ref;
                               compute_error)
work_precision_fixed!(wp, prob_ode, algs2, labels2, dts, alg_ref;
                     compute_error)

plot(wp, [labels1; labels2]; title = "Diffusion benchmark", legend = :topright,
     color = permutedims([1, 3, 13, repeat([4], 3)..., repeat([5], 4)..., repeat([6], 4)...]),
     xlims = (10^-3, 10^0), xticks = 10.0 .^ (-12:1:6),
     ylims = (10^-4, 10^0), yticks = 10.0 .^ (-5:1:0), minorticks = 10)
```

```@example DIFFU
wp = work_precision_fixed(prob_pds, algs1, labels1, dts, alg_ref;
                               compute_error)
work_precision_fixed!(wp, prob_ode, algs3, labels3, dts, alg_ref;
                     compute_error)

plot(wp, [labels1; labels3]; title = "Diffusion benchmark", legend = :topright,
     color = permutedims([1, 3, 13, repeat([4], 5)..., 5, repeat([7], 3)...]),
     xlims = (10^-3, 10^-1), xticks = 10.0 .^ (-3:1:-1),
     ylims = (10^-4, 10^1), yticks = 10.0 .^ (-4:1:1), minorticks = 10)
```

## Package versions

These results were obtained using the following versions.
```@example DIFFU
using InteractiveUtils
versioninfo()
println()

using Pkg
Pkg.status(["PositiveIntegrators", "StaticArrays", "LinearSolve",
            "OrdinaryDiffEqFIRK", "OrdinaryDiffEqRosenbrock",
            "OrdinaryDiffEqSDIRK"],
           mode=PKGMODE_MANIFEST)
nothing # hide
```
