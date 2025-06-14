using Test
using LinearAlgebra
using SparseArrays
using Statistics: mean, median

using DoubleFloats: Double64
using StaticArrays: @SMatrix, SMatrix, MVector, @SVector, SVector, SA

using Unitful: @u_str, ustrip

using ADTypes
using OrdinaryDiffEqLowOrderRK: Euler
using OrdinaryDiffEqRosenbrock: Rosenbrock23, Rodas4P
using OrdinaryDiffEqSDIRK: ImplicitEuler, SDIRK2, TRBDF2
using OrdinaryDiffEqTsit5: Tsit5
using OrdinaryDiffEqVerner: Vern7, Vern9
using PositiveIntegrators

# load RecursiveFactorization to get RFLUFactorization
using RecursiveFactorization: RecursiveFactorization
using LinearSolve: RFLUFactorization, LUFactorization, KLUFactorization, KrylovJL_GMRES

using Aqua: Aqua
using RecipesBase: RecipesBase # only for Aqua tests
using ExplicitImports: check_no_implicit_imports, check_no_stale_explicit_imports

"""
    experimental_orders_of_convergence(prob, alg, dts;
                                      test_time = nothing,
                                      only_first_index = false,
                                      ref_alg = TRBDF2(autodiff = AutoFiniteDiff()))

Solve `prob` with `alg` and fixed time steps taken from `dts`, and compute
the errors at `test_time`. If`test_time` is not specified the error is computed
at the final time.
Return the associated experimental orders of convergence.

If `only_first_index == true`, only the first solution component is used
to compute the error. If no analytic solution is available, a reference
solution is computed using `ref_alg`.
"""
function experimental_orders_of_convergence(prob, alg, dts; test_time = nothing,
                                            only_first_index = false,
                                            ref_alg = Vern7())
    @assert length(dts) > 1
    errors = zeros(eltype(dts), length(dts))

    if !(isnothing(prob.f.analytic))
        # there is an analytic solution
        if isnothing(test_time)
            # we compare the results at the final time
            reference_solution = prob.f.analytic(prob.u0, prob.p, last(prob.tspan))
        else
            # we compare the results at the given time
            reference_solution = prob.f.analytic(prob.u0, prob.p, test_time)
        end
    else
        # we compute a reference solution numerically
        if isnothing(test_time)
            # we compare the results at the final time
            tspan = prob.tspan
        else
            # we compare the results at the given time
            tspan = (first(prob.tspan), test_time)
        end
        dt0 = (tspan[end] - tspan[begin]) / 1e5
        refsol = solve(prob, ref_alg;
                       dt = dt0, adaptive = false,
                       save_everystep = false)
        reference_solution = refsol.u[end]
    end

    for (i, dt) in enumerate(dts)
        if isnothing(test_time)
            sol = solve(prob, alg; dt = dt, adaptive = false, save_everystep = false)
            numerical_solution = sol.u[end]
        else
            sol = solve(prob, alg; dt = dt, adaptive = false)
            numerical_solution = sol(test_time)
        end

        if only_first_index
            errors[i] = norm(first(numerical_solution) - first(reference_solution))
        else
            errors[i] = norm(numerical_solution - reference_solution)
        end
    end

    return experimental_orders_of_convergence(errors, dts)
end

"""
    experimental_orders_of_convergence(errors, dts)

Compute the experimental orders of convergence for given `errors` and
time step sizes `dts`.
"""
function experimental_orders_of_convergence(errors, dts)
    Base.require_one_based_indexing(errors, dts)
    @assert length(errors) == length(dts)
    orders = zeros(eltype(errors), length(errors) - 1)

    for i in eachindex(orders)
        orders[i] = log(errors[i] / errors[i + 1]) / log(dts[i] / dts[i + 1])
    end

    return orders
end

"""
    check_order(orders, alg_order; N = 3, atol = 0.1)

Check if `alg_order` can be found approximately in `N` consecutive elements of `orders`.
"""
function check_order(orders, alg_order; N = 3, atol = 0.1)
    check = false
    # The calculation of the indices of the permissible orders by
    #
    #    indices = findall(x -> isapprox(x, alg_order; atol = atol), orders)
    #
    # sometimes failed because the experimental order was too good. To accept
    # such cases, but avoid to decrease atol, we now use the following asymmetric criterion:
    indices = findall(x -> -atol ≤ x - alg_order ≤ 2 * atol, orders)

    for (i, idx) in enumerate(indices)
        if i + N - 1 ≤ length(indices) && indices[i:(i + N - 1)] == idx:(idx + N - 1)
            check = true
            break
        end
    end
    return check
end

const prob_pds_linmod_array = ConservativePDSProblem(prob_pds_linmod.f,
                                                     Array(prob_pds_linmod.u0),
                                                     prob_pds_linmod.tspan)
const prob_pds_linmod_mvector = ConservativePDSProblem(prob_pds_linmod_inplace.f,
                                                       MVector(prob_pds_linmod.u0),
                                                       prob_pds_linmod.tspan)

# analytic solution of linear model problem
function f_analytic(u0, p, t)
    u₁⁰, u₂⁰ = u0
    a, b = p
    c = a + b
    return ((u₁⁰ + u₂⁰) * [b; a] +
            exp(-c * t) * (a * u₁⁰ - b * u₂⁰) * [1; -1]) / c
end

# This is the usual conservative linear model problem, rewritten as
# u₁' = -3 u₁ + 0.5 u₂ - 2 u₁ + 0.5 u₂ (= -5 u₁ + u₂)
# u₂' =  3 u₁ - 0.5 u₂ - 0.5 u₂ + 2 u₁ (= 5 u₁ - u₂)
# linear model problem - nonconservative - out-of-place
linmodP(u, p, t) = [0.5*u[2] 0.5*u[2]; 3*u[1] 2*u[1]]
linmodD(u, p, t) = [2 * u[1]; 0.5 * u[2]]
const prob_pds_linmod_nonconservative = PDSProblem(linmodP, linmodD, [0.9, 0.1], (0.0, 2.0),
                                                   [5.0, 1.0];
                                                   analytic = f_analytic)

# linear model problem - nonconservative -  in-place
function linmodP!(P, u, p, t)
    P[1, 1] = 0.5 * u[2]
    P[1, 2] = 0.5 * u[2]
    P[2, 1] = 3 * u[1]
    P[2, 2] = 2 * u[1]
    return nothing
end
function linmodD!(D, u, p, t)
    D[1] = 2 * u[1]
    D[2] = 0.5 * u[2]
    return nothing
end
const prob_pds_linmod_nonconservative_inplace = PDSProblem(linmodP!, linmodD!, [0.9, 0.1],
                                                           (0.0, 2.0), [5.0, 1.0];
                                                           analytic = f_analytic)

# Function definitions for testset "Linear advection"
# Functions are defined outside the testset to avoid failing allocation tests,
# see https://github.com/NumericalMathematics/PositiveIntegrators.jl/pull/89.
#
# in-place syntax for f
function linear_advection_fd_upwind_f!(du, u, p, t)
    N = length(u)
    dx = 1 / N
    du[1] = -(u[1] - u[N]) / dx
    for i in 2:N
        du[i] = -(u[i] - u[i - 1]) / dx
    end
end
# in-place sytanx for PDS
function linear_advection_fd_upwind_P!(P, u, p, t)
    P .= 0.0
    N = length(u)
    dx = 1 / N
    P[1, N] = u[N] / dx
    for i in 2:N
        P[i, i - 1] = u[i - 1] / dx
    end
    return nothing
end
function linear_advection_fd_upwind_P!(P::SparseMatrixCSC, u, p, t)
    N = length(u)
    dx = 1 / N
    values = nonzeros(P)
    for col in axes(P, 2)
        for idx in nzrange(P, col)
            values[idx] = u[col] / dx
        end
    end
    return nothing
end
function linear_advection_fd_upwind_D!(D, u, p, t)
    D .= 0.0
    return nothing
end

@testset "PositiveIntegrators.jl tests" begin
    @testset "Aqua.jl" begin
        # The Aqua.jl tests fails in the Downstream CI action
        # of OrdinaryDiffEq.jl but not in our regular CI - we just
        # skip it there.
        if !isempty(get(ENV, "GROUP", ""))
            @info "Skipping tests from Aqua.jl"
        else
            @info "Running tests from Aqua.jl"
            # We do not test ambiguities since we get a lot of
            # false positives from dependencies.
            ambiguities = false
            # The persistent_tasks test fails in the Downgrade CI
            # action but not in regular CI - we just skip it there.
            if !isempty(get(ENV, "POSITIVEINTEGRATORS_DOWNGRADE_CI", ""))
                persistent_tasks = false
            else
                persistent_tasks = true
            end
            Aqua.test_all(PositiveIntegrators;
                          ambiguities = ambiguities,
                          piracies = (; treat_as_own = [RecipesBase.apply_recipe],),
                          persistent_tasks = persistent_tasks)
        end
    end

    @testset "ExplicitImports.jl" begin
        @test isnothing(check_no_implicit_imports(PositiveIntegrators))
        @test isnothing(check_no_stale_explicit_imports(PositiveIntegrators))
    end

    @testset "ODE RHS" begin
        let counter_p = Ref(1), counter_d = Ref(1), counter_rhs = Ref(1)
            # out-of-place
            prod1 = (u, p, t) -> begin
                counter_p[] += 1
                return [0 u[2]; u[1] 0]
            end
            dest1 = (u, p, t) -> begin
                counter_d[] += 1
                return zero(u)
            end
            rhs1 = (u, p, t) -> begin
                counter_rhs[] += 1
                return [-u[1] + u[2], u[1] - u[2]]
            end
            u0 = [1.0, 0.0]
            tspan = (0.0, 1.0)
            prob_default = PDSProblem(prod1, dest1, u0, tspan)
            prob_special = PDSProblem(prod1, dest1, u0, tspan; std_rhs = rhs1)

            counter_p[] = 0
            counter_d[] = 0
            counter_rhs[] = 0
            @inferred prob_default.f(u0, nothing, 0.0)
            @test counter_p[] == 1
            @test counter_d[] == 1
            @test counter_rhs[] == 0

            counter_p[] = 0
            counter_d[] = 0
            counter_rhs[] = 0
            @inferred prob_special.f(u0, nothing, 0.0)
            @test counter_p[] == 0
            @test counter_d[] == 0
            @test counter_rhs[] == 1

            # in-place
            prod1! = (P, u, p, t) -> begin
                counter_p[] += 1
                P[1, 1] = 0
                P[1, 2] = u[2]
                P[2, 1] = u[1]
                P[2, 2] = 0
                return nothing
            end
            dest1! = (D, u, p, t) -> begin
                counter_d[] += 1
                fill!(D, 0)
                return nothing
            end
            rhs1! = (du, u, p, t) -> begin
                counter_rhs[] += 1
                du[1] = -u[1] + u[2]
                du[2] = u[1] - u[2]
                return nothing
            end
            u0 = [1.0, 0.0]
            tspan = (0.0, 1.0)
            prob_default = PDSProblem(prod1!, dest1!, u0, tspan)
            prob_special = PDSProblem(prod1!, dest1!, u0, tspan; std_rhs = rhs1!)

            du = similar(u0)
            counter_p[] = 0
            counter_d[] = 0
            counter_rhs[] = 0
            @inferred prob_default.f(du, u0, nothing, 0.0)
            @test counter_p[] == 1
            @test counter_d[] == 1
            @test counter_rhs[] == 0

            counter_p[] = 0
            counter_d[] = 0
            counter_rhs[] = 0
            @inferred prob_special.f(du, u0, nothing, 0.0)
            @test counter_p[] == 0
            @test counter_d[] == 0
            @test counter_rhs[] == 1

            counter_p[] = 0
            counter_d[] = 0
            counter_rhs[] = 0
            @inferred solve(prob_default, MPE(); dt = 0.1)
            @test 10 <= counter_p[] <= 11
            @test 10 <= counter_d[] <= 11
            @test counter_d[] == counter_p[]
            @test counter_rhs[] == 0

            counter_p[] = 0
            counter_d[] = 0
            counter_rhs[] = 0
            @inferred solve(prob_default, Euler(); dt = 0.1)
            @test 10 <= counter_p[] <= 11
            @test 10 <= counter_d[] <= 11
            @test counter_d[] == counter_p[]
            @test counter_rhs[] == 0

            counter_p[] = 0
            counter_d[] = 0
            counter_rhs[] = 0
            @inferred solve(prob_special, Euler(); dt = 0.1)
            @test counter_p[] == 0
            @test counter_d[] == 0
            @test 10 <= counter_rhs[] <= 11
        end
    end

    @testset "ConservativePDSFunction" begin
        prod_1! = (P, u, p, t) -> begin
            fill!(P, zero(eltype(P)))
            for i in 1:(length(u) - 1)
                P[i, i + 1] = i * u[i + 1]
            end
            return nothing
        end
        prod_2! = (P, u, p, t) -> begin
            fill!(P, zero(eltype(P)))
            for i in 1:(length(u) - 1)
                P[i + 1, i] = i * u[i]
            end
            return nothing
        end
        prod_3! = (P, u, p, t) -> begin
            fill!(P, zero(eltype(P)))
            for i in 1:(length(u) - 1)
                P[i, i + 1] = i * u[i + 1]
                P[i + 1, i] = i * u[i]
            end
            return nothing
        end

        n = 10
        P_tridiagonal = Tridiagonal(rand(n - 1), zeros(n), rand(n - 1))
        P_dense = Matrix(P_tridiagonal)
        P_sparse = sparse(P_tridiagonal)
        u0 = rand(n)
        tspan = (0.0, 1.0)

        du_tridiagonal = similar(u0)
        du_dense = similar(u0)
        du_sparse = similar(u0)

        for prod! in (prod_1!, prod_2!, prod_3!)
            prob_tridiagonal = ConservativePDSProblem(prod!, u0, tspan;
                                                      p_prototype = P_tridiagonal)
            prob_dense = ConservativePDSProblem(prod!, u0, tspan;
                                                p_prototype = P_dense)
            prob_sparse = ConservativePDSProblem(prod!, u0, tspan;
                                                 p_prototype = P_sparse)

            prob_tridiagonal.f(du_tridiagonal, u0, nothing, 0.0)
            prob_dense.f(du_dense, u0, nothing, 0.0)
            prob_sparse.f(du_sparse, u0, nothing, 0.0)

            @test du_tridiagonal ≈ du_dense
            @test du_tridiagonal ≈ du_sparse
            @test du_dense ≈ du_sparse
        end
    end

    @testset "PDSProblem" begin
        @testset "Linear model problem" begin
            # This is an example of a conservative PDS
            #
            # ODE system: u₁' = -5 u₁ + 1 u₂, u₂' = 5 u₁ - 1 u₂
            #
            # Standard: f(t,u) = [-5 1; 5 -1]*u
            #
            # PDS: P = [0 u₂; 5 u₁ 0], D = [0; 0]
            #
            # We implement four variations (in-place and out-of-place
            # for both standard and PDS) of this ODE system and check
            # that we obtain equivalent solutions with a standard
            # solver of OrdinaryDiffEq.jl.

            # initial values
            u0 = [0.9, 0.1]
            # time domain
            tspan = (0.0, 2.0)

            # out-of-place syntax for standard ODE
            A = [-5.0 1.0; 5.0 -1.0]
            linmod(u, p, t) = A * u
            linmod_ODE_op = ODEProblem(linmod, u0, tspan)

            # in-place syntax for standard ODE
            function linmod!(du, u, p, t)
                u1, u2 = u
                du[1] = -5 * u1 + u2
                du[2] = 5 * u1 - u2
                return nothing
            end
            linmod_ODE_ip = ODEProblem(linmod!, u0, tspan)

            # out-of-place syntax for PDS
            linmodP(u, p, t) = [0 u[2]; 5*u[1] 0]
            linmodD(u, p, t) = [0.0; 0.0]
            linmod_PDS_op = PDSProblem(linmodP, linmodD, u0, tspan)

            # in-place sytanx for PDS
            function linmodP!(P, u, p, t)
                fill!(P, zero(eltype(P)))
                P[1, 2] = u[2]
                P[2, 1] = 5 * u[1]
                return nothing
            end
            function linmodD!(D, u, p, t)
                fill!(D, zero(eltype(D)))
                return nothing
            end
            linmod_PDS_ip = PDSProblem(linmodP!, linmodD!, u0, tspan)

            # solutions
            sol_linmod_ODE_op = solve(linmod_ODE_op, Tsit5())
            sol_linmod_ODE_ip = solve(linmod_ODE_ip, Tsit5())
            sol_linmod_PDS_op = solve(linmod_PDS_op, Tsit5())
            sol_linmod_PDS_ip = solve(linmod_PDS_ip, Tsit5())

            # check equality of solutions
            @test sol_linmod_ODE_op.t ≈ sol_linmod_ODE_ip.t
            @test sol_linmod_ODE_op.t ≈ sol_linmod_PDS_op.t
            @test sol_linmod_ODE_op.t ≈ sol_linmod_PDS_ip.t
            @test sol_linmod_ODE_op.u ≈ sol_linmod_ODE_ip.u
            @test sol_linmod_ODE_op.u ≈ sol_linmod_PDS_op.u
            @test sol_linmod_ODE_op.u ≈ sol_linmod_PDS_ip.u

            # check that we really do not use too many additional
            # allocations for in-place implementations
            alloc1 = @allocated(solve(linmod_ODE_ip, Tsit5()))
            alloc2 = @allocated(solve(linmod_PDS_ip, Tsit5()))
            @test 0.9 < alloc1 / alloc2 < 1.1
        end

        @testset "PDSProblem error handling" begin
            P(u, p, t) = 0.0
            D(du, u, p, t) = 0.0
            if VERSION >= v"1.8"
                @test_throws "in-place and out-of-place" PDSProblem(P, D, 0.0, (0.0, 1.0))
            end
        end
    end

    # Here we check that solutions of equivalent ODEProblems, PDSProblems or
    # ConservativePDS Problems are approximately equal.
    # We also check that solvers from OrdinaryDiffEq can solve PDSProblems and
    # ConservativePDSProblems.
    @testset "Check compatibility of PositiveIntegrators and OrdinaryDiffEq" begin
        @testset "Linear model" begin
            # Linear model (conservative)
            u0 = [0.9, 0.1]
            tspan = (0.0, 2.0)
            A = [-5.0 1.0; 5.0 -1.0]
            linmod(u, p, t) = A * u
            linmod_f_op = ODEProblem(linmod, u0, tspan)
            # in-place syntax for f
            function linmod!(du, u, p, t)
                u₁, u₂ = u
                du[1] = -5.0 * u₁ + u₂
                du[2] = 5.0 * u₁ - u₂
            end
            linmod_f_ip = ODEProblem(linmod!, u0, tspan)
            # out-of-place syntax for PDS
            linmodP(u, p, t) = [0.0 u[2]; 5.0*u[1] 0.0]
            linmodD(u, p, t) = [0.0; 0.0]
            linmod_PDS_op = PDSProblem(linmodP, linmodD, u0, tspan)
            linmod_PDS_op_2 = PDSProblem{false}(linmodP, linmodD, u0, tspan)
            linmod_ConsPDS_op = ConservativePDSProblem(linmodP, u0, tspan)
            linmod_ConsPDS_op_2 = ConservativePDSProblem{false}(linmodP, u0, tspan)
            # in-place sytanx for PDS
            function linmodP!(P, u, p, t)
                P .= 0.0
                P[1, 2] = u[2]
                P[2, 1] = 5.0 * u[1]
                return nothing
            end
            function linmodD!(D, u, p, t)
                D .= 0.0
                return nothing
            end
            linmod_PDS_ip = PDSProblem(linmodP!, linmodD!, u0, tspan)
            linmod_PDS_ip_2 = PDSProblem{true}(linmodP!, linmodD!, u0, tspan)
            linmod_ConsPDS_ip = ConservativePDSProblem(linmodP!, u0, tspan)
            linmod_ConsPDS_ip_2 = ConservativePDSProblem{true}(linmodP!, u0, tspan)

            # solutions
            sol_linmod_f_op = solve(linmod_f_op, Tsit5())
            sol_linmod_f_ip = solve(linmod_f_ip, Tsit5())
            sol_linmod_PDS_op = solve(linmod_PDS_op, Tsit5())
            sol_linmod_PDS_op_2 = solve(linmod_PDS_op_2, Tsit5())
            sol_linmod_PDS_ip = solve(linmod_PDS_ip, Tsit5())
            sol_linmod_PDS_ip_2 = solve(linmod_PDS_ip_2, Tsit5())
            sol_linmod_ConsPDS_op = solve(linmod_ConsPDS_op, Tsit5())
            sol_linmod_ConsPDS_op_2 = solve(linmod_ConsPDS_op_2, Tsit5())
            sol_linmod_ConsPDS_ip = solve(linmod_ConsPDS_ip, Tsit5())
            sol_linmod_ConsPDS_ip_2 = solve(linmod_ConsPDS_ip_2, Tsit5())

            # check equality of solutions
            @test sol_linmod_f_op.t ≈ sol_linmod_f_ip.t ≈
                  sol_linmod_PDS_op.t ≈ sol_linmod_PDS_ip.t ≈
                  sol_linmod_PDS_op_2.t ≈ sol_linmod_PDS_ip_2.t ≈
                  sol_linmod_ConsPDS_op.t ≈ sol_linmod_ConsPDS_ip.t ≈
                  sol_linmod_ConsPDS_op_2.t ≈ sol_linmod_ConsPDS_ip_2.t
            @test sol_linmod_f_op.u ≈ sol_linmod_f_ip.u ≈
                  sol_linmod_PDS_op.u ≈ sol_linmod_PDS_ip.u ≈
                  sol_linmod_PDS_op_2.u ≈ sol_linmod_PDS_ip_2.u ≈
                  sol_linmod_ConsPDS_op.u ≈ sol_linmod_ConsPDS_ip.u ≈
                  sol_linmod_ConsPDS_op_2.u ≈ sol_linmod_ConsPDS_ip_2.u

            # check that we really do not use too many additional allocations for in-place implementations
            alloc1 = @allocated(solve(linmod_f_ip, Tsit5()))
            alloc2 = @allocated(solve(linmod_PDS_ip, Tsit5()))
            alloc3 = @allocated(solve(linmod_PDS_ip_2, Tsit5()))
            alloc4 = @allocated(solve(linmod_ConsPDS_ip, Tsit5()))
            alloc5 = @allocated(solve(linmod_ConsPDS_ip_2, Tsit5()))
            @test 0.9 < alloc1 / alloc2 < 1.1
            @test 0.9 < alloc1 / alloc3 < 1.1
            @test 0.9 < alloc1 / alloc4 < 1.1
            @test 0.9 < alloc1 / alloc5 < 1.1
        end

        @testset "Lotka-Volterra" begin
            # Lotka-Volterra (nonconservative)
            u0 = [0.9, 0.1]
            tspan = (0.0, 20.0)
            # out-of-place syntax for f
            lotvol(u, p, t) = [u[1] - u[1] * u[2]; u[1] * u[2] - u[2]]
            lotvol_f_op = ODEProblem(lotvol, u0, tspan)
            # in-place syntax for f
            function lotvol!(du, u, p, t)
                u₁, u₂ = u
                du[1] = u₁ - u₁ * u₂
                du[2] = u₁ * u₂ - u₂
            end
            lotvol_f_ip = ODEProblem(lotvol!, u0, tspan)
            # out-of-place syntax for PDS
            lotvolP(u, p, t) = [u[1] 0.0; u[1]*u[2] 0.0]
            lotvolD(u, p, t) = [0.0; u[2]]
            lotvol_PDS_op = PDSProblem(lotvolP, lotvolD, u0, tspan)
            lotvol_PDS_op_2 = PDSProblem{false}(lotvolP, lotvolD, u0, tspan)
            # in-place sytanx for PDS
            function lotvolP!(P, u, p, t)
                P .= 0.0
                P[1, 1] = u[1]
                P[2, 1] = u[2] * u[1]
                return nothing
            end
            function lotvolD!(D, u, p, t)
                D .= 0.0
                D[2] = u[2]
                return nothing
            end
            lotvol_PDS_ip = PDSProblem(lotvolP!, lotvolD!, u0, tspan)
            lotvol_PDS_ip_2 = PDSProblem{true}(lotvolP!, lotvolD!, u0, tspan)

            # solutions
            sol_lotvol_f_op = solve(lotvol_f_op, Tsit5())
            sol_lotvol_f_ip = solve(lotvol_f_ip, Tsit5())
            sol_lotvol_PDS_op = solve(lotvol_PDS_op, Tsit5())
            sol_lotvol_PDS_op_2 = solve(lotvol_PDS_op_2, Tsit5())
            sol_lotvol_PDS_ip = solve(lotvol_PDS_ip, Tsit5())
            sol_lotvol_PDS_ip_2 = solve(lotvol_PDS_ip_2, Tsit5())

            # check equality of solutions
            @test sol_lotvol_f_op.t ≈ sol_lotvol_f_ip.t ≈
                  sol_lotvol_PDS_op.t ≈ sol_lotvol_PDS_op_2.t ≈
                  sol_lotvol_PDS_ip.t ≈ sol_lotvol_PDS_ip_2.t
            @test sol_lotvol_f_op.u ≈ sol_lotvol_f_ip.u ≈
                  sol_lotvol_PDS_op.u ≈ sol_lotvol_PDS_op_2.u ≈
                  sol_lotvol_PDS_ip.u ≈ sol_lotvol_PDS_ip_2.u

            # check that we really do not use too many additional allocations for in-place implementations
            alloc1 = @allocated(solve(lotvol_f_ip, Tsit5()))
            alloc2 = @allocated(solve(lotvol_PDS_ip, Tsit5()))
            alloc3 = @allocated(solve(lotvol_PDS_ip_2, Tsit5()))
            @test 0.9 < alloc1 / alloc2 < 1.1
            @test 0.9 < alloc1 / alloc3 < 1.1
        end

        @testset "Linear advection" begin
            N = 1000 # number of nodes
            u0 = sin.(π * LinRange(0.0, 1.0, N + 1))[2:end] # initial values
            tspan = (0.0, 1.0)
            # Linear advection discretized with finite differences and upwind, periodic boundary conditions
            linear_advection_fd_upwind_ODE = ODEProblem(linear_advection_fd_upwind_f!, u0,
                                                        tspan)
            # problem with dense matrices
            linear_advection_fd_upwind_PDS_dense = PDSProblem(linear_advection_fd_upwind_P!,
                                                              linear_advection_fd_upwind_D!,
                                                              u0, tspan)
            # problem with sparse matrices
            p_prototype = spdiagm(-1 => ones(eltype(u0), N - 1),
                                  N - 1 => ones(eltype(u0), 1))
            linear_advection_fd_upwind_PDS_sparse = PDSProblem(linear_advection_fd_upwind_P!,
                                                               linear_advection_fd_upwind_D!,
                                                               u0, tspan;
                                                               p_prototype = p_prototype)
            linear_advection_fd_upwind_PDS_sparse_2 = PDSProblem{true}(linear_advection_fd_upwind_P!,
                                                                       linear_advection_fd_upwind_D!,
                                                                       u0, tspan;
                                                                       p_prototype = p_prototype)
            linear_advection_fd_upwind_ConsPDS_sparse = ConservativePDSProblem(linear_advection_fd_upwind_P!,
                                                                               u0, tspan;
                                                                               p_prototype = p_prototype)
            linear_advection_fd_upwind_ConsPDS_sparse_2 = ConservativePDSProblem{true}(linear_advection_fd_upwind_P!,
                                                                                       u0,
                                                                                       tspan;
                                                                                       p_prototype = p_prototype)

            # solutions
            sol_linear_advection_fd_upwind_ODE = solve(linear_advection_fd_upwind_ODE,
                                                       Tsit5())
            sol_linear_advection_fd_upwind_PDS_dense = solve(linear_advection_fd_upwind_PDS_dense,
                                                             Tsit5())
            sol_linear_advection_fd_upwind_PDS_sparse = solve(linear_advection_fd_upwind_PDS_sparse,
                                                              Tsit5())
            sol_linear_advection_fd_upwind_PDS_sparse_2 = solve(linear_advection_fd_upwind_PDS_sparse_2,
                                                                Tsit5())
            sol_linear_advection_fd_upwind_ConsPDS_sparse = solve(linear_advection_fd_upwind_ConsPDS_sparse,
                                                                  Tsit5())
            sol_linear_advection_fd_upwind_ConsPDS_sparse_2 = solve(linear_advection_fd_upwind_ConsPDS_sparse_2,
                                                                    Tsit5())

            # check equality of solutions
            @test sol_linear_advection_fd_upwind_ODE.t ≈
                  sol_linear_advection_fd_upwind_PDS_dense.t ≈
                  sol_linear_advection_fd_upwind_PDS_sparse.t ≈
                  sol_linear_advection_fd_upwind_PDS_sparse_2.t ≈
                  sol_linear_advection_fd_upwind_ConsPDS_sparse.t ≈
                  sol_linear_advection_fd_upwind_ConsPDS_sparse_2.t
            @test sol_linear_advection_fd_upwind_ODE.u ≈
                  sol_linear_advection_fd_upwind_PDS_dense.u ≈
                  sol_linear_advection_fd_upwind_PDS_sparse.u ≈
                  sol_linear_advection_fd_upwind_PDS_sparse_2.u ≈
                  sol_linear_advection_fd_upwind_ConsPDS_sparse.u ≈
                  sol_linear_advection_fd_upwind_ConsPDS_sparse_2.u

            # Check that we really do not use too many additional allocations
            alloc1 = @allocated(solve(linear_advection_fd_upwind_ODE, Tsit5()))
            alloc2 = @allocated(solve(linear_advection_fd_upwind_PDS_dense, Tsit5()))
            alloc3 = @allocated(solve(linear_advection_fd_upwind_PDS_sparse, Tsit5()))
            alloc4 = @allocated(solve(linear_advection_fd_upwind_PDS_sparse_2, Tsit5()))
            alloc5 = @allocated(solve(linear_advection_fd_upwind_ConsPDS_sparse, Tsit5()))
            alloc6 = @allocated(solve(linear_advection_fd_upwind_ConsPDS_sparse_2, Tsit5()))
            @test 0.9 < alloc1 / alloc2 < 1.1
            @test 0.9 < alloc1 / alloc3 < 1.1
            @test 0.9 < alloc1 / alloc4 < 1.1
            @test 0.9 < alloc1 / alloc5 < 1.1
            @test 0.9 < alloc1 / alloc6 < 1.1
        end

        # Here we check that PDSFunctions and ConservativePDSFunctions can be evaluated
        # at (u, p, t) also in cases where u and t carry units.
        # This is required when solvers form OrdinarDiffEq are used to solve a PDS
        # and std_rhs is not specified.
        @testset "(Conservative)PDSFunction and units" begin
            @testset "(Conservative)PDSFunction and units (out-of-place)" begin
                u0 = [0.9u"N", 0.1u"N"]
                u0_static = @SVector [0.9u"N", 0.1u"N"]
                t0 = 0.0u"s"
                tspan = (t0, 2.0u"s")

                f(u, p, t) = [u[2] - 5 * u[1]; -u[2] + 5 * u[1]] / u"s"
                P(u, p, t) = [0u"N/s" u[2]/u"s"; 5 * u[1]/u"s" 0u"N/s"]
                D(u, p, t) = [0.0; 0.0]u"N/s"
                g = PositiveIntegrators.ConservativePDSFunction{false}(P)
                h = PositiveIntegrators.PDSFunction{false}(P, D)

                f_static(u, p, t) = SA[(u[2] - 5 * u[1]) / u"s"; (-u[2] + 5 * u[1]) / u"s"]
                P_static(u, p, t) = SA[0u"N/s" u[2]/u"s"; 5 * u[1]/u"s" 0u"N/s"]
                D_static(u, p, t) = SA[0.0u"N/s"; 0.0u"N/s"]
                g_static = PositiveIntegrators.ConservativePDSFunction{false}(P_static)
                h_static = PositiveIntegrators.PDSFunction{false}(P_static, D_static)

                f1 = f(u0, nothing, t0)
                g1 = g(u0, nothing, t0)
                h1 = h(u0, nothing, t0)
                f2 = f_static(u0_static, nothing, t0)
                g2 = g_static(u0_static, nothing, t0)
                h2 = h_static(u0_static, nothing, t0)

                @test f1 == g1 == h1 == f2 == g2 == h2
            end
            @testset "(Conservative)PDSFunction and units (in-place)" begin
                u0 = [1.0, 1.5, 2.0, 2.5]u"N"
                t0 = 0.0u"s"

                function f!(du, u, p, t)
                    fill!(du, zero(eltype(du)))

                    du[1] = (u[1] - u[2]) / u"s"
                    du[2] = (3 * u[2] - 2 * u[3] - u[1]) / u"s"
                    du[3] = (-2 * u[2] + 5 * u[3] - 3 * u[4]) / u"s"
                    du[4] = (3 * u[4] - 3 * u[3]) / u"s"

                    return nothing
                end
                function P!(P, u, p, t)
                    fill!(P, zero(eltype(P)))
                    for i in 1:(length(u) - 1)
                        P[i, i + 1] = i * u[i] / u"s"
                        P[i + 1, i] = i * u[i + 1] / u"s"
                    end
                    return nothing
                end
                function D!(D, u, p, t)
                    fill!(D, zero(eltype(D)))
                    return nothing
                end

                P_tridiagonal = Tridiagonal([0.1, 0.2, 0.3], zeros(4),
                                            [0.4, 0.5, 0.6])u"N/s"
                P_dense = Matrix(P_tridiagonal)
                P_sparse = sparse(P_tridiagonal)
                D_dense = zeros(4)u"N/s"

                f_dense! = PositiveIntegrators.ConservativePDSFunction{true}(P!,
                                                                             p_prototype = P_dense)
                f_tridiagonal! = PositiveIntegrators.ConservativePDSFunction{true}(P!,
                                                                                   p_prototype = P_tridiagonal)
                f_sparse! = PositiveIntegrators.ConservativePDSFunction{true}(P!,
                                                                              p_prototype = P_sparse)
                g_dense! = PositiveIntegrators.PDSFunction{true}(P!, D!,
                                                                 p_prototype = P_dense,
                                                                 d_prototype = D_dense)
                g_tridiagonal! = PositiveIntegrators.PDSFunction{true}(P!, D!,
                                                                       p_prototype = P_tridiagonal,
                                                                       d_prototype = D_dense)
                g_sparse! = PositiveIntegrators.PDSFunction{true}(P!, D!,
                                                                  p_prototype = P_sparse,
                                                                  d_prototype = D_dense)

                du1 = du2 = du3 = du4 = du5 = du6 = du7 = zeros(4)u"N/s"

                f!(du1, u0, nothing, t0)
                f_dense!(du2, u0, nothing, t0)
                f_tridiagonal!(du3, u0, nothing, t0)
                f_sparse!(du4, u0, nothing, t0)
                g_dense!(du5, u0, nothing, t0)
                g_tridiagonal!(du6, u0, nothing, t0)
                g_sparse!(du7, u0, nothing, t0)

                @test du1 == du2 == du3 == du4 == du5 == du6 == du7
            end
        end

        @testset "Compatibility of OrdinarDiffEq solvers with (Conservative)PDSProblems and units" begin
            @testset "out-of-place" begin
                u0 = [0.9u"N", 0.1u"N"]
                u0_static = @SVector [0.9u"N", 0.1u"N"]
                tspan = (0.0u"s", 2.0u"s")

                f(u, p, t) = [u[2] - 5 * u[1]; -u[2] + 5 * u[1]] / u"s"
                P(u, p, t) = [0u"N/s" u[2]/u"s"; 5 * u[1]/u"s" 0u"N/s"]
                D(u, p, t) = [0.0; 0.0]u"N/s"
                f_static(u, p, t) = SA[(u[2] - 5 * u[1]) / u"s"; (-u[2] + 5 * u[1]) / u"s"]
                P_static(u, p, t) = SA[0u"N/s" u[2]/u"s"; 5 * u[1]/u"s" 0u"N/s"]
                D_static(u, p, t) = SA[0.0u"N/s"; 0.0u"N/s"]

                prob_ode = ODEProblem(f, u0, tspan)
                prob_cpds = ConservativePDSProblem(P, u0, tspan)
                prob_pds = PDSProblem(P, D, u0, tspan)
                prob_cpds_rhs = ConservativePDSProblem(P, u0, tspan; std_rhs = f)
                prob_pds_rhs = PDSProblem(P, D, u0, tspan; std_rhs = f)
                prob_ode_static = ODEProblem(f, u0_static, tspan)
                prob_cpds_static = ConservativePDSProblem(P, u0_static, tspan)
                prob_pds_static = PDSProblem(P, D, u0_static, tspan)
                prob_cpds_static_rhs = ConservativePDSProblem(P, u0_static, tspan;
                                                              std_rhs = f_static)
                prob_pds_static_rhs = PDSProblem(P, D, u0_static, tspan; std_rhs = f_static)

                # implicit solvers and units don't work
                # algs = (Euler(), ImplicitEuler(), Tsit5(), Rosenbrock23(), SDIRK2(),TRBDF2())
                algs = (Euler(), Tsit5(), Vern9())

                for alg in algs
                    sol_ode = solve(prob_ode, alg; dt = 0.1u"s")
                    sol_cpds = solve(prob_cpds, alg; dt = 0.1u"s")
                    sol_pds = solve(prob_pds, alg; dt = 0.1u"s")
                    sol_cpds_rhs = solve(prob_cpds_rhs, alg; dt = 0.1u"s")
                    sol_pds_rhs = solve(prob_pds_rhs, alg; dt = 0.1u"s")
                    sol_ode_static = solve(prob_ode, alg; dt = 0.1u"s")
                    sol_cpds_static = solve(prob_cpds, alg; dt = 0.1u"s")
                    sol_pds_static = solve(prob_pds, alg; dt = 0.1u"s")
                    sol_cpds_static_rhs = solve(prob_cpds_rhs, alg; dt = 0.1u"s")
                    sol_pds_static_rhs = solve(prob_pds_rhs, alg; dt = 0.1u"s")

                    @test sol_ode.t ≈ sol_cpds.t
                    @test sol_ode.t ≈ sol_pds.t
                    @test sol_ode.t ≈ sol_cpds_rhs.t
                    @test sol_ode.t ≈ sol_pds_rhs.t
                    @test sol_ode.t ≈ sol_ode_static.t
                    @test sol_ode.t ≈ sol_cpds_static.t
                    @test sol_ode.t ≈ sol_pds_static.t
                    @test sol_ode.t ≈ sol_cpds_static_rhs.t
                    @test sol_ode.t ≈ sol_pds_static_rhs.t

                    # isapprox(sol0.u,sol.u) tries to check norm(sol0.u-sol.u) < 0.0,
                    # which fails, because we have units on the left, but not on the right
                    # In addition, abstol with units is not allowed.
                    # Hence, we strip the units.
                    # For Vector{Quantity{}} like sol.t there is an isapprox method in quantities.jl, which
                    # takes care of the stripping.
                    @test ustrip.(sol_ode.u) ≈ ustrip.(sol_cpds.u)
                    @test ustrip.(sol_ode.u) ≈ ustrip.(sol_pds.u)
                    @test ustrip.(sol_ode.u) ≈ ustrip.(sol_cpds_rhs.u)
                    @test ustrip.(sol_ode.u) ≈ ustrip.(sol_pds_rhs.u)
                    @test ustrip.(sol_ode.u) ≈ ustrip.(sol_ode_static.u)
                    @test ustrip.(sol_ode.u) ≈ ustrip.(sol_cpds_static.u)
                    @test ustrip.(sol_ode.u) ≈ ustrip.(sol_pds_static.u)
                    @test ustrip.(sol_ode.u) ≈ ustrip.(sol_cpds_static_rhs.u)
                    @test ustrip.(sol_ode.u) ≈ ustrip.(sol_pds_static_rhs.u)
                end
            end
            @testset "in-place" begin
                u0 = [1.0, 1.5, 2.0, 2.5]u"N"
                tspan = (0.0u"s", 1.0u"s")

                function f!(du, u, p, t)
                    fill!(du, zero(eltype(du)))

                    du[1] = (u[1] - u[2]) / u"s"
                    du[2] = (3 * u[2] - 2 * u[3] - u[1]) / u"s"
                    du[3] = (-2 * u[2] + 5 * u[3] - 3 * u[4]) / u"s"
                    du[4] = (3 * u[4] - 3 * u[3]) / u"s"

                    return nothing
                end
                function P!(P, u, p, t)
                    fill!(P, zero(eltype(P)))
                    for i in 1:(length(u) - 1)
                        P[i, i + 1] = i * u[i] / u"s"
                        P[i + 1, i] = i * u[i + 1] / u"s"
                    end
                    return nothing
                end
                function D!(D, u, p, t)
                    fill!(D, zero(eltype(D)))
                    return nothing
                end

                P_tridiagonal = Tridiagonal([0.1, 0.2, 0.3], zeros(4),
                                            [0.4, 0.5, 0.6])u"N/s"
                P_dense = Matrix(P_tridiagonal)
                P_sparse = sparse(P_tridiagonal)
                D_dense = zeros(4)u"N/s"

                # implicit solvers and units don't work
                #algs = (Euler(), ImplicitEuler(), Tsit5(), Rosenbrock23(), SDIRK2(), TRBDF2())
                algs = (Euler(), Tsit5(), Vern9())

                prob0 = ODEProblem(f!, u0, tspan)

                probs = Array{Any}(undef, 14)
                probs[1] = ConservativePDSProblem(P!, u0, tspan)
                probs[2] = ConservativePDSProblem(P!, u0, tspan; p_prototype = P_dense)
                probs[3] = ConservativePDSProblem(P!, u0, tspan; p_prototype = P_dense,
                                                  std_rhs = f!)
                probs[4] = ConservativePDSProblem(P!, u0, tspan;
                                                  p_prototype = P_tridiagonal)
                probs[5] = ConservativePDSProblem(P!, u0, tspan;
                                                  p_prototype = P_tridiagonal, std_rhs = f!)
                probs[6] = ConservativePDSProblem(P!, u0, tspan; p_prototype = P_sparse)
                probs[7] = ConservativePDSProblem(P!, u0, tspan; p_prototype = P_sparse,
                                                  std_rhs = f!)
                probs[8] = PDSProblem(P!, D!, u0, tspan)
                probs[9] = PDSProblem(P!, D!, u0, tspan; p_prototype = P_dense)
                probs[10] = PDSProblem(P!, D!, u0, tspan; p_prototype = P_dense,
                                       std_rhs = f!)
                probs[11] = PDSProblem(P!, D!, u0, tspan; p_prototype = P_tridiagonal)
                probs[12] = PDSProblem(P!, D!, u0, tspan; p_prototype = P_tridiagonal,
                                       std_rhs = f!)
                probs[13] = PDSProblem(P!, D!, u0, tspan; p_prototype = P_sparse)
                probs[14] = PDSProblem(P!, D!, u0, tspan; p_prototype = P_sparse,
                                       std_rhs = f!)

                for alg in algs
                    sol0 = solve(prob0, alg; dt = 0.2u"s")
                    for prob in probs
                        sol = solve(prob, alg; dt = 0.2u"s")

                        @test sol0.t ≈ sol.t
                        # isapprox(sol0.u,sol.u) tries to check norm(sol0.u-sol.u) < 0.0,
                        # which fails, because we have units on the left, but not on the right
                        # In addition, abstol with units is not allowed.
                        # Hence, we strip the units.
                        # For Vector{Quantity{}} like sol.t there is an isapprox method in quantities.jl, which
                        # takes care of the stripping.
                        @test ustrip.(sol0.u) ≈ ustrip.(sol.u)
                    end
                end
            end
        end

        # Here we check that production-destruction form and standard ODE form fit together in predefinded problems,
        # i.e. standard solvers using std_rhs should generate results that are equal to those without specifying std_rhs
        @testset "Check that production-destruction form and standard ODE fit together in predefinded problems" begin
            # non-stiff conservative problems (out-of-place)
            probs = (prob_pds_linmod, prob_pds_nonlinmod, prob_pds_brusselator,
                     prob_pds_sir, prob_pds_npzd)
            algs = (Euler(), ImplicitEuler(), Tsit5(), Rosenbrock23(), SDIRK2(), TRBDF2())
            @testset "$alg" for prob in probs, alg in algs
                dt = (last(prob.tspan) - first(prob.tspan)) / 1e4
                sol = solve(prob, alg; dt, isoutofdomain = isnegative) # use explicit f
                sol2 = solve(ConservativePDSProblem(prob.f.p, prob.u0, prob.tspan), alg; dt,
                             isoutofdomain = isnegative) # use p and d to compute f
                sol3 = solve(ODEProblem(prob.f.std_rhs, prob.u0, prob.tspan), alg; dt,
                             isoutofdomain = isnegative) # use f to create ODEProblem
                @test sol.t ≈ sol2.t ≈ sol3.t
                @test sol.u ≈ sol2.u ≈ sol3.u
            end

            # non-stiff conservative problems (in-place)
            # Requires autodiff = AutoFiniteDiff()
            probs = (prob_pds_linmod_inplace,)
            algs = (Euler(), ImplicitEuler(autodiff = AutoFiniteDiff()), Tsit5(),
                    Rosenbrock23(autodiff = AutoFiniteDiff()),
                    SDIRK2(autodiff = AutoFiniteDiff()),
                    TRBDF2(autodiff = AutoFiniteDiff()))
            @testset "$alg" for prob in probs, alg in algs
                dt = (last(prob.tspan) - first(prob.tspan)) / 1e4
                sol = solve(prob, alg; dt, isoutofdomain = isnegative) # use explicit f
                sol2 = solve(ConservativePDSProblem(prob.f.p, prob.u0, prob.tspan), alg; dt,
                             isoutofdomain = isnegative) # use p and d to compute f
                sol3 = solve(ODEProblem(prob.f.std_rhs, prob.u0, prob.tspan), alg; dt,
                             isoutofdomain = isnegative) # use f to create ODEProblem
                @test sol.t ≈ sol2.t ≈ sol3.t
                @test sol.u ≈ sol2.u ≈ sol3.u
            end

            # non-stiff non-conservative problems (out-of-place)
            probs = (prob_pds_minmapk,)
            algs = (Euler(), ImplicitEuler(), Tsit5(), Rosenbrock23(), SDIRK2(), TRBDF2())
            @testset "$alg" for prob in probs, alg in algs
                dt = (last(prob.tspan) - first(prob.tspan)) / 1e4
                sol = solve(prob, alg; dt, isoutofdomain = isnegative) # use explicit f
                sol2 = solve(PDSProblem(prob.f.p, prob.f.d, prob.u0, prob.tspan), alg; dt,
                             isoutofdomain = isnegative) # use p and d to compute f
                sol3 = solve(ODEProblem(prob.f.std_rhs, prob.u0, prob.tspan), alg; dt,
                             isoutofdomain = isnegative) # use f to create ODEProblem
                @test sol.t ≈ sol2.t ≈ sol3.t
                @test sol.u ≈ sol2.u ≈ sol3.u
            end

            # Robertson problem
            prob = prob_pds_robertson
            algs = (ImplicitEuler(), Rosenbrock23(), SDIRK2(), TRBDF2())
            @testset "$alg" for alg in algs
                dt = 1e-6
                sol = solve(prob, alg; dt) # use explicit f
                sol2 = solve(ConservativePDSProblem(prob.f.p, prob.u0, prob.tspan), alg; dt) # use p and d to compute f
                sol3 = solve(ODEProblem(prob.f.std_rhs, prob.u0, prob.tspan), alg; dt) # use f to create ODEProblem
                @test sol.t ≈ sol2.t ≈ sol3.t
                @test sol.u ≈ sol2.u ≈ sol3.u
            end

            # Bertolazzi problem
            # Did not find any solver configuration to compute a reasonable solution and pass tests.
            # - constant time stepping requires very small dt
            # - adaptive time stepping generates solutions with different number of time steps
            #
            # Nevertheless, the following code shows that the same problem is solved in each case
            # prob = prob_pds_bertolazzi
            # alg = ImplicitEuler()
            # dt = 1e-3
            # sol = solve(prob, alg; dt) # use explicit f
            # sol2 = solve(ConservativePDSProblem(prob.f.p, prob.u0, prob.tspan), alg; dt) # use p and d to compute f
            # sol3 = solve(ODEProblem(prob.f.std_rhs, prob.u0, prob.tspan), alg; dt) # use f to create ODEProblem
            # plot(plot(sol),plot(sol2),plot(sol3))

            # Stratospheric reaction problem
            prob = prob_pds_stratreac
            algs = (ImplicitEuler(), Rosenbrock23(), TRBDF2())
            @testset "$alg" for alg in algs
                dt = 1.0
                sol = solve(prob, alg; dt) # use explicit f
                sol2 = solve(PDSProblem(prob.f.p, prob.f.d, prob.u0, prob.tspan), alg; dt) # use p and d to compute f
                sol3 = solve(ODEProblem(prob.f.std_rhs, prob.u0, prob.tspan), alg; dt) # use f to create ODEProblem
                @test sol.t ≈ sol2.t ≈ sol3.t
                @test sol.u ≈ sol2.u ≈ sol3.u
            end
        end
    end

    @testset "PDS Solvers" begin
        # Here we check that MPRK schemes require a PDSProblem or ConservativePDSProblem.
        # We also check that only permissible parameters can be used.
        @testset "Error handling" begin
            f_oop = (u, p, t) -> u
            f_ip = (du, u, p, t) -> du .= u
            prob_oop = ODEProblem(f_oop, [1.0; 2.0], (0.0, 1.0))
            prob_ip = ODEProblem(f_ip, [1.0; 2.0], (0.0, 1.0))
            @test_throws "MPE can only be applied to production-destruction systems" solve(prob_oop,
                                                                                           MPE(),
                                                                                           dt = 0.25)
            @test_throws "MPE can only be applied to production-destruction systems" solve(prob_ip,
                                                                                           MPE(),
                                                                                           dt = 0.25)
            @test_throws "MPRK22 can only be applied to production-destruction systems" solve(prob_oop,
                                                                                              MPRK22(1.0))
            @test_throws "MPRK22 can only be applied to production-destruction systems" solve(prob_ip,
                                                                                              MPRK22(1.0))
            @test_throws "MPRK22 requires α ≥ 1/2." solve(prob_pds_linmod, MPRK22(0.25))
            @test_throws "MPRK43 can only be applied to production-destruction systems" solve(prob_oop,
                                                                                              MPRK43I(1.0,
                                                                                                      0.5))
            @test_throws "MPRK43 can only be applied to production-destruction systems" solve(prob_ip,
                                                                                              MPRK43I(1.0,
                                                                                                      0.5))
            @test_throws "MPRK43 can only be applied to production-destruction systems" solve(prob_oop,
                                                                                              MPRK43II(0.5))
            @test_throws "MPRK43 can only be applied to production-destruction systems" solve(prob_ip,
                                                                                              MPRK43II(0.5))
            @test_throws "MPRK43I requires α ≥ 1/3 and α ≠ 2/3." solve(prob_pds_linmod,
                                                                       MPRK43I(0.0, 0.5))
            @test_throws "MPRK43I requires α ≥ 1/3 and α ≠ 2/3." solve(prob_pds_linmod,
                                                                       MPRK43I(2.0 / 3.0,
                                                                               0.5))
            @test_throws "For this choice of α MPRK43I requires 2/3 ≤ β ≤ 3α(1-α)." solve(prob_pds_linmod,
                                                                                          MPRK43I(0.5,
                                                                                                  0.5))
            @test_throws "For this choice of α MPRK43I requires 2/3 ≤ β ≤ 3α(1-α)." solve(prob_pds_linmod,
                                                                                          MPRK43I(0.5,
                                                                                                  0.8))
            @test_throws "For this choice of α MPRK43I requires 3α(1-α) ≤ β ≤ 2/3." solve(prob_pds_linmod,
                                                                                          MPRK43I(0.7,
                                                                                                  0.7))
            @test_throws "For this choice of α MPRK43I requires 3α(1-α) ≤ β ≤ 2/3." solve(prob_pds_linmod,
                                                                                          MPRK43I(0.7,
                                                                                                  0.1))
            @test_throws "For this choice of α MPRK43I requires (3α-2)/(6α-3) ≤ β ≤ 2/3." solve(prob_pds_linmod,
                                                                                                MPRK43I(1.0,
                                                                                                        0.7))
            @test_throws "For this choice of α MPRK43I requires (3α-2)/(6α-3) ≤ β ≤ 2/3." solve(prob_pds_linmod,
                                                                                                MPRK43I(1.0,
                                                                                                        0.1))
            @test_throws "MPRK43II requires 3/8 ≤ γ ≤ 3/4." solve(prob_pds_linmod,
                                                                  MPRK43II(1.0))
            @test_throws "MPRK43II requires 3/8 ≤ γ ≤ 3/4." solve(prob_pds_linmod,
                                                                  MPRK43II(0.0))
            @test_throws "SSPMPRK22 can only be applied to production-destruction systems" solve(prob_oop,
                                                                                                 SSPMPRK22(0.5,
                                                                                                           1.0))
            @test_throws "SSPMPRK22 can only be applied to production-destruction systems" solve(prob_ip,
                                                                                                 SSPMPRK22(0.5,
                                                                                                           1.0))
            @test_throws "SSPMPRK22 requires 0 ≤ α ≤ 1, β ≥ 0 and αβ + 1/(2β) ≤ 1." solve(prob_pds_linmod,
                                                                                          SSPMPRK22(-1.0,
                                                                                                    1.0))
            @test_throws "SSPMPRK22 requires 0 ≤ α ≤ 1, β ≥ 0 and αβ + 1/(2β) ≤ 1." solve(prob_pds_linmod,
                                                                                          SSPMPRK22(0.0,
                                                                                                    -1.0))
            @test_throws "SSPMPRK22 requires 0 ≤ α ≤ 1, β ≥ 0 and αβ + 1/(2β) ≤ 1." solve(prob_pds_linmod,
                                                                                          SSPMPRK22(1.0,
                                                                                                    10.0))
            @test_throws "SSPMPRK43 can only be applied to production-destruction systems" solve(prob_oop,
                                                                                                 SSPMPRK43(),
                                                                                                 dt = 0.1)
            @test_throws "SSPMPRK43 can only be applied to production-destruction systems" solve(prob_ip,
                                                                                                 SSPMPRK43(),
                                                                                                 dt = 0.1)
            @test_throws "MPDeC can only be applied to production-destruction systems" solve(prob_ip,
                                                                                             MPDeC(2),
                                                                                             dt = 0.1)
            @test_throws "MPDeC can only be applied to production-destruction systems" solve(prob_oop,
                                                                                             MPDeC(2),
                                                                                             dt = 0.1)
            @test_throws "MPDeC requires 2 ≤ K ≤ 10." solve(prob_pds_linmod, MPDeC(123))
            @test_throws "MPDeC requires 2 ≤ K ≤ 10." solve(prob_pds_linmod,
                                                            MPDeC(123, nodes = :lagrange))
            P = spdiagm(1 => [1.0])
            function prod!(P, u, p, t)
                P[2, 1] = one(eltype(P))
            end
            prob = ConservativePDSProblem(prod!, ones(2, 1), (0.0, 1.0); p_prototype = P)
            @test_throws "must not alter the sparsity pattern" solve(prob, MPDeC(2))
        end

        # Here we check that algorithms which accept input parameters return constants
        # of the same type as the inputs
        @testset "Constant types" begin
            algs = (MPRK22(0.5f0), MPRK22(1.0f0), MPRK22(2.0f0), MPRK43I(1.0f0, 0.5f0),
                    MPRK43I(0.5f0, 0.75f0), MPRK43II(0.5f0), MPRK43II(2.0f0 / 3.0f0),
                    SSPMPRK22(0.5f0, 1.0f0))
            @testset "$i" for (i, alg) in enumerate(algs)
                @test eltype(PositiveIntegrators.get_constant_parameters(alg)) == Float32
            end

            algs = (MPRK22(0.5), MPRK22(1.0), MPRK22(2.0), MPRK43I(1.0, 0.5),
                    MPRK43I(0.5, 0.75), MPRK43II(0.5), MPRK43II(2.0 / 3.0),
                    SSPMPRK22(0.5, 1.0))
            @testset "$i" for (i, alg) in enumerate(algs)
                @test eltype(PositiveIntegrators.get_constant_parameters(alg)) == Float64
            end
        end

        # We check that Float32 inputs result in Float32 outputs
        # TODO: Use a tool such as LIKWID to check that no double-precision calculations are performed internally
        @testset "Float32 input data" begin
            P_linmod(u, p, t) = [0 u[2]; 5*u[1] 0]
            u0 = [0.9f0, 0.1f0]
            prob = ConservativePDSProblem(P_linmod, u0, (0.0f0, 2.0f0))

            algs = (MPRK22(0.5f0), MPRK22(1.0f0), MPRK22(2.0f0), MPRK43I(1.0f0, 0.5f0),
                    MPRK43I(0.5f0, 0.75f0), MPRK43II(0.5f0), MPRK43II(2.0f0 / 3.0f0),
                    SSPMPRK22(0.5f0, 1.0f0), SSPMPRK43(), MPDeC(2),
                    MPDeC(2, nodes = :lagrange))
            for alg in algs
                sol = solve(prob, alg; dt = 0.1f0, save_everystep = false)
                @test sol.t isa Vector{Float32}
                @test sol.u isa Vector{Vector{Float32}}
            end
        end

        # Here we check that MPE equals implicit Euler (IE) for a linear PDS
        @testset "Linear model problem: MPE = IE?" begin
            # problem data
            u0 = [0.9, 0.1]
            tspan = (0.0, 2.0)
            p = [5.0, 1.0]

            # analytic solution
            function f_analytic(u0, p, t)
                u₁⁰, u₂⁰ = u0
                a, b = p
                c = a + b
                return ((u₁⁰ + u₂⁰) * [b; a] +
                        exp(-c * t) * (a * u₁⁰ - b * u₂⁰) * [1; -1]) / c
            end

            # linear model problem - out-of-place
            linmodP(u, p, t) = [0 p[2]*u[2]; p[1]*u[1] 0]
            linmodD(u, p, t) = [0.0; 0.0]
            prob_op = PDSProblem(linmodP, linmodD, u0, tspan, p; analytic = f_analytic)
            prob_op_2 = ConservativePDSProblem(linmodP, u0, tspan, p; analytic = f_analytic)

            dt = 0.25
            sol_MPE_op = solve(prob_op, MPE(); dt)
            sol_MPE_op_2 = solve(prob_op_2, MPE(); dt)
            sol_IE_op = solve(prob_op, ImplicitEuler(autodiff = AutoFiniteDiff());
                              dt, adaptive = false)
            @test sol_MPE_op.t ≈ sol_MPE_op_2.t ≈ sol_IE_op.t
            @test sol_MPE_op.u ≈ sol_MPE_op_2.u ≈ sol_IE_op.u

            # linear model problem - in-place
            function linmodP!(P, u, p, t)
                fill!(P, zero(eltype(P)))
                P[1, 2] = p[2] * u[2]
                P[2, 1] = p[1] * u[1]
                return nothing
            end
            function linmodD!(D, u, p, t)
                fill!(D, zero(eltype(D)))
                return nothing
            end
            prob_ip = PDSProblem(linmodP!, linmodD!, u0, tspan, p; analytic = f_analytic)
            prob_ip_2 = ConservativePDSProblem(linmodP!, u0, tspan, p;
                                               analytic = f_analytic)

            dt = 0.25
            sol_MPE_ip = solve(prob_ip, MPE(); dt)
            sol_MPE_ip_2 = solve(prob_ip_2, MPE(); dt)
            sol_IE_ip = solve(prob_ip, ImplicitEuler(autodiff = AutoFiniteDiff());
                              dt, adaptive = false)
            @test sol_MPE_ip.t ≈ sol_MPE_ip_2.t ≈ sol_IE_ip.t
            @test sol_MPE_ip.u ≈ sol_MPE_ip_2.u ≈ sol_IE_ip.u
        end

        # Here we check that MPRK22(α) = SSPMPRK22(0,α)
        @testset "MPRK22(α) = SSPMPRK22(0, α)" begin
            for α in (0.5, 2.0 / 3.0, 1.0, 2.0)
                # conservative PDS
                sol1 = solve(prob_pds_linmod, MPRK22(α))
                sol2 = solve(prob_pds_linmod, SSPMPRK22(0.0, α))
                sol3 = solve(prob_pds_linmod_inplace, MPRK22(α))
                sol4 = solve(prob_pds_linmod_inplace, SSPMPRK22(0.0, α))
                @test sol1.u ≈ sol2.u ≈ sol3.u ≈ sol4.u

                # nonconservative PDS
                sol1 = solve(prob_pds_linmod_nonconservative, MPRK22(α))
                sol2 = solve(prob_pds_linmod_nonconservative, SSPMPRK22(0.0, α))
                sol3 = solve(prob_pds_linmod_nonconservative_inplace, MPRK22(α))
                sol4 = solve(prob_pds_linmod_nonconservative_inplace, SSPMPRK22(0.0, α))
                @test sol1.u ≈ sol2.u ≈ sol3.u ≈ sol4.u
            end
        end

        # Here we check that MPRK22(1.0) = MPDeC(2)
        @testset "MPRK22(1.0) = MPDeC(2)" begin
            # conservative PDS
            sol1 = solve(prob_pds_linmod, MPRK22(1.0))
            sol2 = solve(prob_pds_linmod, MPDeC(2))
            sol3 = solve(prob_pds_linmod, MPDeC(2, nodes = :lagrange))
            sol4 = solve(prob_pds_linmod_inplace, MPRK22(1.0))
            sol5 = solve(prob_pds_linmod_inplace, MPDeC(2))
            sol6 = solve(prob_pds_linmod_inplace, MPDeC(2, nodes = :lagrange))
            @test sol1.u ≈ sol2.u ≈ sol3.u ≈ sol4.u ≈ sol5.u ≈ sol6.u

            # nonconservative PDS
            sol1 = solve(prob_pds_linmod_nonconservative, MPRK22(1.0))
            sol2 = solve(prob_pds_linmod_nonconservative, MPDeC(2))
            sol3 = solve(prob_pds_linmod_nonconservative, MPDeC(2, nodes = :lagrange))
            sol4 = solve(prob_pds_linmod_nonconservative_inplace, MPRK22(1.0))
            sol5 = solve(prob_pds_linmod_nonconservative_inplace, MPDeC(2))
            sol6 = solve(prob_pds_linmod_nonconservative_inplace,
                         MPDeC(2, nodes = :lagrange))
            @test sol1.u ≈ sol2.u ≈ sol3.u ≈ sol4.u ≈ sol5.u ≈ sol6.u
        end

        # Here we check that different linear solvers can be used
        @testset "Different linear solvers" begin
            # problem data
            u0 = [0.9, 0.1]
            tspan = (0.0, 2.0)
            p = [5.0, 1.0]

            # analytic solution
            function f_analytic(u0, p, t)
                u₁⁰, u₂⁰ = u0
                a, b = p
                c = a + b
                return ((u₁⁰ + u₂⁰) * [b; a] +
                        exp(-c * t) * (a * u₁⁰ - b * u₂⁰) * [1; -1]) / c
            end

            # linear model problem - in-place
            function linmodP!(P, u, p, t)
                fill!(P, zero(eltype(P)))
                P[1, 2] = p[2] * u[2]
                P[2, 1] = p[1] * u[1]
                return nothing
            end
            function linmodD!(D, u, p, t)
                fill!(D, zero(eltype(D)))
                return nothing
            end
            prob_ip = PDSProblem(linmodP!, linmodD!, u0, tspan, p; analytic = f_analytic)
            prob_ip_2 = ConservativePDSProblem(linmodP!, u0, tspan, p;
                                               analytic = f_analytic)

            algs = [MPE, (; kwargs...) -> MPRK22(1.0; kwargs...),
                (; kwargs...) -> MPRK22(0.5; kwargs...),
                (; kwargs...) -> MPRK22(2.0; kwargs...),
                (; kwargs...) -> MPRK43I(1.0, 0.5; kwargs...),
                (; kwargs...) -> MPRK43I(0.5, 0.75; kwargs...),
                (; kwargs...) -> MPRK43II(0.5; kwargs...),
                (; kwargs...) -> MPRK43II(2.0 / 3.0; kwargs...),
                (; kwargs...) -> SSPMPRK22(0.5, 1.0; kwargs...),
                (; kwargs...) -> SSPMPRK43(; kwargs...)]
            for k in 2:10
                push!(algs, (; kwargs...) -> MPDeC(k; kwargs...),
                      (; kwargs...) -> MPDeC(k; nodes = :lagrange, kwargs...))
            end

            for alg in algs
                # Check different linear solvers
                dt = 0.25
                sol1 = solve(prob_ip, alg(); dt)
                sol1_2 = solve(prob_ip_2, alg(); dt)
                sol2 = solve(prob_ip, alg(linsolve = RFLUFactorization()); dt)
                sol2_2 = solve(prob_ip_2, alg(linsolve = RFLUFactorization()); dt)
                sol3 = solve(prob_ip, alg(linsolve = LUFactorization()); dt)
                sol3_2 = solve(prob_ip_2, alg(linsolve = LUFactorization()); dt)
                sol4 = solve(prob_ip, alg(linsolve = KrylovJL_GMRES()); dt)
                sol4_2 = solve(prob_ip_2, alg(linsolve = KrylovJL_GMRES()); dt)
                @test sol1.t ≈ sol2.t ≈ sol3.t ≈ sol4.t ≈ sol1_2.t ≈ sol2_2.t ≈ sol3_2.t ≈
                      sol4_2.t
                @test sol1.u ≈ sol2.u ≈ sol3.u ≈ sol4.u ≈ sol1_2.u ≈ sol2_2.u ≈ sol3_2.u ≈
                      sol4_2.u
            end
        end

        # Here we check that in-place and out-of-place implementations
        # deliver the same results
        @testset "Different matrix types (conservative)" begin
            prod_1! = (P, u, p, t) -> begin
                fill!(P, zero(eltype(P)))
                for i in 1:(length(u) - 1)
                    P[i, i + 1] = i * u[i + 1]
                end
                return nothing
            end

            prod_2! = (P, u, p, t) -> begin
                fill!(P, zero(eltype(P)))
                for i in 1:(length(u) - 1)
                    P[i + 1, i] = i * u[i]
                end
                return nothing
            end

            prod_3! = (P, u, p, t) -> begin
                fill!(P, zero(eltype(P)))
                for i in 1:(length(u) - 1)
                    P[i, i + 1] = i * u[i + 1]
                    P[i + 1, i] = i * u[i]
                end
                return nothing
            end

            n = 4
            P_tridiagonal = Tridiagonal([0.1, 0.2, 0.3],
                                        zeros(n),
                                        [0.4, 0.5, 0.6])
            P_dense = Matrix(P_tridiagonal)
            P_sparse = sparse(P_tridiagonal)
            u0 = [1.0, 1.5, 2.0, 2.5]
            tspan = (0.0, 1.0)
            dt = 0.25

            algs = [
                MPE(),
                MPRK22(0.5),
                MPRK22(1.0),
                MPRK43I(1.0, 0.5),
                MPRK43I(0.5, 0.75),
                MPRK43II(2.0 / 3.0),
                MPRK43II(0.5),
                SSPMPRK22(0.5, 1.0),
                SSPMPRK43()
            ]
            for k in 2:10
                push!(algs, MPDeC(k), MPDeC(k; nodes = :lagrange))
            end
            @testset "$alg" for alg in algs
                for prod! in (prod_1!, prod_2!, prod_3!)
                    prod = (u, p, t) -> begin
                        P = similar(u, (length(u), length(u)))
                        prod!(P, u, p, t)
                        return P
                    end
                    prob_tridiagonal_ip = ConservativePDSProblem(prod!, u0, tspan;
                                                                 p_prototype = P_tridiagonal)
                    prob_tridiagonal_op = ConservativePDSProblem(prod, u0, tspan;
                                                                 p_prototype = P_tridiagonal)
                    prob_dense_ip = ConservativePDSProblem(prod!, u0, tspan;
                                                           p_prototype = P_dense)
                    prob_dense_op = ConservativePDSProblem(prod, u0, tspan;
                                                           p_prototype = P_dense)
                    prob_sparse_ip = ConservativePDSProblem(prod!, u0, tspan;
                                                            p_prototype = P_sparse)
                    prob_sparse_op = ConservativePDSProblem(prod, u0, tspan;
                                                            p_prototype = P_sparse)

                    sol_tridiagonal_ip = solve(prob_tridiagonal_ip, alg; dt,
                                               adaptive = false)
                    sol_tridiagonal_op = solve(prob_tridiagonal_op, alg; dt,
                                               adaptive = false)
                    sol_dense_ip = solve(prob_dense_ip, alg; dt, adaptive = false)
                    sol_dense_op = solve(prob_dense_op, alg; dt, adaptive = false)
                    sol_sparse_ip = solve(prob_sparse_ip, alg; dt, adaptive = false)
                    sol_sparse_op = solve(prob_sparse_op, alg; dt, adaptive = false)

                    @test sol_tridiagonal_ip.t ≈ sol_tridiagonal_op.t
                    @test sol_dense_ip.t ≈ sol_dense_op.t
                    @test sol_sparse_ip.t ≈ sol_sparse_op.t
                    @test sol_tridiagonal_ip.t ≈ sol_dense_ip.t
                    @test sol_tridiagonal_ip.t ≈ sol_sparse_ip.t

                    @test sol_tridiagonal_ip.u ≈ sol_tridiagonal_op.u
                    @test sol_dense_ip.u ≈ sol_dense_op.u
                    @test sol_sparse_ip.u ≈ sol_sparse_op.u
                    @test sol_tridiagonal_ip.u ≈ sol_dense_ip.u
                    @test sol_tridiagonal_ip.u ≈ sol_sparse_ip.u
                end
            end
        end

        @testset "Different matrix types (conservative, adaptive)" begin
            prod_1! = (P, u, p, t) -> begin
                fill!(P, zero(eltype(P)))
                for i in 1:(length(u) - 1)
                    P[i, i + 1] = i * u[i + 1]
                end
                return nothing
            end

            prod_2! = (P, u, p, t) -> begin
                fill!(P, zero(eltype(P)))
                for i in 1:(length(u) - 1)
                    P[i + 1, i] = i * u[i]
                end
                return nothing
            end

            prod_3! = (P, u, p, t) -> begin
                fill!(P, zero(eltype(P)))
                for i in 1:(length(u) - 1)
                    P[i, i + 1] = i * u[i + 1]
                    P[i + 1, i] = i * u[i]
                end
                return nothing
            end

            n = 4
            P_tridiagonal = Tridiagonal([0.1, 0.2, 0.3],
                                        zeros(n),
                                        [0.4, 0.5, 0.6])
            P_dense = Matrix(P_tridiagonal)
            P_sparse = sparse(P_tridiagonal)
            u0 = [1.0, 1.5, 2.0, 2.5]
            tspan = (0.0, 1.0)
            dt = 0.25

            rtol = sqrt(eps(Float32))

            algs = [MPRK22(0.5), MPRK22(1.0),
                MPRK43I(1.0, 0.5), MPRK43I(0.5, 0.75),
                MPRK43II(2.0 / 3.0), MPRK43II(0.5),
                SSPMPRK22(0.5, 1.0)]
            for k in 2:10
                push!(algs, MPDeC(k), MPDeC(k; nodes = :lagrange))
            end
            @testset "$alg" for alg in algs
                for prod! in (prod_1!, prod_2!, prod_3!)
                    prod = (u, p, t) -> begin
                        P = similar(u, (length(u), length(u)))
                        prod!(P, u, p, t)
                        return P
                    end
                    prob_tridiagonal_ip = ConservativePDSProblem(prod!, u0, tspan;
                                                                 p_prototype = P_tridiagonal)
                    prob_tridiagonal_op = ConservativePDSProblem(prod, u0, tspan;
                                                                 p_prototype = P_tridiagonal)
                    prob_dense_ip = ConservativePDSProblem(prod!, u0, tspan;
                                                           p_prototype = P_dense)
                    prob_dense_op = ConservativePDSProblem(prod, u0, tspan;
                                                           p_prototype = P_dense)
                    prob_sparse_ip = ConservativePDSProblem(prod!, u0, tspan;
                                                            p_prototype = P_sparse)
                    prob_sparse_op = ConservativePDSProblem(prod, u0, tspan;
                                                            p_prototype = P_sparse)

                    sol_tridiagonal_ip = solve(prob_tridiagonal_ip, alg)
                    sol_tridiagonal_op = solve(prob_tridiagonal_op, alg)
                    sol_dense_ip = solve(prob_dense_ip, alg)
                    sol_dense_op = solve(prob_dense_op, alg)
                    sol_sparse_ip = solve(prob_sparse_ip, alg)
                    sol_sparse_op = solve(prob_sparse_op, alg)

                    @test isapprox(sol_tridiagonal_ip.t, sol_tridiagonal_op.t; rtol)
                    @test isapprox(sol_dense_ip.t, sol_dense_op.t; rtol)
                    @test isapprox(sol_sparse_ip.t, sol_sparse_op.t; rtol)
                    @test isapprox(sol_tridiagonal_ip.t, sol_dense_ip.t; rtol)
                    @test isapprox(sol_tridiagonal_ip.t, sol_sparse_ip.t; rtol)

                    @test isapprox(sol_tridiagonal_ip.u, sol_tridiagonal_op.u; rtol)
                    @test isapprox(sol_dense_ip.u, sol_dense_op.u; rtol)
                    @test isapprox(sol_sparse_ip.u, sol_sparse_op.u; rtol)
                    @test isapprox(sol_tridiagonal_ip.u, sol_dense_ip.u; rtol)
                    @test isapprox(sol_tridiagonal_ip.u, sol_sparse_ip.u; rtol)
                end
            end
        end

        # Here we check that in-place and out-of-place implementations
        # deliver the same results
        @testset "Different matrix types (nonconservative)" begin
            prod_1! = (P, u, p, t) -> begin
                fill!(P, zero(eltype(P)))
                for i in 1:(length(u) - 1)
                    P[i, i + 1] = i * u[i + 1]
                end
                for i in 1:length(u)
                    P[i, i] = i * u[i]
                end
                return nothing
            end
            dest_1! = (D, u, p, t) -> begin
                fill!(D, zero(eltype(D)))
                for i in 1:length(u)
                    D[i] = (i + 1) * u[i]
                end
                return nothing
            end

            prod_2! = (P, u, p, t) -> begin
                fill!(P, zero(eltype(P)))
                for i in 1:(length(u) - 1)
                    P[i + 1, i] = i * u[i]
                end
                for i in 1:length(u)
                    P[i, i] = (i - 1) * u[i]
                end
                return nothing
            end
            dest_2! = (D, u, p, t) -> begin
                fill!(D, zero(eltype(D)))
                for i in 1:length(u)
                    D[i] = i * u[i]
                end
                return nothing
            end

            prod_3! = (P, u, p, t) -> begin
                fill!(P, zero(eltype(P)))
                for i in 1:(length(u) - 1)
                    P[i, i + 1] = i * u[i + 1]
                    P[i + 1, i] = i * u[i]
                end
                for i in 1:length(u)
                    P[i, i] = (i + 1) * u[i]
                end
                return nothing
            end
            dest_3! = (D, u, p, t) -> begin
                fill!(D, zero(eltype(D)))
                for i in 1:length(u)
                    D[i] = (i - 1) * u[i]
                end
                return nothing
            end

            n = 4
            P_tridiagonal = Tridiagonal([0.1, 0.2, 0.3],
                                        zeros(n),
                                        [0.4, 0.5, 0.6])
            P_dense = Matrix(P_tridiagonal)
            P_sparse = sparse(P_tridiagonal)
            u0 = [1.0, 1.5, 2.0, 2.5]
            D = u0
            tspan = (0.0, 1.0)
            dt = 0.25

            algs = [MPE(),
                MPRK22(0.5), MPRK22(1.0),
                MPRK43I(1.0, 0.5), MPRK43I(0.5, 0.75),
                MPRK43II(2.0 / 3.0), MPRK43II(0.5),
                SSPMPRK22(0.5, 1.0), SSPMPRK43()]
            for k in 2:10
                push!(algs, MPDeC(k), MPDeC(k; nodes = :lagrange))
            end

            @testset "$alg" for alg in algs
                for (prod!, dest!) in zip((prod_1!, prod_2!, prod_3!),
                                          (dest_1!, dest_2!, dest_3!))
                    prod = (u, p, t) -> begin
                        P = similar(u, (length(u), length(u)))
                        prod!(P, u, p, t)
                        return P
                    end
                    dest = (u, p, t) -> begin
                        D = similar(u)
                        dest!(D, u, p, t)
                        return D
                    end
                    prob_tridiagonal_ip = PDSProblem(prod!, dest!, u0, tspan;
                                                     p_prototype = P_tridiagonal)
                    prob_tridiagonal_op = PDSProblem(prod, dest, u0, tspan;
                                                     p_prototype = P_tridiagonal)
                    prob_dense_ip = PDSProblem(prod!, dest!, u0, tspan;
                                               p_prototype = P_dense)
                    prob_dense_op = PDSProblem(prod, dest, u0, tspan;
                                               p_prototype = P_dense)
                    prob_sparse_ip = PDSProblem(prod!, dest!, u0, tspan;
                                                p_prototype = P_sparse)
                    prob_sparse_op = PDSProblem(prod, dest, u0, tspan;
                                                p_prototype = P_sparse)

                    sol_tridiagonal_ip = solve(prob_tridiagonal_ip, alg;
                                               dt, adaptive = false)
                    sol_tridiagonal_op = solve(prob_tridiagonal_op, alg;
                                               dt, adaptive = false)
                    sol_dense_ip = solve(prob_dense_ip, alg;
                                         dt, adaptive = false)
                    sol_dense_op = solve(prob_dense_op, alg;
                                         dt, adaptive = false)
                    sol_sparse_ip = solve(prob_sparse_ip, alg;
                                          dt, adaptive = false)
                    sol_sparse_op = solve(prob_sparse_op, alg;
                                          dt, adaptive = false)

                    @test sol_tridiagonal_ip.t ≈ sol_tridiagonal_op.t
                    @test sol_dense_ip.t ≈ sol_dense_op.t
                    @test sol_sparse_ip.t ≈ sol_sparse_op.t
                    @test sol_tridiagonal_ip.t ≈ sol_dense_ip.t
                    @test sol_tridiagonal_ip.t ≈ sol_sparse_ip.t

                    @test sol_tridiagonal_ip.u ≈ sol_tridiagonal_op.u
                    @test sol_dense_ip.u ≈ sol_dense_op.u
                    @test sol_sparse_ip.u ≈ sol_sparse_op.u
                    @test sol_tridiagonal_ip.u ≈ sol_dense_ip.u
                    @test sol_tridiagonal_ip.u ≈ sol_sparse_ip.u
                end
            end
        end

        @testset "Different matrix types (nonconservative, adaptive)" begin
            prod_1! = (P, u, p, t) -> begin
                fill!(P, zero(eltype(P)))
                for i in 1:(length(u) - 1)
                    P[i, i + 1] = i * u[i + 1]
                end
                for i in 1:length(u)
                    P[i, i] = i * u[i]
                end
                return nothing
            end
            dest_1! = (D, u, p, t) -> begin
                fill!(D, zero(eltype(D)))
                for i in 1:length(u)
                    D[i] = (i + 1) * u[i]
                end
                return nothing
            end

            prod_2! = (P, u, p, t) -> begin
                fill!(P, zero(eltype(P)))
                for i in 1:(length(u) - 1)
                    P[i + 1, i] = i * u[i]
                end
                for i in 1:length(u)
                    P[i, i] = (i - 1) * u[i]
                end
                return nothing
            end
            dest_2! = (D, u, p, t) -> begin
                fill!(D, zero(eltype(D)))
                for i in 1:length(u)
                    D[i] = i * u[i]
                end
                return nothing
            end

            prod_3! = (P, u, p, t) -> begin
                fill!(P, zero(eltype(P)))
                for i in 1:(length(u) - 1)
                    P[i, i + 1] = i * u[i + 1]
                    P[i + 1, i] = i * u[i]
                end
                for i in 1:length(u)
                    P[i, i] = (i + 1) * u[i]
                end
                return nothing
            end
            dest_3! = (D, u, p, t) -> begin
                fill!(D, zero(eltype(D)))
                for i in 1:length(u)
                    D[i] = (i - 1) * u[i]
                end
                return nothing
            end

            n = 4
            P_tridiagonal = Tridiagonal([0.1, 0.2, 0.3],
                                        zeros(n),
                                        [0.4, 0.5, 0.6])
            P_dense = Matrix(P_tridiagonal)
            P_sparse = sparse(P_tridiagonal)
            u0 = [1.0, 1.5, 2.0, 2.5]
            D = u0
            tspan = (0.0, 1.0)
            dt = 0.25

            algs = [MPRK22(0.5), MPRK22(1.0),
                MPRK43I(1.0, 0.5), MPRK43I(0.5, 0.75),
                MPRK43II(2.0 / 3.0), MPRK43II(0.5),
                SSPMPRK22(0.5, 1.0)]
            for k in 2:10
                push!(algs, MPDeC(k), MPDeC(k; nodes = :lagrange))
            end

            rtol = sqrt(eps(Float32))
            @testset "$alg" for alg in algs
                for (prod!, dest!) in zip((prod_1!, prod_2!, prod_3!),
                                          (dest_1!, dest_2!, dest_3!))
                    prod! = prod_3!
                    dest! = dest_3!
                    prod = (u, p, t) -> begin
                        P = similar(u, (length(u), length(u)))
                        prod!(P, u, p, t)
                        return P
                    end
                    dest = (u, p, t) -> begin
                        D = similar(u)
                        dest!(D, u, p, t)
                        return D
                    end
                    prob_tridiagonal_ip = PDSProblem(prod!, dest!, u0, tspan;
                                                     p_prototype = P_tridiagonal)
                    prob_tridiagonal_op = PDSProblem(prod, dest, u0, tspan;
                                                     p_prototype = P_tridiagonal)
                    prob_dense_ip = PDSProblem(prod!, dest!, u0, tspan;
                                               p_prototype = P_dense)
                    prob_dense_op = PDSProblem(prod, dest, u0, tspan;
                                               p_prototype = P_dense)
                    prob_sparse_ip = PDSProblem(prod!, dest!, u0, tspan;
                                                p_prototype = P_sparse)
                    prob_sparse_op = PDSProblem(prod, dest, u0, tspan;
                                                p_prototype = P_sparse)

                    sol_tridiagonal_ip = solve(prob_tridiagonal_ip, alg)
                    sol_tridiagonal_op = solve(prob_tridiagonal_op, alg)
                    sol_dense_ip = solve(prob_dense_ip, alg)
                    sol_dense_op = solve(prob_dense_op, alg)
                    sol_sparse_ip = solve(prob_sparse_ip, alg)
                    sol_sparse_op = solve(prob_sparse_op, alg)

                    @test isapprox(sol_tridiagonal_ip.t, sol_tridiagonal_op.t; rtol)
                    @test isapprox(sol_dense_ip.t, sol_dense_op.t; rtol)
                    @test isapprox(sol_sparse_ip.t, sol_sparse_op.t; rtol)
                    @test isapprox(sol_tridiagonal_ip.t, sol_dense_ip.t; rtol)
                    @test isapprox(sol_tridiagonal_ip.t, sol_sparse_ip.t; rtol)

                    @test isapprox(sol_tridiagonal_ip.u, sol_tridiagonal_op.u; rtol)
                    @test isapprox(sol_dense_ip.u, sol_dense_op.u; rtol)
                    @test isapprox(sol_sparse_ip.u, sol_sparse_op.u; rtol)
                    @test isapprox(sol_tridiagonal_ip.u, sol_dense_ip.u; rtol)
                    @test isapprox(sol_tridiagonal_ip.u, sol_sparse_ip.u; rtol)
                end
            end
        end

        # Here we check that the type of p_prototype actually
        # defines the types of the Ps inside the algorithm caches.
        # We test sparse, tridiagonal, and dense matrices.
        @testset "Prototype type check" begin
            #prod and dest functions
            prod_inner! = (P, u, p, t) -> begin
                fill!(P, zero(eltype(P)))
                for i in 1:(length(u) - 1)
                    P[i, i + 1] = i * u[i + 1]
                end
                return nothing
            end
            prod_sparse! = (P, u, p, t) -> begin
                @test P isa SparseMatrixCSC
                prod_inner!(P, u, p, t)
                return nothing
            end
            prod_tridiagonal! = (P, u, p, t) -> begin
                @test P isa Tridiagonal
                prod_inner!(P, u, p, t)
                return nothing
            end
            prod_dense! = (P, u, p, t) -> begin
                @test P isa Matrix
                prod_inner!(P, u, p, t)
                return nothing
            end
            dest! = (D, u, p, t) -> begin
                fill!(D, zero(eltype(D)))
            end
            #prototypes
            P_tridiagonal = Tridiagonal([0.1, 0.2, 0.3],
                                        [0.0, 0.0, 0.0, 0.0],
                                        [0.4, 0.5, 0.6])
            P_dense = Matrix(P_tridiagonal)
            P_sparse = sparse(P_tridiagonal)
            # problem definition
            u0 = [1.0, 1.5, 2.0, 2.5]
            tspan = (0.0, 1.0)
            dt = 0.5
            ## conservative PDS
            prob_default = ConservativePDSProblem(prod_dense!, u0, tspan)
            prob_tridiagonal = ConservativePDSProblem(prod_tridiagonal!, u0, tspan;
                                                      p_prototype = P_tridiagonal)
            prob_dense = ConservativePDSProblem(prod_dense!, u0, tspan;
                                                p_prototype = P_dense)
            prob_sparse = ConservativePDSProblem(prod_sparse!, u0, tspan;
                                                 p_prototype = P_sparse)
            ## nonconservative PDS
            prob_default2 = PDSProblem(prod_dense!, dest!, u0, tspan)
            prob_tridiagonal2 = PDSProblem(prod_tridiagonal!, dest!, u0, tspan;
                                           p_prototype = P_tridiagonal)
            prob_dense2 = PDSProblem(prod_dense!, dest!, u0, tspan;
                                     p_prototype = P_dense)
            prob_sparse2 = PDSProblem(prod_sparse!, dest!, u0, tspan;
                                      p_prototype = P_sparse)

            algs = [MPE(),
                MPRK22(0.5), MPRK22(1.0),
                MPRK43I(1.0, 0.5), MPRK43I(0.5, 0.75),
                MPRK43II(2.0 / 3.0), MPRK43II(0.5),
                SSPMPRK22(0.5, 1.0), SSPMPRK43()]
            for k in 2:10
                push!(algs, MPDeC(k), MPDeC(k; nodes = :lagrange))
            end
            #solve and test
            for alg in algs
                for prob in (prob_default, prob_tridiagonal, prob_dense, prob_sparse,
                             prob_default2,
                             prob_tridiagonal2, prob_dense2, prob_sparse2)
                    sol1 = solve(prob, alg; dt, adaptive = false)

                    # test get_tmp_cache and integrator interface - modifying
                    # values from the cache should not change the final results
                    integrator = init(prob, alg; dt, adaptive = false)
                    PositiveIntegrators.OrdinaryDiffEqCore.step!(integrator)
                    cache = @inferred PositiveIntegrators.get_tmp_cache(integrator)
                    @test !isempty(cache)
                    tmp = first(cache)
                    fill!(tmp, NaN)
                    sol2 = PositiveIntegrators.solve!(integrator)
                    @test sol1.t ≈ sol2.t
                    @test sol1.u ≈ sol2.u
                end
            end
        end

        # Here we check the convergence order of pth-order schemes for which
        # also an interpolation of order p is available
        @testset "Convergence tests (conservative)" begin
            algs = (MPE(), MPRK22(0.5), MPRK22(1.0), MPRK22(2.0), SSPMPRK22(0.5, 1.0),
                    MPDeC(2), MPDeC(2, nodes = :lagrange))
            dts = 0.5 .^ (4:15)
            problems = (prob_pds_linmod, prob_pds_linmod_array,
                        prob_pds_linmod_mvector, prob_pds_linmod_inplace)

            @testset "$alg" for alg in algs
                alg = MPRK22(1.0)
                for prob in problems
                    prob = problems[1]
                    orders = experimental_orders_of_convergence(prob, alg, dts)
                    @test check_order(orders, PositiveIntegrators.alg_order(alg))

                    test_times = [
                        0.123456789, 1 / pi, exp(-1),
                        1.23456789, 1 + 1 / pi, 1 + exp(-1)
                    ]
                    for test_time in test_times
                        orders = experimental_orders_of_convergence(prob, alg,
                                                                    dts;
                                                                    test_time)
                        @test check_order(orders, PositiveIntegrators.alg_order(alg),
                                          atol = 0.2)
                        orders = experimental_orders_of_convergence(prob, alg,
                                                                    dts;
                                                                    test_time,
                                                                    only_first_index = true)
                        @test check_order(orders, PositiveIntegrators.alg_order(alg),
                                          atol = 0.2)
                    end
                end
            end
        end

        # Here we check the convergence order of pth-order schemes for which
        # also an interpolation of order p is available
        @testset "Convergence tests (nonconservative)" begin
            algs = (MPE(), MPRK22(0.5), MPRK22(1.0), MPRK22(2.0), SSPMPRK22(0.5, 1.0),
                    MPDeC(2), MPDeC(2; nodes = :lagrange))
            dts = 0.5 .^ (4:15)
            problems = (prob_pds_linmod_nonconservative,
                        prob_pds_linmod_nonconservative_inplace)
            @testset "$alg" for alg in algs
                for prob in problems
                    orders = experimental_orders_of_convergence(prob, alg, dts)
                    @test check_order(orders, PositiveIntegrators.alg_order(alg))

                    test_times = [
                        0.123456789, 1 / pi, exp(-1),
                        1.23456789, 1 + 1 / pi, 1 + exp(-1)
                    ]
                    for test_time in test_times
                        orders = experimental_orders_of_convergence(prob, alg,
                                                                    dts;
                                                                    test_time)
                        @test check_order(orders, PositiveIntegrators.alg_order(alg),
                                          atol = 0.3)
                        orders = experimental_orders_of_convergence(prob, alg,
                                                                    dts;
                                                                    test_time,
                                                                    only_first_index = true)
                        @test check_order(orders, PositiveIntegrators.alg_order(alg),
                                          N = 2, atol = 0.3)
                    end
                end
            end
        end

        # Here we check the convergence order of pth-order schemes for which
        # no interpolation of order p is available
        @testset "Convergence tests (conservative)" begin
            dts = 0.5 .^ (4:12)
            problems = (prob_pds_linmod, prob_pds_linmod_array,
                        prob_pds_linmod_mvector, prob_pds_linmod_inplace)
            algs = (MPRK43I(1.0, 0.5), MPRK43I(0.5, 0.75),
                    MPRK43II(0.5), MPRK43II(2.0 / 3.0), SSPMPRK43(), MPDeC(3),
                    MPDeC(3, nodes = :lagrange))
            for alg in algs, prob in problems
                orders = experimental_orders_of_convergence(prob, alg, dts)
                @test check_order(orders, PositiveIntegrators.alg_order(alg), atol = 0.2)
            end
        end

        # Here we check the convergence order of pth-order schemes for which
        # no interpolation of order p is available
        @testset "Convergence tests (nonconservative)" begin
            dts = 0.5 .^ (4:12)
            problems = (prob_pds_linmod_nonconservative,
                        prob_pds_linmod_nonconservative_inplace)
            algs = (MPRK43I(1.0, 0.5), MPRK43I(0.5, 0.75),
                    MPRK43II(0.5), MPRK43II(2.0 / 3.0), SSPMPRK43(),
                    MPDeC(3), MPDeC(3; nodes = :lagrange),
                    MPDeC(4), MPDeC(4; nodes = :lagrange))
            for alg in algs, prob in problems
                orders = experimental_orders_of_convergence(prob, alg, dts)
                @test check_order(orders, PositiveIntegrators.alg_order(alg), atol = 0.2)
            end
        end

        @testset "Interpolation tests (conservative)" begin
            algs = (MPE(), MPRK22(0.5), MPRK22(1.0), MPRK22(2.0), MPRK43I(1.0, 0.5),
                    MPRK43I(0.5, 0.75), MPRK43II(0.5), MPRK43II(2.0 / 3.0),
                    SSPMPRK22(0.5, 1.0), SSPMPRK43(),
                    MPDeC(2), MPDeC(2, nodes = :gausslobatto),
                    MPDeC(3), MPDeC(3, nodes = :gausslobatto),
                    MPDeC(4), MPDeC(4, nodes = :gausslobatto),
                    MPDeC(5), MPDeC(5, nodes = :gausslobatto),
                    MPDeC(6), MPDeC(6, nodes = :gausslobatto),
                    MPDeC(7), MPDeC(7, nodes = :gausslobatto))
            dt = 0.5^6
            problems = (prob_pds_linmod, prob_pds_linmod_array,
                        prob_pds_linmod_mvector, prob_pds_linmod_inplace)
            deriv = -22.0 / (5.0 * exp(3.0))
            for alg in algs
                for prob in problems
                    sol = solve(prob, alg; dt, adaptive = false)
                    # check derivative of interpolation
                    @test isapprox(sol(0.5, Val{1}), [deriv; -deriv], atol = 5e-2)
                    @test isapprox(sol(0.5, Val{1}; idxs = 1), deriv, atol = 5e-2)
                end
            end
        end

        @testset "Check convergence order (nonautonomous conservative PDS)" begin
            prod! = (P, u, p, t) -> begin
                fill!(P, zero(eltype(P)))
                P[1, 2] = sin(t)^2 * u[2]
                P[2, 1] = cos(t)^2 * u[1]
                return nothing
            end
            prod = (u, p, t) -> begin
                P = similar(u, (length(u), length(u)))
                prod!(P, u, p, t)
                return P
            end
            u0 = [1.1; 0.9] # values close to zero may decrease the order
            tspan = (0.0, 1.0)
            prob_oop = ConservativePDSProblem(prod, u0, tspan) #out-of-place
            prob_ip = ConservativePDSProblem(prod!, u0, tspan) #in-place

            dts = 0.5 .^ (4:15)
            algs = (MPE(), MPRK22(0.5), MPRK22(1.0), MPRK43I(1.0, 0.5), MPRK43I(0.5, 0.75),
                    MPRK43II(2.0 / 3.0), MPRK43II(0.5), SSPMPRK22(0.5, 1.0), SSPMPRK43(),
                    MPDeC(2), MPDeC(2, nodes = :lagrange), MPDeC(3),
                    MPDeC(3, nodes = :lagrange))
            @testset "$alg" for alg in algs
                orders = experimental_orders_of_convergence(prob_oop, alg, dts)
                @test check_order(orders, PositiveIntegrators.alg_order(alg), atol = 0.2)
                orders = experimental_orders_of_convergence(prob_ip, alg, dts)
                @test check_order(orders, PositiveIntegrators.alg_order(alg), atol = 0.2)
            end
        end

        @testset "Check convergence order (nonautonomous nonconservative PDS)" begin
            prod! = (P, u, p, t) -> begin
                fill!(P, zero(eltype(P)))
                P[1, 2] = sin(t)^2 * u[2]
                P[2, 1] = cos(2 * t)^2 * u[1]
                P[2, 2] = cos(t)^2 * u[2]
                return nothing
            end
            dest! = (d, u, p, t) -> begin
                fill!(d, zero(eltype(d)))
                d[1] = sin(2 * t)^2 * u[1]
                d[2] = sin(0.5 * t)^2 * u[2]
                return nothing
            end
            prod = (u, p, t) -> begin
                P = similar(u, (length(u), length(u)))
                prod!(P, u, p, t)
                return P
            end
            dest = (u, p, t) -> begin
                d = similar(u, (length(u),))
                dest!(d, u, p, t)
                return d
            end
            u0 = [1.1; 0.9] # stay away from zero
            tspan = (0.0, 1.0)
            prob_oop = PDSProblem(prod, dest, u0, tspan) #out-of-place
            prob_ip = PDSProblem(prod!, dest!, u0, tspan) #in-place

            dts = 0.5 .^ (4:10)
            algs = (MPE(), MPRK22(0.5), MPRK22(1.0), MPRK43I(1.0, 0.5), MPRK43I(0.5, 0.75),
                    MPRK43II(2.0 / 3.0), MPRK43II(0.5), SSPMPRK22(0.5, 1.0), SSPMPRK43(),
                    MPDeC(2), MPDeC(2; nodes = :lagrange), MPDeC(3),
                    MPDeC(3; nodes = :lagrange),
                    MPDeC(4), MPDeC(4; nodes = :lagrange), MPDeC(5),
                    MPDeC(5; nodes = :lagrange))
            @testset "$alg" for alg in algs
                orders = experimental_orders_of_convergence(prob_oop, alg, dts)
                @test check_order(orders, PositiveIntegrators.alg_order(alg), atol = 0.2)
                orders = experimental_orders_of_convergence(prob_ip, alg, dts)
                @test check_order(orders, PositiveIntegrators.alg_order(alg), atol = 0.2)
            end
        end

        @testset "Convergence of higher order schemes (autonomous)" begin
            function f_analytic(u0, p, t)
                u₁⁰, u₂⁰ = u0
                a, b = p
                c = a + b
                return ((u₁⁰ + u₂⁰) * [b; a] +
                        exp(-c * t) * (a * u₁⁰ - b * u₂⁰) * [1; -1]) / c
            end

            # linear model problem - conservative
            function linmodP!(P, u, p, t)
                P[1, 2] = u[2]
                P[2, 1] = 5u[1]
                return nothing
            end
            function linmodP(u, p, t)
                P = zeros(eltype(u), 2, 2)
                linmodP!(P, u, p, t)
                return SMatrix{2, 2}(P)
            end

            # linear model problem - nonconservative -  in-place
            function linmodP_noncons!(P, u, p, t)
                P[1, 1] = u[2] / 2
                P[1, 2] = u[2] / 2
                P[2, 1] = 3u[1]
                P[2, 2] = 2u[1]
                return nothing
            end
            function linmodD_noncons!(D, u, p, t)
                D[1] = 2u[1]
                D[2] = u[2] / 2
                return nothing
            end
            function linmodP_noncons(u, p, t)
                P = zeros(eltype(u), 2, 2)
                linmodP_noncons!(P, u, p, t)
                return SMatrix{2, 2}(P)
            end
            function linmodD_noncons(u, p, t)
                D = zeros(eltype(u), 2, 1)
                linmodD_noncons!(D, u, p, t)
                return SVector{2}(D)
            end

            u0 = [Double64(9) / 10; Double64(1) / 10]
            p = [Double64(5); Double64(1)]
            tspan = (Double64(0), Double64(1))

            prob_ip = ConservativePDSProblem(linmodP!, u0, tspan, p; analytic = f_analytic)
            prob_op = ConservativePDSProblem(linmodP, SVector{2}(u0), tspan, SVector{2}(p);
                                             analytic = f_analytic)
            prob_ip_noncons = PDSProblem(linmodP_noncons!, linmodD_noncons!, u0, tspan, p;
                                         analytic = f_analytic)
            prob_op_noncons = PDSProblem(linmodP_noncons, linmodD_noncons, SVector{2}(u0),
                                         tspan,
                                         SVector{2}(p); analytic = f_analytic)

            dts = 0.5 .^ (8:13)

            algs = MPDeC[]
            for K in 4:10
                push!(algs, MPDeC(K))
                push!(algs, MPDeC(K, nodes = :lagrange))
            end

            atol = 0.3
            @testset "$alg" for alg in algs
                orders = experimental_orders_of_convergence(prob_ip, alg, dts)
                @test check_order(orders, PositiveIntegrators.alg_order(alg); atol)

                orders = experimental_orders_of_convergence(prob_op, alg, dts)
                @test check_order(orders, PositiveIntegrators.alg_order(alg); atol)

                orders = experimental_orders_of_convergence(prob_ip_noncons, alg, dts)
                @test check_order(orders, PositiveIntegrators.alg_order(alg); atol)

                orders = experimental_orders_of_convergence(prob_op_noncons, alg, dts)
                @test check_order(orders, PositiveIntegrators.alg_order(alg); atol)
            end
        end

        @testset "Convergence of higher order schemes (nonautonomous)" begin
            prod! = (P, u, p, t) -> begin
                fill!(P, zero(eltype(P)))
                P[1, 2] = sin(t)^2 * u[2]
                P[2, 1] = cos(t)^2 * u[1]
                return nothing
            end
            prod = (u, p, t) -> begin
                P = similar(u, (length(u), length(u)))
                prod!(P, u, p, t)
                return SMatrix{2, 2}(P)
            end
            function fct!(du, u, p, t)
                du[1] = sin(t)^2 * u[2] - cos(t)^2 * u[1]
                du[2] = -sin(t)^2 * u[2] + cos(t)^2 * u[1]
                return nothing
            end
            fct = (u, p, t) -> begin
                f = similar(u)
                fct!(f, u, p, t)
                return SVector{2}(f)
            end

            prod_noncons! = (P, u, p, t) -> begin
                fill!(P, zero(eltype(P)))
                P[1, 2] = sin(t)^2 * u[2]
                P[2, 1] = cos(2 * t)^2 * u[1]
                P[2, 2] = cos(t)^2 * u[2]
                return nothing
            end
            dest_noncons! = (d, u, p, t) -> begin
                fill!(d, zero(eltype(d)))
                d[1] = sin(2 * t)^2 * u[1]
                d[2] = sin(0.5 * t)^2 * u[2]
                return nothing
            end
            prod_noncons = (u, p, t) -> begin
                P = similar(u, (length(u), length(u)))
                prod_noncons!(P, u, p, t)
                return SMatrix{2, 2}(P)
            end
            dest_noncons = (u, p, t) -> begin
                d = similar(u, (length(u),))
                dest_noncons!(d, u, p, t)
                return SVector{2}(d)
            end
            fct_noncons! = (du, u, p, t) -> begin
                du[1] = sin(t)^2 * u[2] - cos(2 * t)^2 * u[1] - sin(2 * t)^2 * u[1]
                du[2] = -sin(t)^2 * u[2] + cos(2 * t)^2 * u[1] + cos(t)^2 * u[2] -
                        sin(0.5 * t)^2 * u[2]
                return nothing
            end
            fct_noncons = (u, p, t) -> begin
                f = similar(u)
                fct_noncons!(f, u, p, t)
                return SVector{2}(f)
            end

            u0 = [Double64(11) / 10; Double64(9) / 10] # values close to zero may decrease the order
            tspan = (Double64(0), Double64(1))
            prob_op = ConservativePDSProblem(prod, SVector{2}(u0), tspan; std_rhs = fct) #out-of-place
            prob_ip = ConservativePDSProblem(prod!, u0, tspan; std_rhs = fct!) #in-place
            prob_op_noncons = PDSProblem(prod_noncons, dest_noncons, SVector{2}(u0), tspan;
                                         std_rhs = fct_noncons) #out-of-place
            prob_ip_noncons = PDSProblem(prod_noncons!, dest_noncons!, u0, tspan;
                                         std_rhs = fct_noncons!) #in-place

            dts = 0.5 .^ (5:10)

            algs = []
            for K in 4:10
                push!(algs, MPDeC(K))
                push!(algs, MPDeC(K, nodes = :lagrange))
            end

            atol = 0.3
            @testset "$alg" for alg in algs
                orders = experimental_orders_of_convergence(prob_ip, alg, dts)
                @test check_order(orders, PositiveIntegrators.alg_order(alg); atol)

                orders = experimental_orders_of_convergence(prob_op, alg, dts)
                @test check_order(orders, PositiveIntegrators.alg_order(alg); atol)

                orders = experimental_orders_of_convergence(prob_ip_noncons, alg, dts)
                @test check_order(orders, PositiveIntegrators.alg_order(alg); atol)

                orders = experimental_orders_of_convergence(prob_op_noncons, alg, dts)
                @test check_order(orders, PositiveIntegrators.alg_order(alg); atol)
            end
        end

        @testset "Interpolation tests (nonconservative)" begin
            algs = (MPE(), MPRK22(0.5), MPRK22(1.0), MPRK22(2.0), MPRK43I(1.0, 0.5),
                    MPRK43I(0.5, 0.75), MPRK43II(0.5), MPRK43II(2.0 / 3.0),
                    SSPMPRK22(0.5, 1.0), SSPMPRK43())
            dt = 0.5^6
            problems = (prob_pds_linmod_nonconservative,
                        prob_pds_linmod_nonconservative_inplace)
            for alg in algs
                for prob in problems
                    sol = solve(prob, alg; dt, adaptive = false)
                    # check derivative of interpolation
                    @test_nowarn sol(0.5, Val{1})
                    @test_nowarn sol(0.5, Val{1}; idxs = 1)
                end
            end
        end

        # Check that the schemes accept zero initial values
        @testset "Zero initial values" begin
            # Do a single step and check that no NaNs occur
            u0 = [1.0, 0.0]
            dt = 1.0
            tspan = (0.0, dt)
            p = 1000.0
            function prod!(P, u, p, t)
                λ = p
                fill!(P, zero(eltype(P)))
                P[2, 1] = λ * u[1]
            end
            function dest!(D, u, p, t)
                fill!(D, zero(eltype(D)))
            end
            function prod(u, p, t)
                P = similar(u, (length(u), length(u)))
                prod!(P, u, p, t)
                return P
            end
            function dest(u, p, t)
                d = similar(u)
                dest!(d, u, p, t)
                return d
            end

            prob_ip = ConservativePDSProblem(prod!, u0, tspan, p)
            prob_ip_2 = PDSProblem(prod!, dest!, u0, tspan, p)
            prob_oop = ConservativePDSProblem(prod, u0, tspan, p)
            prob_oop_2 = PDSProblem(prod, dest, u0, tspan, p)

            algs = [MPE(), MPRK22(0.5), MPRK22(1.0), MPRK22(2.0),
                MPRK43I(1.0, 0.5), MPRK43I(0.5, 0.75), MPRK43II(0.5),
                MPRK43II(2.0 / 3.0), SSPMPRK22(0.5, 1.0), SSPMPRK43()]
            for k in 2:10
                push!(algs, MPDeC(k), MPDeC(k; nodes = :lagrange))
            end

            for alg in algs
                sol = solve(prob_ip, alg; dt = dt, adaptive = false)
                @test !any(isnan, sol.u[end])
                sol = solve(prob_ip_2, alg; dt = dt, adaptive = false)
                @test !any(isnan, sol.u[end])
                sol = solve(prob_oop, alg; dt = dt, adaptive = false)
                @test !any(isnan, sol.u[end])
                sol = solve(prob_oop_2, alg; dt = dt, adaptive = false)
                @test !any(isnan, sol.u[end])
            end
        end

        # Check that approximations, and thus the Patankar weights,
        # remain positive to avoid division by zero.
        @testset "Positivity check" begin
            # For this problem u[1] decreases montonically to 0 very fast.
            # We perform 10^5 steps and check that u[end] does not contain any NaNs
            u0 = [0.9, 0.1]
            tspan = (0.0, 100.0)
            p = 1000.0
            function prod!(P, u, p, t)
                λ = p
                fill!(P, zero(eltype(P)))
                P[2, 1] = λ * u[1]
            end
            function dest!(D, u, p, t)
                fill!(D, zero(eltype(D)))
            end
            function prod(u, p, t)
                P = similar(u, (length(u), length(u)))
                prod!(P, u, p, t)
                return P
            end
            function dest(u, p, t)
                d = similar(u)
                dest!(d, u, p, t)
                return d
            end

            prob_ip = ConservativePDSProblem(prod!, u0, tspan, p)
            prob_ip_2 = PDSProblem(prod!, dest!, u0, tspan, p)
            prob_oop = ConservativePDSProblem(prod, u0, tspan, p)
            prob_oop_2 = PDSProblem(prod, dest, u0, tspan, p)

            algs = [MPE(), MPRK22(0.5), MPRK22(1.0), MPRK22(2.0),
                MPRK43I(1.0, 0.5), MPRK43I(0.5, 0.75), MPRK43II(0.5),
                MPRK43II(2.0 / 3.0), SSPMPRK22(0.5, 1.0), SSPMPRK43()]
            for k in 2:10
                push!(algs, MPDeC(k), MPDeC(k; nodes = :lagrange))
            end

            dt = 1e-3
            for alg in algs
                sol1 = solve(prob_ip, alg; dt = dt, adaptive = false)
                @test !any(isnan, sol1.u[end])
                sol2 = solve(prob_ip_2, alg; dt = dt, adaptive = false)
                @test !any(isnan, sol2.u[end])
                sol3 = solve(prob_oop, alg; dt = dt, adaptive = false)
                @test !any(isnan, sol3.u[end])
                sol4 = solve(prob_oop_2, alg; dt = dt, adaptive = false)
                @test !any(isnan, sol4.u[end])
                @test sol1.u ≈ sol2.u ≈ sol3.u ≈ sol4.u
            end
        end

        # Here we check that the implemented schemes can solve the predefined PDS
        # (at least for specific parameters)
        @testset "PDS problem library (adaptive schemes)" begin
            algs = [MPRK22(0.5), MPRK22(1.0), MPRK43I(1.0, 0.5), MPRK43I(0.5, 0.75),
                MPRK43II(2.0 / 3.0), MPRK43II(0.5), SSPMPRK22(0.5, 1.0)]
            for k in 2:10
                push!(algs, MPDeC(k), MPDeC(k; nodes = :lagrange))
            end
            probs = (prob_pds_linmod, prob_pds_linmod_inplace, prob_pds_nonlinmod,
                     prob_pds_robertson, prob_pds_bertolazzi, prob_pds_brusselator,
                     prob_pds_npzd,
                     prob_pds_sir, prob_pds_stratreac)
            @testset "$alg" for alg in algs
                @testset "$i" for (i, prob) in enumerate(probs)
                    if prob == prob_pds_stratreac && alg == SSPMPRK22(0.5, 1.0)
                        #TODO: SSPMPRK22(0.5, 1.0) is unstable for prob_pds_stratreac.
                        #Need to figure out if this is a problem of the algorithm or not.
                        break
                    elseif prob == prob_pds_stratreac && alg == MPRK43I(0.5, 0.75)
                        # Not successful on Julia 1.9
                        break
                    elseif prob == prob_pds_stratreac && alg == MPDeC(9; nodes = :lagrange)
                        # unstable
                        break
                    end
                    # later versions of OrdinaryDiffEq.jl use dtmin = 0 by default,
                    # see https://github.com/SciML/OrdinaryDiffEq.jl/pull/2098
                    sol = solve(prob, alg; dtmin = 0.0)
                    @test Int(sol.retcode) == 1
                end
            end
        end

        # Here we check that the implemented schemes can solve the predefined PDS.
        @testset "PDS problem library (non-adaptive schemes)" begin
            algs = (MPE(), SSPMPRK43())
            #prob_pds_robertson not included
            probs = (prob_pds_linmod, prob_pds_linmod_inplace, prob_pds_nonlinmod,
                     prob_pds_bertolazzi, prob_pds_brusselator,
                     prob_pds_npzd, prob_pds_sir, prob_pds_stratreac)
            @testset "$alg" for alg in algs
                @testset "$prob" for prob in probs
                    tspan = prob.tspan
                    dt = (tspan[2] - tspan[1]) / 10
                    sol = solve(prob, alg; dt = dt)
                    @test Int(sol.retcode) == 1
                end
            end
        end

        #Here we check the different possibilities to define small_constant.
        @testset "Different possibilities to set small_constant" begin
            # For this problem u[1] decreases montonically to 0 very fast.
            u0 = [0.9, 0.1]
            tspan = (0.0, 100.0)
            p = 1000.0
            function prod!(P, u, p, t)
                λ = p
                fill!(P, zero(eltype(P)))
                P[2, 1] = λ * u[1]
            end
            function dest!(D, u, p, t)
                fill!(D, zero(eltype(D)))
            end
            function prod(u, p, t)
                P = similar(u, (length(u), length(u)))
                prod!(P, u, p, t)
                return P
            end
            function dest(u, p, t)
                d = similar(u)
                dest!(d, u, p, t)
                return d
            end
            prob_ip = ConservativePDSProblem(prod!, u0, tspan, p)
            prob_ip_2 = PDSProblem(prod!, dest!, u0, tspan, p)
            prob_oop = ConservativePDSProblem(prod, u0, tspan, p)
            prob_oop_2 = PDSProblem(prod, dest, u0, tspan, p)

            probs = (prob_ip, prob_ip_2, prob_oop, prob_oop_2,
                     prob_pds_linmod, prob_pds_linmod_inplace, prob_pds_nonlinmod,
                     prob_pds_robertson, prob_pds_bertolazzi, prob_pds_brusselator,
                     prob_pds_npzd,
                     prob_pds_sir, prob_pds_stratreac)

            algs = [MPE, (; kwargs...) -> MPRK22(1.0; kwargs...),
                (; kwargs...) -> MPRK43I(1.0, 0.5; kwargs...),
                (; kwargs...) -> MPRK43II(0.5; kwargs...),
                (; kwargs...) -> SSPMPRK22(0.5, 1.0; kwargs...),
                (; kwargs...) -> SSPMPRK43(; kwargs...)]
            for k in 2:10
                push!(algs, (; kwargs...) -> MPDeC(k; kwargs...),
                      (; kwargs...) -> MPDeC(k; nodes = :lagrange, kwargs...))
            end
            for alg in algs
                for prob in probs
                    sol1 = solve(prob_pds_linmod, alg(), dt = 0.1)
                    sol2 = solve(prob_pds_linmod, alg(small_constant = floatmin(Float64)),
                                 dt = 0.1)
                    sol3 = solve(prob_pds_linmod, alg(small_constant = floatmin), dt = 0.1)
                    @test sol1.t ≈ sol2.t ≈ sol3.t
                    @test sol1.u ≈ sol2.u ≈ sol3.u
                end
            end
        end

        # Here we check if the RK methods on which the MPRK schemes are based integrate
        # u'(t) = q * t^(q-1) exactly for q from 1 to the order of the method.
        # This is also true for MPDeC as long as the theta matrix is nonnegative, i.e. K = 2.
        # Nevertheless, the results of most MPDeC schemes are good enough to pass this test
        @testset "Exact solutions (RK)" begin
            algs = [MPE(), MPRK22(0.5), MPRK22(1.0), MPRK22(2.0),
                MPRK43I(1.0, 0.5), MPRK43I(0.5, 0.75), MPRK43II(0.5),
                MPRK43II(2.0 / 3.0), SSPMPRK22(0.5, 1.0), SSPMPRK43()]
            for k in 2:10
                push!(algs, MPDeC(k))
                if k != 9
                    push!(algs, MPDeC(k, nodes = :lagrange))
                end
            end
            @testset "$alg, $q" for alg in algs, q in 1:PositiveIntegrators.alg_order(alg)
                f(t) = q * t^(q - 1)
                function prod!(P, u, p, t)
                    fill!(P, zero(eltype(P)))
                    P[1, 1] = f(t)
                end
                function dest!(D, u, p, t)
                    fill!(D, zero(eltype(D)))
                end
                function prod(u, p, t)
                    P = similar(u, (length(u), length(u)))
                    prod!(P, u, p, t)
                    return P
                end
                function dest(u, p, t)
                    d = similar(u)
                    dest!(d, u, p, t)
                    return d
                end
                u0 = [0.0; 2.0]
                prob_oop = PDSProblem(prod, dest, u0, (0.0, 1.0))
                sol_oop = solve(prob_oop, alg, dt = 0.1; adaptive = false)
                @test first(last(sol_oop.u)) ≈ 1.0
                prob_ip = PDSProblem(prod!, dest!, u0, (0.0, 1.0))
                sol_ip = solve(prob_ip, alg, dt = 0.1; adaptive = false)
                @test first(last(sol_ip.u)) ≈ 1.0
            end
        end

        # Here we check that MPE integrates u'(t) = -(a0 * b1)/(b0 + b1  * t)^2 exactly.
        # The solution is u(t) = a0/(b0 + b1 * t)
        @testset "Exact solutions (MPE)" begin
            alg = MPE()
            u_exact(t) = 1 / (2 + 3 * t)
            f(t) = -3 / (2 + 3 * t)^2
            function prod!(P, u, p, t)
                fill!(P, zero(eltype(P)))
            end
            function dest!(D, u, p, t)
                fill!(D, zero(eltype(D)))
                D[1, 1] = -f(t)
            end
            function prod(u, p, t)
                P = similar(u, (length(u), length(u)))
                prod!(P, u, p, t)
                return P
            end
            function dest(u, p, t)
                d = similar(u)
                dest!(d, u, p, t)
                return d
            end
            u0 = [0.5; 0.0]
            prob_oop = PDSProblem(prod, dest, u0, (0.0, 1.0))
            sol_oop = solve(prob_oop, alg, dt = 0.5; adaptive = false)
            @test first(last(sol_oop.u)) ≈ u_exact(last(prob_oop.tspan))
            prob_ip = PDSProblem(prod!, dest!, u0, (0.0, 1.0))
            sol_ip = solve(prob_ip, alg, dt = 0.1; adaptive = false)
            @test first(last(sol_ip.u)) ≈ u_exact(last(prob_ip.tspan))
        end
    end

    @testset "plot" begin
        using Plots

        # plot ode solution
        sol = solve(prob_pds_linmod, MPRK22(1.0))
        @test_nowarn plot(sol)

        # plot work-precision diagram
        prob = prob_pds_nonlinmod
        algs = [MPE(); MPRK22(1.0)]
        labels = ["MPE"; "MPRK"]
        dts = [0.1; 0.01; 0.001]
        alg_ref = Vern7()
        wp = work_precision_fixed(prob, algs, labels, dts, alg_ref)
        @test_nowarn plot(wp, labels; minorticks = 10)
    end

    @testset "utilities" begin
        @testset "is(non)negative" begin
            sol_Euler = solve(prob_pds_linmod, Euler(), dt = 0.25)
            @test isnegative(sol_Euler)
            @test isnegative(sol_Euler.u)
            @test isnegative(sol_Euler.u[2])

            sol_MPE = solve(prob_pds_linmod, MPE(), dt = 0.25)
            @test isnonnegative(sol_MPE)
            @test isnonnegative(sol_MPE.u)

            sol_bruss = solve(prob_pds_brusselator, Euler(), dt = 1.0)
            @test isnegative(sol_bruss)
            @test isnegative(sol_bruss.u)
            @test isnegative(last(sol_bruss.u))
        end

        @testset "errors" begin
            v = [[1.0; 2.0; 3.0], [4.0; 5.0; 6.0]]
            w = [[8.0; 10.0; 12.0], [2.0; 4.0; 6.0]]

            @test rel_max_error_tend(v, w) == 1.0
            @test rel_l1_error_tend(v, w) == 1.25 / 3
            @test rel_l2_error_tend(v, w) == sqrt(1.0625 / 3)
            @test rel_max_error_overall(w, v) == 7.0
        end

        # Here we run a single scheme multiple times and check
        # that the errors are identical and that the computing times
        # differ only slightly.
        @testset "work-precision fixed" begin
            prob = prob_pds_nonlinmod
            alg = MPRK22(1.0)
            algs = [alg; alg; alg; alg; alg]
            labelsA = ["A_1"; "A_2"; "A_3"; "A_4"; "A_5"]
            labelsB = ["B_1"; "B_2"; "B_3"; "B_4"; "B_5"]
            dts = (last(prob.tspan) - first(prob.tspan)) / 10.0 * 0.5 .^ (4:13)
            alg_ref = Vern7()
            wp = work_precision_fixed(prob, algs, labelsA, dts,
                                      alg_ref)
            work_precision_fixed!(wp, prob, algs, labelsB, dts,
                                  alg_ref)

            # check that errors agree
            for (i, _) in enumerate(dts)
                v = [value[i][1] for (key, value) in wp]
                @test all(y -> y == v[1], v)
            end

            # check that computing times are close enough
            for (i, _) in enumerate(dts)
                v = [value[i][2] for (key, value) in wp]
                m1 = mean(v)
                # This test allows computing times that are
                # 2.5 times the mean value. In a loglog plot these
                # differences won't be significant.
                allowed = 1.5
                if Sys.isapple()
                    # On macOS, the time differences are larger
                    # than on other platforms in CI sometimes.
                    allowed = 10.0
                end
                @test maximum((v .- m1) ./ m1) < allowed
            end
        end

        # Here we run a single scheme multiple times and check
        # that the errors are identical and that the computing times
        # differ only slightly.
        @testset "work-precision adaptive" begin
            @testset "adaptive_ref = false" begin
                prob = prob_pds_nonlinmod
                alg = MPRK22(1.0)
                algs = [alg; alg; alg; alg; alg]
                labelsA = ["A_1"; "A_2"; "A_3"; "A_4"; "A_5"]
                labelsB = ["B_1"; "B_2"; "B_3"; "B_4"; "B_5"]
                abstols = 1 ./ 10 .^ (4:8)
                reltols = 1 ./ 10 .^ (3:7)
                alg_ref = Vern7()
                wp = work_precision_adaptive(prob, algs, labelsA,
                                             abstols, reltols, alg_ref)
                work_precision_adaptive!(wp, prob, algs, labelsB,
                                         abstols, reltols, alg_ref)

                # check that errors agree
                for (i, _) in enumerate(abstols)
                    v = [value[i][1] for (key, value) in wp]
                    @test all(y -> y == v[1], v)
                end

                # check that computing times are close enough
                for (i, _) in enumerate(abstols)
                    v = [value[i][2] for (key, value) in wp]
                    m1 = mean(v)
                    # This test allows computing times that are
                    # 2.5 times the mean value. In a loglog plot these
                    # differences won't be significant.
                    allowed = 1.5
                    if Sys.isapple()
                        # On macOS, the time differences are larger
                        # than on other platforms in CI sometimes.
                        allowed = 10.0
                    end
                    @test maximum((v .- m1) ./ m1) < allowed
                end
            end

            @testset "adaptive_ref = true" begin
                prob = prob_pds_robertson
                alg = MPRK22(1.0)
                algs = [alg; alg; alg; alg; alg]
                labelsA = ["A_1"; "A_2"; "A_3"; "A_4"; "A_5"]
                labelsB = ["B_1"; "B_2"; "B_3"; "B_4"; "B_5"]
                abstols = 1 ./ 10 .^ (4:8)
                reltols = 1 ./ 10 .^ (3:7)
                alg_ref = Rodas4P()
                wp = work_precision_adaptive(prob, algs, labelsA,
                                             abstols, reltols, alg_ref; adaptive_ref = true)
                work_precision_adaptive!(wp, prob, algs, labelsB,
                                         abstols, reltols, alg_ref; adaptive_ref = true)

                # check that errors agree
                for (i, _) in enumerate(abstols)
                    v = [value[i][1] for (key, value) in wp]
                    @test all(y -> y == v[1], v)
                end

                # check that computing times are close enough
                for (i, _) in enumerate(abstols)
                    v = [value[i][2] for (key, value) in wp]
                    m1 = mean(v)
                    # This test allows computing times that are
                    # 2.5 times the mean value. In a loglog plot these
                    # differences won't be significant.
                    allowed = 1.5
                    if Sys.isapple()
                        # On macOS, the time differences are larger
                        # than on other platforms in CI sometimes.
                        allowed = 10.0
                    end
                    @test maximum((v .- m1) ./ m1) < allowed
                end
            end
        end
    end
    @testset "linear invariant matrices" begin
        @testset "Check that the structure of the linear invariant field is well implemented" begin
            u0 = [0.9, 0.1]
            tspan = (0.0, 2.0)
            AT = @SMatrix[1 1]

            linmodP(u, p, t) = [0 u[2]; 5*u[1] 0]
            linmodD(u, p, t) = [0.0; 0.0]
            linmod_PDS_op = PDSProblem(linmodP, linmodD, u0, tspan, linear_invariants = AT)

            linmod_PDS_op_cons = ConservativePDSProblem(linmodP, linmodD, u0, tspan,
                                                        linear_invariants = AT)

            @test linmod_PDS_op.f.linear_invariants == AT
            @test linmod_PDS_op_cons.f.linear_invariants == AT
        end

        @testset "Check that the predefined linear invariant matrices fit in predefined problems" begin
            # non-stiff conservative problems (out-of-place)
            probs = (prob_pds_linmod, prob_pds_nonlinmod, prob_pds_brusselator,
                     prob_pds_sir, prob_pds_npzd)
            @testset "$prob" for prob in probs
                dt = (last(prob.tspan) - first(prob.tspan)) / 1e4
                sol = solve(prob, Tsit5(); dt, isoutofdomain = isnegative)
                @test prob.f.linear_invariants * prob.u0 ≈
                      prob.f.linear_invariants * sol.u[end]
            end

            # non-stiff conservative problems (in-place)
            probs = (prob_pds_linmod_inplace,)
            @testset "$prob" for prob in probs
                dt = (last(prob.tspan) - first(prob.tspan)) / 1e4
                sol = solve(prob, Tsit5(); dt, isoutofdomain = isnegative)
                @test prob.f.linear_invariants * prob.u0 ≈
                      prob.f.linear_invariants * sol.u[end]
            end

            # non-stiff non-conservative problems (out-of-place)
            probs = (prob_pds_minmapk,)
            @testset "$prob" for prob in probs
                dt = (last(prob.tspan) - first(prob.tspan)) / 1e4
                sol = solve(prob, Tsit5(); dt, isoutofdomain = isnegative)
                @test prob.f.linear_invariants * prob.u0 ≈
                      prob.f.linear_invariants * sol.u[end]
            end

            # Robertson problem
            prob = prob_pds_robertson
            dt = 1e-6
            sol = solve(prob, Rosenbrock23(); dt)
            @test prob.f.linear_invariants * prob.u0 ≈ prob.f.linear_invariants * sol.u[end]

            # Bertolazzi problem
            prob = prob_pds_bertolazzi
            alg = ImplicitEuler()
            dt = 1e-3
            sol = solve(prob, alg; dt)
            @test prob.f.linear_invariants * prob.u0 ≈ prob.f.linear_invariants * sol.u[end]

            # Stratospheric reaction problem
            prob = prob_pds_stratreac
            dt = 1.0
            sol = solve(prob, Rosenbrock23(); dt)
            @test prob.f.linear_invariants * prob.u0 ≈ prob.f.linear_invariants * sol.u[end]
        end
    end
end;
