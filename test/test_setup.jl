# Packages and names used by the test items. `Test` and `PositiveIntegrators`
# are loaded in every test item by default.
@testsnippet Imports begin
    using LinearAlgebra
    using SparseArrays
    using Statistics: mean, median

    using DoubleFloats: Double64
    using StaticArrays: @SMatrix, SMatrix, MVector, @SVector, SVector, SA

    using Unitful: @u_str, ustrip

    using ADTypes
    using OrdinaryDiffEqLowOrderRK: Euler
    using OrdinaryDiffEqRosenbrock: Rosenbrock23, Rodas4P, ROS2
    using OrdinaryDiffEqSDIRK: ImplicitEuler, SDIRK2, TRBDF2
    using OrdinaryDiffEqTsit5: Tsit5
    using OrdinaryDiffEqVerner: Vern7, Vern9
    import SciMLBase

    using Clarabel
    using JuMP

    # load RecursiveFactorization to get RFLUFactorization
    using RecursiveFactorization: RecursiveFactorization
    using LinearSolve: RFLUFactorization, LUFactorization, KLUFactorization, KrylovJL_GMRES
end

# Helper functions and test problems shared by several test items. They are defined
# in a module, which is evaluated only once, so that the solvers do not have to be
# compiled again for these problems in every test item.
@testmodule TestHelpers begin
    using LinearAlgebra: norm
    using SparseArrays: SparseMatrixCSC, nonzeros, nzrange
    using StaticArrays: MVector
    using OrdinaryDiffEqVerner: Vern7
    using PositiveIntegrators

    export experimental_orders_of_convergence, check_order, check_order_leastsquares
    export prob_pds_linmod_array, prob_pds_linmod_mvector,
           prob_pds_linmod_nonconservative, prob_pds_linmod_nonconservative_inplace
    export linear_advection_fd_upwind_f!, linear_advection_fd_upwind_P!,
           linear_advection_fd_upwind_D!

    """
        experimental_orders_of_convergence(prob, alg, dts;
                                          test_time = nothing,
                                          only_first_index = false,
                                          ref_alg = Vern7(),
                                          use_least_squares = false)

    Solve `prob` with `alg` and fixed time steps taken from `dts`, and compute
    the errors at `test_time`. If`test_time` is not specified the error is computed
    at the final time.
    Return the associated experimental orders of convergence.

    If `only_first_index == true`, only the first solution component is used
    to compute the error. If no analytic solution `prob.f.analytic` is available, a reference
    solution is computed using `ref_alg`.

    If `use_least_squares = true`, the experimental order of convergence is computed by a least-squares fit 
    of the model `log(error) = c + order * log(dt)` to the data `(log(dts), log(errors))`. 
    In this case, the function returns the fitted order instead of the pointwise orders. 
    """
    function experimental_orders_of_convergence(prob, alg, dts; test_time = nothing,
                                                only_first_index = false,
                                                ref_alg = Vern7(),
                                                use_least_squares = false)
        @assert length(dts) > 1
        errors = zeros(eltype(dts), length(dts))

        if !(isnothing(prob.f.analytic))
            # There is an analytic solution
            if isnothing(test_time)
                # Compare results at final time
                reference_solution = prob.f.analytic(prob.u0, prob.p, last(prob.tspan))
            else
                # Compare results at given test_time
                reference_solution = prob.f.analytic(prob.u0, prob.p, test_time)
            end
        else
            # Compute a reference solution numerically
            if isnothing(test_time)
                # Compare results at final time
                tspan = prob.tspan
            else
                # Compare results at given test_time
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

        # --- Mode Selection: Least-Squares vs Pointwise Orders ---
        if use_least_squares
            # Fit model: log(error) = c + order * log(dt)
            X = hcat(ones(eltype(dts), length(dts)), log.(dts))
            fit = X \ log.(errors)
            fitted_order = fit[2]
            return fitted_order
        else
            orders = experimental_orders_of_convergence(errors, dts)
            return orders
        end
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

    """
        check_order_leastsquares(fitted_order, expected_order; atol = 0.2)

    Returns `true` if `fitted_order` is at least `expected_order - atol`.
    """
    function check_order_leastsquares(fitted_order, expected_order; atol = 0.2)
        return expected_order - atol <= fitted_order <= expected_order + 2.5 * atol
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
    const prob_pds_linmod_nonconservative = PDSProblem(linmodP, linmodD, [0.9, 0.1],
                                                       (0.0, 2.0),
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
    const prob_pds_linmod_nonconservative_inplace = PDSProblem(linmodP!, linmodD!,
                                                               [0.9, 0.1],
                                                               (0.0, 2.0), [5.0, 1.0];
                                                               analytic = f_analytic)

    # Function definitions for the test item
    # "OrdinaryDiffEq compatibility: linear advection".
    # Functions are defined outside the test item to avoid failing allocation tests,
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
    # in-place syntax for PDS
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
end
