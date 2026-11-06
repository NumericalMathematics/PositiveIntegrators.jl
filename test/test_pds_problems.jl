@testitem "ODE RHS" setup=[Imports, TestHelpers] begin
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

@testitem "ConservativePDSFunction" setup=[Imports, TestHelpers] begin
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

@testitem "Linear model problem" setup=[Imports, TestHelpers] begin
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

    # in-place syntax for PDS
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

@testitem "PDSProblem error handling" setup=[Imports, TestHelpers] begin
    P(u, p, t) = 0.0
    D(du, u, p, t) = 0.0
    @test_throws "in-place and out-of-place" PDSProblem(P, D, 0.0, (0.0, 1.0))
end

@testitem "Check that the structure of the linear invariant field is well implemented" setup=[
    Imports,
    TestHelpers
] begin
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

@testitem "Check that the predefined linear invariant matrices fit in predefined problems" setup=[
    Imports,
    TestHelpers
] begin
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

    # HIRES problem
    prob = prob_pds_hires
    dt = 1e-3
    sol = solve(prob, Rosenbrock23(); dt)
    @test prob.f.linear_invariants * prob.u0 ≈ prob.f.linear_invariants * sol.u[end]
end
