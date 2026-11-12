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

@testitem "PDSProblem: linear model problem" setup=[Imports, TestHelpers] begin
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

@testitem "remake" setup=[Imports, TestHelpers] begin
    P_op(u, p, t) = [0.0 u[2]; p*u[1] 0.0]
    D_op(u, p, t) = zero(u)
    function P_ip!(P, u, p, t)
        P .= 0
        P[1, 2] = u[2]
        P[2, 1] = p * u[1]
        return nothing
    end
    D_ip!(D, u, p, t) = fill!(D, 0)

    u0 = [0.9, 0.1]
    tspan = (0.0, 1.0)
    u0_new = [0.5, 0.5]
    tspan_new = (0.0, 2.0)
    probs = (ConservativePDSProblem(P_op, u0, tspan, 5.0),
             ConservativePDSProblem(P_ip!, u0, tspan, 5.0),
             PDSProblem(P_op, D_op, u0, tspan, 5.0),
             PDSProblem(P_ip!, D_ip!, u0, tspan, 5.0))
    probs_new = (ConservativePDSProblem(P_op, u0_new, tspan_new, 2.0),
                 ConservativePDSProblem(P_ip!, u0_new, tspan_new, 2.0),
                 PDSProblem(P_op, D_op, u0_new, tspan_new, 2.0),
                 PDSProblem(P_ip!, D_ip!, u0_new, tspan_new, 2.0))
    @testset "$i" for (i, (prob, prob_new)) in enumerate(zip(probs, probs_new))
        prob_remade = SciMLBase.remake(prob; u0 = u0_new, p = 2.0, tspan = tspan_new)
        @test prob_remade isa typeof(prob)
        @test prob_remade.f === prob.f
        @test prob_remade.u0 == u0_new
        @test prob_remade.p == 2.0
        @test prob_remade.tspan == tspan_new
        sol_remade = solve(prob_remade, MPRK22(1.0); dt = 0.1, adaptive = false)
        sol_new = solve(prob_new, MPRK22(1.0); dt = 0.1, adaptive = false)
        @test sol_remade.u == sol_new.u

        @test SciMLBase.remake(prob; f = prob_new.f).f === prob_new.f
    end

    @test_throws "can only be replaced by a ConservativePDSFunction" SciMLBase.remake(probs[1];
                                                                                      f = probs[2].f)
    @test_throws "can only be replaced by a ConservativePDSFunction" SciMLBase.remake(probs[1];
                                                                                      f = probs[3].f)
    @test_throws "can only be replaced by a PDSFunction" SciMLBase.remake(probs[3];
                                                                          f = probs[1].f)
    @test_throws "can only be replaced by a PDSFunction" SciMLBase.remake(probs[3];
                                                                          f = probs[4].f)
    @test_throws "ConservativePDSFunction does not support initialization data" SciMLBase.remake(probs[1].f;
                                                                                                 initialization_data = 1)
    @test_throws "PDSFunction does not support initialization data" SciMLBase.remake(probs[3].f;
                                                                                     initialization_data = 1)

    # Changing the element type of `u0` rebuilds the prototypes
    u0_d64 = Double64.(u0_new)
    tspan_d64 = Double64.(tspan_new)
    alg = MPRK22(Double64(1))
    probs_d64 = (ConservativePDSProblem(P_op, u0_d64, tspan_d64, 2.0),
                 ConservativePDSProblem(P_ip!, u0_d64, tspan_d64, 2.0),
                 PDSProblem(P_op, D_op, u0_d64, tspan_d64, 2.0),
                 PDSProblem(P_ip!, D_ip!, u0_d64, tspan_d64, 2.0))
    @testset "Double64 $i" for (i, (prob, prob_new)) in enumerate(zip(probs, probs_d64))
        prob_remade = SciMLBase.remake(prob; u0 = u0_d64, p = 2.0, tspan = tspan_d64)
        @test typeof(prob_remade) == typeof(prob_new)
        if SciMLBase.isinplace(prob)
            @test eltype(prob_remade.f.p_prototype) == Double64
        end
        sol_remade = solve(prob_remade, alg; dt = Double64(1) / 10, adaptive = false)
        sol_new = solve(prob_new, alg; dt = Double64(1) / 10, adaptive = false)
        @test sol_remade.u == sol_new.u

        # Passing the unchanged `f` is not a replacement of `f`
        prob_remade = SciMLBase.remake(prob; f = prob.f, u0 = u0_d64, p = 2.0,
                                       tspan = tspan_d64)
        @test typeof(prob_remade) == typeof(prob_new)

        # `solve` calls `remake` with the unchanged `f` to apply these keywords
        sol = solve(prob, alg; u0 = u0_d64, p = 2.0, tspan = tspan_d64,
                    dt = Double64(1) / 10, adaptive = false)
        @test sol.u == sol_new.u
    end

    # Sparse prototypes keep their sparsity pattern
    p_prototype = sparse([0.0 1.0; 1.0 0.0])
    p_prototype_d64 = Double64.(p_prototype)
    probs_sparse = (ConservativePDSProblem(P_ip!, u0, tspan, 5.0; p_prototype),
                    PDSProblem(P_ip!, D_ip!, u0, tspan, 5.0; p_prototype))
    probs_sparse_d64 = (ConservativePDSProblem(P_ip!, u0_d64, tspan_d64, 5.0;
                                               p_prototype = p_prototype_d64),
                        PDSProblem(P_ip!, D_ip!, u0_d64, tspan_d64, 5.0;
                                   p_prototype = p_prototype_d64))
    @testset "sparse $i" for (i, (prob, prob_new)) in enumerate(zip(probs_sparse,
                                                                    probs_sparse_d64))
        prob_remade = SciMLBase.remake(prob; u0 = u0_d64, tspan = tspan_d64)
        @test prob_remade.f.p_prototype isa SparseMatrixCSC{Double64}
        @test findnz(prob_remade.f.p_prototype)[1:2] == findnz(p_prototype)[1:2]
        @test typeof(prob_remade) == typeof(prob_new)
        # `Euler` evaluates the default `std_rhs` with its sparse caches
        sol_remade = solve(prob_remade, Euler(); dt = Double64(1) / 10)
        sol_new = solve(prob_new, Euler(); dt = Double64(1) / 10)
        @test sol_remade.u == sol_new.u
        @test eltype(sol_remade.u[end]) == Double64
    end

    # Rebuilding the prototypes keeps `analytic`, `linear_invariants`, and
    # a user-supplied `std_rhs`
    function std_rhs_ip!(du, u, p, t)
        du[1] = u[2] - p * u[1]
        du[2] = -du[1]
        return nothing
    end
    analytic(u0, p, t) = u0
    linear_invariants = [1.0 1.0]
    @testset "keep fields, std_rhs = $std_rhs" for std_rhs in (nothing, std_rhs_ip!)
        kwargs = (; analytic, linear_invariants, std_rhs)
        probs_kw = (ConservativePDSProblem(P_ip!, u0, tspan, 5.0; kwargs...),
                    PDSProblem(P_ip!, D_ip!, u0, tspan, 5.0; kwargs...))
        probs_kw_d64 = (ConservativePDSProblem(P_ip!, u0_d64, tspan_d64, 5.0;
                                               kwargs...),
                        PDSProblem(P_ip!, D_ip!, u0_d64, tspan_d64, 5.0; kwargs...))
        for (prob, prob_new) in zip(probs_kw, probs_kw_d64)
            prob_remade = SciMLBase.remake(prob; u0 = u0_d64, tspan = tspan_d64)
            @test typeof(prob_remade) == typeof(prob_new)
            @test eltype(prob_remade.f.p_prototype) == Double64
            @test prob_remade.f.analytic === prob.f.analytic
            @test prob_remade.f.linear_invariants === prob.f.linear_invariants
            if std_rhs === nothing
                @test prob_remade.f.std_rhs !== prob.f.std_rhs
            else
                @test prob_remade.f.std_rhs === std_rhs
            end
            # `Euler` evaluates `std_rhs`
            sol_remade = solve(prob_remade, Euler(); dt = Double64(1) / 10)
            sol_new = solve(prob_new, Euler(); dt = Double64(1) / 10)
            @test sol_remade.u == sol_new.u
        end
    end

    # `solve` promotes an integer `tspan` and calls `remake` internally
    probs_int = (ConservativePDSProblem(P_op, u0, (0, 1), 5.0),
                 ConservativePDSProblem(P_ip!, u0, (0, 1), 5.0),
                 PDSProblem(P_op, D_op, u0, (0, 1), 5.0),
                 PDSProblem(P_ip!, D_ip!, u0, (0, 1), 5.0))
    @testset "integer tspan $i" for (i, (prob_int, prob)) in enumerate(zip(probs_int,
                                                                           probs))
        sol_int = solve(prob_int, MPRK22(1.0); dt = 0.1, adaptive = false)
        sol = solve(prob, MPRK22(1.0); dt = 0.1, adaptive = false)
        @test sol_int.u == sol.u
    end

    # Changing the units of `tspan` changes the units of the prototypes
    prob = PDSProblem(P_ip!, D_ip!, u0 * u"mol", tspan .* u"s", 5.0 / u"s")
    prob_remade = SciMLBase.remake(prob; tspan = (0.0, 1.0) .* u"minute")
    @test eltype(prob_remade.f.p_prototype) == typeof(1.0u"mol" / 1.0u"minute")
    @test eltype(prob_remade.f.d_prototype) == typeof(1.0u"mol" / 1.0u"minute")

    # The size of `u0` cannot change if the prototypes depend on it
    u0_long = [0.5, 0.3, 0.2]
    @test SciMLBase.remake(probs[1]; u0 = u0_long).u0 == u0_long
    for prob in probs[2:end]
        @test_throws "cannot change the size of `u0`" SciMLBase.remake(prob;
                                                                       u0 = u0_long)
        @test_throws "cannot change the size of `u0`" SciMLBase.remake(prob;
                                                                       f = prob.f,
                                                                       u0 = u0_long)
        @test_throws "cannot change the size of `u0`" solve(prob, MPRK22(1.0);
                                                            u0 = u0_long, dt = 0.1)
    end
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
