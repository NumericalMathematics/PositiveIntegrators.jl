# Here we check that MPRK schemes require a PDSProblem or ConservativePDSProblem.
# We also check that only permissible parameters can be used.
@testitem "PDS solvers: error handling" setup=[Imports, TestHelpers] begin
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
    @test_throws "MPRKO22 can only be applied to production-destruction systems" solve(prob_oop,
                                                                                       MPRKO22(1.0,
                                                                                               0.75))
    @test_throws "MPRKO22 can only be applied to production-destruction systems" solve(prob_ip,
                                                                                       MPRKO22(1.0,
                                                                                               0.75))
    @test_throws "MPRKO22 requires α ≥ 1/2." solve(prob_pds_linmod,
                                                   MPRKO22(0.25, 0.0))
    @test_throws "For α=0.5 MPRKO22 requires 0 ≤ β ≤ 1." solve(prob_pds_linmod,
                                                               MPRKO22(0.5, 2.0))
    @test_throws "For α=1.5 MPRKO22 requires 0.25 ≤ β ≤ 0.75." solve(prob_pds_linmod,
                                                                     MPRKO22(1.5,
                                                                             0.0))
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
# of the requested type, independent of the type of the input parameters
@testitem "Constant types" setup=[Imports, TestHelpers] begin
    algs = (MPRK22(0.5f0), MPRK22(1.0f0), MPRK22(2.0f0), MPRKO22(1.0f0, 0.75f0),
            MPRK43I(1.0f0, 0.5f0),
            MPRK43I(0.5f0, 0.75f0), MPRK43II(0.5f0), MPRK43II(2.0f0 / 3.0f0),
            SSPMPRK22(0.5f0, 1.0f0),
            MPRK22(0.5), MPRK22(1.0), MPRK22(2.0), MPRKO22(1.0, 0.75),
            MPRK43I(1.0, 0.5),
            MPRK43I(0.5, 0.75), MPRK43II(0.5), MPRK43II(2.0 / 3.0),
            SSPMPRK22(0.5, 1.0),
            MPRK22(1), MPRKO22(1, 1 // 2), MPRK43I(1, 1 // 2), MPRK43II(1 // 2),
            SSPMPRK22(1 // 2, 1))
    @testset "$i" for (i, alg) in enumerate(algs)
        for T in (Float32, Float64, BigFloat)
            @test eltype(PositiveIntegrators.get_constant_parameters(alg, T)) == T
        end
    end

    # The parameters are converted before the coefficients are computed, so the
    # coefficients are not rounded to the type of the parameters.
    algs = ((MPRK22(0.6), MPRK22(big(0.6))),
            (MPRKO22(1.0, 0.75), MPRKO22(big(1.0), big(0.75))),
            (MPRK43I(1.0, 0.5), MPRK43I(big(1.0), big(0.5))),
            (MPRK43II(0.5), MPRK43II(big(0.5))),
            (SSPMPRK22(0.5, 1.0), SSPMPRK22(big(0.5), big(1.0))))
    @testset "$i" for (i, (alg, alg_big)) in enumerate(algs)
        @test PositiveIntegrators.get_constant_parameters(alg, BigFloat) ==
              PositiveIntegrators.get_constant_parameters(alg_big, BigFloat)
    end
end

# Algorithms with parameters of one floating point type can be applied to
# problems of another one
@testitem "Parameter types differing from problem types" setup=[Imports, TestHelpers] begin
    P_linmod(u, p, t) = [0 u[2]; 5*u[1] 0]
    function P_linmod!(P, u, p, t)
        fill!(P, zero(eltype(P)))
        P[1, 2] = u[2]
        P[2, 1] = 5 * u[1]
        return nothing
    end
    algs = ((MPRK22(1.0), MPRK22(big(1.0))),
            (MPRK22(1), MPRK22(big(1.0))),
            (MPRKO22(1.0, 0.75), MPRKO22(big(1.0), big(0.75))),
            (MPRK43I(1.0, 0.5), MPRK43I(big(1.0), big(0.5))),
            (MPRK43I(1, 1 // 2), MPRK43I(big(1.0), big(0.5))),
            (MPRK43II(0.5), MPRK43II(big(0.5))),
            (SSPMPRK22(0.5, 1.0), SSPMPRK22(big(0.5), big(1.0))),
            (SSPMPRK22(0, 0.75), SSPMPRK22(big(0.0), big(0.75))))
    @testset "$T" for T in (Float32, BigFloat)
        u0 = T[0.9, 0.1]
        tspan = (T(0), T(2))
        probs = (ConservativePDSProblem(P_linmod, u0, tspan),
                 ConservativePDSProblem(P_linmod!, u0, tspan))
        @testset "$i" for (i, (alg, alg_ref)) in enumerate(algs), prob in probs
            sol = solve(prob, alg; dt = T(0.1), adaptive = false)
            sol_ref = solve(prob, alg_ref; dt = T(0.1), adaptive = false)
            @test SciMLBase.successful_retcode(sol)
            @test sol.t isa Vector{T}
            @test sol.u isa Vector{Vector{T}}
            @test sol.u == sol_ref.u
        end
    end
end

# The parameters of the algorithms are checked when the algorithms are constructed,
# independent of the type of the problem
@testitem "Parameter checks" setup=[Imports, TestHelpers] begin
    @test_throws "MPRK22 requires α ≥ 1/2." MPRK22(0.25)
    @test_throws "MPRKO22 requires α ≥ 1/2." MPRKO22(0.25, 0.5)
    @test_throws "MPRK43I requires α ≥ 1/3 and α ≠ 2/3." MPRK43I(0.25, 0.5)
    @test_throws "MPRK43II requires 3/8 ≤ γ ≤ 3/4." MPRK43II(1.0)
    @test_throws "SSPMPRK22 requires 0 ≤ α ≤ 1, β ≥ 0 and αβ + 1/(2β) ≤ 1." SSPMPRK22(-1.0,
                                                                                      1.0)

    # 2/3 rounded to the type of the parameters
    @testset "$(typeof(alpha))" for alpha in (2 / 3, 2.0f0 / 3, big(2) / 3,
                                              Double64(2) / 3, 2 // 3)
        @test_throws "MPRK43I requires α ≥ 1/3 and α ≠ 2/3." MPRK43I(alpha, one(alpha) / 2)
    end

    # The lower bound of β for MPRK43I changes from 3α(1-α) to (3α-2)/(6α-3) at
    # α0 ≈ 0.8926
    @test MPRK43I(0.85, 0.4) isa MPRK43I
    @test_throws "For this choice of α MPRK43I requires 3α(1-α) ≤ β ≤ 2/3." MPRK43I(0.85,
                                                                                    0.3)
    @test MPRK43I(0.95, 0.4) isa MPRK43I
    @test_throws "For this choice of α MPRK43I requires (3α-2)/(6α-3) ≤ β ≤ 2/3." MPRK43I(0.95,
                                                                                          0.2)

    # Parameters on the boundary of the admissible set are accepted for all
    # problem types, although the bound computed in another type can differ
    alpha = 1.2
    alg = MPRKO22(alpha, (alpha - 1) / (2 * alpha - 1))
    @testset "$T" for T in (Float32, Float64, BigFloat, Double64)
        @test eltype(PositiveIntegrators.get_constant_parameters(alg, T)) == T
    end
end

# We check that Float32 inputs result in Float32 outputs
# TODO: Use a tool such as LIKWID to check that no double-precision calculations are performed internally
@testitem "Float32 input data" setup=[Imports, TestHelpers] begin
    P_linmod(u, p, t) = [0 u[2]; 5*u[1] 0]
    u0 = [0.9f0, 0.1f0]
    prob = ConservativePDSProblem(P_linmod, u0, (0.0f0, 2.0f0))

    algs = (MPRK22(0.5f0), MPRK22(1.0f0), MPRK22(2.0f0), MPRKO22(1.0f0, 0.75f0),
            MPRK43I(1.0f0, 0.5f0),
            MPRK43I(0.5f0, 0.75f0), MPRK43II(0.5f0), MPRK43II(2.0f0 / 3.0f0),
            SSPMPRK22(0.5f0, 1.0f0), SSPMPRK43(), MPDeC(2),
            MPDeC(2, nodes = :lagrange))
    @testset "$(nameof(typeof(alg)))" for alg in algs
        sol = solve(prob, alg; dt = 0.1f0, save_everystep = false)
        @test sol.t isa Vector{Float32}
        @test sol.u isa Vector{Vector{Float32}}
    end

    # in-place
    prob_ip = ConservativePDSProblem((P, u, p, t) -> (P .= P_linmod(u, p, t); nothing),
                                     u0, (0.0f0, 2.0f0))
    @testset "$(nameof(typeof(alg))) (in-place)" for alg in algs
        sol = solve(prob_ip, alg; dt = 0.1f0, save_everystep = false)
        @test SciMLBase.successful_retcode(sol)
        @test sol.t isa Vector{Float32}
        @test sol.u isa Vector{Vector{Float32}}
    end
end

# Here we check that MPE equals implicit Euler (IE) for a linear PDS
@testitem "Linear model problem: MPE = IE?" setup=[Imports, TestHelpers] begin
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
@testitem "MPRK22(α) = SSPMPRK22(0, α) = MPRK22(α, 0)" setup=[Imports, TestHelpers] begin
    for α in (0.5, 2.0 / 3.0, 1.0, 2.0)
        # conservative PDS
        sol1 = solve(prob_pds_linmod, MPRK22(α))
        sol2 = solve(prob_pds_linmod, SSPMPRK22(0.0, α))
        sol3 = solve(prob_pds_linmod_inplace, MPRK22(α))
        sol4 = solve(prob_pds_linmod_inplace, SSPMPRK22(0.0, α))
        @test sol1.u ≈ sol2.u ≈ sol3.u ≈ sol4.u
        if α ≤ 1
            sol5 = solve(prob_pds_linmod, MPRKO22(α, 0.0))
            sol6 = solve(prob_pds_linmod_inplace, MPRKO22(α, 0.0))
            @test sol1.u ≈ sol5.u ≈ sol6.u
        end

        # nonconservative PDS
        sol1 = solve(prob_pds_linmod_nonconservative, MPRK22(α))
        sol2 = solve(prob_pds_linmod_nonconservative, SSPMPRK22(0.0, α))
        sol3 = solve(prob_pds_linmod_nonconservative_inplace, MPRK22(α))
        sol4 = solve(prob_pds_linmod_nonconservative_inplace, SSPMPRK22(0.0, α))
        @test sol1.u ≈ sol2.u ≈ sol3.u ≈ sol4.u
        if α ≤ 1
            sol5 = solve(prob_pds_linmod_nonconservative, MPRKO22(α, 0.0))
            sol6 = solve(prob_pds_linmod_nonconservative_inplace, MPRKO22(α, 0.0))
            @test sol1.u ≈ sol5.u ≈ sol6.u
        end
    end
end

# Here we check that MPRK22(1.0) = MPDeC(2)
@testitem "MPRK22(1.0) = MPDeC(2)" setup=[Imports, TestHelpers] begin
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
@testitem "Different linear solvers" setup=[Imports, TestHelpers] begin
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
        (; kwargs...) -> MPRKO22(1.0, 0.75; kwargs...),
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

# A failed linear solve rejects the step: a fixed-step solve stops with
# ReturnCode.ConvergenceFailure, an adaptive solve retries with a smaller dt.
# For the linear model problems, the system matrices have an off-diagonal entry
# larger than 0.25 times the diagonal entry of its column for dt = 1 and only
# smaller ones for dt = 0.01.
@testitem "Failed linear solves" setup=[Imports, TestHelpers] begin
    linsolve = ThresholdLinsolve(0.25)
    algs = [MPE(; linsolve), MPRK22(1.0; linsolve), MPRK22(0.5; linsolve),
        MPRKO22(1.0, 0.75; linsolve), MPRK43I(1.0, 0.5; linsolve),
        MPRK43II(0.5; linsolve), SSPMPRK22(0.5, 1.0; linsolve),
        SSPMPRK43(; linsolve), MPDeC(3; linsolve)]
    probs = [prob_pds_linmod, prob_pds_linmod_inplace,
        prob_pds_linmod_nonconservative, prob_pds_linmod_nonconservative_inplace]
    increase_dt = SciMLBase.DiscreteCallback((u, t, integrator) -> t >= 0.5,
                                             integrator -> SciMLBase.set_proposed_dt!(integrator,
                                                                                      1.0);
                                             save_positions = (false, false))
    @testset "$alg" for alg in algs
        @testset "$i" for (i, prob) in enumerate(probs)
            sol = solve(prob, alg; dt = 1.0, adaptive = false)
            @test sol.retcode == SciMLBase.ReturnCode.ConvergenceFailure
            @test sol.t == [first(prob.tspan)]
            @test sol.stats.nreject == 1

            sol = solve(prob, alg; dt = 0.01, adaptive = false)
            @test SciMLBase.successful_retcode(sol)
            @test sol.t[end] == last(prob.tspan)

            # A failure after accepted steps keeps the last accepted value
            # as the final value, also if it is saved only at the end.
            sol_all = solve(prob, alg; dt = 0.01, adaptive = false,
                            callback = increase_dt)
            sol_end = solve(prob, alg; dt = 0.01, adaptive = false,
                            callback = increase_dt, save_everystep = false)
            @test sol_end.retcode == SciMLBase.ReturnCode.ConvergenceFailure
            @test sol_end.t[end] == sol_all.t[end] > 0.5
            @test sol_end.u[end] == sol_all.u[end]
            @test sol_end.stats.nreject == 1

            if SciMLBase.isadaptive(alg)
                sol = solve(prob, alg; dt = 1.0)
                @test sol.stats.nreject > 0
                @test SciMLBase.successful_retcode(sol)
                @test sol.t[end] == last(prob.tspan)
                u_analytic = prob.f.analytic(prob.u0, prob.p, last(prob.tspan))
                @test isapprox(sol.u[end], u_analytic; rtol = 1.0e-2)
            end
        end
    end
end

# The diagonal of the production matrix of a ConservativePDSProblem is not used:
# the results with a nonzero diagonal equal those with a zero diagonal for
# out-of-place problems and in-place problems with all matrix types, and the sum of
# the components is conserved.
@testitem "Conservative PDS with nonzero diagonal of P" setup=[Imports, TestHelpers] begin
    P_diag(u, p, t) = [0.3*u[1] u[2]; 5*u[1] 0.2*u[2]]
    P_zero(u, p, t) = [0 u[2]; 5*u[1] 0]
    function P_diag!(P, u, p, t)
        P[1, 1] = 0.3 * u[1]
        P[1, 2] = u[2]
        P[2, 1] = 5 * u[1]
        P[2, 2] = 0.2 * u[2]
        return nothing
    end
    function P_zero!(P, u, p, t)
        P[1, 1] = 0
        P[1, 2] = u[2]
        P[2, 1] = 5 * u[1]
        P[2, 2] = 0
        return nothing
    end
    u0 = [0.9, 0.1]
    tspan = (0.0, 2.0)
    probs = Any[(ConservativePDSProblem(P_diag, u0, tspan),
                 ConservativePDSProblem(P_zero, u0, tspan))]
    for prototype in (zeros(2, 2), Tridiagonal(zeros(2, 2)), sparse(ones(2, 2)))
        push!(probs,
              (ConservativePDSProblem(P_diag!, u0, tspan; p_prototype = prototype),
               ConservativePDSProblem(P_zero!, u0, tspan; p_prototype = prototype)))
    end
    algs = [MPE(), MPRK22(0.5), MPRK22(1.0), MPRKO22(1.0, 0.75), MPRK43I(1.0, 0.5),
        MPRK43II(0.5), SSPMPRK22(0.5, 1.0), SSPMPRK43()]
    for K in 2:10, nodes in (:gausslobatto, :lagrange)
        push!(algs, MPDeC(K; nodes))
    end
    @testset "$alg" for alg in algs
        @testset "$i" for (i, (prob_diag, prob_zero)) in enumerate(probs)
            sol_diag = solve(prob_diag, alg; dt = 0.1, adaptive = false)
            sol_zero = solve(prob_zero, alg; dt = 0.1, adaptive = false)
            @test SciMLBase.successful_retcode(sol_diag)
            @test sol_diag.u == sol_zero.u
            @test all(u -> isapprox(sum(u), sum(u0); rtol = 1.0e-14), sol_diag.u)
        end
    end
end

# The Patankar systems are solved in weight form M w = b with M = diag(σ) + B and
# x = σ .* w. The matrices of all storage types are M_classical * diag(σ) with the
# classical matrix M_classical = I + B diag(1 ./ σ), and the solutions agree with
# the solutions of the classical systems.
@testitem "Patankar systems in weight form" setup=[Imports, TestHelpers] begin
    n = 6
    # tridiagonal production terms, so that all storage types can be used
    P = Matrix(Tridiagonal([0.1 + i / n for i in 1:(n - 1)],
                           [0.3 * i for i in 1:n],
                           [1.5 - i / n for i in 1:(n - 1)]))
    σ = [0.2 + cos(i)^2 for i in 1:n]
    b = [0.5 + sin(i)^2 for i in 1:n]
    dt = 0.7
    @testset "d = $(d === nothing ? "nothing" : "vector")" for d in (nothing,
                                                                     [0.1 * i for i in 1:n])
        dvec = d === nothing ? zeros(n) : d
        B = [i == j ? dt * (dvec[j] + sum(P[k, j] for k in 1:n if k != j)) : -dt * P[i, j]
             for i in 1:n, j in 1:n]
        M_classical = I + B * Diagonal(1 ./ σ)
        M_ref = M_classical * Diagonal(σ)
        rhs = d === nothing ? b : b + dt * diag(P)
        x_ref = M_classical \ rhs

        M_dense = PositiveIntegrators.build_mprk_matrix(P, σ, dt, d)
        M_static = PositiveIntegrators.build_mprk_matrix(SMatrix{n, n}(P), SVector{n}(σ),
                                                         dt,
                                                         d === nothing ? nothing :
                                                         SVector{n}(d))
        M_tri = Tridiagonal(zeros(n, n))
        PositiveIntegrators.build_mprk_matrix!(M_tri, Tridiagonal(P), σ, dt, d)
        M_sparse = PositiveIntegrators.build_mprk_matrix(sparse(P), σ, dt, d)
        for M in (M_dense, M_static, M_tri, M_sparse)
            @test isapprox(Matrix(M), M_ref; rtol = 1.0e-14)
        end

        # out-of-place steps
        for (Pk, σk, bk, dk) in ((P, σ, b, d), (sparse(P), σ, b, d),
                                 (SMatrix{n, n}(P), SVector{n}(σ), SVector{n}(b),
                                  d === nothing ? nothing : SVector{n}(d)))
            x, solved = PositiveIntegrators.basic_patankar_step(bk, Pk, σk, dt,
                                                                LUFactorization(), dk)
            @test solved
            @test isapprox(x, x_ref; rtol = 1.0e-13)
            @test all(>(0), x)
            d === nothing && @test isapprox(sum(x), sum(b); rtol = 1.0e-14)
        end

        # in-place steps with the storage types of in-place problems
        for (Pk, linsolve) in ((P, LUFactorization()), (Tridiagonal(P), LUFactorization()),
                               (sparse(P), KLUFactorization()),
                               (sparse(P), UMFPACKFactorization()), (sparse(P), nothing))
            Pk = copy(Pk)
            rhs_ip = zeros(n)
            cache = init(LinearSolve.LinearProblem(Pk, rhs_ip), linsolve;
                         alias = LinearSolve.LinearAliasSpecifier(; alias_A = true,
                                                                  alias_b = true),
                         assumptions = LinearSolve.OperatorAssumptions(true))
            x = zeros(n)
            @test PositiveIntegrators.basic_patankar_step!(x, b, Pk, d, σ, dt, cache)
            @test isapprox(x, x_ref; rtol = 1.0e-13)
            # the output may alias the Patankar weight denominators
            Pk = copy(cache.A)
            copyto!(Pk isa SparseMatrixCSC ? nonzeros(Pk) : Pk,
                    Pk isa SparseMatrixCSC ? nonzeros(sparse(P)) : P)
            σ_alias = copy(σ)
            @test PositiveIntegrators.basic_patankar_step!(σ_alias, b, Pk, d, σ_alias, dt,
                                                           cache)
            @test isapprox(σ_alias, x_ref; rtol = 1.0e-13)
        end
    end
end

# A Patankar weight denominator σ[j] = Inf gives column j = e_j of the classical
# matrix. In weight form, column j is oneunit(σ[j]) * e_j and x[j] = w[j].
@testitem "MPRK matrix with infinite denominators" setup=[Imports, TestHelpers] begin
    P = [0.0 2.0 0.0; 1.5 0.0 3.0; 0.0 0.7 0.0]
    dt = 0.3
    σ = [0.5, Inf, 2.0]
    b = [0.4, 1.3, 0.8]
    @testset "d = $d" for d in (nothing, [0.2, 0.0, 1.1])
        s = vec(sum(P; dims = 1)) .+ (d === nothing ? 0.0 : d)
        M_classical = [i == j ? 1 + dt * s[i] / σ[i] : -dt * P[i, j] / σ[j]
                       for i in 1:3, j in 1:3]
        @test M_classical[2, 2] == 1
        @test all(iszero, M_classical[[1, 3], 2])
        M_ref = [j == 2 ? float(i == j) :
                 (i == j ? σ[i] + dt * s[i] : -dt * P[i, j])
                 for i in 1:3, j in 1:3]
        rhs = d === nothing ? b : b + dt * diag(P)
        x_ref = M_classical \ rhs

        M_dense = PositiveIntegrators.build_mprk_matrix(P, σ, dt, d)
        M_static = PositiveIntegrators.build_mprk_matrix(SMatrix{3, 3}(P),
                                                         SVector{3}(σ), dt,
                                                         d === nothing ?
                                                         nothing :
                                                         SVector{3}(d))
        M_tri = Tridiagonal(zeros(3, 3))
        PositiveIntegrators.build_mprk_matrix!(M_tri, Tridiagonal(P), σ, dt, d)
        M_sparse = sparse(P + I)
        PositiveIntegrators.build_mprk_matrix!(M_sparse, sparse(P), σ, dt, d)
        M_alias = sparse(P + I)
        for i in 1:3
            M_alias[i, i] = P[i, i]
        end
        PositiveIntegrators.build_mprk_matrix!(M_alias, M_alias, σ, dt, d)

        for M in (M_dense, M_static, M_tri, M_sparse, M_alias)
            @test Matrix(M) ≈ M_ref
            x, solved = PositiveIntegrators.patankar_solve(copy(M), rhs, σ,
                                                           LUFactorization())
            @test solved
            @test x ≈ x_ref
        end
        @test nonzeros(M_alias) == nonzeros(M_sparse)

        # unit columns from σ[j] = Inf are not extreme
        for (Mk, Pk) in ((zeros(3, 3), P), (Tridiagonal(zeros(3, 3)), Tridiagonal(P)))
            @test !PositiveIntegrators.build_mprk_matrix!(Mk, Pk, σ, dt, d)
            @test !PositiveIntegrators.has_extreme_diagonal(Mk)
        end

        # in-place steps, also with the output aliasing σ
        P_sparse = sparse(P + I)
        for i in 1:3
            P_sparse[i, i] = P[i, i]
        end
        for (Pk, linsolve) in ((P, LUFactorization()), (Tridiagonal(P), LUFactorization()),
                               (P_sparse, KLUFactorization()))
            for alias in (false, true)
                Pc = copy(Pk)
                cache = init(LinearSolve.LinearProblem(Pc, zeros(3)), linsolve;
                             alias = LinearSolve.LinearAliasSpecifier(; alias_A = true,
                                                                      alias_b = true),
                             assumptions = LinearSolve.OperatorAssumptions(true))
                σc = copy(σ)
                x = alias ? σc : zeros(3)
                @test PositiveIntegrators.basic_patankar_step!(x, b, Pc, d, σc, dt, cache)
                @test x ≈ x_ref
            end
        end
    end

    # the number of weights must match the size of the matrix
    for M in (zeros(3, 3), Tridiagonal(zeros(3, 3)), sparse(1.0I, 3, 3))
        @test_throws DimensionMismatch PositiveIntegrators.set_unit_columns!(M, [1.0, Inf])
    end
end

# Columns of dense, `Tridiagonal`, and static Patankar systems with tiny diagonal
# entries are scaled by powers of two. Without the scaling, the Patankar weight
# x[c] / σ[c] of a component with σ[c] = floatmin and no outflow that receives
# x[c] > 4 overflows.
@testitem "Column scaling of dense Patankar systems" setup=[Imports, TestHelpers] begin
    bound = 2.0^511
    @test PositiveIntegrators.extreme_column_scale(1.0) == 1.0
    @test PositiveIntegrators.extreme_column_scale(3.0e10) == 1.0
    @test PositiveIntegrators.extreme_column_scale(bound) == 1.0
    @test PositiveIntegrators.extreme_column_scale(inv(bound)) == 1.0
    @test PositiveIntegrators.extreme_column_scale(2 * bound) == 2.0^-512
    @test PositiveIntegrators.extreme_column_scale(floatmin()) == 2.0^1022
    @test PositiveIntegrators.extreme_column_scale(1.0e-310) == 2.0^1023
    @test PositiveIntegrators.extreme_column_scale(0.0) == 1.0
    @test PositiveIntegrators.extreme_column_scale(Inf) == 1.0
    @test PositiveIntegrators.extreme_column_scale(NaN) == 1.0
    @test PositiveIntegrators.extreme_column_scale(floatmin(Float32)) === 2.0f0^126

    # the builders of dense and Tridiagonal matrices report extreme diagonal entries
    P_small = [0.0 1.0; 2.0 0.0]
    for (σ_small, extreme) in (([1.0, 2.0], false), ([floatmin(), 0.0], true))
        for M in (zeros(2, 2), Tridiagonal(zeros(2, 2)))
            dt_small = extreme ? 0.0 : 0.1
            @test PositiveIntegrators.build_mprk_matrix!(M,
                                                         M isa Tridiagonal ?
                                                         Tridiagonal(P_small) :
                                                         P_small,
                                                         σ_small, dt_small) == extreme
            @test PositiveIntegrators.has_extreme_diagonal(M) == extreme
        end
    end

    # chain of 5 components exchanging mass with their neighbors; component 3
    # receives mass but has no outflow
    n = 5
    c = 3
    P = Matrix(Tridiagonal([2.0, 30.0, 0.0, 1.5], zeros(n), [0.5, 0.0, 20.0, 1.0]))
    σ = [1.0, 0.7, floatmin(), 1.3, 0.9]
    b = [2.0, 4.0, 0.0, 3.0, 1.0]
    dt = 0.5
    M_classical = [i == j ? 1 + dt * sum(P[k, j] for k in 1:n if k != j) / σ[j] :
                   -dt * P[i, j] / σ[j] for i in 1:n, j in 1:n]
    x_ref = M_classical \ b
    @test x_ref[c] / σ[c] > floatmax()

    for (Pk, σk, bk) in ((P, σ, b), (Tridiagonal(P), σ, b),
                         (SMatrix{n, n}(P), SVector{n}(σ), SVector{n}(b)))
        x, solved = PositiveIntegrators.basic_patankar_step(bk, Pk, σk, dt,
                                                            LUFactorization(), nothing)
        @test solved
        @test isapprox(x, x_ref; rtol = 1.0e-14)
    end
    for Pk in (copy(P), Tridiagonal(P))
        cache = init(LinearSolve.LinearProblem(Pk, zeros(n)), LUFactorization();
                     alias = LinearSolve.LinearAliasSpecifier(; alias_A = true,
                                                              alias_b = true),
                     assumptions = LinearSolve.OperatorAssumptions(true))
        x = zeros(n)
        @test PositiveIntegrators.basic_patankar_step!(x, b, Pk, nothing, σ, dt, cache)
        @test isapprox(x, x_ref; rtol = 1.0e-14)
    end

    # Component 2 has no outflow and a denominator σ[2] whose diagonal entry is not
    # extreme, but the Patankar weight x[2] / σ[2] overflows. The bound ‖b‖₁ / σ[2] of
    # the weight exceeds 2^e, so column 2 is scaled.
    @test PositiveIntegrators.extreme_column_scale(1.0, true) == 1.0
    @test PositiveIntegrators.extreme_column_scale(0.75, true) == 2.0
    @testset "Weight bound $T" for (T, σ2, b1) in ((Float16, Float16(2)^-7, Float16(600)),
                                                   (Float32, 1.0f-15, 4.0f23),
                                                   (Float64, 2.0^-500, 1.0e160))
        P3 = T[0 0 0; 1 0 1; 0 0 0]
        σ3 = T[1, σ2, 1]
        b3 = T[b1, 0, b1]
        dt3 = one(T)
        M_classical = [i == j ? 1 + dt3 * sum(P3[k, j] for k in 1:3 if k != j) / σ3[j] :
                       -dt3 * P3[i, j] / σ3[j] for i in 1:3, j in 1:3]
        x_ref3 = M_classical \ b3
        @test !isfinite(x_ref3[2] / σ3[2])
        M = PositiveIntegrators.build_mprk_matrix(P3, σ3, dt3)
        @test !PositiveIntegrators.has_extreme_diagonal(M)
        @test PositiveIntegrators.scale_columns(copy(M), σ3, b3)[2] !== nothing
        for (Pk, σk, bk) in ((P3, σ3, b3), (Tridiagonal(P3), σ3, b3),
                             (SMatrix{3, 3}(P3), SVector{3}(σ3), SVector{3}(b3)))
            x, solved = PositiveIntegrators.basic_patankar_step(bk, Pk, σk, dt3,
                                                                LUFactorization(), nothing)
            @test solved
            @test isapprox(x, x_ref3; rtol = 4 * eps(T))
        end
        for Pk in (copy(P3), Tridiagonal(P3))
            cache = init(LinearSolve.LinearProblem(Pk, zeros(T, 3)), LUFactorization();
                         alias = LinearSolve.LinearAliasSpecifier(; alias_A = true,
                                                                  alias_b = true),
                         assumptions = LinearSolve.OperatorAssumptions(true))
            x = zeros(T, 3)
            @test PositiveIntegrators.basic_patankar_step!(x, b3, Pk, nothing, σ3, dt3,
                                                           cache)
            @test isapprox(x, x_ref3; rtol = 4 * eps(T))
        end
    end
end

# The scaling of columns with extreme diagonal entries leaves the other components of
# the solution bitwise unchanged; without extreme entries, the solution is bitwise the
# solution of the unscaled system.
@testitem "Column scaling of dense Patankar systems (bitwise)" setup=[Imports, TestHelpers] begin
    n = 4
    P = Matrix(Tridiagonal([1.0, 3.0, 0.5], zeros(n), [2.0, 0.7, 1.5]))
    dt = 0.3
    b = [1.0, 0.5, 2.0, 1.5]
    # σ[3] = 1e-200 gives the extreme diagonal entry M[3, 3] = σ[3] for a component
    # without outflow, but no overflow of the weight x[3] / σ[3]
    P_extreme = copy(P)
    P_extreme[2, 3] = 0.0
    P_extreme[4, 3] = 0.0
    @testset "$(c === nothing ? "normal" : "extreme")" for (c, Pc, σ) in ((nothing, P,
                                                                           [
                                                                               0.5,
                                                                               1.0,
                                                                               2.0,
                                                                               0.7
                                                                           ]),
                                                                          (3, P_extreme,
                                                                           [
                                                                               0.5,
                                                                               1.0,
                                                                               1.0e-200,
                                                                               0.7
                                                                           ]))
        M = PositiveIntegrators.build_mprk_matrix(Pc, σ, dt)
        @test PositiveIntegrators.has_extreme_diagonal(M) == (c !== nothing)
        # out-of-place solves scale only systems with extreme diagonal entries
        for Mk in (copy(M), Tridiagonal(M), SMatrix{n, n}(M))
            _, factors = PositiveIntegrators.scale_columns(Mk, σ, b)
            @test (factors === nothing) == (c === nothing)
        end
        x_unscaled = σ .* (lu(M) \ b)
        x_tri_unscaled = σ .* (lu(Tridiagonal(M)) \ b)
        x_static_unscaled = σ .* (lu(SMatrix{n, n}(M)) \ SVector{n}(b))
        results = Any[]
        x, _ = PositiveIntegrators.basic_patankar_step(b, Pc, σ, dt, LUFactorization(),
                                                       nothing)
        push!(results, (x, x_unscaled))
        x, _ = PositiveIntegrators.basic_patankar_step(SVector{n}(b), SMatrix{n, n}(Pc),
                                                       SVector{n}(σ), dt,
                                                       LUFactorization(), nothing)
        push!(results, (x, x_static_unscaled))
        for (Pk, x_ref) in ((copy(Pc), x_unscaled), (Tridiagonal(Pc), x_tri_unscaled))
            cache = init(LinearSolve.LinearProblem(Pk, zeros(n)), LUFactorization();
                         alias = LinearSolve.LinearAliasSpecifier(; alias_A = true,
                                                                  alias_b = true),
                         assumptions = LinearSolve.OperatorAssumptions(true))
            x = zeros(n)
            @test PositiveIntegrators.basic_patankar_step!(x, b, Pk, nothing, σ, dt, cache)
            push!(results, (x, x_ref))
        end
        for (x, x_ref) in results
            for i in 1:n
                if i == c
                    @test isapprox(x[i], x_ref[i]; rtol = 1.0e-14)
                else
                    @test x[i] === x_ref[i]
                end
            end
        end
    end
end

# The columns of sparse Patankar systems are scaled by powers of two. Without the
# scaling, the following system with a component that receives mass but has a
# Patankar weight denominator σ[c] = floatmin and no outflow (as at the front of the
# porous medium equation) makes KLU report a failure, UMFPACK return non-finite values
# with a successful return code, and LinearSolve's default sparse solver return a wrong
# finite solution with a successful return code.
@testitem "Column scaling of sparse Patankar systems" setup=[Imports, TestHelpers] begin
    @test PositiveIntegrators.column_scale(3.0) == 0.5
    @test PositiveIntegrators.column_scale(0.75) == 2.0
    @test PositiveIntegrators.column_scale(1.0) == 1.0
    @test PositiveIntegrators.column_scale(floatmax()) == 2.0^-1023
    @test PositiveIntegrators.column_scale(1.0e-310) == 2.0^1023
    @test PositiveIntegrators.column_scale(nextfloat(0.0)) == 2.0^1023
    @test PositiveIntegrators.column_scale(0.0) == 1.0
    @test PositiveIntegrators.column_scale(Inf) == 1.0
    @test PositiveIntegrators.column_scale(NaN) == 1.0
    @test PositiveIntegrators.column_scale(3.0f0) === 0.5f0
    @test PositiveIntegrators.column_scale(1.0f-40) === 2.0f0^127
    @test PositiveIntegrators.column_scale(prevfloat(2.0)) == 1.0
    @test PositiveIntegrators.column_scale(floatmin()) == 2.0^1022
    @test PositiveIntegrators.column_scale(floatmax(Float16)) === Float16(2)^-15

    # 4 x 4 grid with exchange between neighboring cells
    m = 4
    n = m^2
    L = LinearIndices((m, m))
    rows = Int[]
    cols = Int[]
    for c in CartesianIndices((m, m)), offset in ((1, 0), (-1, 0), (0, 1), (0, -1))
        c2 = c + CartesianIndex(offset)
        checkbounds(Bool, L, c2) || continue
        push!(rows, L[c2])
        push!(cols, L[c])
    end
    rates = [10.0^(4 * mod(0.37 * k^2 + 0.11 * k, 1) - 2) for k in eachindex(rows)]
    P = sparse(vcat(rows, 1:n), vcat(cols, 1:n), vcat(rates, zeros(n)), n, n)
    σ = [10.0^(2 * mod(0.61 * k^2 + 0.3 * k, 1) - 1) for k in 1:n]
    b = [10.0^(2 * mod(0.47 * k^2 + 0.7 * k, 1) - 1) for k in 1:n]
    dt = 0.1
    # The component c receives mass but has no outflow.
    c = 6
    for idx in nzrange(P, c)
        rowvals(P)[idx] != c && (nonzeros(P)[idx] = 0.0)
    end
    σ[c] = floatmin()
    b[c] = 0.0

    Pd = Matrix(P)
    M_classical = [i == j ?
                   1 + dt * sum(Pd[k, j] for k in 1:n if k != j) / σ[j] :
                   -dt * Pd[i, j] / σ[j] for i in 1:n, j in 1:n]
    x_ref = M_classical \ b
    @test x_ref[c] > 0.1

    @testset "$(nameof(typeof(linsolve)))" for linsolve in (KLUFactorization(),
                                                            UMFPACKFactorization(),
                                                            LUFactorization(), nothing)
        x, solved = PositiveIntegrators.basic_patankar_step(b, P, σ, dt, linsolve,
                                                            nothing)
        @test solved
        @test isapprox(x, x_ref; rtol = 1.0e-13)

        Pk = copy(P)
        cache = init(LinearSolve.LinearProblem(Pk, zeros(n)), linsolve;
                     alias = LinearSolve.LinearAliasSpecifier(; alias_A = true,
                                                              alias_b = true),
                     assumptions = LinearSolve.OperatorAssumptions(true))
        x = zeros(n)
        @test PositiveIntegrators.basic_patankar_step!(x, b, Pk, nothing, σ, dt, cache)
        @test isapprox(x, x_ref; rtol = 1.0e-13)
    end
end

# A tiny Patankar weight denominator σ[1] = floatmin makes the entries
# dt * P[i, 1] / σ[1] of the classical matrix overflow to Inf. The matrix of the
# weight form stays finite, and the solution is accurate.
@testitem "Patankar systems with tiny denominators" setup=[Imports, TestHelpers] begin
    P = [0.0 2.0 0.0 1.0; 30.0 0.0 1.0 0.0; 10.0 4.0 0.0 2.0; 20.0 0.0 3.0 0.0]
    σ = [floatmin(), 0.5, 1.5, 2.0]
    b = [0.0, 1.0, 0.7, 0.3]
    dt = 1.0
    d = [0.1, 0.2, 0.0, 0.3]
    @test !all(isfinite, -dt * P[:, 1] / σ[1])

    @testset "d = $(dd === nothing ? "nothing" : "vector")" for dd in (nothing, d)
        dvec = dd === nothing ? zeros(4) : dd
        rhs = dd === nothing ? b : b + dt * diag(P)
        # reference: classical form in high precision
        x_ref = setprecision(BigFloat, 2048) do
            Pb = big.(P)
            σb = big.(σ)
            dtb = big(dt)
            Mb = [i == j ?
                  1 + dtb * (big(dvec[j]) + sum(Pb[k, j] for k in 1:4 if k != j)) / σb[j] :
                  -dtb * Pb[i, j] / σb[j] for i in 1:4, j in 1:4]
            Float64.(Mb \ big.(rhs))
        end
        @test all(>(0), x_ref)

        for (Pk, σk, bk, dk) in ((P, σ, b, dd), (sparse(P), σ, b, dd),
                                 (SMatrix{4, 4}(P), SVector{4}(σ), SVector{4}(b),
                                  dd === nothing ? nothing : SVector{4}(dd)))
            M = PositiveIntegrators.build_mprk_matrix(Pk, σk, dt, dk)
            @test all(isfinite, M)
            x, solved = PositiveIntegrators.basic_patankar_step(bk, Pk, σk, dt,
                                                                LUFactorization(), dk)
            @test solved
            @test isapprox(x, x_ref; rtol = 1.0e-12)
        end
    end
end

# The column scaling of dual numbers is determined by their primal values. Without it,
# the Patankar weight of the second component overflows in the first step. All schemes
# use coefficients and small constants of the primal type.
@testitem "Patankar systems with dual numbers" setup=[Imports, TestHelpers] begin
    Dual = ForwardDiff.Dual
    @test PositiveIntegrators.column_scale(Dual(3.0, 1.0)) === 0.5
    @test PositiveIntegrators.column_scale(Dual(Dual(3.0f0, 1.0f0), 1.0f0)) === 0.5f0
    @test PositiveIntegrators.extreme_bounds(typeof(Dual(1.0, 1.0))) ==
          PositiveIntegrators.extreme_bounds(Float64)
    @test PositiveIntegrators.extreme_column_scale(Dual(floatmin(), 1.0)) === 2.0^1022
    @test PositiveIntegrators.extreme_column_scale(Dual(1.0, 1.0)) === 1.0
    @test PositiveIntegrators.extreme_column_scale(Dual(0.75, 1.0), true) === 2.0

    # u[1] is converted into u[2] with rate p; u[2] has no outflow
    u0 = [100.0, 0.0]
    dt = 0.25
    function final_value(alg, p; static)
        u0p = convert.(typeof(p), u0)
        if static
            prob = ConservativePDSProblem((u, p, t) -> @SMatrix([0 0; p*u[1] 0]),
                                          SVector{2}(u0p), (0.0, 1.0), p)
        else
            prob = ConservativePDSProblem((u, p, t) -> [0 0; p*u[1] 0], u0p, (0.0, 1.0),
                                          p)
        end
        sol = solve(prob, alg; dt, adaptive = false)
        @test SciMLBase.successful_retcode(sol)
        return sol.u[end][2]
    end

    algs = (MPE(), MPRK22(0.5), MPRK22(1.0), MPRKO22(1.0, 0.75), MPRK43I(1.0, 0.5),
            MPRK43II(0.5), SSPMPRK22(0.5, 1.0), SSPMPRK43(), MPDeC(2), MPDeC(4))
    @testset "$alg, static = $static" for alg in algs, static in (false, true)
        # The coefficients and small constants of the schemes have the type of the primal
        # values, so that the primal values agree with a Float64 solve.
        dual_value = final_value(alg, Dual(1.0, 1.0); static)
        @test ForwardDiff.value(dual_value) === final_value(alg, 1.0; static)

        derivative = ForwardDiff.derivative(p -> final_value(alg, p; static), 1.0)
        if alg isa MPE
            # MPE is the implicit Euler method for this linear problem, so that
            # u[2] = 100 - 100 / (1 + p * dt)^4 at t = 1.
            @test derivative ≈ 400 * dt / (1 + dt)^5
        else
            h = 1.0e-6
            fd = (final_value(alg, 1 + h; static) - final_value(alg, 1 - h; static)) /
                 (2 * h)
            @test isapprox(derivative, fd; rtol = 1.0e-6)
        end
    end
end

# LinearSolve's default solver throws an exception for some singular `Tridiagonal`
# matrices instead of reporting a failure. With `linsolve = nothing`, the schemes
# use `LUFactorization` for `Tridiagonal` matrices, so that the step is rejected.
@testitem "Tridiagonal matrices with linsolve = nothing" setup=[Imports, TestHelpers] begin
    @test PositiveIntegrators.linsolve_algorithm(nothing,
                                                 Tridiagonal(ones(2), ones(3),
                                                             ones(2))) isa
          LUFactorization
    @test PositiveIntegrators.linsolve_algorithm(nothing, ones(3, 3)) === nothing
    @test PositiveIntegrators.linsolve_algorithm(KLUFactorization(),
                                                 Tridiagonal(ones(2), ones(3),
                                                             ones(2))) isa
          KLUFactorization

    # singular matrix: the first two rows are equal
    M = Tridiagonal([1.0, 1.0], [1.0, 1.0, 1.0], [1.0, 0.0])
    x, solved = PositiveIntegrators.patankar_solve(M, ones(3), ones(3), nothing)
    @test !solved

    # production terms that are NaN for t ≥ 0.5
    n = 4
    function prod!(P, u, p, t)
        fill!(P, zero(eltype(P)))
        for i in 1:(n - 1)
            P[i, i + 1] = (t >= 0.5 ? NaN : 1.0) * u[i + 1]
            P[i + 1, i] = 2 * u[i]
        end
        return nothing
    end
    u0 = [1.0, 1.5, 2.0, 2.5]
    P_tri = Tridiagonal(ones(n - 1), zeros(n), ones(n - 1))
    prob = ConservativePDSProblem(prod!, u0, (0.0, 1.0); p_prototype = P_tri)
    @testset "$alg" for alg in (MPE(; linsolve = nothing),
                                MPRK22(1.0; linsolve = nothing),
                                SSPMPRK43(; linsolve = nothing),
                                MPDeC(3; linsolve = nothing))
        sol = solve(prob, alg; dt = 0.25, adaptive = false)
        @test sol.retcode == SciMLBase.ReturnCode.ConvergenceFailure
        @test sol.t[end] <= 0.5
        @test sol.stats.nreject == 1
        @test all(isfinite, sol.u[end])
    end
end

# A step attempt with a failed linear solve or with non-finite values is rejected.
# A fixed-step solve stops with ReturnCode.ConvergenceFailure; an adaptive solve
# retries with a smaller dt and stops with ReturnCode.DtLessThanMin if dt falls below
# dtmin or with ReturnCode.Unstable if dt becomes too small to advance the time. In all
# cases, the last accepted value is the final value.
@testitem "Rejected step attempts" setup=[Imports, TestHelpers] begin
    # Linear model problem whose production terms are NaN for t > 0, so that
    # every step attempt starting at t = 0 produces non-finite values.
    nan_factor(t) = t > 0 ? NaN : 1.0
    nanP(u, p, t) = nan_factor(t) * [0 u[2]; 5*u[1] 0]
    function nanP!(P, u, p, t)
        P[1, 1] = 0
        P[1, 2] = nan_factor(t) * u[2]
        P[2, 1] = nan_factor(t) * 5 * u[1]
        P[2, 2] = 0
        return nothing
    end
    nanD(u, p, t) = zero(u)
    nanD!(D, u, p, t) = fill!(D, 0)
    u0 = [0.9, 0.1]
    tspan = (0.0, 2.0)
    probs = (ConservativePDSProblem(nanP, u0, tspan),
             ConservativePDSProblem(nanP!, u0, tspan),
             PDSProblem(nanP, nanD, u0, tspan),
             PDSProblem(nanP!, nanD!, u0, tspan),
             ConservativePDSProblem(nanP, SVector{2}(u0), tspan))
    algs = (MPRK22(0.5), MPRK22(1.0), MPRKO22(1.0, 0.75), MPRK43I(1.0, 0.5),
            MPRK43II(0.5), SSPMPRK22(0.5, 1.0), MPDeC(3), MPDeC(6; nodes = :lagrange))
    @testset "$alg" for alg in algs
        @testset "$i" for (i, prob) in enumerate(probs)
            # 7 rejections reduce dt = 0.1 to 0.1 / 2^7 < dtmin
            sol = @test_logs (:warn, r"dtmin") solve(prob, alg; dt = 0.1, dtmin = 1.0e-3)
            @test sol.retcode == SciMLBase.ReturnCode.DtLessThanMin
            @test sol.t == [0.0]
            @test sol.u == [u0]
            @test sol.stats.nreject == 7
            @test sol.stats.naccept == 0

            sol = solve(prob, alg; dt = 0.1, adaptive = false)
            @test sol.retcode == SciMLBase.ReturnCode.ConvergenceFailure
            @test sol.t == [0.0]
            @test sol.u == [u0]
            @test sol.stats.nreject == 1
        end
    end
    # MPE evaluates the production terms only at the beginning of the step.
    @testset "$alg" for (alg, t_end) in ((MPE(), 0.1), (SSPMPRK43(), 0.0))
        @testset "$i" for (i, prob) in enumerate(probs)
            sol = solve(prob, alg; dt = 0.1)
            @test sol.retcode == SciMLBase.ReturnCode.ConvergenceFailure
            @test sol.t[end] == t_end
            @test all(isfinite, sol.u[end])
            @test sol.stats.nreject == 1
        end
    end

    # With the default dtmin = 0, the rejections reduce dt until it is too small to
    # advance the time, which takes more than 1000 rejections at t = 0.
    @testset "Default dtmin $i" for (i, prob) in enumerate(probs[1:2])
        sol = @test_logs (:warn, r"dt was forced below floating point epsilon") solve(prob,
                                                                                      MPRK22(1.0);
                                                                                      dt = 0.1)
        @test sol.retcode == SciMLBase.ReturnCode.Unstable
        @test sol.t == [0.0]
        @test sol.u == [u0]
        @test sol.stats.nreject > 1000
    end

    ODECore = PositiveIntegrators.OrdinaryDiffEqCore
    verbose = ODECore.DEVerbosity(ODECore.SciMLLogging.None())

    # After a rejection, the final value is the last accepted value, also if an
    # earlier attempt was rejected by the error test and only the final value is
    # saved. Solves 1-2 are exact, solve 3 is scaled by 3 (large but finite
    # error estimate), solve 4 is exact, and all later solves give NaN.
    action(k, u) = k == 3 ? (u .*= 3) : (k >= 5 ? fill!(u, NaN) : u)
    u_end = map((prob_pds_linmod_nonconservative,
                 prob_pds_linmod_nonconservative_inplace)) do prob
        sol = solve(prob, MPRK22(1.0; linsolve = CountingLinsolve(action));
                    dt = 0.1, abstol = 0.1, reltol = 0.1, save_everystep = false,
                    verbose)
        @test sol.retcode == SciMLBase.ReturnCode.Unstable
        @test sol.t == [0.0, 0.1]
        sol.u[end]
    end
    sol_ref = solve(prob_pds_linmod_nonconservative, MPRK22(1.0); dt = 0.1,
                    adaptive = false)
    @test all(u -> isapprox(u, sol_ref.u[2]; rtol = 1.0e-12), u_end)

    # The integration continues after rejected attempts: every fifth solve gives NaN,
    # and the integration succeeds.
    @testset "Recovery after rejections $i" for (i, prob) in enumerate((prob_pds_linmod,
                                                                        prob_pds_linmod_inplace))
        linsolve = CountingLinsolve((k, u) -> k % 5 == 0 ? fill!(u, NaN) : u)
        sol = solve(prob, MPRK22(1.0; linsolve); dt = 0.1)
        @test SciMLBase.successful_retcode(sol)
        @test sol.t[end] == last(prob.tspan)
        @test sol.stats.nreject > 0
    end

    # Non-finite scaled error components give an infinite error estimate, which
    # the step size controller rejects; the attempt is not rejected with
    # `force_stepfail`.
    get_EEst = PositiveIntegrators.OrdinaryDiffEqCore.get_EEst
    for prob in (prob_pds_linmod, prob_pds_linmod_inplace)
        integrator = init(prob, MPRK22(1.0); dt = 0.1, verbose)
        for atmp in ([NaN, 0.0], [-Inf, 0.0], SVector(0.0, Inf))
            PositiveIntegrators.set_error_estimate!(integrator, atmp, 0.0)
            @test get_EEst(integrator) == Inf
            @test !integrator.force_stepfail
            @test integrator.stats.nreject == 0
        end
        PositiveIntegrators.set_error_estimate!(integrator, [0.3, 0.4], 0.0)
        @test get_EEst(integrator) ≈ sqrt(0.125)
    end

    # A rejection resets the state to the last accepted value.
    for alg in (MPE(), SSPMPRK43(), MPRK22(1.0))
        for prob in (prob_pds_linmod, prob_pds_linmod_inplace)
            integrator = init(prob, alg; dt = 0.1)
            if SciMLBase.isinplace(prob)
                integrator.u .*= 2
            else
                integrator.u = 2 * integrator.uprev
            end
            @test integrator.u != integrator.uprev
            PositiveIntegrators.reject_step!(integrator)
            @test integrator.force_stepfail
            @test integrator.stats.nreject == 1
            @test integrator.u == integrator.uprev
        end
    end

    # A Patankar weight denominator σ = Inf of the embedded solution (here the stage
    # value u[2]^2 / (0 + floatmin) of MPRK22(0.5) overflows) gives an infinite error
    # estimate. Fixed time steps ignore the error estimate, and adaptive time steps
    # reduce the time step until σ is finite.
    infP(u, p, t) = [0 0; u[1] 0]
    function infP!(P, u, p, t)
        fill!(P, 0)
        P[2, 1] = u[1]
        return nothing
    end
    u0_inf = [1.0e6, 0.0]
    u_exact = [1.0e6 * exp(-1.0), 1.0e6 * (1 - exp(-1.0))]
    probs_inf = (ConservativePDSProblem(infP, u0_inf, (0.0, 1.0)),
                 ConservativePDSProblem(infP!, u0_inf, (0.0, 1.0)))
    @testset "Infinite denominators $i" for (i, prob) in enumerate(probs_inf)
        for alg in (MPRK22(0.5), MPRKO22(0.5, 1.0))
            sol = solve(prob, alg; dt = 0.1, adaptive = false)
            @test sol.retcode == SciMLBase.ReturnCode.Success
            @test sol.stats.nreject == 0
            @test isapprox(sol.u[end], u_exact; rtol = 1.0e-2)

            sol = solve(prob, alg; dt = 0.1)
            @test sol.retcode == SciMLBase.ReturnCode.Success
            @test sol.t[end] == 1.0
            @test isapprox(sol.u[end], u_exact; rtol = 1.0e-3)
        end
    end
end

# Building the MPRK matrix in place of the production matrix (M === P)
# must give the same values as building it in a separate matrix.
@testitem "Sparse MPRK matrix in place of P" setup=[Imports, TestHelpers] begin
    n = 20
    P = sparse([mod(i * j + i, 7) == 0 ? 1 + i / j : 0.0 for i in 1:n, j in 1:n]) +
        spdiagm(0 => ones(n))
    # Nonzero diagonal: the diagonal of the MPRK matrix must not include P[j, j].
    for i in 1:n
        P[i, i] = 0.5 + i / n
    end
    dt = 0.7
    σ = [0.1 + sin(i)^2 for i in 1:n]
    d_vec = [cos(i)^2 for i in 1:n]
    @testset "d = $(d === nothing ? "nothing" : "vector")" for d in (nothing, d_vec)
        M = copy(P)
        PositiveIntegrators.build_mprk_matrix!(M, P, σ, dt, d)
        M_alias = copy(P)
        PositiveIntegrators.build_mprk_matrix!(M_alias, M_alias, σ, dt, d)
        @test rowvals(M_alias) == rowvals(M)
        @test nonzeros(M_alias) == nonzeros(M)
    end

    # Sparse matrices that are not stored in compressed sparse column format use the
    # generic builder, which reports extreme diagonal entries.
    struct WrappedSparse{T} <: SparseArrays.AbstractSparseMatrix{T, Int}
        A::SparseMatrixCSC{T, Int}
    end
    Base.size(W::WrappedSparse) = size(W.A)
    Base.getindex(W::WrappedSparse, i::Int, j::Int) = W.A[i, j]
    Base.setindex!(W::WrappedSparse, v, i::Int, j::Int) = (W.A[i, j] = v; W)
    @testset "Wrapped sparse matrix, d = $(d === nothing ? "nothing" : "vector")" for d in (nothing,
                                                                                            d_vec)
        M_ref = Matrix(P)
        PositiveIntegrators.build_mprk_matrix!(M_ref, Matrix(P), σ, dt, d)
        W = WrappedSparse(copy(P))
        @test PositiveIntegrators.build_mprk_matrix!(W, W, σ, dt, d) === false
        @test Matrix(W.A) ≈ M_ref
    end

    P_nodiag = sparse([0.0 1.0; 2.0 0.0])
    build! = PositiveIntegrators.build_mprk_matrix!
    @test_throws "does not store the diagonal entry in column 1" build!(P_nodiag,
                                                                        P_nodiag,
                                                                        ones(2),
                                                                        0.1)
end

# Here we check that in-place and out-of-place implementations
# deliver the same results
@testitem "Different matrix types (conservative)" setup=[Imports, TestHelpers] begin
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
        MPRKO22(1.0, 0.75),
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

@testitem "Different matrix types (conservative, adaptive)" setup=[Imports, TestHelpers] begin
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

    algs = [MPRK22(0.5), MPRK22(1.0), MPRKO22(1.0, 0.75),
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
@testitem "Different matrix types (nonconservative)" setup=[Imports, TestHelpers] begin
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
        MPRK22(0.5), MPRK22(1.0), MPRKO22(1.0, 0.75),
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

@testitem "Different matrix types (nonconservative, adaptive)" setup=[Imports, TestHelpers] begin
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

    algs = [MPRK22(0.5), MPRK22(1.0), MPRKO22(1.0, 0.75),
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

# Out-of-place production functions returning sparse matrices with
# and without stored diagonal entries must give the same results as
# dense production matrices.
@testitem "Sparse matrices (out-of-place)" setup=[Imports, TestHelpers] begin
    prod_c! = (P, u, p, t) -> begin
        fill!(P, zero(eltype(P)))
        for i in 1:(length(u) - 1)
            P[i, i + 1] = i * u[i + 1]
            P[i + 1, i] = i * u[i]
        end
        return nothing
    end

    # P[1, 1] is zero, the other diagonal entries are positive
    prod_nc! = (P, u, p, t) -> begin
        prod_c!(P, u, p, t)
        for i in 1:length(u)
            P[i, i] = (i - 1) * u[i]
        end
        return nothing
    end
    dest_nc = (u, p, t) -> [i * u[i] for i in 1:length(u)]

    n = 4
    P_sparse = sparse(Tridiagonal([0.1, 0.2, 0.3],
                                  zeros(n),
                                  [0.4, 0.5, 0.6]))
    u0 = [1.0, 1.5, 2.0, 2.5]
    tspan = (0.0, 1.0)
    dt = 0.25

    algs = [MPE(),
        MPRK22(0.5), MPRK22(1.0), MPRKO22(1.0, 0.75),
        MPRK43I(1.0, 0.5), MPRK43I(0.5, 0.75),
        MPRK43II(2.0 / 3.0), MPRK43II(0.5),
        SSPMPRK22(0.5, 1.0), SSPMPRK43(),
        MPDeC(3), MPDeC(3; nodes = :lagrange)]

    @testset "$alg" for alg in algs
        for prod! in (prod_c!, prod_nc!)
            # stored diagonal entries, possibly zero
            prod_stored = (u, p, t) -> begin
                P = copy(P_sparse)
                prod!(P, u, p, t)
                return P
            end
            # no stored zero entries
            prod_dropped = (u, p, t) -> dropzeros!(prod_stored(u, p, t))
            prod_dense = (u, p, t) -> Matrix(prod_stored(u, p, t))

            @test prod_stored(u0, nothing, 0.0) isa SparseMatrixCSC
            @test nnz(diag(prod_stored(u0, nothing, 0.0))) == n
            @test prod_dropped(u0, nothing, 0.0) isa SparseMatrixCSC
            @test nnz(diag(prod_dropped(u0, nothing, 0.0))) < n

            if prod! === prod_c!
                prob_dense = ConservativePDSProblem(prod_dense, u0, tspan)
                prob_stored = ConservativePDSProblem(prod_stored, u0, tspan;
                                                     p_prototype = P_sparse)
                prob_dropped = ConservativePDSProblem(prod_dropped, u0, tspan;
                                                      p_prototype = P_sparse)
            else
                prob_dense = PDSProblem(prod_dense, dest_nc, u0, tspan)
                prob_stored = PDSProblem(prod_stored, dest_nc, u0, tspan;
                                         p_prototype = P_sparse)
                prob_dropped = PDSProblem(prod_dropped, dest_nc, u0, tspan;
                                          p_prototype = P_sparse)
            end

            sol_dense = solve(prob_dense, alg; dt, adaptive = false)
            sol_stored = solve(prob_stored, alg; dt, adaptive = false)
            sol_dropped = solve(prob_dropped, alg; dt, adaptive = false)

            @test SciMLBase.successful_retcode(sol_dense)
            @test SciMLBase.successful_retcode(sol_stored)
            @test SciMLBase.successful_retcode(sol_dropped)
            @test sol_stored.u ≈ sol_dense.u
            @test sol_dropped.u ≈ sol_dense.u
        end
    end
end

# Here we check that the type of p_prototype actually
# defines the types of the Ps inside the algorithm caches.
# We test sparse, tridiagonal, and dense matrices.
@testitem "Prototype type check" setup=[Imports, TestHelpers] begin
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
        MPRK22(0.5), MPRK22(1.0), MPRKO22(1.0, 0.75),
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
@testitem "Convergence tests (conservative, with interpolation)" setup=[
    Imports,
    TestHelpers
] begin
    algs = (MPE(), MPRK22(0.5), MPRK22(1.0), MPRK22(2.0), MPRKO22(1.0, 0.75),
            SSPMPRK22(0.5, 1.0),
            MPDeC(2), MPDeC(2, nodes = :lagrange))
    dts = 0.2 .* 0.5 .^ (6:10)
    problems = (prob_pds_linmod, prob_pds_linmod_array,
                prob_pds_linmod_mvector, prob_pds_linmod_inplace)

    @testset "Problem $i" for (i, prob) in enumerate(problems)
        @testset "Algorithm $j" for (j, alg) in enumerate(algs)
            orders = experimental_orders_of_convergence(prob, alg, dts)
            @test check_order(orders, PositiveIntegrators.alg_order(alg))

            test_times = [
                0.123456789, 1 / pi, exp(-1),
                1.23456789, 1 + 1 / pi, 1 + exp(-1)
            ]
            @testset "Time $k" for (k, test_time) in enumerate(test_times)
                orders = experimental_orders_of_convergence(prob, alg,
                                                            dts;
                                                            test_time,
                                                            use_least_squares = true)
                @test check_order_leastsquares(orders,
                                               PositiveIntegrators.alg_order(alg),
                                               atol = 0.2)
                orders = experimental_orders_of_convergence(prob, alg,
                                                            dts;
                                                            test_time,
                                                            only_first_index = true,
                                                            use_least_squares = true)
                @test check_order_leastsquares(orders,
                                               PositiveIntegrators.alg_order(alg),
                                               atol = 0.2)
            end
        end
    end
end

# Here we check the convergence order of pth-order schemes for which
# also an interpolation of order p is available
@testitem "Convergence tests (nonconservative, with interpolation)" setup=[
    Imports,
    TestHelpers
] begin
    algs = (MPE(), MPRK22(0.5), MPRK22(1.0), MPRK22(2.0), MPRKO22(1.0, 0.75),
            SSPMPRK22(0.5, 1.0),
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
@testitem "Convergence tests (conservative, without interpolation)" setup=[
    Imports,
    TestHelpers
] begin
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
@testitem "Convergence tests (nonconservative, without interpolation)" setup=[
    Imports,
    TestHelpers
] begin
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

@testitem "Interpolation tests (conservative)" setup=[Imports, TestHelpers] begin
    algs = (MPE(), MPRK22(0.5), MPRK22(1.0), MPRK22(2.0), MPRKO22(1.0, 0.75),
            MPRK43I(1.0, 0.5),
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

@testitem "Check convergence order (nonautonomous conservative PDS)" setup=[
    Imports,
    TestHelpers
] begin
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
    algs = (MPE(), MPRK22(0.5), MPRK22(1.0), MPRKO22(1.0, 0.75), MPRK43I(1.0, 0.5),
            MPRK43I(0.5, 0.75),
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

@testitem "Check convergence order (nonautonomous nonconservative PDS)" setup=[
    Imports,
    TestHelpers
] begin
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
    algs = (MPE(), MPRK22(0.5), MPRK22(1.0), MPRKO22(1.0, 0.75), MPRK43I(1.0, 0.5),
            MPRK43I(0.5, 0.75),
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

@testitem "Convergence of higher order schemes (autonomous)" setup=[Imports, TestHelpers] begin
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

@testitem "Convergence of higher order schemes (nonautonomous)" setup=[Imports, TestHelpers] begin
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

@testitem "Interpolation tests (nonconservative)" setup=[Imports, TestHelpers] begin
    algs = (MPE(), MPRK22(0.5), MPRK22(1.0), MPRK22(2.0), MPRKO22(1.0, 0.75),
            MPRK43I(1.0, 0.5),
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
@testitem "Zero initial values" setup=[Imports, TestHelpers] begin
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

    algs = [MPE(), MPRK22(0.5), MPRK22(1.0), MPRK22(2.0), MPRKO22(1.0, 0.75),
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
@testitem "Positivity check" setup=[Imports, TestHelpers] begin
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

    algs = [MPE(), MPRK22(0.5), MPRK22(1.0), MPRK22(2.0), MPRKO22(1.0, 0.75),
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
@testitem "PDS problem library (adaptive schemes)" setup=[Imports, TestHelpers] begin
    algs = [MPRK22(0.5), MPRK22(1.0), MPRKO22(1.0, 0.75), MPRK43I(1.0, 0.5),
        MPRK43I(0.5, 0.75),
        MPRK43II(2.0 / 3.0), MPRK43II(0.5), SSPMPRK22(0.5, 1.0)]
    for k in 2:10
        push!(algs, MPDeC(k), MPDeC(k; nodes = :lagrange))
    end
    probs = (prob_pds_linmod, prob_pds_linmod_inplace, prob_pds_nonlinmod,
             prob_pds_robertson, prob_pds_bertolazzi, prob_pds_brusselator,
             prob_pds_npzd,
             prob_pds_sir, prob_pds_stratreac, prob_pds_jakstat,
             prob_pds_hires)
    @testset "$alg" for alg in algs
        @testset "$i" for (i, prob) in enumerate(probs)
            if prob == prob_pds_stratreac && alg == SSPMPRK22(0.5, 1.0)
                #TODO: SSPMPRK22(0.5, 1.0) is unstable for prob_pds_stratreac.
                #Need to figure out if this is a problem of the algorithm or not.
                continue
            elseif prob == prob_pds_stratreac && alg == MPRK43I(0.5, 0.75)
                # Not successful on Julia 1.9
                continue
            elseif prob == prob_pds_stratreac && alg == MPDeC(9; nodes = :lagrange)
                # unstable
                continue
            elseif prob == prob_pds_jakstat && alg == MPDeC(10; nodes = :lagrange)
                # unstable
                continue
            end
            # later versions of OrdinaryDiffEq.jl use dtmin = 0 by default,
            # see https://github.com/SciML/OrdinaryDiffEq.jl/pull/2098
            sol = solve(prob, alg; dtmin = 0.0)
            @test SciMLBase.successful_retcode(sol)
        end
    end
end

# Here we check that the implemented schemes can solve the predefined PDS.
@testitem "PDS problem library (non-adaptive schemes)" setup=[Imports, TestHelpers] begin
    algs = (MPE(), SSPMPRK43())
    #prob_pds_robertson not included
    probs = (prob_pds_linmod, prob_pds_linmod_inplace, prob_pds_nonlinmod,
             prob_pds_bertolazzi, prob_pds_brusselator,
             prob_pds_npzd, prob_pds_sir, prob_pds_stratreac, prob_pds_jakstat,
             prob_pds_hires)
    @testset "$alg" for alg in algs
        @testset "$prob" for prob in probs
            tspan = prob.tspan
            dt = (tspan[2] - tspan[1]) / 10
            sol = solve(prob, alg; dt = dt)
            @test SciMLBase.successful_retcode(sol)
        end
    end
end

@testitem "HIRES reference solution" setup=[Imports, TestHelpers] begin
    # Reference values at t = 321.8122 from the Test Set for IVP Solvers,
    # https://archimede.uniba.it/~testset/
    u_ref = [0.7371312573325668e-3, 0.1442485726316185e-3,
        0.5888729740967575e-4, 0.1175651343283149e-2,
        0.2386356198831331e-2, 0.6238968252742796e-2,
        0.2849998395185769e-2, 0.2850001604814231e-2]
    prob = prob_pds_hires
    sol = solve(prob, Rodas4P(); abstol = 1e-14, reltol = 1e-12)
    @test sol.u[end]≈u_ref rtol=1e-10
    sol = solve(prob, MPRK43I(1.0, 0.5); abstol = 1e-12, reltol = 1e-10)
    @test sol.u[end]≈u_ref rtol=1e-5
end

#Here we check the different possibilities to define small_constant.
@testitem "Different possibilities to set small_constant" setup=[Imports, TestHelpers] begin
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
        (; kwargs...) -> MPRKO22(1.0, 0.75; kwargs...),
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
@testitem "Exact solutions (RK)" setup=[Imports, TestHelpers] begin
    algs = [MPE(), MPRK22(0.5), MPRK22(1.0), MPRK22(2.0), MPRKO22(1.0, 0.75),
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
@testitem "Exact solutions (MPE)" setup=[Imports, TestHelpers] begin
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
