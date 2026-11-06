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
# of the same type as the inputs
@testitem "Constant types" setup=[Imports, TestHelpers] begin
    algs = (MPRK22(0.5f0), MPRK22(1.0f0), MPRK22(2.0f0), MPRKO22(1.0f0, 0.75f0),
            MPRK43I(1.0f0, 0.5f0),
            MPRK43I(0.5f0, 0.75f0), MPRK43II(0.5f0), MPRK43II(2.0f0 / 3.0f0),
            SSPMPRK22(0.5f0, 1.0f0))
    @testset "$i" for (i, alg) in enumerate(algs)
        @test eltype(PositiveIntegrators.get_constant_parameters(alg)) == Float32
    end

    algs = (MPRK22(0.5), MPRK22(1.0), MPRK22(2.0), MPRKO22(1.0, 0.75),
            MPRK43I(1.0, 0.5),
            MPRK43I(0.5, 0.75), MPRK43II(0.5), MPRK43II(2.0 / 3.0),
            SSPMPRK22(0.5, 1.0))
    @testset "$i" for (i, alg) in enumerate(algs)
        @test eltype(PositiveIntegrators.get_constant_parameters(alg)) == Float64
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
    for alg in algs
        sol = solve(prob, alg; dt = 0.1f0, save_everystep = false)
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

# A Patankar weight denominator sigma[j] = Inf must give the finite
# limit: M[j, j] = 1 and M[i, j] = 0 for i != j.
@testitem "MPRK matrix with infinite denominators" setup=[Imports, TestHelpers] begin
    P = [0.0 2.0 0.0; 1.5 0.0 3.0; 0.0 0.7 0.0]
    dt = 0.3
    σ = [0.5, Inf, 2.0]
    @testset "d = $d" for d in (nothing, [0.2, 0.0, 1.1])
        s = vec(sum(P; dims = 1)) .+ (d === nothing ? 0.0 : d)
        M_ref = [i == j ? 1 + dt * s[i] / σ[i] : -dt * P[i, j] / σ[j]
                 for i in 1:3, j in 1:3]
        @test M_ref[2, 2] == 1
        @test all(iszero, M_ref[[1, 3], 2])

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

        for M in (M_dense, M_static, M_tri, M_sparse)
            @test Matrix(M) ≈ M_ref
        end
    end
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
