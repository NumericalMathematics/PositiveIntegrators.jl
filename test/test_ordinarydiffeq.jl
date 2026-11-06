# Here we check that solutions of equivalent ODEProblems, PDSProblems or
# ConservativePDS Problems are approximately equal.
# We also check that solvers from OrdinaryDiffEq can solve PDSProblems and
# ConservativePDSProblems.

@testitem "Linear model" setup=[Imports, TestHelpers] begin
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
    # in-place syntax for PDS
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

@testitem "Lotka-Volterra" setup=[Imports, TestHelpers] begin
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
    # in-place syntax for PDS
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

@testitem "Linear advection" setup=[Imports, TestHelpers] begin
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
# This is required when solvers from OrdinaryDiffEq are used to solve a PDS
# and std_rhs is not specified.
@testitem "(Conservative)PDSFunction and units" setup=[Imports, TestHelpers] begin
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

@testitem "Compatibility of OrdinaryDiffEq solvers with (Conservative)PDSProblems and units" setup=[
    Imports,
    TestHelpers
] begin
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

# Here we check that production-destruction form and standard ODE form fit together in predefined problems,
# i.e. standard solvers using std_rhs should generate results that are equal to those without specifying std_rhs
@testitem "Check that production-destruction form and standard ODE fit together in predefined problems" setup=[
    Imports,
    TestHelpers
] begin
    # non-stiff conservative problems (out-of-place)
    probs = (prob_pds_linmod, prob_pds_nonlinmod, prob_pds_brusselator,
             prob_pds_sir, prob_pds_npzd)
    algs = (Euler(), ImplicitEuler(), Tsit5(), Rosenbrock23(), SDIRK2(), TRBDF2())
    abstol = 1e-5
    reltol = 1e-5
    @testset "$alg, $prob" for prob in probs, alg in algs
        dt = (last(prob.tspan) - first(prob.tspan)) / 1e4
        sol = solve(prob, alg; dt, abstol, reltol,
                    isoutofdomain = isnegative) # use explicit f
        sol2 = solve(ConservativePDSProblem(prob.f.p, prob.u0, prob.tspan), alg; dt,
                     abstol, reltol, isoutofdomain = isnegative) # use p and d to compute f
        sol3 = solve(ODEProblem(prob.f.std_rhs, prob.u0, prob.tspan), alg; dt,
                     abstol, reltol, isoutofdomain = isnegative) # use f to create ODEProblem
        t = range(prob.tspan..., length = 10_000)
        @test sol.(t) ≈ sol2.(t) ≈ sol3.(t)
    end

    # non-stiff conservative problems (in-place)
    # Requires autodiff = AutoFiniteDiff()
    probs = (prob_pds_linmod_inplace,)
    algs = (Euler(), ImplicitEuler(autodiff = AutoFiniteDiff()), Tsit5(),
            Rosenbrock23(autodiff = AutoFiniteDiff()),
            SDIRK2(autodiff = AutoFiniteDiff()),
            TRBDF2(autodiff = AutoFiniteDiff()))
    @testset "$alg, $prob" for prob in probs, alg in algs
        dt = (last(prob.tspan) - first(prob.tspan)) / 1e4
        sol = solve(prob, alg; dt, abstol, reltol,
                    isoutofdomain = isnegative) # use explicit f
        sol2 = solve(ConservativePDSProblem(prob.f.p, prob.u0, prob.tspan), alg; dt,
                     abstol, reltol, isoutofdomain = isnegative) # use p and d to compute f
        sol3 = solve(ODEProblem(prob.f.std_rhs, prob.u0, prob.tspan), alg; dt,
                     abstol, reltol, isoutofdomain = isnegative) # use f to create ODEProblem
        t = range(prob.tspan..., length = 10_000)
        @test sol.(t) ≈ sol2.(t) ≈ sol3.(t)
    end

    # non-stiff non-conservative problems (out-of-place)
    probs = (prob_pds_minmapk,)
    algs = (Euler(), ImplicitEuler(), Tsit5(), Rosenbrock23(), SDIRK2(), TRBDF2())
    @testset "$alg, $prob" for prob in probs, alg in algs
        dt = (last(prob.tspan) - first(prob.tspan)) / 1e4
        sol = solve(prob, alg; dt, abstol, reltol, isoutofdomain = isnegative) # use explicit f
        sol2 = solve(PDSProblem(prob.f.p, prob.f.d, prob.u0, prob.tspan), alg; dt,
                     abstol, reltol, isoutofdomain = isnegative) # use p and d to compute f
        sol3 = solve(ODEProblem(prob.f.std_rhs, prob.u0, prob.tspan), alg; dt,
                     abstol, reltol, isoutofdomain = isnegative) # use f to create ODEProblem
        t = range(prob.tspan..., length = 10_000)
        @test sol.(t) ≈ sol2.(t) ≈ sol3.(t)
    end

    # Robertson problem
    prob = prob_pds_robertson
    algs = (ImplicitEuler(), Rosenbrock23(), SDIRK2(), TRBDF2())
    @testset "$alg" for alg in algs
        dt = 1e-6
        sol = solve(prob, alg; dt, abstol, reltol) # use explicit f
        sol2 = solve(ConservativePDSProblem(prob.f.p, prob.u0, prob.tspan), alg; dt,
                     abstol, reltol) # use p and d to compute f
        sol3 = solve(ODEProblem(prob.f.std_rhs, prob.u0, prob.tspan), alg; dt,
                     abstol, reltol) # use f to create ODEProblem
        t = range(prob.tspan..., length = 10_000)
        @test sol.(t) ≈ sol2.(t) ≈ sol3.(t)
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
        sol = solve(prob, alg; dt, abstol, reltol) # use explicit f
        sol2 = solve(PDSProblem(prob.f.p, prob.f.d, prob.u0, prob.tspan), alg; dt,
                     abstol, reltol) # use p and d to compute f
        sol3 = solve(ODEProblem(prob.f.std_rhs, prob.u0, prob.tspan), alg; dt,
                     abstol, reltol) # use f to create ODEProblem
        t = range(prob.tspan..., length = 10_000)
        @test sol.(t) ≈ sol2.(t) ≈ sol3.(t)
    end

    # HIRES problem
    prob = prob_pds_hires
    algs = (ImplicitEuler(), Rosenbrock23(), TRBDF2())
    @testset "$alg" for alg in algs
        dt = 1e-3
        sol = solve(prob, alg; dt, abstol, reltol) # use explicit f
        sol2 = solve(PDSProblem(prob.f.p, prob.f.d, prob.u0, prob.tspan), alg; dt,
                     abstol, reltol) # use p and d to compute f
        sol3 = solve(ODEProblem(prob.f.std_rhs, prob.u0, prob.tspan), alg; dt,
                     abstol, reltol) # use f to create ODEProblem
        t = range(prob.tspan..., length = 10_000)
        @test sol.(t) ≈ sol2.(t) ≈ sol3.(t)
    end
end

@testitem "getproperty and SciML metadata delegation" setup=[Imports, TestHelpers] begin
    # Valid dummy PDS matrix functions: P(P, u, p, t) and D(d, u, p, t)
    dummy_P(P, u, p, t) = (P .= 0.0)
    dummy_D(d, u, p, t) = (d .= 0.0)
    dummy_rhs(du, u, p, t) = (du .= u)

    dummy_u0 = [1.0, 2.0]
    dummy_tspan = (0.0, 1.0)
    dummy_p = nothing
    dummy_proto = [1.0 0.0; 0.0 1.0]
    dummy_sparsity = [true false; false true]
    dummy_color = [1, 2]

    # -----------------------------------------------------------------
    # 1. Tests for PDSFunction
    # -----------------------------------------------------------------
    @testset "PDSFunction" begin
        # Plain function
        prob_plain = PDSProblem(dummy_P, dummy_D, dummy_u0, dummy_tspan, dummy_p;
                                std_rhs = dummy_rhs)
        f_plain = prob_plain.f

        @test f_plain.mass_matrix === I
        @test f_plain.std_rhs === dummy_rhs
        @test f_plain.jac_prototype === nothing
        @test f_plain.sparsity === nothing
        @test_throws Exception f_plain.non_existent_field_xyz

        # ODEFunction with jac_prototype
        ode_proto = SciMLBase.ODEFunction(dummy_rhs; jac_prototype = dummy_proto)
        prob_proto = PDSProblem(dummy_P, dummy_D, dummy_u0, dummy_tspan, dummy_p;
                                std_rhs = ode_proto)
        f_proto = prob_proto.f

        @test f_proto.mass_matrix === I
        @test f_proto.jac_prototype === dummy_proto
        @test f_proto.sparsity === dummy_proto  # Cascades from jac_prototype

        # ODEFunction with explicit sparsity
        ode_sparsity = SciMLBase.ODEFunction(dummy_rhs; jac_prototype = dummy_proto,
                                             sparsity = dummy_sparsity)
        prob_sparsity = PDSProblem(dummy_P, dummy_D, dummy_u0, dummy_tspan, dummy_p;
                                   std_rhs = ode_sparsity)
        f_sparsity = prob_sparsity.f

        @test f_sparsity.sparsity === dummy_sparsity

        # ODEFunction with general metadata (colorvec)
        ode_color = SciMLBase.ODEFunction(dummy_rhs; colorvec = dummy_color)
        prob_color = PDSProblem(dummy_P, dummy_D, dummy_u0, dummy_tspan, dummy_p;
                                std_rhs = ode_color)
        f_color = prob_color.f

        @test f_color.colorvec === dummy_color
    end

    # -----------------------------------------------------------------
    # 2. Tests for ConservativePDSFunction
    # -----------------------------------------------------------------
    @testset "ConservativePDSFunction" begin
        # Plain function
        prob_plain = ConservativePDSProblem(dummy_P, dummy_u0, dummy_tspan, dummy_p;
                                            std_rhs = dummy_rhs)
        f_plain = prob_plain.f

        @test f_plain.mass_matrix === I
        @test f_plain.std_rhs === dummy_rhs
        @test f_plain.jac_prototype === nothing
        @test f_plain.sparsity === nothing
        @test_throws Exception f_plain.non_existent_field_xyz

        # ODEFunction with jac_prototype
        ode_proto = SciMLBase.ODEFunction(dummy_rhs; jac_prototype = dummy_proto)
        prob_proto = ConservativePDSProblem(dummy_P, dummy_u0, dummy_tspan, dummy_p;
                                            std_rhs = ode_proto)
        f_proto = prob_proto.f

        @test f_proto.mass_matrix === I
        @test f_proto.jac_prototype === dummy_proto
        @test f_proto.sparsity === dummy_proto  # Cascades from jac_prototype

        # ODEFunction with explicit sparsity
        ode_sparsity = SciMLBase.ODEFunction(dummy_rhs; jac_prototype = dummy_proto,
                                             sparsity = dummy_sparsity)
        prob_sparsity = ConservativePDSProblem(dummy_P, dummy_u0, dummy_tspan,
                                               dummy_p;
                                               std_rhs = ode_sparsity)
        f_sparsity = prob_sparsity.f

        @test f_sparsity.sparsity === dummy_sparsity

        # ODEFunction with general metadata (colorvec)
        ode_color = SciMLBase.ODEFunction(dummy_rhs; colorvec = dummy_color)
        prob_color = ConservativePDSProblem(dummy_P, dummy_u0, dummy_tspan, dummy_p;
                                            std_rhs = ode_color)
        f_color = prob_color.f

        @test f_color.colorvec === dummy_color
    end

    # -----------------------------------------------------------------
    # 3. Explicit coverage for sparsity fallback branch (val === nothing)
    # -----------------------------------------------------------------
    @testset "Explicit sparsity fallback branch" begin
        # Dummy subtype of AbstractODEFunction where getproperty(f, :sparsity) returns nothing directly
        struct DummyODEFunc{F} <: SciMLBase.AbstractODEFunction{false}
            f::F
            jac_prototype::Any
            sparsity::Nothing
        end

        dummy_ode = DummyODEFunc(dummy_rhs, dummy_proto, nothing)

        f_pds = PDSProblem(dummy_P, dummy_D, dummy_u0, dummy_tspan, dummy_p;
                           std_rhs = dummy_ode).f
        @test f_pds.sparsity === dummy_proto

        f_cpds = ConservativePDSProblem(dummy_P, dummy_u0, dummy_tspan, dummy_p;
                                        std_rhs = dummy_ode).f
        @test f_cpds.sparsity === dummy_proto
    end
end
