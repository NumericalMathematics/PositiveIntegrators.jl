@testitem "Sandu projection is positive" setup=[Imports, TestHelpers] begin
    # standard arrays
    sol = solve(prob_pds_linmod_array, Euler(); dt = 0.25)
    @test isnegative(sol)

    cb = SanduProjection(Model(Clarabel.Optimizer), [1.0 1.0], [1.0])
    sol_cb = solve(prob_pds_linmod_array, Euler(); dt = 0.25,
                   save_everystep = false,
                   callback = cb)
    @test isnonnegative(sol_cb)

    # static arrays
    sol = solve(prob_ode_stratreac_scaled, ROS2())
    @test isnegative(sol)

    AT = linear_invariants_stratreac_scaled()
    b = AT * prob_ode_stratreac_scaled.u0

    cb = SanduProjection(Model(Clarabel.Optimizer), AT, b)
    sol_cb = solve(prob_ode_stratreac_scaled, ROS2(); save_everystep = false,
                   callback = cb)
    @test isnonnegative(sol_cb)

    cb = SanduProjection(Model(Clarabel.Optimizer), AT, b; verbose = true)
    sol_cb = solve(prob_ode_stratreac_scaled, ROS2(); save_everystep = false,
                   callback = cb)
    @test isnonnegative(sol_cb)

    cb2 = SanduProjection(Model(Clarabel.Optimizer), AT, b, 0.1)
    sol_cb2 = solve(prob_ode_stratreac_scaled, ROS2(); save_everystep = false,
                    callback = cb2)
    @test isnonnegative(sol_cb2)

    cb3 = SanduProjection(Model(Clarabel.Optimizer), AT, b, 0.1 * ones(6))
    sol_cb3 = solve(prob_ode_stratreac_scaled, ROS2(); save_everystep = false,
                    callback = cb3)
    @test isnonnegative(sol_cb3)
end

@testitem "Sandu projection save" setup=[Imports, TestHelpers] begin
    AT = linear_invariants_stratreac_scaled()
    b = AT * prob_ode_stratreac_scaled.u0
    cb = SanduProjection(Model(Clarabel.Optimizer), AT, b; save = false)
    sol_cb = solve(prob_ode_stratreac_scaled, ROS2(); save_everystep = false,
                   callback = cb)

    @test length(sol_cb.t) == 2
end

@testitem "Sandu projection get_numsteps_SanduProjection" setup=[Imports, TestHelpers] begin
    AT = linear_invariants_stratreac_scaled()
    b = AT * prob_ode_stratreac_scaled.u0
    cb = SanduProjection(Model(Clarabel.Optimizer), AT, b; save = false)
    sol_cb = solve(prob_ode_stratreac_scaled, ROS2(); save_everystep = false,
                   callback = cb)

    @test get_numsteps_SanduProjection(cb) > 0
end
