@testitem "plot" setup=[Imports, TestHelpers] begin
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

@testitem "is(non)negative" setup=[Imports, TestHelpers] begin
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

@testitem "errors" setup=[Imports, TestHelpers] begin
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
@testitem "work-precision fixed" setup=[Imports, TestHelpers] begin
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
@testitem "work-precision adaptive" setup=[Imports, TestHelpers] begin
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
