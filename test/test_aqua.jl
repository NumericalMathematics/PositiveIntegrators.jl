@testitem "Aqua.jl" begin
    using Aqua: Aqua
    using RecipesBase: RecipesBase

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
        # `isnegative` is part of Base since Julia 1.13, but
        # PositiveIntegrators extends it with methods for vectors
        # and ODE solutions.
        treat_as_own = Function[RecipesBase.apply_recipe]
        if isdefined(Base, :isnegative)
            push!(treat_as_own, Base.isnegative)
        end
        Aqua.test_all(PositiveIntegrators;
                      ambiguities = ambiguities,
                      piracies = (; treat_as_own = treat_as_own),
                      persistent_tasks = persistent_tasks)
    end
end

@testitem "ExplicitImports.jl" begin
    using ExplicitImports: check_no_implicit_imports, check_no_stale_explicit_imports

    @test isnothing(check_no_implicit_imports(PositiveIntegrators))
    @test isnothing(check_no_stale_explicit_imports(PositiveIntegrators))
end
