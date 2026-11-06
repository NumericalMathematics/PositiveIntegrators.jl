# Run selected `@testitem`s of the test suite in `test/` by name, e.g., in a persistent
# Julia session using the project `run_agents/` (see `AGENTS.md`):
#   include("/path/to/PositiveIntegrators.jl/utils/run_testitems.jl")
#   run_testitems("Float32 input data", "HIRES reference solution")
# The names must be exactly the names of `@testitem`s in `test/test_*.jl`. Nested
# `@testset`s cannot be selected; select the enclosing test item instead. Each call
# evaluates the test setups (`test/test_setup.jl`) again.

using TestItemRunner: TestItemRunner

const REPO_DIR = dirname(@__DIR__)

# Names of all test items, read from the lines starting with `@testitem "..."`
function testitem_names()
    names = String[]
    for file in readdir(joinpath(REPO_DIR, "test"); join = true)
        endswith(file, ".jl") || continue
        for line in eachline(file)
            m = match(r"^@testitem \"([^\"]*)\"", line)
            isnothing(m) || push!(names, m.captures[1])
        end
    end
    return names
end

function run_testitems(names::AbstractString...; verbose = true)
    isempty(names) && throw(ArgumentError("pass the names of the test items to run"))
    missing_names = setdiff(names, testitem_names())
    isempty(missing_names) ||
        throw(ArgumentError("test items not found in $(joinpath(REPO_DIR, "test")): " *
                            "$missing_names"))

    println("Running ", length(names), " test item(s): ", collect(names))
    return TestItemRunner.run_tests(REPO_DIR; filter = ti -> ti.name in names, verbose)
end
