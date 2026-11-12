# Run selected `@testset`s of `test/runtests.jl` by name, e.g., in a persistent Julia
# session using the project `run_agents/` (see `AGENTS.md`):
#   include("/path/to/PositiveIntegrators.jl/utils/run_testsets.jl")
#   run_testsets("Float32 input data", "HIRES reference solution")
# Each call parses `test/runtests.jl` again and evaluates the setup code before the
# outermost `@testset` (imports and helper functions) and all `@testset`s whose name is
# exactly one of the given strings (wherever they are nested) in a new module. Testsets
# with interpolated names such as `"$alg"` cannot be selected. A selected testset must not
# depend on variables defined in an enclosing testset; select the enclosing one instead.

using Test: Test, @testset

const RUNTESTS_FILE = joinpath(dirname(@__DIR__), "test", "runtests.jl")

is_testset(ex) = Meta.isexpr(ex, :macrocall) && ex.args[1] === Symbol("@testset")

# The parser wraps a top-level expression ending with `;` in an `Expr(:toplevel, ...)`
function contains_testset(ex)
    is_testset(ex) ||
        (Meta.isexpr(ex, :toplevel) && any(contains_testset, ex.args))
end

function find_testsets!(found, ex, names)
    if is_testset(ex) && any(arg -> arg isa String && arg in names, ex.args)
        push!(found, ex)
    elseif ex isa Expr
        foreach(arg -> find_testsets!(found, arg, names), ex.args)
    end
    return found
end

function run_testsets(names::AbstractString...)
    isempty(names) && throw(ArgumentError("pass the names of the testsets to run"))
    code = Meta.parseall(read(RUNTESTS_FILE, String); filename = RUNTESTS_FILE)
    index_main = findfirst(contains_testset, code.args)
    setup = code.args[1:(index_main - 1)]
    found = find_testsets!(Expr[], code.args[index_main], names)
    found_names = [only(filter(arg -> arg isa String, ex.args)) for ex in found]
    missing_names = setdiff(collect(names), found_names)
    isempty(missing_names) ||
        throw(ArgumentError("testsets not found in $RUNTESTS_FILE: $missing_names"))

    mod = Module(:SelectedTestsets)
    foreach(ex -> Core.eval(mod, ex), setup)
    Core.eval(mod, :(using Test: @testset))
    println("Running ", length(found), " testset(s): ", found_names)
    return Core.eval(mod, :(@testset "Selected testsets" begin
                                $(found...)
                            end))
end
