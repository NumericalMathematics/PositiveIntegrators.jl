module SanduProjectionExt

using StaticArrays: StaticArray, SVector # Why do we need this here? 
using JuMP: @variable, @objective, @constraint, print, set_silent,
            optimize!, is_solved_and_feasible, value, set_string_names_on_creation,
            set_objective_coefficient, @expression
using SciMLBase: DiscreteCallback
using PositiveIntegrators

mutable struct SanduProjection{M} <: PositiveIntegrators.SanduProjection
    model::M
    cnt::Int
end

"""
    SanduProjection(model, AT, b, eps = nothing; [save = true, verbose = false])

A projection method which ensures conservation of prescribed linear invariants and positivity.
If the current approximation ``\\mathbf{u}`` has negative components then a projection ``\\mathbf{z}`` is computed such that
```math
\\min \\lVert \\mathbf{z} - \\mathbf{u} \\rVert_G,\\quad \\mathbf{A}^T\\mathbf{z}=\\mathbf{b},\\quad \\mathbf{z}≥ \\mathbf{0},
```
is satisfied, where the matrix ``\\mathbf{A^T}`` and the vector ``\\mathbf{b}`` define the linear invariants. 
See Sandu (2001) for details.

To use this callback one must also specify `save_everystep = false`.

## Arguments

  - `model`: A [`JuMP Model`](https://jump.dev/JuMP.jl/stable/api/JuMP/#JuMP.Model) to solve the minimization problem.
  - `AT`: The matrix ``\\mathbf{A}^T`` defining the linear invariants.
  - `b`: The vector ``\\mathbf{b}`` defining the linear invariants.
  - `eps`: It may be helpful for the optimization solver that feasible solutions are bounded away from 0. To achieve this one can specify the optional parameter `eps`. The positivity constraint is then replaced by ``\\mathbf{z}≥```eps`, where `eps` can either be a scalar or a vector. 

## Keyword Arguments

  - `save`: If the keyword argument `save` is set to `false` only the initial value and the last approximation will be saved.
            The default value is `true`.
  - `verbose`: Enables additional output of the optimization solver. The default value is `false`.

## References

- Adrian Sandu.
  "Positive numerical integration methods for chemical kinetic systems."
  Journal of Computational Physics 170 (2001): 589-602.
  [DOI: 10.1006/jcph.2001.6750](https://doi.org/10.1006/jcph.2001.6750)
"""
function PositiveIntegrators.SanduProjection(args...; kwargs...)
    SanduProjection(args...; kwargs...)
end

function SanduProjection(model, AT, b, eps = nothing; save = true, verbose = false)
    if isnothing(eps) || eps isa Number
        epsv = zeros(eltype(AT), size(AT, 2))
        if eps isa Number
            fill!(epsv, eps)
        end
    else
        epsv = eps
    end

    # Set up optimization problem
    s = size(AT, 2)
    if !verbose
        set_silent(model)
    end
    set_string_names_on_creation(model, false)

    @variable(model, z[i = 1:s]>=epsv[i])
    @constraint(model, AT * z.==b)
    # This just initializes the objective. The correct coefficients will be set later
    @expression(model, obj_exp, sum(z .^ 2)+sum(z))
    @objective(model, Min, obj_exp)

    if verbose
        print(model)
    end

    affect! = SanduProjection(model, 0)

    return DiscreteCallback(Returns(true), affect!; save_positions = (false, save),
                            finalize = finalize_sandu_projection,
                            initialize = initialize_sandu_projection)
end

function initialize_sandu_projection(c, u, t, integrator)
    return initialize_sandu_projection(c.affect!)
end

function initialize_sandu_projection(proj::SanduProjection)
    proj.cnt = 0
end

function finalize_sandu_projection(c, u, t, integrator)
    return finalize_sandu_projection(c.affect!)
end

function finalize_sandu_projection(proj::SanduProjection)
end

function (proj::SanduProjection)(integrator)
    u = integrator.u

    if isnegative(u)
        proj.cnt += 1

        rtol = integrator.opts.reltol
        atol = integrator.opts.abstol
        model = proj.model

        s = length(u)
        g = @. 1 / (s * (atol + rtol * abs(u))^2)

        # set coefficients of quadratic terms
        set_objective_coefficient(model, model[:z], model[:z], Vector(1 / 2 .* g))
        # set coefficients of linear terms
        set_objective_coefficient(model, model[:z], Vector(-g .* u))

        # solve optimization problem
        optimize!(model)
        if !is_solved_and_feasible(model)
            error("Solver did not find an optimal solution")
        end

        if integrator.u isa StaticArray
            integrator.u = SVector{length(integrator.u)}(value.(model[:z]))
        else
            integrator.u = value.(model[:z])
        end
    end
    return nothing
end

"""
    get_SanduProjection_steps(proj)

If `proj` is a `SanduProjection`, this function returns the number of the performed projection steps.

"""
function PositiveIntegrators.get_numsteps_SanduProjection(proj)
    return proj.affect!.cnt
end

end
