module SanduProjectionExt

using StaticArrays: StaticArray, SVector # Why do we need this here? 
using JuMP: @variable, @objective, @constraint, print, set_silent,
            optimize!, is_solved_and_feasible, value, set_string_names_on_creation
using SciMLBase: DiscreteCallback
using PositiveIntegrators

#import PositiveIntegrators: SanduProjection
#Base.retry_load_extensions()

mutable struct SanduProjection{M} <: PositiveIntegrators.SanduProjection
    model::M
    cnt::Int
end

"""
    SanduProjection(model, AT, b, eps = nothing; [save = true])

A projection method which ensures conservation of prescribed linear invariants and positivity.

Given an approximation ``\\mathbf{u}`` the projection ``\\mathbf{z}`` is computed as 
```math
\\min \\lVert \\mathbf{z} - \\mathbf{u} \\rVert_G,\\quad \\mathbf{A}^T\\mathbf{z}=\\mathbf{b},\\quad \\mathbf{z}≥ \\mathbf{0},
```
where `model` is used to solve the optimization problem, `AT` refers to ``\\math{A}^T``, and `b` refers to ``\\mathbf{b}``.
See Sandu (2001) for details.

The projection is implemented as a [`DiscreteCallback`](https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions/#SciMLBase.DiscreteCallback).
To use this callback one must also specify `save_everystep = false`.

To avoid negative elements of ``\\mathbf{z}`` due to roundoff one can specify the optional parameter `eps`.
The positivity constraint is then replaced by ``\\mathbf{z}≥\\mathrm{eps}``, where `eps` can either be a scalar or a vector.

If the keyword argument `save` is set to `false` only the initial value and the last approximation will be saved.
The default value is `true`.

## References

- Adrian Sandu
  "Positive numerical integration methods for chemical kinetic systems."
  Journal of Computational Physics 170 (2001): 589-602.
  [DOI: 10.1006/jcph.2001.6750](https://doi.org/10.1006/jcph.2001.6750)
"""
function PositiveIntegrators.SanduProjection(args...; kwargs...)
    SanduProjection(args...; kwargs...)
end

function SanduProjection(model, AT, b, eps = nothing; save = true)
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
    set_silent(model)
    set_string_names_on_creation(model, false)

    #TODO: This does not work as intended. z may be less than epsv!
    @variable(model, z[i = 1:s]>=epsv[i])

    @constraint(model, AT * z.==b)
    #print(model)

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
    print("\nNumber of Sandu-Projection steps: $(proj.cnt)\n\n")
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

        # update of minimization problem
        # TODO: Instead of changing the objective we could just replace the coefficients
        # See https://jump.dev/JuMP.jl/stable/api/JuMP/#set_normalized_coefficient
        # and use (5.1) in Sandu's Paper
        @objective(model, Min, 1 / 2*sum(g .* (model[:z] - u) .^ 2))
        #@objective(model, Min, 1/2 * sum(g .* model[:z].^2) - sum(g .* u .* model[:z]) )
        #print(model)

        # solve optimization problem
        optimize!(model)
        if !is_solved_and_feasible(model)
            error("Solver did not find an optimal solution")
        end

        if integrator.u isa StaticArray
            # Tried to use @SVector, but failed :(. WHY?
            integrator.u = SVector{length(integrator.u)}(value.(model[:z]))
        else
            integrator.u = value.(model[:z])
        end
    end
    return nothing
end

end
