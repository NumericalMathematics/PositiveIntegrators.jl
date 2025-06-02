module SanduProjectionExt

using PositiveIntegrators
using JuMP: Model, @variable, @objective, @constraint, print, objective_value, set_silent,
            optimize!, is_solved_and_feasible, value, set_string_names_on_creation
using SciMLBase: DiscreteCallback

mutable struct SanduProjection{M} <: PositiveIntegrators.SanduProjection
    model::M
    cnt::Int
end

PositiveIntegrators.SanduProjection(args...) = SanduProjection(args...)

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
    @variable(model, z[i = 1:s]>=epsv[i])
    @constraint(model, AT * z.==b)

    affect! = SanduProjection(model, 0)

    return DiscreteCallback(Returns(true), affect!; save_positions = (false, save),
                            finalize = finalize_sandu_projection)
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

        # update minimization problem
        # TODO: Instead of changing the objective we could just replace the coefficients
        # See https://jump.dev/JuMP.jl/stable/api/JuMP/#set_normalized_coefficient
        # and use (5.1) in Sandu's Paper
        @objective(model, Min, 1 / 2*sum(g .* (model[:z] - u) .^ 2))

        # solve optimization problem
        optimize!(model)
        if !is_solved_and_feasible(model)
            error("Solver did not find an optimal solution")
        end

        integrator.u = value.(model[:z])
    end
    return nothing
end

end
