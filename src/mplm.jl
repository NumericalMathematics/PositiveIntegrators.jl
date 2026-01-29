struct MPLM22{F, T} <: OrdinaryDiffEqAlgorithm
    linsolve::F
    small_constant_function::T
end

alg_order(alg::MPLM22) = 2
alg_extrapolates(alg::MPLM22) = true

@cache mutable struct MPLM22ConstantCache{uType, rateType, T} <: OrdinaryDiffEqConstantCache
    uprevprev::uType
    k2::rateType
    step::Int
    small_constant::T
end

function MPLM22(; linsolve = LUFactorization(), small_constant = nothing)
    if isnothing(small_constant)
        small_constant_function = floatmin
    elseif small_constant isa Number
        small_constant_function = Returns(small_constant)
    else # assume small_constant isa Function
        small_constant_function = small_constant
    end
    MPLM22(linsolve, small_constant_function)
end

function alg_cache(alg::MPLM22, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    k2 = rate_prototype
    MPLM22ConstantCache(u, k2, 1, alg.small_constant_function(uEltypeNoUnits))
end

function initialize!(integrator,
                     cache::Union{MPLM22ConstantCache})
    #=
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal

    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    =#
end

@muladd function perform_step!(integrator, cache::MPLM22ConstantCache, repeat_step = false)
    (; alg, t, dt, uprev, uprev2, u, f, p) = integrator
    (; uprevprev, k2, small_constant) = cache
    k1 = integrator.fsalfirst
    if integrator.u_modified
        cache.step = 1
    end

    if cache.step <= 1
        cache.step += 1

        # evaluate production matrix
        P, d = evaluate_pds(f, uprev, p, t)
        integrator.stats.nf += 1

        # avoid division by zero due to zero Patankar weights
        σ = add_small_constant(uprev, small_constant)

        u = basic_patankar_step(uprev, P, σ, dt, alg.linsolve, d)
        integrator.stats.nsolve += 1
    else
        # evaluate production matrix
        P, d = evaluate_pds(f, uprev, p, t)
        integrator.stats.nf += 1

        # avoid division by zero due to zero Patankar weights
        σ = add_small_constant(uprev, small_constant)

        σ = basic_patankar_step(uprev, P, σ, dt, alg.linsolve, d)
        integrator.stats.nsolve += 1

        # avoid division by zero due to zero Patankar weights
        σ = add_small_constant(σ, small_constant)

        u = basic_patankar_step(uprevprev, P, σ, 2 * dt, alg.linsolve, d)
        integrator.stats.nsolve += 1
    end

    #TODO: Should be possible to use uprev2. But uprev2 is currently not updated.
    #TODO: ConstantCache contains non-constant uprevprev. This is confusing.
    cache.uprevprev = uprev

    #=
    #TODO: Don't need k1 and k2 and initialize! must be changed accordingly.
    k2 = k1
    cache.k2 = k2
    integrator.fsallast = f(u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    =#
    integrator.u = u
end
