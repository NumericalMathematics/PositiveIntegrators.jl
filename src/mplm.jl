### Structs and caches ##########################################################################
struct MPLM22{F, T} <: OrdinaryDiffEqAlgorithm
    linsolve::F
    small_constant_function::T
end

alg_order(alg::MPLM22) = 2
alg_extrapolates(alg::MPLM22) = true # TODO: Should probably be false

@cache mutable struct MPLM22ConstantCache{uType, T} <: OrdinaryDiffEqConstantCache
    uprevprev::uType
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
    MPLM22ConstantCache(u, 1, alg.small_constant_function(uEltypeNoUnits))
end

struct MPLM33{F, T} <: OrdinaryDiffEqAlgorithm
    linsolve::F
    small_constant_function::T
end

alg_order(alg::MPLM33) = 3
alg_extrapolates(alg::MPLM33) = true # TODO: Should probably be false

@cache mutable struct MPLM33ConstantCache{uType, PType, dType, T, T2} <:
                      OrdinaryDiffEqConstantCache
    uprevprev::uType
    uprev3::uType
    P2::PType
    P3::PType
    d2::dType
    d3::dType
    α1::T
    α2::T
    α3::T
    β1::T
    β2::T
    β3::T
    step::Int
    small_constant::T2
end

function MPLM33(; linsolve = LUFactorization(), small_constant = nothing)
    if isnothing(small_constant)
        small_constant_function = floatmin
    elseif small_constant isa Number
        small_constant_function = Returns(small_constant)
    else # assume small_constant isa Function
        small_constant_function = small_constant
    end
    MPLM33(linsolve, small_constant_function)
end

function alg_cache(alg::MPLM33, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    α1 = zero(uEltypeNoUnits)
    α2 = zero(uEltypeNoUnits)
    α3 = one(uEltypeNoUnits)
    β1 = 9 / 4 * one(uEltypeNoUnits)
    β2 = zero(uEltypeNoUnits)
    β3 = 3 / 4 * one(uEltypeNoUnits)

    # TODO: This is currently necessary to get the correct type of P (d is of type rateType)
    P, d = evaluate_pds(f, u, p, t)
    # TODO: integrator_stats_nf = 1

    MPLM33ConstantCache(u, u, P, P, d, d, α1, α2, α3, β1, β2, β3, 1,
                        alg.small_constant_function(uEltypeNoUnits))
end

function initialize!(integrator,
                     cache::Union{MPLM22ConstantCache,
                                  MPLM33ConstantCache})
end

### perform_step! ##########################################################################
@muladd function perform_step!(integrator, cache::MPLM22ConstantCache, repeat_step = false)
    (; alg, t, dt, uprev, uprev2, u, f, p) = integrator
    (; uprevprev, small_constant) = cache
    if integrator.u_modified
        cache.step = 1
    end

    if cache.step <= 1
        cache.step += 1
        # Here we perform one step of MPE

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

    integrator.u = u
end

@muladd function perform_step!(integrator, cache::MPLM33ConstantCache, repeat_step = false)
    (; alg, t, dt, uprev, uprev2, f, p) = integrator
    (; uprevprev, uprev3, P2, P3, d2, d3, α1, α2, α3, β1, β2, β3, small_constant) = cache

    #TODO: is this necessary?
    if integrator.u_modified
        cache.step = 1
    end

    if cache.step <= 2 # We perform two steps of MPLM22 to initialize MPLM33
        mplm22_cache = MPLM22ConstantCache(uprevprev, cache.step, small_constant)
        perform_step!(integrator, mplm22_cache, repeat_step)

        #a21, a31, a32, b1, b2, b3, c2, c3, beta1, beta2, q1, q2 = get_constant_parameters(MPRK43I(1.0, 0.5))
        #mprk43_cache = MPRK43ConstantCache(a21, a31, a32, b1, b2, b3, c2, c3,
        #                      beta1, beta2, q1, q2, small_constant)
        #perform_step!(integrator, mprk43_cache, repeat_step)

        #a21, b1, b2 = get_constant_parameters(MPRK22(1.0))
        #mprk22_cache = MPRK22ConstantCache(a21, b1, b2, small_constant)
        #perform_step!(integrator, mprk22_cache, repeat_step) 

        # increase step count
        cache.step += 1

        # This is currently necessary to obtain the correct history of P evaluations
        # Otherwise, MPLM22 would have to cache P and d
        P, d = evaluate_pds(f, uprev, p, t)
        integrator.stats.nf += 1

        u = integrator.u
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

        σ = basic_patankar_step(uprevprev, P, σ, 2 * dt, alg.linsolve, d)
        integrator.stats.nsolve += 1

        # avoid division by zero due to zero Patankar weights
        σ = add_small_constant(σ, small_constant)

        Ptmp, dtmp = lincomb(β1, P, d, β2, P2, d2, β3, P3, d3)
        v = α1 * uprev + α2 * uprevprev + α3 * uprev3
        u = basic_patankar_step(v, Ptmp, σ, dt, alg.linsolve, dtmp)
        integrator.stats.nsolve += 1
    end

    integrator.u = u

    #TODO: Should be possible to use uprev2. But uprev2 is currently not updated.
    #TODO: ConstantCache contains non-constant uprevprev. This is confusing.
    cache.uprev3 = uprevprev
    cache.uprevprev = uprev
    cache.P3 = P2
    cache.P2 = P
    cache.d3 = d2
    cache.d2 = d
end
