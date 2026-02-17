### Structs and caches ##########################################################################
struct MPLM22{F, T} <: OrdinaryDiffEqAlgorithm
    linsolve::F
    small_constant_function::T
end

alg_order(alg::MPLM22) = 2
alg_extrapolates(alg::MPLM22) = true # TODO: Should probably be false

@cache mutable struct MPLM22oopCache{uType, T} <: OrdinaryDiffEqConstantCache
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
    MPLM22oopCache(u, 1, alg.small_constant_function(uEltypeNoUnits))
end

struct MPLM33{F, T} <: OrdinaryDiffEqAlgorithm
    linsolve::F
    small_constant_function::T
end

alg_order(alg::MPLM33) = 3
alg_extrapolates(alg::MPLM33) = true # TODO: Should probably be false

@cache mutable struct MPLM33oopCache{uType, PType, dType, T, T2} <:
                      OrdinaryDiffEqConstantCache
    uprevprev::uType
    uprev3::uType
    P2::PType
    P3::PType
    d2::dType
    d3::dType
    αβ::NTuple{6, T}
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
    αβ = (α1, α2, α3, β1, β2, β3)

    # TODO: This is currently necessary to get the correct type of P (d is of type rateType)
    P, d = evaluate_pds(f, u, p, t)
    # TODO: integrator_stats_nf = 1

    MPLM33oopCache(u, u, P, P, d, d, αβ, 1, alg.small_constant_function(uEltypeNoUnits))
end

struct MPLM43{F, T} <: OrdinaryDiffEqAlgorithm
    linsolve::F
    small_constant_function::T
end

alg_order(alg::MPLM43) = 3
alg_extrapolates(alg::MPLM43) = true # TODO: Should probably be false

@cache mutable struct MPLM43oopCache{uType, PType, dType, T, T2} <:
                      OrdinaryDiffEqConstantCache
    uprevprev::uType
    uprev3::uType
    uprev4::uType
    P2::PType
    P3::PType
    P4::PType
    d2::dType
    d3::dType
    d4::dType
    αβ::NTuple{8, T}
    step::Int
    small_constant::T2
end

function MPLM43(; linsolve = LUFactorization(), small_constant = nothing)
    if isnothing(small_constant)
        small_constant_function = floatmin
    elseif small_constant isa Number
        small_constant_function = Returns(small_constant)
    else # assume small_constant isa Function
        small_constant_function = small_constant
    end
    MPLM43(linsolve, small_constant_function)
end

function alg_cache(alg::MPLM43, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    # TODO: This is currently necessary to get the correct type of P (d is of type rateType)
    P, d = evaluate_pds(f, u, p, t)
    # TODO: integrator_stats_nf = 1

    #=
    α1 = zero(uEltypeNoUnits)
    α2 = zero(uEltypeNoUnits)
    α3 = one(uEltypeNoUnits)
    β1 = 9 / 4 * one(uEltypeNoUnits)
    β2 = zero(uEltypeNoUnits)
    β3 = 3 / 4 * one(uEltypeNoUnits)
    cache = MPLM33oopCache(u, u, P, P, d, d, α1, α2, α3, β1, β2, β3, 1,
                           alg.small_constant_function(uEltypeNoUnits))
    =#

    α1 = 1 / 4 * one(uEltypeNoUnits)
    α2 = zero(uEltypeNoUnits)
    α3 = 3 / 4 * one(uEltypeNoUnits)
    α4 = zero(uEltypeNoUnits)
    β1 = 35 / 18 * one(uEltypeNoUnits)
    β2 = 1 / 3 * one(uEltypeNoUnits)
    β3 = zero(uEltypeNoUnits)
    β4 = 2 / 9 * one(uEltypeNoUnits)
    αβ = (α1, α2, α3, α4, β1, β2, β3, β4)
    MPLM43oopCache(u, u, u, P, P, P, d, d, d, αβ, 1,
                   alg.small_constant_function(uEltypeNoUnits))
end

struct MPLM54{F, T} <: OrdinaryDiffEqAlgorithm
    linsolve::F
    small_constant_function::T
end

alg_order(alg::MPLM54) = 4
alg_extrapolates(alg::MPLM54) = true # TODO: Should probably be false

@cache mutable struct MPLM54oopCache{uType, PType, dType, T, T2} <:
                      OrdinaryDiffEqConstantCache
    uprevprev::uType
    uprev3::uType
    uprev4::uType
    P2::PType
    P3::PType
    P4::PType
    P5::PType
    d2::dType
    d3::dType
    d4::dType
    d5::dType
    αβ::NTuple{10, T}
    step::Int
    small_constant::T2
end

function MPLM54(; linsolve = LUFactorization(), small_constant = nothing)
    if isnothing(small_constant)
        small_constant_function = floatmin
    elseif small_constant isa Number
        small_constant_function = Returns(small_constant)
    else # assume small_constant isa Function
        small_constant_function = small_constant
    end
    MPLM54(linsolve, small_constant_function)
end

function alg_cache(alg::MPLM54, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    # TODO: This is currently necessary to get the correct type of P (d is of type rateType)
    P, d = evaluate_pds(f, u, p, t)
    # TODO: integrator_stats_nf = 1

    α1 = zero(uEltypeNoUnits)
    α2 = zero(uEltypeNoUnits)
    α3 = zero(uEltypeNoUnits)
    α4 = zero(uEltypeNoUnits)
    α5 = one(uEltypeNoUnits)
    β1 = 225 / 96 * one(uEltypeNoUnits)
    β2 = zero(uEltypeNoUnits)
    β3 = 50 / 96 * one(uEltypeNoUnits)
    β4 = 200 / 96 * one(uEltypeNoUnits)
    β5 = 5 / 96 * one(uEltypeNoUnits)
    αβ = (α1, α2, α3, α4, α5, β1, β2, β3, β4, β5)
    MPLM54oopCache(u, u, u, P, P, P, P, d, d, d, d, αβ, 1,
                   alg.small_constant_function(uEltypeNoUnits))
end

function initialize!(integrator,
                     cache::Union{MPLM22oopCache,
                                  MPLM33oopCache, MPLM43oopCache})
end

### perform_step! ##########################################################################
@muladd function perform_step_MPLM22_oop(t, dt, uprev, uprev2, f, p, uprevprev,
                                         small_constant, linsolve)

    # evaluate production matrix
    P, d = evaluate_pds(f, uprev, p, t)

    # avoid division by zero due to zero Patankar weights
    σ = add_small_constant(uprev, small_constant)

    σ = basic_patankar_step(uprev, P, σ, dt, linsolve, d)

    # avoid division by zero due to zero Patankar weights
    σ = add_small_constant(σ, small_constant)

    u = basic_patankar_step(uprevprev, P, σ, 2 * dt, linsolve, d)

    # statistics: 1 nf, 2 nsolve

    return u
end

@muladd function perform_step!(integrator, cache::MPLM22oopCache, repeat_step = false)
    (; alg, t, dt, uprev, uprev2, f, p) = integrator
    (; uprevprev, small_constant) = cache

    if integrator.u_modified
        cache.step = 1
    end

    if cache.step <= 1

        # increase step counter
        cache.step += 1

        # One macro step of MPE with reduced time step size
        v, nf, ns = perform_substeps_MPE_oop(t, dt, 4, 1, uprev, f, p, small_constant,
                                             alg.linsolve)
        u = v[4]

        integrator.stats.nf += nf
        integrator.stats.nsolve += ns
    else
        u = perform_step_MPLM22_oop(t, dt, uprev, uprev2, f, p, uprevprev, small_constant,
                                    alg.linsolve)

        integrator.stats.nf += 1
        integrator.stats.nsolve += 2
    end

    #TODO: Should be possible to use uprev2. But uprev2 is currently not updated.
    cache.uprevprev = uprev

    integrator.u = u
end

@muladd function perform_substeps_MPLM22_oop(t, dt, num_sub_steps, num_macro_steps, uprev,
                                             uprev2, f, p,
                                             uprevprev, small_constant, linsolve)
    nfunc = 0
    nsolve = 0

    # we use v to store u at the end of each macro step
    v = Vector{typeof(uprev)}()

    dt_sub = dt / num_sub_steps

    # one step of size dt_sub
    v_mpe, nf, ns = perform_substeps_MPE_oop(t, dt_sub, 4, 1, uprev, f, p, small_constant,
                                             linsolve)

    u = v_mpe[4]
    push!(v, u)

    t += dt_sub

    nfunc += nf
    nsolve += ns

    # substeps 2-4 / macro step 1
    for _ in 1:(num_sub_steps - 1)
        uprevprev = uprev
        uprev = u

        u = perform_step_MPLM22_oop(t, dt_sub, uprev, uprev2, f, p, uprevprev,
                                    small_constant,
                                    linsolve)
        push!(v, u)

        t += dt_sub

        nfunc += 1
        nsolve += 2
    end

    # remaining macro steps 
    if num_macro_steps > 1
        for _ in 1:(num_macro_steps - 1)
            for _ in 1:num_sub_steps
                uprevprev = uprev
                uprev = u

                u = perform_step_MPLM22_oop(t, dt_sub, uprev, uprev2, f, p, uprevprev,
                                            small_constant, linsolve)
                push!(v, u)

                t += dt_sub
                nfunc += 1
                nsolve += 2
            end
        end
    end

    return v, nfunc, nsolve
end

@muladd function perform_step_MPLM33_oop(P_tup, d_tup, dt, u_tup, linsolve, αβ,
                                         small_constant)
    #   (; α1, α2, α3, β1, β2, β3, small_constant) = cache
    #(; P2, P3, d2, d3) = cache

    P, P2, P3 = P_tup
    d, d2, d3 = d_tup
    uprev, uprevprev, uprev3 = u_tup
    α1, α2, α3, β1, β2, β3 = αβ

    # avoid division by zero due to zero Patankar weights
    σ = add_small_constant(uprev, small_constant)

    σ = basic_patankar_step(uprev, P, σ, dt, linsolve, d)

    # avoid division by zero due to zero Patankar weights
    σ = add_small_constant(σ, small_constant)

    σ = basic_patankar_step(uprevprev, P, σ, 2 * dt, linsolve, d)

    # avoid division by zero due to zero Patankar weights
    σ = add_small_constant(σ, small_constant)

    Ptmp, dtmp = lincomb(β1, P, d, β2, P2, d2, β3, P3, d3)
    v = α1 * uprev + α2 * uprevprev + α3 * uprev3
    u = basic_patankar_step(v, Ptmp, σ, dt, linsolve, dtmp)

    # statistics: 3 nsolve

    return u
end

@muladd function perform_step!(integrator, cache::MPLM33oopCache, repeat_step = false)
    (; alg, t, dt, uprev, uprev2, f, p) = integrator
    (; uprevprev, uprev3, P2, P3, d2, d3, αβ, small_constant) = cache

    #TODO: is this necessary?
    if integrator.u_modified
        cache.step = 1
    end

    if cache.step == 1
        # increase step count
        cache.step += 1

        # evaluate production matrix at tspan[1]
        P, d = evaluate_pds(f, uprev, p, t)
        integrator.stats.nf += 1

        # compute starting values using MPLM22 with reduced time step size 
        v, nf, ns = perform_substeps_MPLM22_oop(t, dt, 4, 2, uprev, uprev2, f, p,
                                                uprevprev, small_constant, alg.linsolve)

        integrator.stats.nf += nf
        integrator.stats.nsolve += ns

        # u at time tspan[1] + dt
        u = v[4]

        cache.uprevprev = uprev

        # we use uprev3 as temporary storage for the value of u needed in step 2.
        cache.uprev3 = v[8]

    elseif cache.step == 2
        # increase step count
        cache.step += 1

        # evaluate production matrix at tspan[1] + dt
        P, d = evaluate_pds(f, uprev, p, t)
        integrator.stats.nf += 1

        # u at time tspan[1] + 2*dt (this was computed in step 1)
        u = cache.uprev3

        cache.uprev3 = uprevprev
        cache.uprevprev = uprev
    else

        # evaluate production matrix
        P, d = evaluate_pds(f, uprev, p, t)
        integrator.stats.nf += 1

        P_tup = (P, P2, P3)
        d_tup = (d, d2, d3)
        u_tup = (uprev, uprevprev, uprev3)

        u = perform_step_MPLM33_oop(P_tup, d_tup, dt, u_tup, alg.linsolve, αβ,
                                    small_constant)
        integrator.stats.nsolve += 3

        cache.uprev3 = uprevprev
        cache.uprevprev = uprev
    end

    integrator.u = u

    cache.P3 = cache.P2
    cache.P2 = P
    cache.d3 = cache.d2
    cache.d2 = d
end

@muladd function perform_substeps_MPLM33_oop(t, dt, num_sub_steps, num_macro_steps, uprev,
                                             uprev2, f, p,
                                             uprevprev, small_constant, linsolve, αβ)
    nfunc = 0
    nsolve = 0

    # we use v to store u at the end of each macro step
    v = Vector{typeof(uprev)}()

    dt_sub = dt / num_sub_steps

    # 1st and 2nd substep of macro step 1
    u22, nf, ns = perform_substeps_MPLM22_oop(t, dt_sub, 4, 2, uprev, uprev2, f, p,
                                              uprevprev, small_constant, linsolve)

    push!(v, u22[4])
    push!(v, u22[8])

    t += 2 * dt_sub

    nfunc += nf
    nsolve += ns

    uprevprev = uprev
    uprev = u22[4]
    u = u22[8]

    P2, d2 = evaluate_pds(f, uprevprev, p, t)
    nfunc += 1

    P, d = evaluate_pds(f, uprev, p, t + dt_sub)
    nfunc += 1

    # substeps 3 and 4 of macro step 1
    for _ in 3:4
        uprev3 = uprevprev
        uprevprev = uprev
        uprev = u

        P3 = P2
        P2 = P
        d3 = d2
        d2 = d

        P, d = evaluate_pds(f, uprev, p, t)
        nfunc += 1

        P_tup = (P, P2, P3)
        d_tup = (d, d2, d3)
        u_tup = (uprev, uprevprev, uprev3)

        u = perform_step_MPLM33_oop(P_tup, d_tup, dt_sub, u_tup, linsolve, αβ,
                                    small_constant)
        push!(v, u)

        nsolve += 3

        t += dt_sub
    end

    # remaining macro steps 
    if num_macro_steps > 1
        for _ in 1:(num_macro_steps - 1)
            for _ in 1:num_sub_steps
                uprev3 = uprevprev
                uprevprev = uprev
                uprev = u

                P3 = P2
                P2 = P
                d3 = d2
                d2 = d

                P, d = evaluate_pds(f, uprev, p, t)
                nfunc += 1

                P_tup = (P, P2, P3)
                d_tup = (d, d2, d3)
                u_tup = (uprev, uprevprev, uprev3)

                u = perform_step_MPLM33_oop(P_tup, d_tup, dt_sub, u_tup, linsolve, αβ,
                                            small_constant)
                push!(v, u)
                nsolve += 3

                t += dt_sub
            end
        end
    end

    return v, nfunc, nsolve
end

@muladd function perform_step!(integrator, cache::MPLM43oopCache, repeat_step = false)
    (; alg, t, dt, uprev, uprev2, f, p) = integrator
    (; uprevprev, uprev3, uprev4, P2, P3, P4, d2, d3, d4, αβ, small_constant) = cache
    α1, α2, α3, α4, β1, β2, β3, β4 = αβ

    #TODO: is this necessary?
    if integrator.u_modified
        cache.step = 1
    end

    if cache.step == 1
        # increase step count
        cache.step += 1

        # evaluate production matrix at tspan[1]
        P, d = evaluate_pds(f, uprev, p, t)
        integrator.stats.nf += 1

        # compute starting values using MPLM33 with reduced time step size
        αβ33 = (0, 0, 1, 9 / 4, 0, 3 / 4)
        v, nf, ns = perform_substeps_MPLM33_oop(t, dt, 4, 3, uprev, uprev2, f, p,
                                                uprevprev, small_constant, alg.linsolve,
                                                αβ33)
        integrator.stats.nf += nf
        integrator.stats.nsolve += ns

        # u at time tspan[1] + dt
        u = v[4]

        cache.uprevprev = uprev

        # we use uprev3 as temporary storage for the value of u needed in step 2.
        cache.uprev3 = v[8]
        # we use uprev4 as temporary storage for the value of u needed in step 3.
        cache.uprev4 = v[12]

    elseif cache.step == 2
        # increase step count
        cache.step += 1

        # evaluate production matrix at tspan[1] + dt
        P, d = evaluate_pds(f, uprev, p, t)
        integrator.stats.nf += 1

        # u at time tspan[1] + 2*dt (this was computed in step 1)
        u = cache.uprev3

        cache.uprev3 = uprevprev
        cache.uprevprev = uprev
    elseif cache.step == 3
        # increase step count
        cache.step += 1

        # evaluate production matrix at tspan[1] + 3*dt
        P, d = evaluate_pds(f, uprev, p, t)
        integrator.stats.nf += 1

        # u at time tspan[1] + 3*dt (this was computed in step 1)
        u = cache.uprev4

        cache.uprev4 = uprev3
        cache.uprev3 = uprevprev
        cache.uprevprev = uprev
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

        Ptmp, dtmp = lincomb(β1, P, d, β2, P2, d2, β3, P3, d3, β4, P4, d4)
        v = α1 * uprev + α2 * uprevprev + α3 * uprev3 + α4 * uprev4
        u = basic_patankar_step(v, Ptmp, σ, dt, alg.linsolve, dtmp)
        integrator.stats.nsolve += 1

        cache.uprev4 = uprev3
        cache.uprev3 = uprevprev
        cache.uprevprev = uprev
    end

    integrator.u = u

    cache.P4 = P3
    cache.P3 = P2
    cache.P2 = P
    cache.d4 = d3
    cache.d3 = d2
    cache.d2 = d
end
