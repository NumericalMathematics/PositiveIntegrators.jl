########################################################################################
### Structs and caches #################################################################
########################################################################################
#### MPLM22 ############################################################################
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

#### MPLM33 ############################################################################
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

#### MPLM43 ############################################################################
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

#### MPLM54 ############################################################################
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
    uprev5::uType
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
    MPLM54oopCache(u, u, u, u, P, P, P, P, d, d, d, d, αβ, 1,
                   alg.small_constant_function(uEltypeNoUnits))
end

#### MPLM75 ############################################################################
struct MPLM75{F, T} <: OrdinaryDiffEqAlgorithm
    linsolve::F
    small_constant_function::T
end

alg_order(alg::MPLM75) = 5
alg_extrapolates(alg::MPLM75) = true # TODO: Should probably be false

@cache mutable struct MPLM75oopCache{uType, PType, dType, T, T2} <:
                      OrdinaryDiffEqConstantCache
    uprevprev::uType
    uprev3::uType
    uprev4::uType
    uprev5::uType
    uprev6::uType
    uprev7::uType
    P2::PType
    P3::PType
    P4::PType
    P5::PType
    P6::PType
    P7::PType
    d2::dType
    d3::dType
    d4::dType
    d5::dType
    d6::dType
    d7::dType
    αβ::NTuple{14, T}
    step::Int
    small_constant::T2
end

function MPLM75(; linsolve = LUFactorization(), small_constant = nothing)
    if isnothing(small_constant)
        small_constant_function = floatmin
    elseif small_constant isa Number
        small_constant_function = Returns(small_constant)
    else # assume small_constant isa Function
        small_constant_function = small_constant
    end
    MPLM75(linsolve, small_constant_function)
end

function alg_cache(alg::MPLM75, u, rate_prototype, ::Type{uEltypeNoUnits},
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
    α5 = zero(uEltypeNoUnits)
    α6 = zero(uEltypeNoUnits)
    α7 = one(uEltypeNoUnits)
    β1 = 12 / 5 * one(uEltypeNoUnits)
    β2 = zero(uEltypeNoUnits)
    β3 = 197 / 720 * one(uEltypeNoUnits)
    β4 = 701 / 360 * one(uEltypeNoUnits)
    β5 = 43 / 30 * one(uEltypeNoUnits)
    β6 = 107 / 360 * one(uEltypeNoUnits)
    β7 = 467 / 720 * one(uEltypeNoUnits)
    αβ = (α1, α2, α3, α4, α5, α6, α7, β1, β2, β3, β4, β5, β6, β7)
    MPLM75oopCache(u, u, u, u, u, u, P, P, P, P, P, P,
                   d, d, d, d, d, d, αβ, 1, alg.small_constant_function(uEltypeNoUnits))
end
#### MPLM106 ###########################################################################
struct MPLM106{F, T} <: OrdinaryDiffEqAlgorithm
    linsolve::F
    small_constant_function::T
end

alg_order(alg::MPLM106) = 6
alg_extrapolates(alg::MPLM106) = true # TODO: Should probably be false

@cache mutable struct MPLM106oopCache{uType, PType, dType, T, T2} <:
                      OrdinaryDiffEqConstantCache
    uprevprev::uType
    uprev3::uType
    uprev4::uType
    uprev5::uType
    uprev6::uType
    uprev7::uType
    uprev8::uType
    uprev9::uType
    uprev10::uType
    P2::PType
    P3::PType
    P4::PType
    P5::PType
    P6::PType
    P7::PType
    P8::PType
    P9::PType
    P10::PType
    d2::dType
    d3::dType
    d4::dType
    d5::dType
    d6::dType
    d7::dType
    d8::dType
    d9::dType
    d10::dType
    αβ::NTuple{20, T}
    step::Int
    small_constant::T2
end

function MPLM106(; linsolve = LUFactorization(), small_constant = nothing)
    if isnothing(small_constant)
        small_constant_function = floatmin
    elseif small_constant isa Number
        small_constant_function = Returns(small_constant)
    else # assume small_constant isa Function
        small_constant_function = small_constant
    end
    MPLM106(linsolve, small_constant_function)
end

function alg_cache(alg::MPLM106, u, rate_prototype, ::Type{uEltypeNoUnits},
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
    α5 = zero(uEltypeNoUnits)
    α6 = zero(uEltypeNoUnits)
    α7 = zero(uEltypeNoUnits)
    α8 = zero(uEltypeNoUnits)
    α9 = zero(uEltypeNoUnits)
    α10 = one(uEltypeNoUnits)

    β1 = 11125 / 4536 * one(uEltypeNoUnits)
    β2 = zero(uEltypeNoUnits)
    β3 = zero(uEltypeNoUnits)
    β4 = 50 / 27 * one(uEltypeNoUnits)
    β5 = 85 / 36 * one(uEltypeNoUnits)
    β6 = zero(uEltypeNoUnits)
    β7 = zero(uEltypeNoUnits)
    β8 = 125 / 63 * one(uEltypeNoUnits)
    β9 = 25 / 24 * one(uEltypeNoUnits)
    β10 = 25 / 81 * one(uEltypeNoUnits)
    αβ = (α1, α2, α3, α4, α5, α6, α7, α8, α9, α10,
          β1, β2, β3, β4, β5, β6, β7, β8, β9, β10)
    MPLM106oopCache(u, u, u, u, u, u, u, u, u,
                    P, P, P, P, P, P, P, P, P,
                    d, d, d, d, d, d, d, d, d,
                    αβ, 1, alg.small_constant_function(uEltypeNoUnits))
end

function initialize!(integrator,
                     cache::Union{MPLM22oopCache, MPLM33oopCache, MPLM43oopCache,
                                  MPLM54oopCache, MPLM75oopCache, MPLM106oopCache})
end

########################################################################################
### perform_substeps ###################################################################
########################################################################################
# This function performs the computation of initial values for the m-step methods.
@muladd function perform_substeps_MPLM_oop(t, dt, uprev, f, p, small_constant, linsolve, αβ,
                                           num_macro_steps, num_sub_steps,
                                           startup_func, startup_steps, step_func,
                                           nsolve_step)
    nfunc = 0
    nsolve = 0

    total_steps = num_sub_steps * num_macro_steps

    # we use v to store u at the end of each macro step
    v = Vector{typeof(uprev)}()

    # we use history_u, history_P, history_d to store the values of u, P, d needed for the step function
    history_u = Vector{typeof(uprev)}(undef, startup_steps + 1)
    history_P = Vector{Any}(undef, startup_steps + 1)
    history_d = Vector{Any}(undef, startup_steps + 1)

    # set substep size
    dt_sub = dt / num_sub_steps

    # Startup phase using startup_func with reduced time step size
    v_start, nf, ns = startup_func(t, dt_sub, 4, startup_steps, uprev, f, p, small_constant,
                                   linsolve)
    for i in 1:startup_steps
        push!(v, v_start[num_sub_steps * i])
    end
    u = v[startup_steps] # == v_start[num_sub_steps * startup_steps]

    # Set t to time at end of startup phase
    t += startup_steps * dt_sub

    nfunc += nf
    nsolve += ns

    # Fill history of u 
    # v = ..., uprevprev, ..., uprev, ..., u  
    # history_u = uprev, uprevprev, uprev3, ...
    for i in 1:(startup_steps - 1)
        history_u[i] = v[startup_steps - i]
    end
    history_u[startup_steps] = uprev

    # Fill history of P and d
    for i in 1:startup_steps
        history_P[i], history_d[i] = evaluate_pds(f, history_u[i], p, t - i * dt_sub)
    end
    nfunc += startup_steps

    for _ in (startup_steps + 1):total_steps
        for i in (startup_steps + 1):-1:2
            history_u[i] = history_u[i - 1]
            history_P[i] = history_P[i - 1]
            history_d[i] = history_d[i - 1]
        end
        history_u[1] = u

        history_P[1], history_d[1] = evaluate_pds(f, history_u[1], p, t)
        nfunc += 1

        P_tup = tuple(history_P...)
        d_tup = tuple(history_d...)
        u_tup = tuple(history_u...)

        u = step_func(P_tup, d_tup, dt_sub, u_tup, linsolve, αβ,
                      small_constant)
        push!(v, u)
        nsolve += nsolve_step

        t += dt_sub
    end

    return v, nfunc, nsolve
end

@muladd function perform_substeps_MPLM22_oop(t, dt, num_sub_steps, num_macro_steps, uprev,
                                             f, p, small_constant, linsolve)
    αβ22 = (0, 1, 2, 0)

    return perform_substeps_MPLM_oop(t, dt, uprev, f, p, small_constant, linsolve, αβ22,
                                     num_macro_steps, # macro steps ofsize dt
                                     num_sub_steps, # substeps per macro step
                                     perform_substeps_MPE_oop, # substep function for startup phase
                                     1, # number of startup steps (number of substeps computed with startup function) 
                                     perform_step_MPLM22_oop, # step function for main phase
                                     2)
end

@muladd function perform_substeps_MPLM33_oop(t, dt, num_sub_steps, num_macro_steps, uprev,
                                             f, p, small_constant, linsolve)
    αβ33 = (0, 0, 1, 9 / 4, 0, 3 / 4)

    return perform_substeps_MPLM_oop(t, dt, uprev, f, p, small_constant, linsolve, αβ33,
                                     num_macro_steps, # macro steps ofsize dt
                                     num_sub_steps, # substeps per macro step
                                     perform_substeps_MPLM22_oop, # substep function for startup phase
                                     2, # number of startup steps (number of substeps computed with startup function) 
                                     perform_step_MPLM33_oop, # step function for main phase
                                     3)
end

@muladd function perform_substeps_MPLM43_oop(t, dt, num_sub_steps, num_macro_steps, uprev,
                                             f, p, small_constant, linsolve)
    αβ43 = (1 / 4, 0, 3 / 4, 0, 35 / 18, 1 / 3, 0, 2 / 9)

    return perform_substeps_MPLM_oop(t, dt, uprev, f, p, small_constant, linsolve, αβ43,
                                     num_macro_steps, # macro steps ofsize dt
                                     num_sub_steps, # substeps per macro step
                                     perform_substeps_MPLM33_oop, # substep function for startup phase
                                     3, # number of startup steps (number of substeps computed with startup function) 
                                     perform_step_MPLM43_oop, # step function for main phase
                                     3)
end

@muladd function perform_substeps_MPLM54_oop(t, dt, num_sub_steps, num_macro_steps, uprev,
                                             f, p, small_constant, linsolve)
    αβ54 = (0, 0, 0, 0, 1, 225 / 96, 0, 50 / 96, 200 / 96, 5 / 96)

    return perform_substeps_MPLM_oop(t, dt, uprev, f, p, small_constant, linsolve, αβ54,
                                     num_macro_steps, # macro steps ofsize dt
                                     num_sub_steps, # substeps per macro step
                                     perform_substeps_MPLM43_oop, # substep function for startup phase
                                     4, # number of startup steps (number of substeps computed with startup function) 
                                     perform_step_MPLM54_oop, # step function for main phase
                                     4)
end

@muladd function perform_substeps_MPLM75_oop(t, dt, num_sub_steps, num_macro_steps, uprev,
                                             f, p, small_constant, linsolve)
    αβ75 = (0, 0, 0, 0, 0, 0, 1, 12 / 5, 0, 197 / 720, 701 / 360, 43 / 30, 107 / 360,
            467 / 720)

    return perform_substeps_MPLM_oop(t, dt, uprev, f, p, small_constant, linsolve, αβ75,
                                     num_macro_steps, # macro steps ofsize dt
                                     num_sub_steps, # substeps per macro step
                                     perform_substeps_MPLM54_oop, # substep function for startup phase
                                     6, # number of startup steps (number of substeps computed with startup function) 
                                     perform_step_MPLM75_oop, # step function for main phase
                                     5)
end

########################################################################################
### perform_step! ######################################################################
########################################################################################
#### MPLM22 ############################################################################
#TODO Use αβ in MPLM22
@muladd function perform_step_MPLM22_oop(P_tup, d_tup, dt, u_tup, linsolve, αβ,
                                         small_constant)
    uprev, uprevprev = u_tup
    P = P_tup[1]
    d = d_tup[1]

    # First σ approximation
    σ = add_small_constant(uprev, small_constant)

    σ = basic_patankar_step(uprev, P, σ, dt, linsolve, d)

    # Main step 
    σ = add_small_constant(σ, small_constant)

    u = basic_patankar_step(uprevprev, P, σ, 2 * dt, linsolve, d)

    # statistics: 2 nsolve

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
        # evaluate production matrix
        P, d = evaluate_pds(f, uprev, p, t)
        integrator.stats.nf += 1

        u = perform_step_MPLM22_oop((P,), (d,), dt, (uprev, uprevprev), alg.linsolve,
                                    nothing, small_constant)
        integrator.stats.nsolve += 2
    end

    #TODO: Should be possible to use uprev2. But uprev2 is currently not updated.
    cache.uprevprev = uprev

    integrator.u = u
end

#### MPLM33 ############################################################################
@muladd function perform_step_MPLM33_oop(P_tup, d_tup, dt, u_tup, linsolve, αβ,
                                         small_constant)
    P, P2, P3 = P_tup
    d, d2, d3 = d_tup
    uprev, uprevprev, uprev3 = u_tup
    α1, α2, α3, β1, β2, β3 = αβ

    # First σ approximation
    σ = add_small_constant(uprev, small_constant)

    σ = basic_patankar_step(uprev, P, σ, dt, linsolve, d)

    # Second σ approximation
    σ = add_small_constant(σ, small_constant)

    σ = basic_patankar_step(uprevprev, P, σ, 2 * dt, linsolve, d)

    # Main step 
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
        v, nf, ns = perform_substeps_MPLM22_oop(t, dt, 4, 2, uprev, f, p,
                                                small_constant, alg.linsolve)

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

#### MPLM43 ############################################################################
@muladd function perform_step_MPLM43_oop(P_tup, d_tup, dt, u_tup, linsolve, αβ,
                                         small_constant)
    P, P2, P3, P4 = P_tup
    d, d2, d3, d4 = d_tup
    uprev, uprevprev, uprev3, uprev4 = u_tup
    α1, α2, α3, α4, β1, β2, β3, β4 = αβ

    # First σ approximation 
    σ = add_small_constant(uprev, small_constant)

    σ = basic_patankar_step(uprev, P, σ, dt, linsolve, d)

    # Second σ approximation
    σ = add_small_constant(σ, small_constant)

    σ = basic_patankar_step(uprevprev, P, σ, 2 * dt, linsolve, d)

    # Main step 
    σ = add_small_constant(σ, small_constant)

    Ptmp, dtmp = lincomb(β1, P, d, β2, P2, d2, β3, P3, d3, β4, P4, d4)
    v = α1 * uprev + α2 * uprevprev + α3 * uprev3 + α4 * uprev4
    u = basic_patankar_step(v, Ptmp, σ, dt, linsolve, dtmp)

    # statistics: 3 nsolve

    return u
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
        #v, nf, ns = perform_substeps_MPLM33_oop(t, dt, 4, 3, uprev, uprev2, f, p,
        #                                        uprevprev, small_constant, alg.linsolve,
        #                                        αβ33)
        #TODO Check that nsolve_step is correct! EVERYWHERE!
        v, nf, ns = perform_substeps_MPLM_oop(t, dt, uprev, f, p, small_constant,
                                              alg.linsolve, αβ33,
                                              3, # macro steps ofsize dt
                                              4, # substeps per macro step
                                              perform_substeps_MPLM22_oop, # substep function for startup phase
                                              2, # number of startup steps (number of substeps computed with startup function) 
                                              perform_step_MPLM33_oop, # step function for main phase
                                              3)
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

        P_tup = (P, P2, P3, P4)
        d_tup = (d, d2, d3, d4)
        u_tup = (uprev, uprevprev, uprev3, uprev4)

        u = perform_step_MPLM43_oop(P_tup, d_tup, dt, u_tup, alg.linsolve, αβ,
                                    small_constant)
        integrator.stats.nsolve += 3

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

#### MPLM54 ############################################################################
@muladd function perform_step_MPLM54_oop(P_tup, d_tup, dt, u_tup, linsolve, αβ,
                                         small_constant)
    P, P2, P3, P4, P5 = P_tup
    d, d2, d3, d4, d5 = d_tup
    uprev, uprevprev, uprev3, uprev4, uprev5 = u_tup
    α1, α2, α3, α4, α5, β1, β2, β3, β4, β5 = αβ

    # First σ approximation 
    σ = add_small_constant(uprev, small_constant)

    σ = basic_patankar_step(uprev, P, σ, dt, linsolve, d)

    # Second σ approximation
    σ = add_small_constant(σ, small_constant)

    σ = basic_patankar_step(uprevprev, P, σ, 2 * dt, linsolve, d)

    # Third σ approximation
    σ = add_small_constant(σ, small_constant)

    Ptmp, dtmp = lincomb(35 / 18, P, d, 1 / 3, P2, d2, 2 / 9, P4, d4)
    v = 1 / 4 * uprev + 3 / 4 * uprev3
    σ = basic_patankar_step(v, Ptmp, σ, dt, linsolve, dtmp)

    # Main step 
    σ = add_small_constant(σ, small_constant)

    Ptmp, dtmp = lincomb(β1, P, d, β2, P2, d2, β3, P3, d3, β4, P4, d4, β5, P5, d5)
    v = α1 * uprev + α2 * uprevprev + α3 * uprev3 + α4 * uprev4 + α5 * uprev5
    u = basic_patankar_step(v, Ptmp, σ, dt, linsolve, dtmp)

    # statistics: 4 nsolve

    return u
end

@muladd function perform_step!(integrator, cache::MPLM54oopCache, repeat_step = false)
    (; alg, t, dt, uprev, uprev2, f, p) = integrator
    (; uprevprev, uprev3, uprev4, uprev5, P2, P3, P4, P5, d2, d3, d4, d5, αβ, small_constant) = cache
    α1, α2, α3, α4, α5, β1, β2, β3, β4, β5 = αβ

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
        αβ43 = (1 / 4, 0, 3 / 4, 0, 35 / 18, 1 / 3, 0, 2 / 9)
        v, nf, ns = perform_substeps_MPLM_oop(t, dt, uprev, f, p, small_constant,
                                              alg.linsolve, αβ43,
                                              4, # macro steps ofsize dt
                                              4, # substeps per macro step
                                              perform_substeps_MPLM33_oop, # substep function for startup phase
                                              3, # number of startup steps (number of substeps computed with startup function) 
                                              perform_step_MPLM43_oop, # step function for main phase
                                              3)
        integrator.stats.nf += nf
        integrator.stats.nsolve += ns

        # u at time tspan[1] + dt
        u = v[4]

        cache.uprevprev = uprev

        # we use uprev3 as temporary storage for the value of u needed in step 2.
        cache.uprev3 = v[8]
        # we use uprev4 as temporary storage for the value of u needed in step 3.
        cache.uprev4 = v[12]
        # we use uprev5 as temporary storage for the value of u needed in step 4.
        cache.uprev5 = v[16]

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
    elseif cache.step == 4
        # increase step count
        cache.step += 1

        # evaluate production matrix at tspan[1] + 4*dt
        P, d = evaluate_pds(f, uprev, p, t)
        integrator.stats.nf += 1

        # u at time tspan[1] + 4*dt (this was computed in step 1)
        u = cache.uprev5

        cache.uprev5 = uprev4
        cache.uprev4 = uprev3
        cache.uprev3 = uprevprev
        cache.uprevprev = uprev
    else
        # increase step count
        cache.step += 1

        # evaluate production matrix
        P, d = evaluate_pds(f, uprev, p, t)
        integrator.stats.nf += 1

        P_tup = (P, P2, P3, P4, P5)
        d_tup = (d, d2, d3, d4, d5)
        u_tup = (uprev, uprevprev, uprev3, uprev4, uprev5)

        u = perform_step_MPLM54_oop(P_tup, d_tup, dt, u_tup, alg.linsolve, αβ,
                                    small_constant)
        integrator.stats.nsolve += 4

        cache.uprev5 = uprev4
        cache.uprev4 = uprev3
        cache.uprev3 = uprevprev
        cache.uprevprev = uprev
    end

    integrator.u = u

    cache.P5 = P4
    cache.P4 = P3
    cache.P3 = P2
    cache.P2 = P
    cache.d5 = d4
    cache.d4 = d3
    cache.d3 = d2
    cache.d2 = d
end
#### MPLM75 ############################################################################
@muladd function perform_step_MPLM75_oop(P_tup, d_tup, dt, u_tup, linsolve, αβ,
                                         small_constant)
    P, P2, P3, P4, P5, P6, P7 = P_tup
    d, d2, d3, d4, d5, d6, d7 = d_tup
    uprev, uprevprev, uprev3, uprev4, uprev5, uprev6, uprev7 = u_tup
    α1, α2, α3, α4, α5, α6, α7, β1, β2, β3, β4, β5, β6, β7 = αβ

    # First σ approximation 
    σ = add_small_constant(uprev, small_constant)

    σ = basic_patankar_step(uprev, P, σ, dt, linsolve, d)

    # Second σ approximation
    σ = add_small_constant(σ, small_constant)

    σ = basic_patankar_step(uprevprev, P, σ, 2 * dt, linsolve, d)

    # Third σ approximation
    σ = add_small_constant(σ, small_constant)

    Ptmp, dtmp = lincomb(35 / 18, P, d, 1 / 3, P2, d2, 2 / 9, P4, d4)
    v = 1 / 4 * uprev + 3 / 4 * uprev3
    σ = basic_patankar_step(v, Ptmp, σ, dt, linsolve, dtmp)

    # Fourth σ approximation
    σ = add_small_constant(σ, small_constant)

    Ptmp, dtmp = lincomb(225 / 96, P, d, 50 / 96, P3, d3, 200 / 96, P4, d4, 5 / 96, P5, d5)
    v = uprev5
    σ = basic_patankar_step(v, Ptmp, σ, dt, linsolve, dtmp)

    # Main step 
    σ = add_small_constant(σ, small_constant)

    Ptmp, dtmp = lincomb(β1, P, d, β2, P2, d2, β3, P3, d3, β4, P4, d4, β5, P5, d5, β6, P6,
                         d6, β7, P7, d7)
    v = α1 * uprev + α2 * uprevprev + α3 * uprev3 + α4 * uprev4 + α5 * uprev5 +
        α6 * uprev6 + α7 * uprev7
    u = basic_patankar_step(v, Ptmp, σ, dt, linsolve, dtmp)

    # statistics: 5 nsolve

    return u
end

@muladd function perform_step!(integrator, cache::MPLM75oopCache, repeat_step = false)
    (; alg, t, dt, uprev, uprev2, f, p) = integrator
    (; uprevprev, uprev3, uprev4, uprev5, uprev6, uprev7, P2, P3, P4, P5, P6, P7, d2, d3, d4, d5, d6, d7, αβ, small_constant) = cache
    α1, α2, α3, α4, α5, α6, α7, β1, β2, β3, β4, β5 = αβ

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
        αβ54 = (0, 0, 0, 0, 1, 225 / 96, 0, 50 / 96, 200 / 96, 5 / 96)
        v, nf, ns = perform_substeps_MPLM_oop(t, dt, uprev, f, p, small_constant,
                                              alg.linsolve, αβ54,
                                              6, # macro steps of size dt
                                              4, # substeps per macro step
                                              perform_substeps_MPLM43_oop, # substep function for startup phase
                                              4, # number of startup steps (number of substeps computed with startup function) 
                                              perform_step_MPLM54_oop, # step function for main phase
                                              4)
        integrator.stats.nf += nf
        integrator.stats.nsolve += ns

        # u at time tspan[1] + dt
        u = v[4]

        cache.uprevprev = uprev

        # we use uprev3 as temporary storage for the value of u needed in step 2.
        cache.uprev3 = v[8]
        # we use uprev4 as temporary storage for the value of u needed in step 3.
        cache.uprev4 = v[12]
        # we use uprev5 as temporary storage for the value of u needed in step 4.
        cache.uprev5 = v[16]
        # we use uprev6 as temporary storage for the value of u needed in step 5.
        cache.uprev6 = v[20]
        # we use uprev7 as temporary storage for the value of u needed in step 6.
        cache.uprev7 = v[24]

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
    elseif cache.step == 4
        # increase step count
        cache.step += 1

        # evaluate production matrix at tspan[1] + 4*dt
        P, d = evaluate_pds(f, uprev, p, t)
        integrator.stats.nf += 1

        # u at time tspan[1] + 4*dt (this was computed in step 1)
        u = cache.uprev5

        cache.uprev5 = uprev4
        cache.uprev4 = uprev3
        cache.uprev3 = uprevprev
        cache.uprevprev = uprev
    elseif cache.step == 5
        # increase step count
        cache.step += 1

        # evaluate production matrix at tspan[1] + 4*dt
        P, d = evaluate_pds(f, uprev, p, t)
        integrator.stats.nf += 1

        # u at time tspan[1] + 5*dt (this was computed in step 1)
        u = cache.uprev6

        cache.uprev6 = uprev5
        cache.uprev5 = uprev4
        cache.uprev4 = uprev3
        cache.uprev3 = uprevprev
        cache.uprevprev = uprev
    elseif cache.step == 6
        # increase step count
        cache.step += 1

        # evaluate production matrix at tspan[1] + 4*dt
        P, d = evaluate_pds(f, uprev, p, t)
        integrator.stats.nf += 1

        # u at time tspan[1] + 6*dt (this was computed in step 1)
        u = cache.uprev7

        cache.uprev7 = uprev6
        cache.uprev6 = uprev5
        cache.uprev5 = uprev4
        cache.uprev4 = uprev3
        cache.uprev3 = uprevprev
        cache.uprevprev = uprev
    else
        # increase step count
        cache.step += 1

        # evaluate production matrix
        P, d = evaluate_pds(f, uprev, p, t)
        integrator.stats.nf += 1

        P_tup = (P, P2, P3, P4, P5, P6, P7)
        d_tup = (d, d2, d3, d4, d5, d6, d7)
        u_tup = (uprev, uprevprev, uprev3, uprev4, uprev5, uprev6, uprev7)

        u = perform_step_MPLM75_oop(P_tup, d_tup, dt, u_tup, alg.linsolve, αβ,
                                    small_constant)
        integrator.stats.nsolve += 5

        cache.uprev7 = uprev6
        cache.uprev6 = uprev5
        cache.uprev5 = uprev4
        cache.uprev4 = uprev3
        cache.uprev3 = uprevprev
        cache.uprevprev = uprev
    end

    integrator.u = u

    cache.P7 = P6
    cache.P6 = P5
    cache.P5 = P4
    cache.P4 = P3
    cache.P3 = P2
    cache.P2 = P
    cache.d7 = d6
    cache.d6 = d5
    cache.d5 = d4
    cache.d4 = d3
    cache.d3 = d2
    cache.d2 = d
end

#### MPLM106 ############################################################################
@muladd function perform_step_MPLM106_oop(P_tup, d_tup, dt, u_tup, linsolve, αβ,
                                          small_constant)
    P, P2, P3, P4, P5, P6, P7, P8, P9, P10 = P_tup
    d, d2, d3, d4, d5, d6, d7, d8, d9, d10 = d_tup
    uprev, uprevprev, uprev3, uprev4, uprev5, uprev6, uprev7, uprev8, uprev9, uprev10 = u_tup
    α1, α2, α3, α4, α5, α6, α7, α8, α9, α10, β1, β2, β3, β4, β5, β6, β7, β8, β9, β10 = αβ

    # First σ approximation 
    σ = add_small_constant(uprev, small_constant)

    σ = basic_patankar_step(uprev, P, σ, dt, linsolve, d)

    # Second σ approximation
    σ = add_small_constant(σ, small_constant)

    σ = basic_patankar_step(uprevprev, P, σ, 2 * dt, linsolve, d)

    # Third σ approximation
    σ = add_small_constant(σ, small_constant)

    Ptmp, dtmp = lincomb(35 / 18, P, d, 1 / 3, P2, d2, 2 / 9, P4, d4)
    v = 1 / 4 * uprev + 3 / 4 * uprev3
    σ = basic_patankar_step(v, Ptmp, σ, dt, linsolve, dtmp)

    # Fourth σ approximation
    σ = add_small_constant(σ, small_constant)

    Ptmp, dtmp = lincomb(225 / 96, P, d, 50 / 96, P3, d3, 200 / 96, P4, d4, 5 / 96, P5, d5)
    v = uprev5
    σ = basic_patankar_step(v, Ptmp, σ, dt, linsolve, dtmp)

    # Fifth σ approximation
    σ = add_small_constant(σ, small_constant)

    Ptmp, dtmp = lincomb(12 / 5, P, d, 197 / 720, P3, d3, 701 / 360, P4, d4, 43 / 30, P5,
                         d5, 107 / 360, P6, d6, 467 / 720, P7, d7)
    v = uprev7
    σ = basic_patankar_step(v, Ptmp, σ, dt, linsolve, dtmp)

    # Main step 
    σ = add_small_constant(σ, small_constant)

    Ptmp, dtmp = lincomb(β1, P, d, β2, P2, d2, β3, P3, d3, β4, P4, d4, β5, P5, d5, β6, P6,
                         d6, β7, P7, d7, β8, P8, d8, β9, P9, d9, β10, P10, d10)
    v = α1 * uprev + α2 * uprevprev + α3 * uprev3 + α4 * uprev4 + α5 * uprev5 +
        α6 * uprev6 + α7 * uprev7 + α8 * uprev8 + α9 * uprev9 + α10 * uprev10
    u = basic_patankar_step(v, Ptmp, σ, dt, linsolve, dtmp)

    # statistics: 6 nsolve

    return u
end

@muladd function perform_step!(integrator, cache::MPLM106oopCache, repeat_step = false)
    (; alg, t, dt, uprev, uprev2, f, p) = integrator
    (; uprevprev, uprev3, uprev4, uprev5, uprev6, uprev7, uprev8, uprev9, uprev10,
    P2, P3, P4, P5, P6, P7, P8, P9, P10, d2, d3, d4, d5, d6, d7, d8, d9, d10,
    αβ, small_constant) = cache
    α1, α2, α3, α4, α5, α6, α7, α8, α9, α10, β1, β2, β3, β4, β5, β6, β7, β8, β9, β10 = αβ

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
        αβ75 = (0, 0, 0, 0, 0, 0, 1, 12 / 5, 0, 197 / 720, 701 / 360, 43 / 30, 107 / 360,
                467 / 720)
        v, nf, ns = perform_substeps_MPLM_oop(t, dt, uprev, f, p, small_constant,
                                              alg.linsolve, αβ75,
                                              9, # macro steps of size dt
                                              4, # substeps per macro step
                                              perform_substeps_MPLM54_oop, # substep function for startup phase
                                              6, # number of startup steps (number of substeps computed with startup function) 
                                              perform_step_MPLM75_oop, # step function for main phase
                                              5)
        integrator.stats.nf += nf
        integrator.stats.nsolve += ns

        # u at time tspan[1] + dt
        u = v[4]

        cache.uprevprev = uprev

        # we use uprev3 as temporary storage for the value of u needed in step 2.
        cache.uprev3 = v[8]
        # we use uprev4 as temporary storage for the value of u needed in step 3.
        cache.uprev4 = v[12]
        # we use uprev5 as temporary storage for the value of u needed in step 4.
        cache.uprev5 = v[16]
        # we use uprev6 as temporary storage for the value of u needed in step 5.
        cache.uprev6 = v[20]
        # we use uprev7 as temporary storage for the value of u needed in step 6.
        cache.uprev7 = v[24]
        # we use uprev8 as temporary storage for the value of u needed in step 7.
        cache.uprev8 = v[28]
        # we use uprev9 as temporary storage for the value of u needed in step 8.
        cache.uprev9 = v[32]
        # we use uprev10 as temporary storage for the value of u needed in step 9.
        cache.uprev10 = v[36]

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
    elseif cache.step == 4
        # increase step count
        cache.step += 1

        # evaluate production matrix at tspan[1] + 4*dt
        P, d = evaluate_pds(f, uprev, p, t)
        integrator.stats.nf += 1

        # u at time tspan[1] + 4*dt (this was computed in step 1)
        u = cache.uprev5

        cache.uprev5 = uprev4
        cache.uprev4 = uprev3
        cache.uprev3 = uprevprev
        cache.uprevprev = uprev
    elseif cache.step == 5
        # increase step count
        cache.step += 1

        # evaluate production matrix at tspan[1] + 4*dt
        P, d = evaluate_pds(f, uprev, p, t)
        integrator.stats.nf += 1

        # u at time tspan[1] + 5*dt (this was computed in step 1)
        u = cache.uprev6

        cache.uprev6 = uprev5
        cache.uprev5 = uprev4
        cache.uprev4 = uprev3
        cache.uprev3 = uprevprev
        cache.uprevprev = uprev
    elseif cache.step == 6
        # increase step count
        cache.step += 1

        # evaluate production matrix at tspan[1] + 4*dt
        P, d = evaluate_pds(f, uprev, p, t)
        integrator.stats.nf += 1

        # u at time tspan[1] + 6*dt (this was computed in step 1)
        u = cache.uprev7

        cache.uprev7 = uprev6
        cache.uprev6 = uprev5
        cache.uprev5 = uprev4
        cache.uprev4 = uprev3
        cache.uprev3 = uprevprev
        cache.uprevprev = uprev
    elseif cache.step == 7
        # increase step count
        cache.step += 1

        # evaluate production matrix at tspan[1] + 4*dt
        P, d = evaluate_pds(f, uprev, p, t)
        integrator.stats.nf += 1

        # u at time tspan[1] + 7*dt (this was computed in step 1)
        u = cache.uprev8

        cache.uprev8 = uprev7
        cache.uprev7 = uprev6
        cache.uprev6 = uprev5
        cache.uprev5 = uprev4
        cache.uprev4 = uprev3
        cache.uprev3 = uprevprev
        cache.uprevprev = uprev
    elseif cache.step == 8
        # increase step count
        cache.step += 1

        # evaluate production matrix at tspan[1] + 4*dt
        P, d = evaluate_pds(f, uprev, p, t)
        integrator.stats.nf += 1

        # u at time tspan[1] + 8*dt (this was computed in step 1)
        u = cache.uprev9

        cache.uprev9 = uprev8
        cache.uprev8 = uprev7
        cache.uprev7 = uprev6
        cache.uprev6 = uprev5
        cache.uprev5 = uprev4
        cache.uprev4 = uprev3
        cache.uprev3 = uprevprev
        cache.uprevprev = uprev
    elseif cache.step == 9
        # increase step count
        cache.step += 1

        # evaluate production matrix at tspan[1] + 4*dt
        P, d = evaluate_pds(f, uprev, p, t)
        integrator.stats.nf += 1

        # u at time tspan[1] + 9*dt (this was computed in step 1)
        u = cache.uprev10

        cache.uprev10 = uprev9
        cache.uprev9 = uprev8
        cache.uprev8 = uprev7
        cache.uprev7 = uprev6
        cache.uprev6 = uprev5
        cache.uprev5 = uprev4
        cache.uprev4 = uprev3
        cache.uprev3 = uprevprev
        cache.uprevprev = uprev
    else

        # evaluate production matrix
        P, d = evaluate_pds(f, uprev, p, t)
        integrator.stats.nf += 1

        P_tup = (P, P2, P3, P4, P5, P6, P7, P8, P9, P10)
        d_tup = (d, d2, d3, d4, d5, d6, d7, d8, d9, d10)
        u_tup = (uprev, uprevprev, uprev3, uprev4, uprev5, uprev6, uprev7, uprev8, uprev9,
                 uprev10)

        u = perform_step_MPLM106_oop(P_tup, d_tup, dt, u_tup, alg.linsolve, αβ,
                                     small_constant)
        integrator.stats.nsolve += 6

        cache.uprev10 = uprev9
        cache.uprev9 = uprev8
        cache.uprev8 = uprev7
        cache.uprev7 = uprev6
        cache.uprev6 = uprev5
        cache.uprev5 = uprev4
        cache.uprev4 = uprev3
        cache.uprev3 = uprevprev
        cache.uprevprev = uprev
    end

    integrator.u = u

    cache.P10 = P9
    cache.P9 = P8
    cache.P8 = P7
    cache.P7 = P6
    cache.P6 = P5
    cache.P5 = P4
    cache.P4 = P3
    cache.P3 = P2
    cache.P2 = P
    cache.d10 = d9
    cache.d9 = d8
    cache.d8 = d7
    cache.d7 = d6
    cache.d6 = d5
    cache.d5 = d4
    cache.d4 = d3
    cache.d3 = d2
    cache.d2 = d
end
