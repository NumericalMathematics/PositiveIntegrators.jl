########################################################################################
### Structs and caches #################################################################
########################################################################################
abstract type MPLMMutableCache <: OrdinaryDiffEqMutableCache end
get_fsalfirstlast(cache::MPLMMutableCache, rate_prototype) = (nothing, nothing)

#### MPLM22 ############################################################################
struct MPLM22{F, T} <: OrdinaryDiffEqAlgorithm
    linsolve::F
    small_constant_function::T
end

alg_order(alg::MPLM22) = 2
isfsal(::MPLM22) = false
alg_extrapolates(alg::MPLM22) = false # TODO: Should probably be false

#TODO: Check if OrdinaryDiffEqConstantCache is correct supertype
@cache mutable struct MPLM22oopCache{uType, T} <: OrdinaryDiffEqConstantCache
    uprevprev::uType
    step::Int
    small_constant::T
end

@cache mutable struct MPLM22Cache{uType, dType, T, PType, F} <: MPLMMutableCache
    uprevprev::uType
    step::Int
    small_constant::T
    tmp::uType
    P::PType
    P2::PType
    D::dType
    D2::dType
    σ::uType
    linsolve::F
end

get_tmp_cache(integrator, ::MPLM22, cache::MPLMMutableCache) = (cache.σ,)

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
                   ::Val{false},
                   verbose) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    MPLM22oopCache(u, 1, alg.small_constant_function(uEltypeNoUnits))
end

# In-place
function alg_cache(alg::MPLM22, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
                   uprev, uprev2, f, t, dt, reltol, p, calck,
                   ::Val{true},
                   verbose) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    uprevprev = zero(u)

    step = 1
    small_constant = alg.small_constant_function(uEltypeNoUnits)
    tmp = zero(u)
    P = p_prototype(u, f)
    # We use P2 to store the last evaluation of the PDS
    # as well as to store the system matrix of the linear system
    P2 = p_prototype(u, f)
    σ = zero(u)

    if f isa ConservativePDSFunction
        # The right hand side of the linear system is always uprev. But using
        # tmp instead of uprev for the rhs we allow `alias_b=true`. uprev must
        # not be altered, since it is needed to compute the adaptive time step
        # size.
        linprob = LinearProblem(P2, _vec(tmp))
        linsolve = init(linprob, alg.linsolve,
                        alias = LinearSolve.LinearAliasSpecifier(; alias_A = true,
                                                                 alias_b = true),
                        assumptions = LinearSolve.OperatorAssumptions(true))

        MPLM22Cache(uprevprev, step, small_constant, tmp, P, P2, nothing, nothing, σ,
                    linsolve)
    elseif f isa PDSFunction
        linprob = LinearProblem(P2, _vec(tmp))
        linsolve = init(linprob, alg.linsolve,
                        alias = LinearSolve.LinearAliasSpecifier(; alias_A = true,
                                                                 alias_b = true),
                        assumptions = LinearSolve.OperatorAssumptions(true))

        MPLM22Cache(uprevprev, step, small_constant, tmp, P, P2,
                    similar(u), # D
                    similar(u), # D2
                    σ,
                    linsolve)
    else
        throw(ArgumentError("MPLM22 can only be applied to production-destruction systems"))
    end
end

#### MPLM33 ############################################################################
struct MPLM33{F, T} <: OrdinaryDiffEqAlgorithm
    linsolve::F
    small_constant_function::T
end

alg_order(alg::MPLM33) = 3
isfsal(::MPLM33) = false
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
                   ::Val{false},
                   verbose) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
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

@cache mutable struct MPLM33Cache{uType, dType, T, PType, F, TabType} <: MPLMMutableCache
    uprevprev::uType
    uprev3::uType
    v::uType
    vprev::uType
    vprev2::uType
    step::Int
    small_constant::T
    b::uType # rhs of the linear system
    P::PType
    P2::PType
    P3::PType
    A::PType # system matrix of the linear system
    d::dType
    d2::dType
    d3::dType
    d_sys::dType
    σ::uType
    linsolve::F
    αβ::TabType
end

function alg_cache(alg::MPLM33, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
                   uprev, uprev2, f, t, dt, reltol, p, calck,
                   ::Val{true},
                   verbose) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    uprevprev = zero(u)
    uprev3 = zero(u)
    v = zero(u)
    vprev = zero(u)
    vprev2 = zero(u)
    step = 1
    small_constant = alg.small_constant_function(uEltypeNoUnits)
    b = zero(u)
    P = p_prototype(u, f)
    P2 = p_prototype(u, f)
    P3 = p_prototype(u, f)
    A = p_prototype(u, f)
    σ = zero(u)

    # MPLM33 coefficients 
    α1 = zero(uEltypeNoUnits)
    α2 = zero(uEltypeNoUnits)
    α3 = one(uEltypeNoUnits)
    β1 = 9 / 4 * one(uEltypeNoUnits)
    β2 = zero(uEltypeNoUnits)
    β3 = 3 / 4 * one(uEltypeNoUnits)
    αβ = (α1, α2, α3, β1, β2, β3)

    if f isa ConservativePDSFunction
        linprob = LinearProblem(A, _vec(b))
        linsolve = init(linprob, alg.linsolve,
                        alias = LinearSolve.LinearAliasSpecifier(; alias_A = true,
                                                                 alias_b = true),
                        assumptions = LinearSolve.OperatorAssumptions(true))

        MPLM33Cache(uprevprev, uprev3, v, vprev, vprev2, step, small_constant, b, P, P2, P3,
                    A, nothing, nothing, nothing, nothing, σ,
                    linsolve, αβ)
    elseif f isa PDSFunction
        linprob = LinearProblem(A, _vec(b))
        linsolve = init(linprob, alg.linsolve,
                        alias = LinearSolve.LinearAliasSpecifier(; alias_A = true,
                                                                 alias_b = true),
                        assumptions = LinearSolve.OperatorAssumptions(true))

        MPLM33Cache(uprevprev, uprev3, v, vprev, vprev2, step, small_constant, b, P, P2, P3,
                    A,
                    similar(u), similar(u), similar(u), similar(u),
                    σ, linsolve, αβ)
    else
        throw(ArgumentError("MPLM33 can only be applied to production-destruction systems"))
    end
end

#### MPLM43 ############################################################################
struct MPLM43{F, T} <: OrdinaryDiffEqAlgorithm
    linsolve::F
    small_constant_function::T
end

alg_order(alg::MPLM43) = 3
isfsal(::MPLM43) = false
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
                   ::Val{false},
                   verbose) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
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
isfsal(::MPLM54) = false
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
                   ::Val{false},
                   verbose) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
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
isfsal(::MPLM75) = false
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
                   ::Val{false},
                   verbose) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
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
isfsal(::MPLM106) = false
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
                   ::Val{false},
                   verbose) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
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
                                  MPLM54oopCache, MPLM75oopCache, MPLM106oopCache,
                                  MPLMMutableCache})
end

########################################################################################
### perform_step! ######################################################################
########################################################################################
#### MPLM22 ############################################################################
@muladd function start_MPLM22_oop(P, d, t, dt, uprev, f, p, small_constant, linsolve)
    dt = dt / 4

    u = uprev
    @inbounds for _ in 1:4
        u = perform_step_MPE_oop(P, d, dt, u, small_constant, linsolve)

        P, d = evaluate_pds(f, u, p, t)

        t += dt
    end

    # 4 function evals and 4 solves
    nf = 4
    ns = 4

    return u, t, nf, ns
end

@muladd function start_MPLM22!(u, P, d, t, dt, uprev, σ, f, p, small_constant, linsolve)
    dt = dt / 4

    u .= uprev
    @inbounds for _ in 1:4
        perform_step_MPE!(u, P, d, dt, u, σ, small_constant, linsolve)

        evaluate_pds!(P, d, f, u, p, t)

        t += dt
    end

    # 4 function evals and 4 solves
    nf = 4
    ns = 4

    return nf, ns
end

#TODO Use αβ in MPLM22
@muladd function perform_step_MPLM22_oop(P, d, dt, uprev, uprevprev, linsolve,
                                         small_constant)

    # First σ approximation
    σ = add_small_constant(uprev, small_constant)

    σ = basic_patankar_step(uprev, P, σ, dt, linsolve, d)

    # Main step 
    σ = add_small_constant(σ, small_constant)

    u = basic_patankar_step(uprevprev, P, σ, 2 * dt, linsolve, d)
    # statistics: 2 nsolve

    return u
end

@muladd function perform_step_MPLM22!(u, P, P2, d, d2, dt, uprev, uprevprev, σ, linsolve,
                                      small_constant)
    # First σ approximation
    @.. broadcast=false σ=uprev + small_constant

    # use lincomb! to handle cases in which d2 is nothing
    #P2 .= P
    #lincomb!(d2, 1, d)

    #basic_patankar_step!(σ, uprev, P2, d2, linsolve.A, σ, dt, linsolve)
    basic_patankar_step!(σ, uprev, P, d, linsolve.A, σ, dt, linsolve)

    # Main step 
    @.. broadcast=false σ=σ + small_constant

    basic_patankar_step!(u, uprevprev, P, d, linsolve.A, σ, 2 * dt, linsolve)

    # statistics: 2 nsolve

    return nothing
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

        P, d = evaluate_pds(f, uprev, p, t)
        integrator.stats.nf += 1

        # compute initial values for MPLM22
        u, _, nf, ns = start_MPLM22_oop(P, d, t, dt, uprev, f, p, small_constant,
                                        alg.linsolve)

        integrator.stats.nf += nf
        integrator.stats.nsolve += ns
    else
        # increase step counter
        cache.step += 1

        # evaluate production matrix
        P, d = evaluate_pds(f, uprev, p, t)
        integrator.stats.nf += 1

        u = perform_step_MPLM22_oop(P, d, dt, uprev, uprevprev, alg.linsolve,
                                    small_constant)
        integrator.stats.nsolve += 2
    end

    #TODO: Should be possible to use uprev2. But uprev2 is currently not updated.
    cache.uprevprev = uprev

    integrator.u = u
end

@muladd function perform_step!(integrator, cache::MPLM22Cache, repeat_step = false)
    (; t, dt, u, uprev, uprev2, f, p) = integrator
    (; uprevprev, small_constant, P, P2, D, D2, σ, linsolve) = cache

    if integrator.u_modified
        cache.step = 1
    end

    if cache.step <= 1

        # increase step counter
        cache.step += 1

        evaluate_pds!(P, D, f, uprev, p, t)
        integrator.stats.nf += 1

        # compute initial values for MPLM22
        nf, ns = start_MPLM22!(u, P, D, t, dt, uprev, σ, f, p, small_constant, linsolve)

        integrator.stats.nf += nf
        integrator.stats.nsolve += ns
    else

        # increase step counter
        cache.step += 1

        # evaluate production matrix    
        evaluate_pds!(P, D, f, uprev, p, t)
        integrator.stats.nf += 1

        perform_step_MPLM22!(u, P, P2, D, D2, dt, uprev, uprevprev, σ, linsolve,
                             small_constant)

        integrator.stats.nsolve += 2
    end

    #TODO: Should be possible to use uprev2. But uprev2 is currently not updated.
    uprevprev .= uprev
end

#### MPLM33 ############################################################################
@muladd function start_MPLM33_oop(P, d, t, dt, uprev, f, p, small_constant, linsolve)

    # 1 macro step consists of 4 substeps                                  
    dts = dt / 4

    ### first macro step ###############################################################
    # substep 1
    u, t, nf, ns = start_MPLM22_oop(P, d, t, dts, uprev, f, p, small_constant, linsolve)

    # substeps 2 - 4                                    
    for _ in 1:3
        uprevprev = uprev
        uprev = u

        P, d = evaluate_pds(f, uprev, p, t)
        nf += 1

        u = perform_step_MPLM22_oop(P, d, dts, uprev, uprevprev, linsolve, small_constant)
        t += dts
        ns += 2
    end

    v = u

    ### second macro step ############################################################
    for _ in 1:4
        uprevprev = uprev
        uprev = u

        P, d = evaluate_pds(f, uprev, p, t)
        nf += 1

        u = perform_step_MPLM22_oop(P, d, dts, uprev, uprevprev, linsolve, small_constant)
        t += dts
        ns += 2
    end

    return (v, u), t, nf, ns
end

@muladd function start_MPLM33!(v, tmp, P, P2, d, d2, t, dt, vprev, vprev2, σ, f, p,
                               small_constant, linsolve)

    # 1 macro step consists of 4 substeps                                  
    dts = dt / 4

    ### first macro step ###############################################################
    # substep 1
    nf, ns = start_MPLM22!(v, P, d, t, dts, vprev, σ, f, p, small_constant, linsolve)

    # substeps 2 - 4                                    
    for _ in 1:3
        vprev2 .= vprev
        vprev .= v

        evaluate_pds!(P, d, f, vprev, p, t)
        nf += 1

        perform_step_MPLM22!(v, P, P2, d, d2, dts, vprev, vprev2, σ, linsolve,
                             small_constant)
        t += dts
        ns += 2
    end

    tmp .= v

    ### second macro step ############################################################
    for _ in 1:4
        vprev2 .= vprev
        vprev .= v

        evaluate_pds!(P, d, f, vprev, p, t)
        nf += 1

        perform_step_MPLM22!(v, P, P2, d, d2, dts, vprev, vprev2, σ, linsolve,
                             small_constant)
        t += dts
        ns += 2
    end

    return nf, ns
end

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

@muladd function perform_step_MPLM33!(u, P_tup, d_tup, dt, u_tup, σ, linsolve, αβ,
                                      small_constant)
    P, P2, P3 = P_tup
    d, d2, d3 = d_tup
    uprev, uprevprev, uprev3 = u_tup
    α1, α2, α3, β1, β2, β3 = αβ

    # First σ approximation
    @.. broadcast=false σ=uprev + small_constant

    basic_patankar_step!(σ, uprev, P, d, linsolve.A, σ, dt, linsolve)

    # Second σ approximation
    @.. broadcast=false σ=σ + small_constant

    basic_patankar_step!(σ, uprevprev, P, d, linsolve.A, σ, 2 * dt, linsolve)

    # Main step 
    @.. broadcast=false σ=σ + small_constant

    lincomb!(P3, β1, P, β2, P2, β3, P3)
    lincomb!(d3, β1, d, β2, d2, β3, d3)

    @.. broadcast=false linsolve.b=α1 * uprev + α2 * uprevprev + α3 * uprev3

    basic_patankar_step!(u, linsolve.b, P3, d3, linsolve.A, σ, dt, linsolve)

    # statistics: 3 nsolve

    return nothing
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

        # compute initial values for MPLM33
        v, t, nf, ns = start_MPLM33_oop(P, d, t, dt, uprev, f, p, small_constant,
                                        alg.linsolve)
        integrator.stats.nf += nf
        integrator.stats.nsolve += ns

        # u at time tspan[1] + dt
        u = v[1]

        cache.uprevprev = uprev

        # we use uprev3 as temporary storage for the value of u needed in step 2.
        cache.uprev3 = v[2]

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
        # increase step count
        cache.step += 1

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

@muladd function perform_step!(integrator, cache::MPLM33Cache, repeat_step = false)
    (; alg, t, dt, uprev, uprev2, u, f, p) = integrator
    (; uprevprev, uprev3, v, vprev, vprev2, P, P2, P3, d, d2, d3, σ, αβ, small_constant, linsolve) = cache

    #TODO: is this necessary?
    if integrator.u_modified
        cache.step = 1
    end

    if cache.step == 1
        # increase step count
        cache.step += 1

        # initilialze vprev
        vprev .= uprev

        # evaluate production matrix at tspan[1]
        evaluate_pds!(P, d, f, uprev, p, t)
        integrator.stats.nf += 1

        # save current P
        P3 .= P

        # compute initial values for MPLM33

        # we use uprev3 as temporary storage for the value of u needed in step 1.
        nf, ns = start_MPLM33!(v, uprev3, P, P2, d, d2, t, dt, vprev, vprev2, σ, f, p,
                               small_constant,
                               linsolve)
        integrator.stats.nf += nf
        integrator.stats.nsolve += ns

        # reset P
        P .= P3

        # u at time tspan[1] + dt
        u .= uprev3

        # we use uprev3 as temporary storage for the value of u needed in step 2.
        uprev3 .= v

        uprevprev .= uprev

    elseif cache.step == 2
        # increase step count
        cache.step += 1

        # evaluate production matrix at tspan[1] + dt
        evaluate_pds!(P, d, f, uprev, p, t)
        integrator.stats.nf += 1

        # u at time tspan[1] + 2*dt (this was computed in step 1)
        u .= uprev3

        uprev3 .= uprevprev
        uprevprev .= uprev

    else
        # increase step count
        cache.step += 1

        # evaluate production matrix
        evaluate_pds!(P, d, f, uprev, p, t)
        integrator.stats.nf += 1

        P_tup = (P, P2, P3)
        d_tup = (d, d2, d3)
        u_tup = (uprev, uprevprev, uprev3)

        perform_step_MPLM33!(u, P_tup, d_tup, dt, u_tup, σ, linsolve, αβ,
                             small_constant)
        integrator.stats.nsolve += 3

        uprev3 .= uprevprev
        uprevprev .= uprev
    end

    P3 .= P2
    P2 .= P
    !isnothing(d2) && (d3 .= d2)
    !isnothing(d) && (d2 .= d)
end

#### MPLM43 ############################################################################
@muladd function start_MPLM43_oop(P, d, t, dt, uprev, f, p, small_constant, linsolve)
    αβ33 = (0, 0, 1, 9 / 4, 0, 3 / 4)

    # 1 macro step consists of  4 substeps 
    dts = dt / 4

    ### first macro step ###############################################################
    # substeps 1 - 2
    v, t, nf, ns = start_MPLM33_oop(P, d, t, dts, uprev, f, p, small_constant, linsolve)

    uprevprev = uprev
    P2 = P
    d2 = d

    uprev = v[1]
    P, d = evaluate_pds(f, uprev, p, t)
    nf += 1

    u = v[2]

    # substeps 3-4
    for _ in 1:2
        uprev3 = uprevprev
        uprevprev = uprev
        uprev = u

        P3 = P2
        P2 = P
        d3 = d2
        d2 = d

        P, d = evaluate_pds(f, uprev, p, t)
        nf += 1

        P_tup = (P, P2, P3)
        d_tup = (d, d2, d3)
        u_tup = (uprev, uprevprev, uprev3)

        u = perform_step_MPLM33_oop(P_tup, d_tup, dts, u_tup, linsolve, αβ33,
                                    small_constant)
        t += dts
        ns += 3
    end

    v1 = u

    ### second macro step ############################################################
    for _ in 1:4
        uprev3 = uprevprev
        uprevprev = uprev
        uprev = u

        P3 = P2
        P2 = P
        d3 = d2
        d2 = d

        P, d = evaluate_pds(f, uprev, p, t)
        nf += 1

        P_tup = (P, P2, P3)
        d_tup = (d, d2, d3)
        u_tup = (uprev, uprevprev, uprev3)

        u = perform_step_MPLM33_oop(P_tup, d_tup, dts, u_tup, linsolve, αβ33,
                                    small_constant)
        t += dts
        ns += 3
    end

    v2 = u

    ### third macro step ############################################################
    for _ in 1:4
        uprev3 = uprevprev
        uprevprev = uprev
        uprev = u

        P3 = P2
        P2 = P
        d3 = d2
        d2 = d

        P, d = evaluate_pds(f, uprev, p, t)
        nf += 1

        P_tup = (P, P2, P3)
        d_tup = (d, d2, d3)
        u_tup = (uprev, uprevprev, uprev3)

        u = perform_step_MPLM33_oop(P_tup, d_tup, dts, u_tup, linsolve, αβ33,
                                    small_constant)
        t += dts
        ns += 3
    end

    return (v1, v2, u), t, nf, ns
end
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

        # compute initial values for MPLM43
        v, _, nf, ns = start_MPLM43_oop(P, d, t, dt, uprev, f, p, small_constant,
                                        alg.linsolve)
        integrator.stats.nf += nf
        integrator.stats.nsolve += ns

        # u at time tspan[1] + dt
        u = v[1]

        cache.uprevprev = uprev

        # we use uprev3 as temporary storage for the value of u needed in step 2.
        cache.uprev3 = v[2]
        # we use uprev4 as temporary storage for the value of u needed in step 3.
        cache.uprev4 = v[3]

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
#TODO Check nf and ns everywhere!
@muladd function start_MPLM54_oop(P, d, t, dt, uprev, f, p, small_constant, linsolve)
    αβ43 = (1 / 4, 0, 3 / 4, 0, 35 / 18, 1 / 3, 0, 2 / 9)

    # 1 macro step consists of 4 substeps
    dts = dt / 4

    ### first macro step ###############################################################
    # substep 1 - 3
    v, t, nf, ns = start_MPLM43_oop(P, d, t, dts, uprev, f, p, small_constant, linsolve)

    uprev3 = uprev
    P3 = P
    d3 = d

    uprevprev = v[1]
    P2, d2 = evaluate_pds(f, uprevprev, p, t)
    nf += 1

    uprev = v[2]
    P, d = evaluate_pds(f, uprev, p, t)
    nf += 1

    u = v[3]

    # substep 4
    uprev4 = uprev3
    uprev3 = uprevprev
    uprevprev = uprev
    uprev = u

    P4 = P3
    P3 = P2
    P2 = P
    d4 = d3
    d3 = d2
    d2 = d

    P, d = evaluate_pds(f, uprev, p, t)
    nf += 1

    P_tup = (P, P2, P3, P4)
    d_tup = (d, d2, d3, d4)
    u_tup = (uprev, uprevprev, uprev3, uprev4)

    u = perform_step_MPLM43_oop(P_tup, d_tup, dts, u_tup, linsolve, αβ43,
                                small_constant)
    t += dts
    ns += 3

    v1 = u

    ### second macro step ############################################################
    for _ in 1:4
        uprev4 = uprev3
        uprev3 = uprevprev
        uprevprev = uprev
        uprev = u

        P4 = P3
        P3 = P2
        P2 = P
        d4 = d3
        d3 = d2
        d2 = d

        P, d = evaluate_pds(f, uprev, p, t)
        nf += 1

        P_tup = (P, P2, P3, P4)
        d_tup = (d, d2, d3, d4)
        u_tup = (uprev, uprevprev, uprev3, uprev4)

        u = perform_step_MPLM43_oop(P_tup, d_tup, dts, u_tup, linsolve, αβ43,
                                    small_constant)
        t += dts
        ns += 3
    end

    v2 = u

    ### third macro step ############################################################
    for _ in 1:4
        uprev4 = uprev3
        uprev3 = uprevprev
        uprevprev = uprev
        uprev = u

        P4 = P3
        P3 = P2
        P2 = P
        d4 = d3
        d3 = d2
        d2 = d

        P, d = evaluate_pds(f, uprev, p, t)
        nf += 1

        P_tup = (P, P2, P3, P4)
        d_tup = (d, d2, d3, d4)
        u_tup = (uprev, uprevprev, uprev3, uprev4)

        u = perform_step_MPLM43_oop(P_tup, d_tup, dts, u_tup, linsolve, αβ43,
                                    small_constant)
        t += dts
        ns += 3
    end

    v3 = u

    # fourth macro step
    for _ in 1:4
        uprev4 = uprev3
        uprev3 = uprevprev
        uprevprev = uprev
        uprev = u

        P4 = P3
        P3 = P2
        P2 = P
        d4 = d3
        d3 = d2
        d2 = d

        P, d = evaluate_pds(f, uprev, p, t)
        nf += 1

        P_tup = (P, P2, P3, P4)
        d_tup = (d, d2, d3, d4)
        u_tup = (uprev, uprevprev, uprev3, uprev4)

        u = perform_step_MPLM43_oop(P_tup, d_tup, dts, u_tup, linsolve, αβ43,
                                    small_constant)
        t += dts
        ns += 3
    end

    return (v1, v2, v3, u), t, nf, ns
end
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

        # compute initial values for MPLM54
        v, _, nf, ns = start_MPLM54_oop(P, d, t, dt, uprev, f, p, small_constant,
                                        alg.linsolve)
        integrator.stats.nf += nf
        integrator.stats.nsolve += ns

        # u at time tspan[1] + dt
        u = v[1]

        cache.uprevprev = uprev

        # we use uprev3 as temporary storage for the value of u needed in step 2.
        cache.uprev3 = v[2]
        # we use uprev4 as temporary storage for the value of u needed in step 3.
        cache.uprev4 = v[3]
        # we use uprev5 as temporary storage for the value of u needed in step 4.
        cache.uprev5 = v[4]

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
@muladd function start_MPLM75_oop(P, d, t, dt, uprev, f, p, small_constant, linsolve)
    αβ54 = (0, 0, 0, 0, 1, 225 / 96, 0, 50 / 96, 200 / 96, 5 / 96)

    # 1 macro steps consists of 4 substeps
    dts = dt / 4

    ### first macro step ###############################################################
    # substeps 1 - 4
    v, t, nf, ns = start_MPLM54_oop(P, d, t, dts, uprev, f, p, small_constant, linsolve)

    uprev4 = uprev
    P4 = P
    d4 = d

    uprev3 = v[1]
    P3, d3 = evaluate_pds(f, uprev3, p, t)
    nf += 1

    uprevprev = v[2]
    P2, d2 = evaluate_pds(f, uprevprev, p, t)
    nf += 1

    uprev = v[3]
    P, d = evaluate_pds(f, uprev, p, t)
    nf += 1

    u = v[4]

    v1 = u

    ### second macro step ############################################################
    for _ in 1:4
        uprev5 = uprev4
        uprev4 = uprev3
        uprev3 = uprevprev
        uprevprev = uprev
        uprev = u

        P5 = P4
        P4 = P3
        P3 = P2
        P2 = P
        d5 = d4
        d4 = d3
        d3 = d2
        d2 = d

        P, d = evaluate_pds(f, uprev, p, t)
        nf += 1

        P_tup = (P, P2, P3, P4, P5)
        d_tup = (d, d2, d3, d4, d5)
        u_tup = (uprev, uprevprev, uprev3, uprev4, uprev5)

        u = perform_step_MPLM54_oop(P_tup, d_tup, dts, u_tup, linsolve, αβ54,
                                    small_constant)
        t += dts
        ns += 4
    end

    v2 = u

    ### third macro step ############################################################
    for _ in 1:4
        uprev5 = uprev4
        uprev4 = uprev3
        uprev3 = uprevprev
        uprevprev = uprev
        uprev = u

        P5 = P4
        P4 = P3
        P3 = P2
        P2 = P
        d5 = d4
        d4 = d3
        d3 = d2
        d2 = d

        P, d = evaluate_pds(f, uprev, p, t)
        nf += 1

        P_tup = (P, P2, P3, P4, P5)
        d_tup = (d, d2, d3, d4, d5)
        u_tup = (uprev, uprevprev, uprev3, uprev4, uprev5)

        u = perform_step_MPLM54_oop(P_tup, d_tup, dts, u_tup, linsolve, αβ54,
                                    small_constant)
        t += dts
        ns += 4
    end

    v3 = u

    ### fourth macro step ############################################################
    for _ in 1:4
        uprev5 = uprev4
        uprev4 = uprev3
        uprev3 = uprevprev
        uprevprev = uprev
        uprev = u

        P5 = P4
        P4 = P3
        P3 = P2
        P2 = P
        d5 = d4
        d4 = d3
        d3 = d2
        d2 = d

        P, d = evaluate_pds(f, uprev, p, t)
        nf += 1

        P_tup = (P, P2, P3, P4, P5)
        d_tup = (d, d2, d3, d4, d5)
        u_tup = (uprev, uprevprev, uprev3, uprev4, uprev5)

        u = perform_step_MPLM54_oop(P_tup, d_tup, dts, u_tup, linsolve, αβ54,
                                    small_constant)
        t += dts
        ns += 4
    end

    v4 = u

    ### fifth macro step ############################################################
    for _ in 1:4
        uprev5 = uprev4
        uprev4 = uprev3
        uprev3 = uprevprev
        uprevprev = uprev
        uprev = u

        P5 = P4
        P4 = P3
        P3 = P2
        P2 = P
        d5 = d4
        d4 = d3
        d3 = d2
        d2 = d

        P, d = evaluate_pds(f, uprev, p, t)
        nf += 1

        P_tup = (P, P2, P3, P4, P5)
        d_tup = (d, d2, d3, d4, d5)
        u_tup = (uprev, uprevprev, uprev3, uprev4, uprev5)

        u = perform_step_MPLM54_oop(P_tup, d_tup, dts, u_tup, linsolve, αβ54,
                                    small_constant)
        t += dts
        ns += 4
    end

    v5 = u

    ### sixth macro step ############################################################
    for _ in 1:4
        uprev5 = uprev4
        uprev4 = uprev3
        uprev3 = uprevprev
        uprevprev = uprev
        uprev = u

        P5 = P4
        P4 = P3
        P3 = P2
        P2 = P
        d5 = d4
        d4 = d3
        d3 = d2
        d2 = d

        P, d = evaluate_pds(f, uprev, p, t)
        nf += 1

        P_tup = (P, P2, P3, P4, P5)
        d_tup = (d, d2, d3, d4, d5)
        u_tup = (uprev, uprevprev, uprev3, uprev4, uprev5)

        u = perform_step_MPLM54_oop(P_tup, d_tup, dts, u_tup, linsolve, αβ54,
                                    small_constant)
        t += dts
        ns += 4
    end

    return (v1, v2, v3, v4, v5, u), t, nf, ns
end

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

        # compute initial values for MPLM75
        v, _, nf, ns = start_MPLM75_oop(P, d, t, dt, uprev, f, p, small_constant,
                                        alg.linsolve)
        integrator.stats.nf += nf
        integrator.stats.nsolve += ns

        # u at time tspan[1] + dt
        u = v[1]

        cache.uprevprev = uprev

        # we use uprev3 as temporary storage for the value of u needed in step 2.
        cache.uprev3 = v[2]
        # we use uprev4 as temporary storage for the value of u needed in step 3.
        cache.uprev4 = v[3]
        # we use uprev5 as temporary storage for the value of u needed in step 4.
        cache.uprev5 = v[4]
        # we use uprev6 as temporary storage for the value of u needed in step 5.
        cache.uprev6 = v[5]
        # we use uprev7 as temporary storage for the value of u needed in step 6.
        cache.uprev7 = v[6]

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
@muladd function start_MPLM106_oop(P, d, t, dt, uprev, f, p, small_constant, linsolve)
    αβ75 = (0, 0, 0, 0, 0, 0, 1, 12 / 5, 0, 197 / 720, 701 / 360, 43 / 30, 107 / 360,
            467 / 720)

    # 1 macro step consists of 4 substeps            
    dts = dt / 4

    ### 1.5 macro steps ###############################################################
    v, t, nf, ns = start_MPLM75_oop(P, d, t, dts, uprev, f, p, small_constant, linsolve)

    uprev6 = uprev
    P6 = P
    d6 = d

    uprev5 = v[1]
    P5, d5 = evaluate_pds(f, uprev5, p, t)
    nf += 1

    uprev4 = v[2]
    P4, d4 = evaluate_pds(f, uprev4, p, t)
    nf += 1

    uprev3 = v[3]
    P3, d3 = evaluate_pds(f, uprev3, p, t)
    nf += 1

    uprevprev = v[4]
    P2, d2 = evaluate_pds(f, uprevprev, p, t)
    nf += 1

    uprev = v[5]
    P, d = evaluate_pds(f, uprev, p, t)
    nf += 1

    u = v[6]

    v1 = v[4]

    ### second macro step ############################################################
    # substeps 3 - 4
    for _ in 1:2
        uprev7 = uprev6
        uprev6 = uprev5
        uprev5 = uprev4
        uprev4 = uprev3
        uprev3 = uprevprev
        uprevprev = uprev
        uprev = u

        P7 = P6
        P6 = P5
        P5 = P4
        P4 = P3
        P3 = P2
        P2 = P
        d7 = d6
        d6 = d5
        d5 = d4
        d4 = d3
        d3 = d2
        d2 = d

        P, d = evaluate_pds(f, uprev, p, t)
        nf += 1

        P_tup = (P, P2, P3, P4, P5, P6, P7)
        d_tup = (d, d2, d3, d4, d5, d6, d7)
        u_tup = (uprev, uprevprev, uprev3, uprev4, uprev5, uprev6, uprev7)

        u = perform_step_MPLM75_oop(P_tup, d_tup, dts, u_tup, linsolve, αβ75,
                                    small_constant)
        t += dts
        ns += 5
    end

    v2 = u

    ### third macro step ############################################################
    for _ in 1:4
        uprev7 = uprev6
        uprev6 = uprev5
        uprev5 = uprev4
        uprev4 = uprev3
        uprev3 = uprevprev
        uprevprev = uprev
        uprev = u

        P7 = P6
        P6 = P5
        P5 = P4
        P4 = P3
        P3 = P2
        P2 = P
        d7 = d6
        d6 = d5
        d5 = d4
        d4 = d3
        d3 = d2
        d2 = d

        P, d = evaluate_pds(f, uprev, p, t)
        nf += 1

        P_tup = (P, P2, P3, P4, P5, P6, P7)
        d_tup = (d, d2, d3, d4, d5, d6, d7)
        u_tup = (uprev, uprevprev, uprev3, uprev4, uprev5, uprev6, uprev7)

        u = perform_step_MPLM75_oop(P_tup, d_tup, dts, u_tup, linsolve, αβ75,
                                    small_constant)
        t += dts
        ns += 5
    end

    v3 = u

    ### fourth macro step ############################################################
    for _ in 1:4
        uprev7 = uprev6
        uprev6 = uprev5
        uprev5 = uprev4
        uprev4 = uprev3
        uprev3 = uprevprev
        uprevprev = uprev
        uprev = u

        P7 = P6
        P6 = P5
        P5 = P4
        P4 = P3
        P3 = P2
        P2 = P
        d7 = d6
        d6 = d5
        d5 = d4
        d4 = d3
        d3 = d2
        d2 = d

        P, d = evaluate_pds(f, uprev, p, t)
        nf += 1

        P_tup = (P, P2, P3, P4, P5, P6, P7)
        d_tup = (d, d2, d3, d4, d5, d6, d7)
        u_tup = (uprev, uprevprev, uprev3, uprev4, uprev5, uprev6, uprev7)

        u = perform_step_MPLM75_oop(P_tup, d_tup, dts, u_tup, linsolve, αβ75,
                                    small_constant)
        t += dts
        ns += 5
    end

    v4 = u

    ### fifth macro step ############################################################
    for _ in 1:4
        uprev7 = uprev6
        uprev6 = uprev5
        uprev5 = uprev4
        uprev4 = uprev3
        uprev3 = uprevprev
        uprevprev = uprev
        uprev = u

        P7 = P6
        P6 = P5
        P5 = P4
        P4 = P3
        P3 = P2
        P2 = P
        d7 = d6
        d6 = d5
        d5 = d4
        d4 = d3
        d3 = d2
        d2 = d

        P, d = evaluate_pds(f, uprev, p, t)
        nf += 1

        P_tup = (P, P2, P3, P4, P5, P6, P7)
        d_tup = (d, d2, d3, d4, d5, d6, d7)
        u_tup = (uprev, uprevprev, uprev3, uprev4, uprev5, uprev6, uprev7)

        u = perform_step_MPLM75_oop(P_tup, d_tup, dts, u_tup, linsolve, αβ75,
                                    small_constant)
        t += dts
        ns += 5
    end

    v5 = u

    ### sixth macro step ############################################################
    for _ in 1:4
        uprev7 = uprev6
        uprev6 = uprev5
        uprev5 = uprev4
        uprev4 = uprev3
        uprev3 = uprevprev
        uprevprev = uprev
        uprev = u

        P7 = P6
        P6 = P5
        P5 = P4
        P4 = P3
        P3 = P2
        P2 = P
        d7 = d6
        d6 = d5
        d5 = d4
        d4 = d3
        d3 = d2
        d2 = d

        P, d = evaluate_pds(f, uprev, p, t)
        nf += 1

        P_tup = (P, P2, P3, P4, P5, P6, P7)
        d_tup = (d, d2, d3, d4, d5, d6, d7)
        u_tup = (uprev, uprevprev, uprev3, uprev4, uprev5, uprev6, uprev7)

        u = perform_step_MPLM75_oop(P_tup, d_tup, dts, u_tup, linsolve, αβ75,
                                    small_constant)
        t += dts
        ns += 5
    end

    v6 = u

    ### seventh macro step ############################################################
    for _ in 1:4
        uprev7 = uprev6
        uprev6 = uprev5
        uprev5 = uprev4
        uprev4 = uprev3
        uprev3 = uprevprev
        uprevprev = uprev
        uprev = u

        P7 = P6
        P6 = P5
        P5 = P4
        P4 = P3
        P3 = P2
        P2 = P
        d7 = d6
        d6 = d5
        d5 = d4
        d4 = d3
        d3 = d2
        d2 = d

        P, d = evaluate_pds(f, uprev, p, t)
        nf += 1

        P_tup = (P, P2, P3, P4, P5, P6, P7)
        d_tup = (d, d2, d3, d4, d5, d6, d7)
        u_tup = (uprev, uprevprev, uprev3, uprev4, uprev5, uprev6, uprev7)

        u = perform_step_MPLM75_oop(P_tup, d_tup, dts, u_tup, linsolve, αβ75,
                                    small_constant)
        t += dts
        ns += 5
    end

    v7 = u

    ### eighth macro step ############################################################
    for _ in 1:4
        uprev7 = uprev6
        uprev6 = uprev5
        uprev5 = uprev4
        uprev4 = uprev3
        uprev3 = uprevprev
        uprevprev = uprev
        uprev = u

        P7 = P6
        P6 = P5
        P5 = P4
        P4 = P3
        P3 = P2
        P2 = P
        d7 = d6
        d6 = d5
        d5 = d4
        d4 = d3
        d3 = d2
        d2 = d

        P, d = evaluate_pds(f, uprev, p, t)
        nf += 1

        P_tup = (P, P2, P3, P4, P5, P6, P7)
        d_tup = (d, d2, d3, d4, d5, d6, d7)
        u_tup = (uprev, uprevprev, uprev3, uprev4, uprev5, uprev6, uprev7)

        u = perform_step_MPLM75_oop(P_tup, d_tup, dts, u_tup, linsolve, αβ75,
                                    small_constant)
        t += dts
        ns += 5
    end

    v8 = u

    ### ninth macro step ############################################################
    for _ in 1:4
        uprev7 = uprev6
        uprev6 = uprev5
        uprev5 = uprev4
        uprev4 = uprev3
        uprev3 = uprevprev
        uprevprev = uprev
        uprev = u

        P7 = P6
        P6 = P5
        P5 = P4
        P4 = P3
        P3 = P2
        P2 = P
        d7 = d6
        d6 = d5
        d5 = d4
        d4 = d3
        d3 = d2
        d2 = d

        P, d = evaluate_pds(f, uprev, p, t)
        nf += 1

        P_tup = (P, P2, P3, P4, P5, P6, P7)
        d_tup = (d, d2, d3, d4, d5, d6, d7)
        u_tup = (uprev, uprevprev, uprev3, uprev4, uprev5, uprev6, uprev7)

        u = perform_step_MPLM75_oop(P_tup, d_tup, dts, u_tup, linsolve, αβ75,
                                    small_constant)
        t += dts
        ns += 5
    end

    return (v1, v2, v3, v4, v5, v6, v7, v8, u), t, nf, ns
end
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

        # compute initial values for MPLM106
        v, _, nf, ns = start_MPLM106_oop(P, d, t, dt, uprev, f, p, small_constant,
                                         alg.linsolve)
        integrator.stats.nf += nf
        integrator.stats.nsolve += ns

        # u at time tspan[1] + dt
        u = v[1]

        cache.uprevprev = uprev

        # we use uprev3 as temporary storage for the value of u needed in step 2.
        cache.uprev3 = v[2]
        # we use uprev4 as temporary storage for the value of u needed in step 3.
        cache.uprev4 = v[3]
        # we use uprev5 as temporary storage for the value of u needed in step 4.
        cache.uprev5 = v[4]
        # we use uprev6 as temporary storage for the value of u needed in step 5.
        cache.uprev6 = v[5]
        # we use uprev7 as temporary storage for the value of u needed in step 6.
        cache.uprev7 = v[6]
        # we use uprev8 as temporary storage for the value of u needed in step 7.
        cache.uprev8 = v[7]
        # we use uprev9 as temporary storage for the value of u needed in step 8.
        cache.uprev9 = v[8]
        # we use uprev10 as temporary storage for the value of u needed in step 9.
        cache.uprev10 = v[9]

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
