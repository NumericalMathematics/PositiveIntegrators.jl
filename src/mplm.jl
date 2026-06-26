########################################################################################
### Structs and caches #################################################################
########################################################################################
abstract type MPLMMutableCache <: OrdinaryDiffEqMutableCache end
get_fsalfirstlast(cache::MPLMMutableCache, rate_prototype) = (nothing, nothing)

#### MPLM22 ############################################################################
struct MPLM22{F, T} <: OrdinaryDiffEqAlgorithm
    substep_level::Int
    linsolve::F
    small_constant_function::T
end

alg_order(alg::MPLM22) = 2
isfsal(::MPLM22) = false

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

function MPLM22(n)
    if !(n isa Integer && n ≥ 1)
        throw(ArgumentError("MPLM22 requires a positive integer for the substep level."))
    end
    MPLM22(; substep_level = n)
end

function MPLM22(; substep_level = 1, linsolve = LUFactorization(), small_constant = nothing)
    if isnothing(small_constant)
        small_constant_function = floatmin
    elseif small_constant isa Number
        small_constant_function = Returns(small_constant)
    else # assume small_constant isa Function
        small_constant_function = small_constant
    end
    MPLM22(substep_level, linsolve, small_constant_function)
end

function alg_cache(alg::MPLM22, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{false},
                   verbose) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    if !(f isa PDSFunction || f isa ConservativePDSFunction)
        throw(ArgumentError("MPLM22 can only be applied to production-destruction systems"))
    end
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
    substep_level::Int
    linsolve::F
    small_constant_function::T
end

alg_order(alg::MPLM33) = 3
isfsal(::MPLM33) = false

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

get_tmp_cache(integrator, ::MPLM33, cache::MPLMMutableCache) = (cache.σ,)

function MPLM33(n)
    if !(n isa Integer && n ≥ 1)
        throw(ArgumentError("MPLM33 requires a positive integer for the substep level."))
    end
    MPLM33(; substep_level = n)
end

function MPLM33(; substep_level = 1, linsolve = LUFactorization(), small_constant = nothing)
    if isnothing(small_constant)
        small_constant_function = floatmin
    elseif small_constant isa Number
        small_constant_function = Returns(small_constant)
    else # assume small_constant isa Function
        small_constant_function = small_constant
    end
    MPLM33(substep_level, linsolve, small_constant_function)
end

function get_constant_parameters(alg::MPLM33, ::Type{uEltypeNoUnits}) where {uEltypeNoUnits}
    α1 = zero(uEltypeNoUnits)
    α2 = zero(uEltypeNoUnits)
    α3 = one(uEltypeNoUnits)
    β1 = (9 * one(uEltypeNoUnits)) / 4
    β2 = zero(uEltypeNoUnits)
    β3 = (3 * one(uEltypeNoUnits)) / 4
    return (α1, α2, α3, β1, β2, β3)
end

function alg_cache(alg::MPLM33, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{false},
                   verbose) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    if !(f isa PDSFunction || f isa ConservativePDSFunction)
        throw(ArgumentError("MPLM33 can only be applied to production-destruction systems"))
    end
    αβ = get_constant_parameters(alg, uEltypeNoUnits)

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
    αβ = get_constant_parameters(alg, uEltypeNoUnits)

    if f isa ConservativePDSFunction
        linprob = LinearProblem(A, _vec(b))
        linsolve = init(linprob, alg.linsolve,
                        alias = LinearSolve.LinearAliasSpecifier(; alias_A = true,
                                                                 alias_b = true),
                        assumptions = LinearSolve.OperatorAssumptions(true))

        MPLM33Cache(uprevprev, uprev3, v, vprev, vprev2, step, small_constant, b, P, P2, P3,
                    A, nothing, nothing, nothing, σ,
                    linsolve, αβ)
    elseif f isa PDSFunction
        linprob = LinearProblem(A, _vec(b))
        linsolve = init(linprob, alg.linsolve,
                        alias = LinearSolve.LinearAliasSpecifier(; alias_A = true,
                                                                 alias_b = true),
                        assumptions = LinearSolve.OperatorAssumptions(true))

        MPLM33Cache(uprevprev, uprev3, v, vprev, vprev2, step, small_constant, b, P, P2, P3,
                    A,
                    similar(u), similar(u), similar(u),
                    σ, linsolve, αβ)
    else
        throw(ArgumentError("MPLM33 can only be applied to production-destruction systems"))
    end
end

#### MPLM43 ############################################################################
struct MPLM43{F, T} <: OrdinaryDiffEqAlgorithm
    substep_level::Int
    linsolve::F
    small_constant_function::T
end

alg_order(alg::MPLM43) = 3
isfsal(::MPLM43) = false

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

@cache mutable struct MPLM43Cache{uType, dType, T, PType, F, TabType} <: MPLMMutableCache
    uprevprev::uType
    uprev3::uType
    uprev4::uType
    v::uType
    vprev::uType
    vprev2::uType
    vprev3::uType
    step::Int
    small_constant::T
    b::uType # rhs of the linear system
    P::PType
    P2::PType
    P3::PType
    P4::PType
    A::PType # system matrix of the linear system
    d::dType
    d2::dType
    d3::dType
    d4::dType
    σ::uType
    linsolve::F
    αβ::TabType
end

get_tmp_cache(integrator, ::MPLM43, cache::MPLMMutableCache) = (cache.σ,)

function MPLM43(n)
    if !(n isa Integer && n ≥ 1)
        throw(ArgumentError("MPLM43 requires a positive integer for the substep level."))
    end
    MPLM43(; substep_level = n)
end

function MPLM43(; substep_level = 1, linsolve = LUFactorization(), small_constant = nothing)
    if isnothing(small_constant)
        small_constant_function = floatmin
    elseif small_constant isa Number
        small_constant_function = Returns(small_constant)
    else # assume small_constant isa Function
        small_constant_function = small_constant
    end
    MPLM43(substep_level, linsolve, small_constant_function)
end

function get_constant_parameters(alg::MPLM43, ::Type{uEltypeNoUnits}) where {uEltypeNoUnits}
    α1 = one(uEltypeNoUnits) / 4
    α2 = zero(uEltypeNoUnits)
    α3 = (3 * one(uEltypeNoUnits)) / 4
    α4 = zero(uEltypeNoUnits)
    β1 = (35 * one(uEltypeNoUnits)) / 18
    β2 = one(uEltypeNoUnits) / 3
    β3 = zero(uEltypeNoUnits)
    β4 = (2 * one(uEltypeNoUnits)) / 9

    return (α1, α2, α3, α4, β1, β2, β3, β4)
end

function alg_cache(alg::MPLM43, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{false},
                   verbose) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    if !(f isa PDSFunction || f isa ConservativePDSFunction)
        throw(ArgumentError("MPLM43 can only be applied to production-destruction systems"))
    end
    # TODO: This is currently necessary to get the correct type of P (d is of type rateType)
    P, d = evaluate_pds(f, u, p, t)
    # TODO: integrator_stats_nf = 1

    αβ = get_constant_parameters(alg, uEltypeNoUnits)
    MPLM43oopCache(u, u, u, P, P, P, d, d, d, αβ, 1,
                   alg.small_constant_function(uEltypeNoUnits))
end

function alg_cache(alg::MPLM43, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
                   uprev, uprev2, f, t, dt, reltol, p, calck,
                   ::Val{true},
                   verbose) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    uprevprev = zero(u)
    uprev3 = zero(u)
    uprev4 = zero(u)
    v = zero(u)
    vprev = zero(u)
    vprev2 = zero(u)
    vprev3 = zero(u)
    step = 1
    small_constant = alg.small_constant_function(uEltypeNoUnits)
    b = zero(u)
    P = p_prototype(u, f)
    P2 = p_prototype(u, f)
    P3 = p_prototype(u, f)
    P4 = p_prototype(u, f)
    A = p_prototype(u, f)
    σ = zero(u)

    # MPLM43 coefficients 
    αβ = get_constant_parameters(alg, uEltypeNoUnits)

    if f isa ConservativePDSFunction
        linprob = LinearProblem(A, _vec(b))
        linsolve = init(linprob, alg.linsolve,
                        alias = LinearSolve.LinearAliasSpecifier(; alias_A = true,
                                                                 alias_b = true),
                        assumptions = LinearSolve.OperatorAssumptions(true))

        MPLM43Cache(uprevprev, uprev3, uprev4, v, vprev, vprev2, vprev3, step,
                    small_constant, b, P, P2, P3, P4,
                    A, nothing, nothing, nothing, nothing, σ,
                    linsolve, αβ)
    elseif f isa PDSFunction
        linprob = LinearProblem(A, _vec(b))
        linsolve = init(linprob, alg.linsolve,
                        alias = LinearSolve.LinearAliasSpecifier(; alias_A = true,
                                                                 alias_b = true),
                        assumptions = LinearSolve.OperatorAssumptions(true))

        MPLM43Cache(uprevprev, uprev3, uprev4, v, vprev, vprev2, vprev3, step,
                    small_constant, b, P, P2, P3, P4,
                    A,
                    similar(u), similar(u), similar(u), similar(u),
                    σ, linsolve, αβ)
    else
        throw(ArgumentError("MPLM43 can only be applied to production-destruction systems"))
    end
end

#### MPLM54 ############################################################################
struct MPLM54{F, T} <: OrdinaryDiffEqAlgorithm
    substep_level::Int
    linsolve::F
    small_constant_function::T
end

alg_order(alg::MPLM54) = 4
isfsal(::MPLM54) = false

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

@cache mutable struct MPLM54Cache{uType, dType, T, PType, F, TabType} <: MPLMMutableCache
    uprevprev::uType
    uprev3::uType
    uprev4::uType
    uprev5::uType
    v::uType
    vprev::uType
    vprev2::uType
    vprev3::uType
    vprev4::uType
    step::Int
    small_constant::T
    b::uType # rhs of the linear system
    P::PType
    P2::PType
    P3::PType
    P4::PType
    P5::PType
    A::PType # system matrix of the linear system
    d::dType
    d2::dType
    d3::dType
    d4::dType
    d5::dType
    d_tmp::dType
    σ::uType
    linsolve::F
    αβ::TabType
end

get_tmp_cache(integrator, ::MPLM54, cache::MPLMMutableCache) = (cache.σ,)

function MPLM54(n)
    if !(n isa Integer && n ≥ 1)
        throw(ArgumentError("MPLM54 requires a positive integer for the substep level."))
    end
    MPLM54(; substep_level = n)
end

function MPLM54(; substep_level = 1, linsolve = LUFactorization(), small_constant = nothing)
    if isnothing(small_constant)
        small_constant_function = floatmin
    elseif small_constant isa Number
        small_constant_function = Returns(small_constant)
    else # assume small_constant isa Function
        small_constant_function = small_constant
    end
    MPLM54(substep_level, linsolve, small_constant_function)
end

function get_constant_parameters(alg::MPLM54, ::Type{uEltypeNoUnits}) where {uEltypeNoUnits}
    α1 = zero(uEltypeNoUnits)
    α2 = zero(uEltypeNoUnits)
    α3 = zero(uEltypeNoUnits)
    α4 = zero(uEltypeNoUnits)
    α5 = one(uEltypeNoUnits)
    β1 = (225 * one(uEltypeNoUnits)) / 96
    β2 = zero(uEltypeNoUnits)
    β3 = (50 * one(uEltypeNoUnits)) / 96
    β4 = (200 * one(uEltypeNoUnits)) / 96
    β5 = (5 * one(uEltypeNoUnits)) / 96

    return (α1, α2, α3, α4, α5, β1, β2, β3, β4, β5)
end

function alg_cache(alg::MPLM54, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{false},
                   verbose) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    if !(f isa PDSFunction || f isa ConservativePDSFunction)
        throw(ArgumentError("MPLM54 can only be applied to production-destruction systems"))
    end
    # TODO: This is currently necessary to get the correct type of P (d is of type rateType)
    P, d = evaluate_pds(f, u, p, t)
    # TODO: integrator_stats_nf = 1

    αβ = get_constant_parameters(alg, uEltypeNoUnits)

    MPLM54oopCache(u, u, u, u, P, P, P, P, d, d, d, d, αβ, 1,
                   alg.small_constant_function(uEltypeNoUnits))
end

function alg_cache(alg::MPLM54, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
                   uprev, uprev2, f, t, dt, reltol, p, calck,
                   ::Val{true},
                   verbose) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    uprevprev = zero(u)
    uprev3 = zero(u)
    uprev4 = zero(u)
    uprev5 = zero(u)
    v = zero(u)
    vprev = zero(u)
    vprev2 = zero(u)
    vprev3 = zero(u)
    vprev4 = zero(u)
    step = 1
    small_constant = alg.small_constant_function(uEltypeNoUnits)
    b = zero(u)
    P = p_prototype(u, f)
    P2 = p_prototype(u, f)
    P3 = p_prototype(u, f)
    P4 = p_prototype(u, f)
    P5 = p_prototype(u, f)
    A = p_prototype(u, f)
    σ = zero(u)

    # MPLM54 coefficients
    αβ = get_constant_parameters(alg, uEltypeNoUnits)

    if f isa ConservativePDSFunction
        linprob = LinearProblem(A, _vec(b))
        linsolve = init(linprob, alg.linsolve,
                        alias = LinearSolve.LinearAliasSpecifier(; alias_A = true,
                                                                 alias_b = true),
                        assumptions = LinearSolve.OperatorAssumptions(true))

        MPLM54Cache(uprevprev, uprev3, uprev4, uprev5, v, vprev, vprev2, vprev3, vprev4,
                    step,
                    small_constant, b, P, P2, P3, P4, P5,
                    A, nothing, nothing, nothing, nothing, nothing, nothing, σ,
                    linsolve, αβ)
    elseif f isa PDSFunction
        linprob = LinearProblem(A, _vec(b))
        linsolve = init(linprob, alg.linsolve,
                        alias = LinearSolve.LinearAliasSpecifier(; alias_A = true,
                                                                 alias_b = true),
                        assumptions = LinearSolve.OperatorAssumptions(true))

        MPLM54Cache(uprevprev, uprev3, uprev4, uprev5, v, vprev, vprev2, vprev3, vprev4,
                    step,
                    small_constant, b, P, P2, P3, P4, P5,
                    A,
                    similar(u), similar(u), similar(u), similar(u), similar(u), similar(u),
                    σ, linsolve, αβ)
    else
        throw(ArgumentError("MPLM54 can only be applied to production-destruction systems"))
    end
end

#### MPLM75 ############################################################################
struct MPLM75{F, T} <: OrdinaryDiffEqAlgorithm
    substep_level::Int
    linsolve::F
    small_constant_function::T
end

alg_order(alg::MPLM75) = 5
isfsal(::MPLM75) = false

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

@cache mutable struct MPLM75Cache{uType, dType, T, PType, F, TabType} <: MPLMMutableCache
    uprevprev::uType
    uprev3::uType
    uprev4::uType
    uprev5::uType
    uprev6::uType
    uprev7::uType
    v::uType
    vprev::uType
    vprev2::uType
    vprev3::uType
    vprev4::uType
    vprev5::uType
    vprev6::uType
    step::Int
    small_constant::T
    b::uType # rhs of the linear system
    P::PType
    P2::PType
    P3::PType
    P4::PType
    P5::PType
    P6::PType
    P7::PType
    A::PType # system matrix of the linear system
    d::dType
    d2::dType
    d3::dType
    d4::dType
    d5::dType
    d6::dType
    d7::dType
    d_tmp::dType
    σ::uType
    linsolve::F
    αβ::TabType
end

get_tmp_cache(integrator, ::MPLM75, cache::MPLMMutableCache) = (cache.σ,)

function MPLM75(n)
    if !(n isa Integer && n ≥ 1)
        throw(ArgumentError("MPLM75 requires a positive integer for the substep level."))
    end
    MPLM75(; substep_level = n)
end

function MPLM75(; substep_level = 1, linsolve = LUFactorization(), small_constant = nothing)
    if isnothing(small_constant)
        small_constant_function = floatmin
    elseif small_constant isa Number
        small_constant_function = Returns(small_constant)
    else # assume small_constant isa Function
        small_constant_function = small_constant
    end
    MPLM75(substep_level, linsolve, small_constant_function)
end

function get_constant_parameters(alg::MPLM75, ::Type{uEltypeNoUnits}) where {uEltypeNoUnits}
    α1 = zero(uEltypeNoUnits)
    α2 = zero(uEltypeNoUnits)
    α3 = zero(uEltypeNoUnits)
    α4 = zero(uEltypeNoUnits)
    α5 = zero(uEltypeNoUnits)
    α6 = zero(uEltypeNoUnits)
    α7 = one(uEltypeNoUnits)
    β1 = (12 * one(uEltypeNoUnits)) / 5
    β2 = zero(uEltypeNoUnits)
    β3 = (197 * one(uEltypeNoUnits)) / 720
    β4 = (701 * one(uEltypeNoUnits)) / 360
    β5 = (43 * one(uEltypeNoUnits)) / 30
    β6 = (107 * one(uEltypeNoUnits)) / 360
    β7 = (467 * one(uEltypeNoUnits)) / 720

    return (α1, α2, α3, α4, α5, α6, α7, β1, β2, β3, β4, β5, β6, β7)
end

function alg_cache(alg::MPLM75, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{false},
                   verbose) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    if !(f isa PDSFunction || f isa ConservativePDSFunction)
        throw(ArgumentError("MPLM75 can only be applied to production-destruction systems"))
    end
    # TODO: This is currently necessary to get the correct type of P (d is of type rateType)
    P, d = evaluate_pds(f, u, p, t)
    # TODO: integrator_stats_nf = 1

    αβ = get_constant_parameters(alg, uEltypeNoUnits)
    MPLM75oopCache(u, u, u, u, u, u, P, P, P, P, P, P,
                   d, d, d, d, d, d, αβ, 1, alg.small_constant_function(uEltypeNoUnits))
end

function alg_cache(alg::MPLM75, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
                   uprev, uprev2, f, t, dt, reltol, p, calck,
                   ::Val{true},
                   verbose) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    uprevprev = zero(u)
    uprev3 = zero(u)
    uprev4 = zero(u)
    uprev5 = zero(u)
    uprev6 = zero(u)
    uprev7 = zero(u)
    v = zero(u)
    vprev = zero(u)
    vprev2 = zero(u)
    vprev3 = zero(u)
    vprev4 = zero(u)
    vprev5 = zero(u)
    vprev6 = zero(u)
    step = 1
    small_constant = alg.small_constant_function(uEltypeNoUnits)
    b = zero(u)
    P = p_prototype(u, f)
    P2 = p_prototype(u, f)
    P3 = p_prototype(u, f)
    P4 = p_prototype(u, f)
    P5 = p_prototype(u, f)
    P6 = p_prototype(u, f)
    P7 = p_prototype(u, f)
    A = p_prototype(u, f)
    σ = zero(u)

    # MPLM75 coefficients
    αβ = get_constant_parameters(alg, uEltypeNoUnits)

    if f isa ConservativePDSFunction
        linprob = LinearProblem(A, _vec(b))
        linsolve = init(linprob, alg.linsolve,
                        alias = LinearSolve.LinearAliasSpecifier(; alias_A = true,
                                                                 alias_b = true),
                        assumptions = LinearSolve.OperatorAssumptions(true))

        MPLM75Cache(uprevprev, uprev3, uprev4, uprev5, uprev6, uprev7, v, vprev, vprev2,
                    vprev3, vprev4, vprev5, vprev6,
                    step,
                    small_constant, b, P, P2, P3, P4, P5, P6, P7,
                    A, nothing, nothing, nothing, nothing, nothing, nothing, nothing,
                    nothing, σ,
                    linsolve, αβ)
    elseif f isa PDSFunction
        linprob = LinearProblem(A, _vec(b))
        linsolve = init(linprob, alg.linsolve,
                        alias = LinearSolve.LinearAliasSpecifier(; alias_A = true,
                                                                 alias_b = true),
                        assumptions = LinearSolve.OperatorAssumptions(true))

        MPLM75Cache(uprevprev, uprev3, uprev4, uprev5, uprev6, uprev7, v, vprev, vprev2,
                    vprev3, vprev4, vprev5, vprev6,
                    step,
                    small_constant, b, P, P2, P3, P4, P5, P6, P7,
                    A,
                    similar(u), similar(u), similar(u), similar(u), similar(u), similar(u),
                    similar(u), similar(u),
                    σ, linsolve, αβ)
    else
        throw(ArgumentError("MPLM75 can only be applied to production-destruction systems"))
    end
end
#### MPLM106 ###########################################################################
struct MPLM106{F, T} <: OrdinaryDiffEqAlgorithm
    substep_level::Int
    linsolve::F
    small_constant_function::T
end

alg_order(alg::MPLM106) = 6
isfsal(::MPLM106) = false

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

@cache mutable struct MPLM106Cache{uType, dType, T, PType, F, TabType} <: MPLMMutableCache
    uprevprev::uType
    uprev3::uType
    uprev4::uType
    uprev5::uType
    uprev6::uType
    uprev7::uType
    uprev8::uType
    uprev9::uType
    uprev10::uType
    v::uType
    vprev::uType
    vprev2::uType
    vprev3::uType
    vprev4::uType
    vprev5::uType
    vprev6::uType
    vprev7::uType
    step::Int
    small_constant::T
    b::uType # rhs of the linear system
    P::PType
    P2::PType
    P3::PType
    P4::PType
    P5::PType
    P6::PType
    P7::PType
    P8::PType
    P9::PType
    P10::PType
    A::PType # system matrix of the linear system
    d::dType
    d2::dType
    d3::dType
    d4::dType
    d5::dType
    d6::dType
    d7::dType
    d8::dType
    d9::dType
    d10::dType
    d_tmp::dType
    σ::uType
    linsolve::F
    αβ::TabType
end

get_tmp_cache(integrator, ::MPLM106, cache::MPLMMutableCache) = (cache.σ,)

function MPLM106(n)
    if !(n isa Integer && n ≥ 1)
        throw(ArgumentError("MPLM106 requires a positive integer for the substep level."))
    end
    MPLM106(; substep_level = n)
end

function MPLM106(; substep_level = 1, linsolve = LUFactorization(),
                 small_constant = nothing)
    if isnothing(small_constant)
        small_constant_function = floatmin
    elseif small_constant isa Number
        small_constant_function = Returns(small_constant)
    else # assume small_constant isa Function
        small_constant_function = small_constant
    end
    MPLM106(substep_level, linsolve, small_constant_function)
end

function get_constant_parameters(alg::MPLM106,
                                 ::Type{uEltypeNoUnits}) where {uEltypeNoUnits}
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

    β1 = (11125 * one(uEltypeNoUnits)) / 4536
    β2 = zero(uEltypeNoUnits)
    β3 = zero(uEltypeNoUnits)
    β4 = (50 * one(uEltypeNoUnits)) / 27
    β5 = (85 * one(uEltypeNoUnits)) / 36
    β6 = zero(uEltypeNoUnits)
    β7 = zero(uEltypeNoUnits)
    β8 = (125 * one(uEltypeNoUnits)) / 63
    β9 = (25 * one(uEltypeNoUnits)) / 24
    β10 = (25 * one(uEltypeNoUnits)) / 81

    return (α1, α2, α3, α4, α5, α6, α7, α8, α9, α10,
            β1, β2, β3, β4, β5, β6, β7, β8, β9, β10)
end

function alg_cache(alg::MPLM106, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{false},
                   verbose) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    if !(f isa PDSFunction || f isa ConservativePDSFunction)
        throw(ArgumentError("MPLM106 can only be applied to production-destruction systems"))
    end
    P = get_inplace_p_prototype(u, uEltypeNoUnits)
    d = get_inplace_d_prototype(u, f)

    αβ = get_constant_parameters(alg, uEltypeNoUnits)
    MPLM106oopCache(u, u, u, u, u, u, u, u, u,
                    P, P, P, P, P, P, P, P, P,
                    d, d, d, d, d, d, d, d, d,
                    αβ, 1, alg.small_constant_function(uEltypeNoUnits))
end

function alg_cache(alg::MPLM106, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
                   uprev, uprev2, f, t, dt, reltol, p, calck,
                   ::Val{true},
                   verbose) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    uprevprev = zero(u)
    uprev3 = zero(u)
    uprev4 = zero(u)
    uprev5 = zero(u)
    uprev6 = zero(u)
    uprev7 = zero(u)
    uprev8 = zero(u)
    uprev9 = zero(u)
    uprev10 = zero(u)
    v = zero(u)
    vprev = zero(u)
    vprev2 = zero(u)
    vprev3 = zero(u)
    vprev4 = zero(u)
    vprev5 = zero(u)
    vprev6 = zero(u)
    vprev7 = zero(u)
    step = 1
    small_constant = alg.small_constant_function(uEltypeNoUnits)
    b = zero(u)
    P = p_prototype(u, f)
    P2 = p_prototype(u, f)
    P3 = p_prototype(u, f)
    P4 = p_prototype(u, f)
    P5 = p_prototype(u, f)
    P6 = p_prototype(u, f)
    P7 = p_prototype(u, f)
    P8 = p_prototype(u, f)
    P9 = p_prototype(u, f)
    P10 = p_prototype(u, f)
    A = p_prototype(u, f)
    σ = zero(u)

    # MPLM106 coefficients
    αβ = get_constant_parameters(alg, uEltypeNoUnits)

    if f isa ConservativePDSFunction
        linprob = LinearProblem(A, _vec(b))
        linsolve = init(linprob, alg.linsolve,
                        alias = LinearSolve.LinearAliasSpecifier(; alias_A = true,
                                                                 alias_b = true),
                        assumptions = LinearSolve.OperatorAssumptions(true))

        MPLM106Cache(uprevprev, uprev3, uprev4, uprev5, uprev6, uprev7, uprev8, uprev9,
                     uprev10, v, vprev, vprev2,
                     vprev3, vprev4, vprev5, vprev6, vprev7,
                     step,
                     small_constant, b, P, P2, P3, P4, P5, P6, P7, P8, P9, P10,
                     A, nothing, nothing, nothing, nothing, nothing, nothing, nothing,
                     nothing, nothing, nothing, nothing,
                     σ, linsolve, αβ)
    elseif f isa PDSFunction
        linprob = LinearProblem(A, _vec(b))
        linsolve = init(linprob, alg.linsolve,
                        alias = LinearSolve.LinearAliasSpecifier(; alias_A = true,
                                                                 alias_b = true),
                        assumptions = LinearSolve.OperatorAssumptions(true))

        MPLM106Cache(uprevprev, uprev3, uprev4, uprev5, uprev6, uprev7, uprev8, uprev9,
                     uprev10, v, vprev, vprev2,
                     vprev3, vprev4, vprev5, vprev6, vprev7,
                     step,
                     small_constant, b, P, P2, P3, P4, P5, P6, P7, P8, P9, P10,
                     A,
                     similar(u), similar(u), similar(u), similar(u), similar(u), similar(u),
                     similar(u), similar(u), similar(u), similar(u), similar(u),
                     σ, linsolve, αβ)
    else
        throw(ArgumentError("MPLM106 can only be applied to production-destruction systems"))
    end
end

function initialize!(integrator,
                     cache::Union{MPLM22oopCache, MPLM33oopCache, MPLM43oopCache,
                                  MPLM54oopCache, MPLM75oopCache, MPLM106oopCache,
                                  MPLMMutableCache})
end

########################################################################################
### sigma approximations ###############################################################
########################################################################################
@inline function sigma_approx_1(uprev, P, d, dt, linsolve, small_constant)
    σ = add_small_constant(uprev, small_constant)
    return basic_patankar_step(uprev, P, σ, dt, linsolve, d)
end

@inline function sigma_approx_2(σ, uprevprev, P, d, dt, linsolve, small_constant)
    σ = add_small_constant(σ, small_constant)
    return basic_patankar_step(uprevprev, P, σ, 2 * dt, linsolve, d)
end

@inline function sigma_approx_3(σ, uprev, uprev3, P_tup, d_tup, dt, linsolve,
                                small_constant)
    T = eltype(σ)
    c1, c2, c3, c4, c5 = T(35) / 18, T(1) / 3, T(2) / 9, T(1) / 4, T(3) / 4

    P, P2, _, P4 = P_tup
    d, d2, _, d4 = d_tup
    σ = add_small_constant(σ, small_constant)
    Ptmp, dtmp = lincomb(c1, P, d, c2, P2, d2, c3, P4, d4)
    v = c4 * uprev + c5 * uprev3
    return basic_patankar_step(v, Ptmp, σ, dt, linsolve, dtmp)
end

@inline function sigma_approx_4(σ, uprev5, P_tup, d_tup, dt, linsolve, small_constant)
    T = eltype(σ)
    c1, c2, c3, c4 = T(225) / 96, T(50) / 96, T(200) / 96, T(5) / 96

    P, _, P3, P4, P5 = P_tup
    d, _, d3, d4, d5 = d_tup
    σ = add_small_constant(σ, small_constant)
    Ptmp, dtmp = lincomb(c1, P, d, c2, P3, d3, c3, P4, d4, c4, P5, d5)
    return basic_patankar_step(uprev5, Ptmp, σ, dt, linsolve, dtmp)
end

@inline function sigma_approx_5(σ, uprev7, P_tup, d_tup, dt, linsolve, small_constant)
    T = eltype(σ)
    c1, c2, c3, c4, c5, c6 = T(12) / 5, T(197) / 720, T(701) / 360, T(43) / 30,
                             T(107) / 360,
                             T(467) / 720

    P, _, P3, P4, P5, P6, P7 = P_tup
    d, _, d3, d4, d5, d6, d7 = d_tup
    σ = add_small_constant(σ, small_constant)
    Ptmp, dtmp = lincomb(c1, P, d, c2, P3, d3, c3, P4, d4, c4, P5,
                         d5, c5, P6, d6, c6, P7, d7)
    return basic_patankar_step(uprev7, Ptmp, σ, dt, linsolve, dtmp)
end

@inline function sigma_approx_1!(σ, uprev, P, d, dt, linsolve, small_constant)
    @.. broadcast=false σ=uprev + small_constant
    basic_patankar_step!(σ, uprev, P, d, linsolve.A, σ, dt, linsolve)
end

@inline function sigma_approx_2!(σ, uprev2, P, d, dt, linsolve, small_constant)
    @.. broadcast=false σ=σ + small_constant
    basic_patankar_step!(σ, uprev2, P, d, linsolve.A, σ, 2 * dt, linsolve)
end

@inline function sigma_approx_3!(σ, uprev, uprev3, P_tup, d_tup, d_tmp, dt, linsolve,
                                 small_constant)
    T = eltype(σ)
    c1, c2, c3, c4, c5 = T(35) / 18, T(1) / 3, T(2) / 9, T(1) / 4, T(3) / 4

    P, P2, _, P4 = P_tup
    d, d2, _, d4 = d_tup
    @.. broadcast=false σ=σ + small_constant
    lincomb!(linsolve.A, c1, P, c2, P2, c3, P4)
    lincomb!(d_tmp, c1, d, c2, d2, c3, d4)
    @.. broadcast=false linsolve.b=c4 * uprev + c5 * uprev3
    basic_patankar_step!(σ, linsolve.b, linsolve.A, d_tmp, linsolve.A, σ, dt, linsolve)
end

@inline function sigma_approx_4!(σ, uprev5, P_tup, d_tup, d_tmp, dt, linsolve,
                                 small_constant)
    T = eltype(σ)
    c1, c2, c3, c4 = T(225) / 96, T(50) / 96, T(200) / 96, T(5) / 96

    P, _, P3, P4, P5 = P_tup
    d, _, d3, d4, d5 = d_tup
    @.. broadcast=false σ=σ + small_constant
    lincomb!(linsolve.A, c1, P, c2, P3, c3, P4, c4, P5)
    lincomb!(d_tmp, c1, d, c2, d3, c3, d4, c4, d5)
    basic_patankar_step!(σ, uprev5, linsolve.A, d_tmp, linsolve.A, σ, dt, linsolve)
end

@inline function sigma_approx_5!(σ, uprev7, P_tup, d_tup, d_tmp, dt, linsolve,
                                 small_constant)
    T = eltype(σ)
    c1, c2, c3, c4, c5, c6 = T(12) / 5, T(197) / 720, T(701) / 360, T(43) / 30,
                             T(107) / 360,
                             T(467) / 720

    P, _, P3, P4, P5, P6, P7 = P_tup
    d, _, d3, d4, d5, d6, d7 = d_tup
    @.. broadcast=false σ=σ + small_constant
    lincomb!(linsolve.A, c1, P, c2, P3, c3, P4, c4, P5, c5,
             P6, c6, P7)
    lincomb!(d_tmp, c1, d, c2, d3, c3, d4, c4, d5, c5, d6,
             c6, d7)
    basic_patankar_step!(σ, uprev7, linsolve.A, d_tmp, linsolve.A, σ, dt, linsolve)
end

########################################################################################
### perform_step! ######################################################################
########################################################################################
#### MPLM22 ############################################################################
@muladd function start_MPLM22(P, d, t, dt, uprev, f, p, small_constant, linsolve;
                              substep_exp = 2)
    substeps = 2^substep_exp
    dt = dt / substeps

    u = uprev
    @inbounds for _ in 1:substeps
        u = perform_step_MPE(P, d, dt, u, small_constant, linsolve)
        t += dt

        P, d = evaluate_pds(f, u, p, t)
    end

    # 4 function evals and 4 solves
    nf = substeps
    ns = substeps

    return u, t, nf, ns
end

@muladd function start_MPLM22!(u, P, d, t, dt, uprev, σ, f, p, small_constant, linsolve;
                               substep_exp = 2)
    substeps = 2^substep_exp
    dt = dt / substeps

    u .= uprev
    @inbounds for _ in 1:substeps
        perform_step_MPE!(u, P, d, dt, u, σ, small_constant, linsolve)
        t += dt

        evaluate_pds!(P, d, f, u, p, t)
    end

    # 4 function evals and 4 solves
    nf = substeps
    ns = substeps

    return t, nf, ns
end

#TODO Use αβ in MPLM22
@muladd function perform_step_MPLM22(P, d, dt, uprev, uprevprev, linsolve,
                                     small_constant)

    # First σ approximation
    σ = sigma_approx_1(uprev, P, d, dt, linsolve, small_constant)

    # Main step 
    σ = add_small_constant(σ, small_constant)
    u = basic_patankar_step(uprevprev, P, σ, 2 * dt, linsolve, d)

    # statistics: 2 nsolve

    return u
end

@muladd function perform_step_MPLM22!(u, P, P2, d, d2, dt, uprev, uprevprev, σ, linsolve,
                                      small_constant)
    # First σ approximation
    sigma_approx_1!(σ, uprev, P, d, dt, linsolve, small_constant)

    # Main step 
    @.. broadcast=false σ=σ + small_constant
    basic_patankar_step!(u, uprevprev, P, d, linsolve.A, σ, 2 * dt, linsolve)

    # statistics: 2 nsolve

    return nothing
end

@muladd function perform_step!(integrator, cache::MPLM22oopCache, repeat_step = false)
    (; alg, t, dt, uprev, uprev2, f, p) = integrator
    (; uprevprev, small_constant) = cache

    # TODO: Check that this actually works
    if integrator.derivative_discontinuity
        cache.step = 1
    end

    if cache.step <= 1

        # increase step counter
        cache.step += 1

        P, d = evaluate_pds(f, uprev, p, t)
        integrator.stats.nf += 1

        # compute initial values
        u, _, nf, ns = start_MPLM22(P, d, t, dt, uprev, f, p, small_constant,
                                    alg.linsolve; substep_exp = alg.substep_level + 1)

        integrator.stats.nf += nf
        integrator.stats.nsolve += ns
    else
        # increase step counter
        cache.step += 1

        # evaluate production matrix
        P, d = evaluate_pds(f, uprev, p, t)
        integrator.stats.nf += 1

        u = perform_step_MPLM22(P, d, dt, uprev, uprevprev, alg.linsolve,
                                small_constant)
        integrator.stats.nsolve += 2
    end

    #TODO: Should be possible to use uprev2. But uprev2 is currently not updated.
    cache.uprevprev = uprev

    integrator.u = u
end

@muladd function perform_step!(integrator, cache::MPLM22Cache, repeat_step = false)
    (; alg, t, dt, u, uprev, uprev2, f, p) = integrator
    (; uprevprev, small_constant, P, P2, D, D2, σ, linsolve) = cache

    # TODO: Check that this actually works
    if integrator.derivative_discontinuity
        cache.step = 1
    end

    if cache.step <= 1

        # increase step counter
        cache.step += 1

        evaluate_pds!(P, D, f, uprev, p, t)
        integrator.stats.nf += 1

        # compute initial values
        t, nf, ns = start_MPLM22!(u, P, D, t, dt, uprev, σ, f, p, small_constant, linsolve;
                                  substep_exp = alg.substep_level + 1)

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
@muladd function start_MPLM33(P, d, t, dt, uprev, f, p, small_constant, linsolve;
                              substep_exp = 2)
    substeps = 2^substep_exp
    dts = dt / substeps

    ### first macro step ###############################################################
    # substep 1
    u, t, nf, ns = start_MPLM22(P, d, t, dts, uprev, f, p, small_constant, linsolve;
                                substep_exp)

    # substeps 2 - substeps                                    
    for _ in 1:(substeps - 1)
        uprevprev, uprev = uprev, u

        P, d = evaluate_pds(f, uprev, p, t)
        nf += 1

        u = perform_step_MPLM22(P, d, dts, uprev, uprevprev, linsolve, small_constant)
        t += dts
        ns += 2
    end

    v = u

    ### second macro step ############################################################
    for _ in 1:substeps
        uprevprev, uprev = uprev, u

        P, d = evaluate_pds(f, uprev, p, t)
        nf += 1

        u = perform_step_MPLM22(P, d, dts, uprev, uprevprev, linsolve, small_constant)
        t += dts
        ns += 2
    end

    return (v, u), t, nf, ns
end

@muladd function start_MPLM33!(v, tmp, P, P2, d, d2, t, dt, vprev, vprev2, σ, f, p,
                               small_constant, linsolve; substep_exp = 2)
    substeps = 2^substep_exp
    dts = dt / substeps

    ### first macro step ###############################################################
    # substep 1
    t, nf, ns = start_MPLM22!(v, P, d, t, dts, vprev, σ, f, p, small_constant, linsolve;
                              substep_exp)

    # substeps 2 - substeps                                    
    for _ in 1:(substeps - 1)
        shift!(v, vprev, vprev2)

        evaluate_pds!(P, d, f, vprev, p, t)
        nf += 1

        perform_step_MPLM22!(v, P, P2, d, d2, dts, vprev, vprev2, σ, linsolve,
                             small_constant)
        t += dts
        ns += 2
    end

    tmp .= v

    ### second macro step ############################################################
    for _ in 1:substeps
        shift!(v, vprev, vprev2)

        evaluate_pds!(P, d, f, vprev, p, t)
        nf += 1

        perform_step_MPLM22!(v, P, P2, d, d2, dts, vprev, vprev2, σ, linsolve,
                             small_constant)
        t += dts
        ns += 2
    end

    return t, nf, ns
end

@muladd function perform_step_MPLM33(P_tup, d_tup, dt, u_tup, linsolve, αβ,
                                     small_constant)
    P, P2, P3 = P_tup
    d, d2, d3 = d_tup
    uprev, uprevprev, uprev3 = u_tup
    α1, α2, α3, β1, β2, β3 = αβ

    # σ approximations 
    σ = sigma_approx_1(uprev, P, d, dt, linsolve, small_constant)
    σ = sigma_approx_2(σ, uprevprev, P, d, dt, linsolve, small_constant)

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

    # σ approximations
    sigma_approx_1!(σ, uprev, P, d, dt, linsolve, small_constant)
    sigma_approx_2!(σ, uprevprev, P, d, dt, linsolve, small_constant)

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

    # TODO: Check that this actually works
    if integrator.derivative_discontinuity
        cache.step = 1
    end

    if cache.step == 1
        # increase step count
        cache.step += 1

        # evaluate production matrix at tspan[1]
        P, d = evaluate_pds(f, uprev, p, t)
        integrator.stats.nf += 1

        # compute initial values 
        v, t, nf, ns = start_MPLM33(P, d, t, dt, uprev, f, p, small_constant,
                                    alg.linsolve; substep_exp = alg.substep_level + 1)
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

        cache.uprev3, cache.uprevprev = (uprevprev, uprev)
    else
        # increase step count
        cache.step += 1

        # evaluate production matrix
        P, d = evaluate_pds(f, uprev, p, t)
        integrator.stats.nf += 1

        P_tup = (P, P2, P3)
        d_tup = (d, d2, d3)
        u_tup = (uprev, uprevprev, uprev3)

        u = perform_step_MPLM33(P_tup, d_tup, dt, u_tup, alg.linsolve, αβ,
                                small_constant)
        integrator.stats.nsolve += 3

        cache.uprev3, cache.uprevprev = (uprevprev, uprev)
    end

    integrator.u = u

    cache.P3, cache.P2 = (P2, P)
    cache.d3, cache.d2 = (d2, d)
end

@muladd function perform_step!(integrator, cache::MPLM33Cache, repeat_step = false)
    (; alg, t, dt, uprev, uprev2, u, f, p) = integrator
    (; uprevprev, uprev3, v, vprev, vprev2, P, P2, P3, d, d2, d3, σ, αβ, small_constant, linsolve) = cache
    #TODO Check if number of v-vectors can be reduced. 
    # vprev2 and vprev3 are only used in the initialization phase.

    # TODO: Check that this actually works
    if integrator.derivative_discontinuity
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

        # save current P and d 
        P3 .= P
        !isnothing(d) && (d3 .= d)

        # compute initial values 
        # we use uprev3 as temporary storage for the value of u needed in step 1.
        _, nf, ns = start_MPLM33!(v, uprev3, P, P2, d, d2, t, dt, vprev, vprev2, σ, f, p,
                                  small_constant,
                                  linsolve; substep_exp = alg.substep_level + 1)
        integrator.stats.nf += nf
        integrator.stats.nsolve += ns

        # reset P and d
        P .= P3
        !isnothing(d) && (d .= d3)

        # u at time tspan[1] + dt
        u .= uprev3

        uprevprev .= uprev
    elseif cache.step == 2
        # increase step count
        cache.step += 1

        # evaluate production matrix at tspan[1] + dt
        evaluate_pds!(P, d, f, uprev, p, t)
        integrator.stats.nf += 1

        # u at time tspan[1] + 2*dt (this was computed in step 1)
        u .= v

        shift!(uprev, uprevprev, uprev3)
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

        shift!(uprev, uprevprev, uprev3)
    end

    shift!(P, P2, P3)
    shift!(d, d2, d3)
end

#### MPLM43 ############################################################################
@muladd function start_MPLM43(P, d, t, dt, uprev, f, p, small_constant, linsolve;
                              substep_exp = 2)
    # αβ33 = (0, 0, 1, 9 / 4, 0, 3 / 4)
    αβ33 = get_constant_parameters(MPLM33(), eltype(uprev))

    substeps = 2^substep_exp
    dts = dt / substeps

    ### first macro step ###############################################################
    # substeps 1 - 2
    v, t, nf, ns = start_MPLM33(P, d, t, dts, uprev, f, p, small_constant, linsolve;
                                substep_exp)

    uprevprev = uprev
    P2 = P
    d2 = d

    uprev = v[1]
    P, d = evaluate_pds(f, uprev, p, t - dts)
    nf += 1

    u = v[2]

    # substeps 3-substeps
    for _ in 1:(substeps - 2)
        uprev3, uprevprev, uprev = (uprevprev, uprev, u)
        P3, P2 = (P2, P)
        d3, d2 = (d2, d)

        P, d = evaluate_pds(f, uprev, p, t)
        nf += 1

        P_tup = (P, P2, P3)
        d_tup = (d, d2, d3)
        u_tup = (uprev, uprevprev, uprev3)

        u = perform_step_MPLM33(P_tup, d_tup, dts, u_tup, linsolve, αβ33,
                                small_constant)
        t += dts
        ns += 3
    end

    v1 = u

    ### second macro step ############################################################
    for _ in 1:substeps
        uprev3, uprevprev, uprev = (uprevprev, uprev, u)
        P3, P2 = (P2, P)
        d3, d2 = (d2, d)

        P, d = evaluate_pds(f, uprev, p, t)
        nf += 1

        P_tup = (P, P2, P3)
        d_tup = (d, d2, d3)
        u_tup = (uprev, uprevprev, uprev3)

        u = perform_step_MPLM33(P_tup, d_tup, dts, u_tup, linsolve, αβ33,
                                small_constant)
        t += dts
        ns += 3
    end

    v2 = u

    ### third macro step ############################################################
    for _ in 1:substeps
        uprev3, uprevprev, uprev = (uprevprev, uprev, u)
        P3, P2 = (P2, P)
        d3, d2 = (d2, d)

        P, d = evaluate_pds(f, uprev, p, t)
        nf += 1

        P_tup = (P, P2, P3)
        d_tup = (d, d2, d3)
        u_tup = (uprev, uprevprev, uprev3)

        u = perform_step_MPLM33(P_tup, d_tup, dts, u_tup, linsolve, αβ33,
                                small_constant)
        t += dts
        ns += 3
    end

    return (v1, v2, u), t, nf, ns
end

@muladd function start_MPLM43!(v, tmp, tmp2, P, P2, P3, d, d2, d3, t, dt, vprev, vprev2,
                               vprev3, σ, f, p,
                               small_constant, linsolve; substep_exp = 2)
    αβ33 = get_constant_parameters(MPLM33(), eltype(vprev))

    tmps = (tmp, tmp2)
    P_tup = (P, P2, P3)
    d_tup = (d, d2, d3)
    v_tup = (vprev, vprev2, vprev3)

    substeps = 2^substep_exp
    dts = dt / substeps

    # save current P and d
    P3 .= P
    !isnothing(d) && (d3 .= d)

    ### first macro step ###############################################################
    # substep 1 - 2
    t, nf, ns = start_MPLM33!(v, tmp, P, P2, d, d2, t, dts, vprev, vprev2, σ, f, p,
                              small_constant, linsolve; substep_exp)

    # vprev3 must be initialized as uprev.
    vprev2 .= vprev3 # == uprev
    P2 .= P3
    !isnothing(d3) && (d2 .= d3)

    vprev .= tmp
    evaluate_pds!(P, d, f, vprev, p, t - dts)
    nf += 1

    sub_steps = (substeps - 2, substeps, substeps)

    for (step_idx, n_iter) in enumerate(sub_steps)
        for _ in 1:n_iter
            shift!(v, v_tup...)
            shift!(P_tup...)
            shift!(d_tup...)

            evaluate_pds!(P, d, f, vprev, p, t)
            nf += 1

            perform_step_MPLM33!(v, P_tup, d_tup, dts, v_tup, σ, linsolve, αβ33,
                                 small_constant)
            t += dts
            ns += 3
        end

        # save initial data
        if step_idx < length(sub_steps) # last initial value is stored in v
            tmps[step_idx] .= v
        end
    end

    return t, nf, ns
end

@muladd function perform_step_MPLM43(P_tup, d_tup, dt, u_tup, linsolve, αβ,
                                     small_constant)
    P, P2, P3, P4 = P_tup
    d, d2, d3, d4 = d_tup
    uprev, uprevprev, uprev3, uprev4 = u_tup
    α1, α2, α3, α4, β1, β2, β3, β4 = αβ

    # σ approximations
    σ = sigma_approx_1(uprev, P, d, dt, linsolve, small_constant)
    σ = sigma_approx_2(σ, uprevprev, P, d, dt, linsolve, small_constant)

    # Main step 
    σ = add_small_constant(σ, small_constant)
    Ptmp, dtmp = lincomb(β1, P, d, β2, P2, d2, β3, P3, d3, β4, P4, d4)
    v = α1 * uprev + α2 * uprevprev + α3 * uprev3 + α4 * uprev4
    u = basic_patankar_step(v, Ptmp, σ, dt, linsolve, dtmp)

    # statistics: 3 nsolve

    return u
end

@muladd function perform_step_MPLM43!(u, P_tup, d_tup, dt, u_tup, σ, linsolve, αβ,
                                      small_constant)
    P, P2, P3, P4 = P_tup
    d, d2, d3, d4 = d_tup
    uprev, uprevprev, uprev3, uprev4 = u_tup
    α1, α2, α3, α4, β1, β2, β3, β4 = αβ

    # σ approximations
    sigma_approx_1!(σ, uprev, P, d, dt, linsolve, small_constant)
    sigma_approx_2!(σ, uprevprev, P, d, dt, linsolve, small_constant)

    # Main step 
    @.. broadcast=false σ=σ + small_constant
    lincomb!(P4, β1, P, β2, P2, β3, P3, β4, P4)
    lincomb!(d4, β1, d, β2, d2, β3, d3, β4, d4)
    @.. broadcast=false linsolve.b=α1 * uprev + α2 * uprevprev + α3 * uprev3 + α4 * uprev4
    basic_patankar_step!(u, linsolve.b, P4, d4, linsolve.A, σ, dt, linsolve)

    # statistics: 3 nsolve

    return nothing
end

@muladd function perform_step!(integrator, cache::MPLM43oopCache, repeat_step = false)
    (; alg, t, dt, uprev, uprev2, f, p) = integrator
    (; uprevprev, uprev3, uprev4, P2, P3, P4, d2, d3, d4, αβ, small_constant) = cache

    # TODO: Check that this actually works
    if integrator.derivative_discontinuity
        cache.step = 1
    end

    if cache.step == 1
        # increase step count
        cache.step += 1

        # evaluate production matrix at tspan[1]
        P, d = evaluate_pds(f, uprev, p, t)
        integrator.stats.nf += 1

        # compute initial values 
        v, _, nf, ns = start_MPLM43(P, d, t, dt, uprev, f, p, small_constant,
                                    alg.linsolve; substep_exp = alg.substep_level + 1)
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

        cache.uprev3, cache.uprevprev = (uprevprev, uprev)
    elseif cache.step == 3
        # increase step count
        cache.step += 1

        # evaluate production matrix at tspan[1] + 2*dt
        P, d = evaluate_pds(f, uprev, p, t)
        integrator.stats.nf += 1

        # u at time tspan[1] + 3*dt (this was computed in step 1)
        u = cache.uprev4

        cache.uprev4, cache.uprev3, cache.uprevprev = (uprev3, uprevprev, uprev)
    else
        # evaluate production matrix
        P, d = evaluate_pds(f, uprev, p, t)
        integrator.stats.nf += 1

        P_tup = (P, P2, P3, P4)
        d_tup = (d, d2, d3, d4)
        u_tup = (uprev, uprevprev, uprev3, uprev4)

        u = perform_step_MPLM43(P_tup, d_tup, dt, u_tup, alg.linsolve, αβ,
                                small_constant)
        integrator.stats.nsolve += 3

        cache.uprev4, cache.uprev3, cache.uprevprev = (uprev3, uprevprev, uprev)
    end

    integrator.u = u

    cache.P4, cache.P3, cache.P2 = (P3, P2, P)
    cache.d4, cache.d3, cache.d2 = (d3, d2, d)
end

@muladd function perform_step!(integrator, cache::MPLM43Cache, repeat_step = false)
    (; alg, t, dt, uprev, uprev2, u, f, p) = integrator
    (; uprevprev, uprev3, uprev4, v, vprev, vprev2, vprev3, P, P2, P3, P4, d, d2, d3, d4, σ, αβ, small_constant, linsolve) = cache
    #TODO Check if number of v-vectors can be reduced. 
    # vprev2 and vprev3 are only used in the initialization phase.

    # TODO: Check that this actually works
    if integrator.derivative_discontinuity
        cache.step = 1
    end

    if cache.step == 1
        # increase step count
        cache.step += 1

        # initilialze v vectors 
        vprev .= uprev
        vprev2 .= uprev
        vprev3 .= uprev

        # evaluate production matrix at tspan[1]
        evaluate_pds!(P, d, f, uprev, p, t)
        integrator.stats.nf += 1

        # save current P and d
        P4 .= P
        !isnothing(d) && (d4 .= d)

        # compute initial values 
        # we use uprevprev as temporary storage for the value of u needed in step 1.
        # we use uprev3 as temporary storage for the value of u needed in step 2.
        # we use v as temporary storage for the value of u needed in step 3.
        _, nf, ns = start_MPLM43!(v, uprevprev, uprev3, P, P2, P3, d, d2, d3, t, dt, vprev,
                                  vprev2, vprev3, σ, f, p,
                                  small_constant,
                                  linsolve; substep_exp = alg.substep_level + 1)
        integrator.stats.nf += nf
        integrator.stats.nsolve += ns

        # reset P
        P .= P4
        !isnothing(d4) && (d .= d4)

        # u at time tspan[1] + dt
        u .= uprevprev

        uprevprev .= uprev
    elseif cache.step == 2
        # increase step count
        cache.step += 1

        # evaluate production matrix at tspan[1] + dt
        evaluate_pds!(P, d, f, uprev, p, t)
        integrator.stats.nf += 1

        # u at time tspan[1] + 2*dt (this was computed in step 1)
        u .= uprev3

        shift!(uprev, uprevprev, uprev3)
    elseif cache.step == 3
        # increase step count
        cache.step += 1

        # evaluate production matrix at tspan[1] + 2*dt
        evaluate_pds!(P, d, f, uprev, p, t)
        integrator.stats.nf += 1

        # u at time tspan[1] + 3*dt (this was computed in step 1)
        u .= v

        shift!(uprev, uprevprev, uprev3, uprev4)
    else
        # increase step count
        cache.step += 1

        # evaluate production matrix
        evaluate_pds!(P, d, f, uprev, p, t)
        integrator.stats.nf += 1

        P_tup = (P, P2, P3, P4)
        d_tup = (d, d2, d3, d4)
        u_tup = (uprev, uprevprev, uprev3, uprev4)

        perform_step_MPLM43!(u, P_tup, d_tup, dt, u_tup, σ, linsolve, αβ,
                             small_constant)
        integrator.stats.nsolve += 3

        shift!(uprev, uprevprev, uprev3, uprev4)
    end

    shift!(P, P2, P3, P4)
    shift!(d, d2, d3, d4)
end

#### MPLM54 ############################################################################
#TODO Check nf and ns everywhere!
@muladd function start_MPLM54(P, d, t, dt, uprev, f, p, small_constant, linsolve;
                              substep_exp = 2)
    αβ43 = get_constant_parameters(MPLM43(), eltype(uprev))

    substeps = 2^substep_exp
    dts = dt / substeps

    ### first macro step ###############################################################
    # substep 1 - 3
    v, t, nf, ns = start_MPLM43(P, d, t, dts, uprev, f, p, small_constant, linsolve;
                                substep_exp)

    uprev3 = uprev
    P3 = P
    d3 = d

    uprevprev = v[1]
    P2, d2 = evaluate_pds(f, uprevprev, p, t - 2 * dts)
    nf += 1

    uprev = v[2]
    P, d = evaluate_pds(f, uprev, p, t - dts)
    nf += 1

    u = v[3]

    # substep 4 - substeps
    for _ in 1:(substeps - 3)
        uprev4, uprev3, uprevprev, uprev = (uprev3, uprevprev, uprev, u)
        P4, P3, P2 = (P3, P2, P)
        d4, d3, d2 = (d3, d2, d)

        P, d = evaluate_pds(f, uprev, p, t)
        nf += 1

        P_tup = (P, P2, P3, P4)
        d_tup = (d, d2, d3, d4)
        u_tup = (uprev, uprevprev, uprev3, uprev4)

        u = perform_step_MPLM43(P_tup, d_tup, dts, u_tup, linsolve, αβ43,
                                small_constant)
        t += dts
        ns += 3
    end

    v1 = u

    ### second macro step ############################################################
    for _ in 1:substeps
        uprev4, uprev3, uprevprev, uprev = (uprev3, uprevprev, uprev, u)
        P4, P3, P2 = (P3, P2, P)
        d4, d3, d2 = (d3, d2, d)

        P, d = evaluate_pds(f, uprev, p, t)
        nf += 1

        P_tup = (P, P2, P3, P4)
        d_tup = (d, d2, d3, d4)
        u_tup = (uprev, uprevprev, uprev3, uprev4)

        u = perform_step_MPLM43(P_tup, d_tup, dts, u_tup, linsolve, αβ43,
                                small_constant)
        t += dts
        ns += 3
    end

    v2 = u

    ### third macro step ############################################################
    for _ in 1:substeps
        uprev4, uprev3, uprevprev, uprev = (uprev3, uprevprev, uprev, u)
        P4, P3, P2 = (P3, P2, P)
        d4, d3, d2 = (d3, d2, d)

        P, d = evaluate_pds(f, uprev, p, t)
        nf += 1

        P_tup = (P, P2, P3, P4)
        d_tup = (d, d2, d3, d4)
        u_tup = (uprev, uprevprev, uprev3, uprev4)

        u = perform_step_MPLM43(P_tup, d_tup, dts, u_tup, linsolve, αβ43,
                                small_constant)
        t += dts
        ns += 3
    end

    v3 = u

    # fourth macro step
    for _ in 1:substeps
        uprev4, uprev3, uprevprev, uprev = (uprev3, uprevprev, uprev, u)
        P4, P3, P2 = (P3, P2, P)
        d4, d3, d2 = (d3, d2, d)

        P, d = evaluate_pds(f, uprev, p, t)
        nf += 1

        P_tup = (P, P2, P3, P4)
        d_tup = (d, d2, d3, d4)
        u_tup = (uprev, uprevprev, uprev3, uprev4)

        u = perform_step_MPLM43(P_tup, d_tup, dts, u_tup, linsolve, αβ43,
                                small_constant)
        t += dts
        ns += 3
    end

    return (v1, v2, v3, u), t, nf, ns
end

@muladd function start_MPLM54!(v, vprev, vprev2, vprev3, vprev4,
                               tmp, tmp2, tmp3, P, P2, P3, P4,
                               d, d2, d3, d4, t, dt, σ, f, p,
                               small_constant, linsolve; substep_exp = 2)
    αβ43 = get_constant_parameters(MPLM43(), eltype(vprev))

    tmps = (tmp, tmp2, tmp3)
    P_tup = (P, P2, P3, P4)
    d_tup = (d, d2, d3, d4)
    v_tup = (vprev, vprev2, vprev3, vprev4)

    substeps = 2^substep_exp
    dts = dt / substeps

    # save current P and d
    P4 .= P
    !isnothing(d) && (d4 .= d)

    ### first macro step ###############################################################
    # substep 1 - 3
    t, nf, ns = start_MPLM43!(v, tmp, tmp2, P, P2, P3, d, d2, d3, t, dts, vprev, vprev2,
                              vprev3, σ, f, p,
                              small_constant, linsolve; substep_exp)

    # vprev4 must be initialized as uprev.
    vprev3 .= vprev4 # == uprev
    P3 .= P4 # == P(uprev)
    !isnothing(d4) && (d3 .= d4)

    vprev2 .= tmp
    evaluate_pds!(P2, d2, f, vprev2, p, t - 2 * dts)
    nf += 1

    vprev .= tmp2
    evaluate_pds!(P, d, f, vprev, p, t - dts)
    nf += 1

    # we have four substeps per macro step, except for the first macro step
    sub_steps = (substeps - 3, substeps, substeps, substeps)

    for (step_idx, n_iter) in enumerate(sub_steps)
        for _ in 1:n_iter
            shift!(v, v_tup...)
            shift!(P_tup...)
            shift!(d_tup...)

            evaluate_pds!(P, d, f, vprev, p, t)
            nf += 1

            perform_step_MPLM43!(v, P_tup, d_tup, dts, v_tup, σ, linsolve, αβ43,
                                 small_constant)
            t += dts
            ns += 3
        end

        # save initial data 
        if step_idx < length(sub_steps) # last initial value is stored in v
            tmps[step_idx] .= v
        end
    end

    return t, nf, ns
end

@muladd function perform_step_MPLM54(P_tup, d_tup, dt, u_tup, linsolve, αβ,
                                     small_constant)
    P, P2, P3, P4, P5 = P_tup
    d, d2, d3, d4, d5 = d_tup
    uprev, uprevprev, uprev3, uprev4, uprev5 = u_tup
    α1, α2, α3, α4, α5, β1, β2, β3, β4, β5 = αβ

    # σ approximations
    σ = sigma_approx_1(uprev, P, d, dt, linsolve, small_constant)
    σ = sigma_approx_2(σ, uprevprev, P, d, dt, linsolve, small_constant)
    σ = sigma_approx_3(σ, uprev, uprev3, P_tup, d_tup, dt, linsolve, small_constant)

    # Main step 
    σ = add_small_constant(σ, small_constant)
    Ptmp, dtmp = lincomb(β1, P, d, β2, P2, d2, β3, P3, d3, β4, P4, d4, β5, P5, d5)
    v = α1 * uprev + α2 * uprevprev + α3 * uprev3 + α4 * uprev4 + α5 * uprev5
    u = basic_patankar_step(v, Ptmp, σ, dt, linsolve, dtmp)

    # statistics: 4 nsolve

    return u
end

@muladd function perform_step_MPLM54!(u, P_tup, d_tup, d_tmp, dt, u_tup, σ, linsolve, αβ,
                                      small_constant)
    P, P2, P3, P4, P5 = P_tup
    d, d2, d3, d4, d5 = d_tup
    uprev, uprevprev, uprev3, uprev4, uprev5 = u_tup
    α1, α2, α3, α4, α5, β1, β2, β3, β4, β5 = αβ

    # σ approximations
    sigma_approx_1!(σ, uprev, P, d, dt, linsolve, small_constant)
    sigma_approx_2!(σ, uprevprev, P, d, dt, linsolve, small_constant)
    sigma_approx_3!(σ, uprev, uprev3, P_tup, d_tup, d_tmp, dt, linsolve, small_constant)

    # Main step 
    @.. broadcast=false σ=σ + small_constant
    lincomb!(P5, β1, P, β2, P2, β3, P3, β4, P4, β5, P5)
    lincomb!(d5, β1, d, β2, d2, β3, d3, β4, d4, β5, d5)
    @.. broadcast=false linsolve.b=α1 * uprev + α2 * uprevprev + α3 * uprev3 + α4 * uprev4 +
                                   α5 * uprev5
    basic_patankar_step!(u, linsolve.b, P5, d5, linsolve.A, σ, dt, linsolve)

    # statistics: 3 nsolve

    return nothing
end

@muladd function perform_step!(integrator, cache::MPLM54oopCache, repeat_step = false)
    (; alg, t, dt, uprev, uprev2, f, p) = integrator
    (; uprevprev, uprev3, uprev4, uprev5, P2, P3, P4, P5, d2, d3, d4, d5, αβ, small_constant) = cache

    # TODO: Check that this actually works
    if integrator.derivative_discontinuity
        cache.step = 1
    end

    if cache.step == 1
        # increase step count
        cache.step += 1

        # evaluate production matrix at tspan[1]
        P, d = evaluate_pds(f, uprev, p, t)
        integrator.stats.nf += 1

        # compute initial values 
        v, _, nf, ns = start_MPLM54(P, d, t, dt, uprev, f, p, small_constant,
                                    alg.linsolve; substep_exp = alg.substep_level + 1)
        integrator.stats.nf += nf
        integrator.stats.nsolve += ns

        # u at time tspan[1] + dt
        u = v[1]
        #u = f.analytic(u,p,t+dt)

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
        #u = f.analytic(u,p,t+dt)

        cache.uprev3, cache.uprevprev = (uprevprev, uprev)
    elseif cache.step == 3
        # increase step count
        cache.step += 1

        # evaluate production matrix at tspan[1] + 3*dt
        P, d = evaluate_pds(f, uprev, p, t)
        integrator.stats.nf += 1

        # u at time tspan[1] + 3*dt (this was computed in step 1)
        u = cache.uprev4
        #u = f.analytic(u,p,t+dt)

        cache.uprev4, cache.uprev3, cache.uprevprev = (uprev3, uprevprev, uprev)
    elseif cache.step == 4
        # increase step count
        cache.step += 1

        # evaluate production matrix at tspan[1] + 4*dt
        P, d = evaluate_pds(f, uprev, p, t)
        integrator.stats.nf += 1

        # u at time tspan[1] + 4*dt (this was computed in step 1)
        u = cache.uprev5
        #u = f.analytic(u,p,t+dt)

        cache.uprev5, cache.uprev4, cache.uprev3, cache.uprevprev = (uprev4, uprev3,
                                                                     uprevprev, uprev)
    else
        # increase step count
        cache.step += 1

        # evaluate production matrix
        P, d = evaluate_pds(f, uprev, p, t)
        integrator.stats.nf += 1

        P_tup = (P, P2, P3, P4, P5)
        d_tup = (d, d2, d3, d4, d5)
        u_tup = (uprev, uprevprev, uprev3, uprev4, uprev5)

        u = perform_step_MPLM54(P_tup, d_tup, dt, u_tup, alg.linsolve, αβ,
                                small_constant)
        integrator.stats.nsolve += 4

        cache.uprev5, cache.uprev4, cache.uprev3, cache.uprevprev = (uprev4, uprev3,
                                                                     uprevprev, uprev)
    end

    integrator.u = u

    cache.P5, cache.P4, cache.P3, cache.P2 = (P4, P3, P2, P)
    cache.d5, cache.d4, cache.d3, cache.d2 = (d4, d3, d2, d)
end

@muladd function perform_step!(integrator, cache::MPLM54Cache, repeat_step = false)
    (; alg, t, dt, uprev, uprev2, u, f, p) = integrator
    (; uprevprev, uprev3, uprev4, uprev5, v, vprev, vprev2, vprev3, vprev4, P, P2, P3, P4, P5, d, d2, d3, d4, d5, d_tmp, σ, αβ, small_constant, linsolve) = cache
    #TODO Check if number of v-vectors can be reduced. 
    # vprev2, vprev3, vprev4 are only used in the initialization phase.

    # TODO: Check that this actually works
    if integrator.derivative_discontinuity
        cache.step = 1
    end

    if cache.step == 1
        # increase step count
        cache.step += 1

        # initilialze v vectors 
        vprev .= uprev
        vprev2 .= uprev
        vprev3 .= uprev
        vprev4 .= uprev

        # evaluate production matrix at tspan[1]
        evaluate_pds!(P, d, f, uprev, p, t)
        integrator.stats.nf += 1

        # save current P and d
        P5 .= P
        !isnothing(d) && (d5 .= d)

        # compute initial values 
        # we use uprevprev as temporary storage for the value of u needed in step 1.
        # we use uprev3 as temporary storage for the value of u needed in step 2.
        # we use uprev4 as temporary storage for the value of u needed in step 3.
        # we use v as temporary storage for the value of u needed in step 4.
        _, nf, ns = start_MPLM54!(v, vprev, vprev2, vprev3, vprev4, uprevprev, uprev3,
                                  uprev4,
                                  P, P2, P3, P4, d, d2, d3, d4, t,
                                  dt, σ, f, p,
                                  small_constant,
                                  linsolve; substep_exp = alg.substep_level + 1)
        integrator.stats.nf += nf
        integrator.stats.nsolve += ns

        # reset P
        P .= P5
        !isnothing(d5) && (d .= d5)

        # u at time tspan[1] + dt
        u .= uprevprev

        uprevprev .= uprev
    elseif cache.step == 2
        # increase step count
        cache.step += 1

        # evaluate production matrix at tspan[1] + dt
        evaluate_pds!(P, d, f, uprev, p, t)
        integrator.stats.nf += 1

        # u at time tspan[1] + 2*dt (this was computed in step 1)
        u .= uprev3

        shift!(uprev, uprevprev, uprev3)
    elseif cache.step == 3
        # increase step count
        cache.step += 1

        # evaluate production matrix at tspan[1] + 2*dt
        evaluate_pds!(P, d, f, uprev, p, t)
        integrator.stats.nf += 1

        # u at time tspan[1] + 3*dt (this was computed in step 1)
        u .= uprev4

        shift!(uprev, uprevprev, uprev3, uprev4)
    elseif cache.step == 4
        # increase step count
        cache.step += 1

        # evaluate production matrix at tspan[1] + 2*dt
        evaluate_pds!(P, d, f, uprev, p, t)
        integrator.stats.nf += 1

        # u at time tspan[1] + 4*dt (this was computed in step 1)
        u .= v

        shift!(uprev, uprevprev, uprev3, uprev4, uprev5)
    else
        # increase step count
        cache.step += 1

        # evaluate production matrix
        evaluate_pds!(P, d, f, uprev, p, t)
        integrator.stats.nf += 1

        P_tup = (P, P2, P3, P4, P5)
        d_tup = (d, d2, d3, d4, d5)
        u_tup = (uprev, uprevprev, uprev3, uprev4, uprev5)

        perform_step_MPLM54!(u, P_tup, d_tup, d_tmp, dt, u_tup, σ, linsolve, αβ,
                             small_constant)
        integrator.stats.nsolve += 4

        shift!(uprev, uprevprev, uprev3, uprev4, uprev5)
    end

    shift!(P, P2, P3, P4, P5)
    shift!(d, d2, d3, d4, d5)
end

#### MPLM75 ############################################################################
@muladd function start_MPLM75(P, d, t, dt, uprev, f, p, small_constant, linsolve;
                              substep_exp = 2)
    αβ54 = get_constant_parameters(MPLM54(), eltype(uprev))

    substeps = 2^substep_exp
    dts = dt / substeps

    ### first macro step ###############################################################
    # substeps 1 - 4
    v, t, nf, ns = start_MPLM54(P, d, t, dts, uprev, f, p, small_constant, linsolve;
                                substep_exp)

    uprev4 = uprev
    P4 = P
    d4 = d

    uprev3 = v[1]
    P3, d3 = evaluate_pds(f, uprev3, p, t - 3 * dts)
    nf += 1

    uprevprev = v[2]
    P2, d2 = evaluate_pds(f, uprevprev, p, t - 2 * dts)
    nf += 1

    uprev = v[3]
    P, d = evaluate_pds(f, uprev, p, t - dts)
    nf += 1

    u = v[4]

    for _ in 1:(substeps - 4)
        uprev5, uprev4, uprev3, uprevprev, uprev = (uprev4, uprev3, uprevprev, uprev, u)
        P5, P4, P3, P2 = (P4, P3, P2, P)
        d5, d4, d3, d2 = (d4, d3, d2, d)

        P, d = evaluate_pds(f, uprev, p, t)
        nf += 1

        P_tup = (P, P2, P3, P4, P5)
        d_tup = (d, d2, d3, d4, d5)
        u_tup = (uprev, uprevprev, uprev3, uprev4, uprev5)

        u = perform_step_MPLM54(P_tup, d_tup, dts, u_tup, linsolve, αβ54,
                                small_constant)
        t += dts
        ns += 4
    end

    v1 = u

    ### second macro step ############################################################
    for _ in 1:substeps
        uprev5, uprev4, uprev3, uprevprev, uprev = (uprev4, uprev3, uprevprev, uprev, u)
        P5, P4, P3, P2 = (P4, P3, P2, P)
        d5, d4, d3, d2 = (d4, d3, d2, d)

        P, d = evaluate_pds(f, uprev, p, t)
        nf += 1

        P_tup = (P, P2, P3, P4, P5)
        d_tup = (d, d2, d3, d4, d5)
        u_tup = (uprev, uprevprev, uprev3, uprev4, uprev5)

        u = perform_step_MPLM54(P_tup, d_tup, dts, u_tup, linsolve, αβ54,
                                small_constant)
        t += dts
        ns += 4
    end

    v2 = u

    ### third macro step ############################################################
    for _ in 1:substeps
        uprev5, uprev4, uprev3, uprevprev, uprev = (uprev4, uprev3, uprevprev, uprev, u)
        P5, P4, P3, P2 = (P4, P3, P2, P)
        d5, d4, d3, d2 = (d4, d3, d2, d)

        P, d = evaluate_pds(f, uprev, p, t)
        nf += 1

        P_tup = (P, P2, P3, P4, P5)
        d_tup = (d, d2, d3, d4, d5)
        u_tup = (uprev, uprevprev, uprev3, uprev4, uprev5)

        u = perform_step_MPLM54(P_tup, d_tup, dts, u_tup, linsolve, αβ54,
                                small_constant)
        t += dts
        ns += 4
    end

    v3 = u

    ### fourth macro step ############################################################
    for _ in 1:substeps
        uprev5, uprev4, uprev3, uprevprev, uprev = (uprev4, uprev3, uprevprev, uprev, u)
        P5, P4, P3, P2 = (P4, P3, P2, P)
        d5, d4, d3, d2 = (d4, d3, d2, d)

        P, d = evaluate_pds(f, uprev, p, t)
        nf += 1

        P_tup = (P, P2, P3, P4, P5)
        d_tup = (d, d2, d3, d4, d5)
        u_tup = (uprev, uprevprev, uprev3, uprev4, uprev5)

        u = perform_step_MPLM54(P_tup, d_tup, dts, u_tup, linsolve, αβ54,
                                small_constant)
        t += dts
        ns += 4
    end

    v4 = u

    ### fifth macro step ############################################################
    for _ in 1:substeps
        uprev5, uprev4, uprev3, uprevprev, uprev = (uprev4, uprev3, uprevprev, uprev, u)
        P5, P4, P3, P2 = (P4, P3, P2, P)
        d5, d4, d3, d2 = (d4, d3, d2, d)

        P, d = evaluate_pds(f, uprev, p, t)
        nf += 1

        P_tup = (P, P2, P3, P4, P5)
        d_tup = (d, d2, d3, d4, d5)
        u_tup = (uprev, uprevprev, uprev3, uprev4, uprev5)

        u = perform_step_MPLM54(P_tup, d_tup, dts, u_tup, linsolve, αβ54,
                                small_constant)
        t += dts
        ns += 4
    end

    v5 = u

    ### sixth macro step ############################################################
    for _ in 1:substeps
        uprev5, uprev4, uprev3, uprevprev, uprev = (uprev4, uprev3, uprevprev, uprev, u)
        P5, P4, P3, P2 = (P4, P3, P2, P)
        d5, d4, d3, d2 = (d4, d3, d2, d)

        P, d = evaluate_pds(f, uprev, p, t)
        nf += 1

        P_tup = (P, P2, P3, P4, P5)
        d_tup = (d, d2, d3, d4, d5)
        u_tup = (uprev, uprevprev, uprev3, uprev4, uprev5)

        u = perform_step_MPLM54(P_tup, d_tup, dts, u_tup, linsolve, αβ54,
                                small_constant)
        t += dts
        ns += 4
    end

    return (v1, v2, v3, v4, v5, u), t, nf, ns
end

@muladd function start_MPLM75!(v, vprev, vprev2, vprev3, vprev4, vprev5,
                               tmp, tmp2, tmp3, tmp4, tmp5,
                               P, P2, P3, P4, P5, d, d2, d3, d4, d5, d_tmp,
                               t, dt, σ, f, p, small_constant, linsolve; substep_exp = 2)
    αβ54 = get_constant_parameters(MPLM54(), eltype(vprev))

    tmps = (tmp, tmp2, tmp3, tmp4, tmp5)
    P_tup = (P, P2, P3, P4, P5)
    d_tup = (d, d2, d3, d4, d5)
    v_tup = (vprev, vprev2, vprev3, vprev4, vprev5)

    substeps = 2^substep_exp
    dts = dt / substeps

    # save current P and d
    P5 .= P
    !isnothing(d) && (d5 .= d)

    ### first macro step ###############################################################
    # substep 1 - 4
    t, nf, ns = start_MPLM54!(v, vprev, vprev2, vprev3, vprev4, tmp, tmp2, tmp3, P, P2, P3,
                              P4,
                              d, d2, d3, d4, t, dts, σ, f, p,
                              small_constant, linsolve; substep_exp)

    # vprev5 must be initialized as uprev.
    vprev4 .= vprev5 # == uprev
    P4 .= P5 # == P(uprev)
    !isnothing(d5) && (d4 .= d5)

    vprev3 .= tmp
    evaluate_pds!(P3, d3, f, vprev3, p, t - 3 * dts)
    nf += 1

    vprev2 .= tmp2
    evaluate_pds!(P2, d2, f, vprev2, p, t - 2 * dts)
    nf += 1

    vprev .= tmp3
    evaluate_pds!(P, d, f, vprev, p, t - dts)
    nf += 1

    tmp .= v

    sub_steps = (substeps - 4, substeps, substeps, substeps, substeps, substeps)

    for (step_idx, n_iter) in enumerate(sub_steps)
        # step_idx 1 corresponds to "second macro step"
        # we need to use tmps[step_idx + 1] because tmp was already handled 

        for _ in 1:n_iter
            shift!(v, v_tup...)
            shift!(P_tup...)
            shift!(d_tup...)

            evaluate_pds!(P, d, f, vprev, p, t)
            nf += 1

            perform_step_MPLM54!(v, P_tup, d_tup, d_tmp, dts, v_tup, σ, linsolve, αβ54,
                                 small_constant)
            t += dts
            ns += 4
        end

        # save initial data
        if step_idx < length(sub_steps) # last initial value is stored in v
            tmps[step_idx] .= v
        end
    end

    return t, nf, ns
end

@muladd function perform_step_MPLM75(P_tup, d_tup, dt, u_tup, linsolve, αβ,
                                     small_constant)
    P, P2, P3, P4, P5, P6, P7 = P_tup
    d, d2, d3, d4, d5, d6, d7 = d_tup
    uprev, uprevprev, uprev3, uprev4, uprev5, uprev6, uprev7 = u_tup
    α1, α2, α3, α4, α5, α6, α7, β1, β2, β3, β4, β5, β6, β7 = αβ

    # σ approximations
    σ = sigma_approx_1(uprev, P, d, dt, linsolve, small_constant)
    σ = sigma_approx_2(σ, uprevprev, P, d, dt, linsolve, small_constant)
    σ = sigma_approx_3(σ, uprev, uprev3, P_tup, d_tup, dt, linsolve, small_constant)
    σ = sigma_approx_4(σ, uprev5, P_tup, d_tup, dt, linsolve, small_constant)

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

@muladd function perform_step_MPLM75!(u, P_tup, d_tup, d_tmp, dt, u_tup, σ, linsolve, αβ,
                                      small_constant)
    P, P2, P3, P4, P5, P6, P7 = P_tup
    d, d2, d3, d4, d5, d6, d7 = d_tup
    uprev, uprevprev, uprev3, uprev4, uprev5, uprev6, uprev7 = u_tup
    α1, α2, α3, α4, α5, α6, α7, β1, β2, β3, β4, β5, β6, β7 = αβ

    # σ approximations
    sigma_approx_1!(σ, uprev, P, d, dt, linsolve, small_constant)
    sigma_approx_2!(σ, uprevprev, P, d, dt, linsolve, small_constant)
    sigma_approx_3!(σ, uprev, uprev3, P_tup, d_tup, d_tmp, dt, linsolve, small_constant)
    sigma_approx_4!(σ, uprev5, P_tup, d_tup, d_tmp, dt, linsolve, small_constant)

    # Main step 
    @.. broadcast=false σ=σ + small_constant
    lincomb!(P7, β1, P, β2, P2, β3, P3, β4, P4, β5, P5, β6, P6, β7, P7)
    lincomb!(d7, β1, d, β2, d2, β3, d3, β4, d4, β5, d5, β6, d6, β7, d7)
    @.. broadcast=false linsolve.b=α1 * uprev + α2 * uprevprev + α3 * uprev3 + α4 * uprev4 +
                                   α5 * uprev5 + α6 * uprev6 + α7 * uprev7
    basic_patankar_step!(u, linsolve.b, P7, d7, linsolve.A, σ, dt, linsolve)

    # statistics: 5 nsolves

    return nothing
end

@muladd function perform_step!(integrator, cache::MPLM75oopCache, repeat_step = false)
    (; alg, t, dt, uprev, uprev2, f, p) = integrator
    (; uprevprev, uprev3, uprev4, uprev5, uprev6, uprev7, P2, P3, P4, P5, P6, P7, d2, d3, d4, d5, d6, d7, αβ, small_constant) = cache

    # TODO: Check that this actually works
    if integrator.derivative_discontinuity
        cache.step = 1
    end

    if cache.step == 1
        # increase step count
        cache.step += 1

        # evaluate production matrix at tspan[1]
        P, d = evaluate_pds(f, uprev, p, t)
        integrator.stats.nf += 1

        # compute initial values 
        v, _, nf, ns = start_MPLM75(P, d, t, dt, uprev, f, p, small_constant,
                                    alg.linsolve; substep_exp = alg.substep_level + 1)
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

        cache.uprev3, cache.uprevprev = (uprevprev, uprev)
    elseif cache.step == 3
        # increase step count
        cache.step += 1

        # evaluate production matrix at tspan[1] + 2*dt
        P, d = evaluate_pds(f, uprev, p, t)
        integrator.stats.nf += 1

        # u at time tspan[1] + 3*dt (this was computed in step 1)
        u = cache.uprev4

        cache.uprev4, cache.uprev3, cache.uprevprev = (uprev3, uprevprev, uprev)
    elseif cache.step == 4
        # increase step count
        cache.step += 1

        # evaluate production matrix at tspan[1] + 3*dt
        P, d = evaluate_pds(f, uprev, p, t)
        integrator.stats.nf += 1

        # u at time tspan[1] + 4*dt (this was computed in step 1)
        u = cache.uprev5

        cache.uprev5, cache.uprev4, cache.uprev3, cache.uprevprev = (uprev4, uprev3,
                                                                     uprevprev, uprev)
    elseif cache.step == 5
        # increase step count
        cache.step += 1

        # evaluate production matrix at tspan[1] + 4*dt
        P, d = evaluate_pds(f, uprev, p, t)
        integrator.stats.nf += 1

        # u at time tspan[1] + 5*dt (this was computed in step 1)
        u = cache.uprev6

        cache.uprev6, cache.uprev5, cache.uprev4, cache.uprev3, cache.uprevprev = (uprev5,
                                                                                   uprev4,
                                                                                   uprev3,
                                                                                   uprevprev,
                                                                                   uprev)
    elseif cache.step == 6
        # increase step count
        cache.step += 1

        # evaluate production matrix at tspan[1] + 5*dt
        P, d = evaluate_pds(f, uprev, p, t)
        integrator.stats.nf += 1

        # u at time tspan[1] + 6*dt (this was computed in step 1)
        u = cache.uprev7

        cache.uprev7, cache.uprev6, cache.uprev5, cache.uprev4, cache.uprev3, cache.uprevprev = (uprev6,
                                                                                                 uprev5,
                                                                                                 uprev4,
                                                                                                 uprev3,
                                                                                                 uprevprev,
                                                                                                 uprev)
    else
        # increase step count
        cache.step += 1

        # evaluate production matrix
        P, d = evaluate_pds(f, uprev, p, t)
        integrator.stats.nf += 1

        P_tup = (P, P2, P3, P4, P5, P6, P7)
        d_tup = (d, d2, d3, d4, d5, d6, d7)
        u_tup = (uprev, uprevprev, uprev3, uprev4, uprev5, uprev6, uprev7)

        u = perform_step_MPLM75(P_tup, d_tup, dt, u_tup, alg.linsolve, αβ,
                                small_constant)
        integrator.stats.nsolve += 5

        cache.uprev7, cache.uprev6, cache.uprev5, cache.uprev4, cache.uprev3,
        cache.uprevprev = (uprev6,
                           uprev5,
                           uprev4,
                           uprev3,
                           uprevprev,
                           uprev)
    end

    integrator.u = u

    cache.P7, cache.P6, cache.P5, cache.P4, cache.P3, cache.P2 = (P6, P5, P4, P3, P2, P)
    cache.d7, cache.d6, cache.d5, cache.d4, cache.d3, cache.d2 = (d6, d5, d4, d3, d2, d)
end

@muladd function perform_step!(integrator, cache::MPLM75Cache, repeat_step = false)
    (; alg, t, dt, uprev, uprev2, u, f, p) = integrator
    (; uprevprev, uprev3, uprev4, uprev5, uprev6, uprev7,
    v, vprev, vprev2, vprev3, vprev4, vprev5, P, P2, P3, P4, P5, P6, P7,
    d, d2, d3, d4, d5, d6, d7, d_tmp, σ, αβ, small_constant, linsolve) = cache
    #TODO Check if number of v-vectors can be reduced. 
    # vprevX are only used in the initialization phase.

    # TODO: Check that this actually works
    if integrator.derivative_discontinuity
        cache.step = 1
    end

    if cache.step == 1
        # increase step count
        cache.step += 1

        # initilialze v vectors 
        vprev .= uprev
        vprev2 .= uprev
        vprev3 .= uprev
        vprev4 .= uprev
        vprev5 .= uprev

        # evaluate production matrix at tspan[1]
        evaluate_pds!(P, d, f, uprev, p, t)
        integrator.stats.nf += 1

        # save current P and d
        P7 .= P
        !isnothing(d) && (d7 .= d)

        # compute initial values 
        # we use uprevprev as temporary storage for the value of u needed in step 1.
        # we use uprev3 as temporary storage for the value of u needed in step 2.
        # we use uprev4 as temporary storage for the value of u needed in step 3.
        # we use uprev5 as temporary storage for the value of u needed in step 4.
        # we use uprev6 as temporary storage for the value of u needed in step 5.
        # we use v as temporary storage for the value of u needed in step 6.
        _, nf, ns = start_MPLM75!(v, vprev, vprev2, vprev3, vprev4, vprev5,
                                  uprevprev, uprev3, uprev4, uprev5, uprev6,
                                  P, P2, P3, P4, P5, d, d2, d3, d4, d5, d_tmp,
                                  t, dt, σ, f, p, small_constant, linsolve;
                                  substep_exp = alg.substep_level + 1)
        integrator.stats.nf += nf
        integrator.stats.nsolve += ns

        # reset P
        P .= P7
        !isnothing(d7) && (d .= d7)

        # u at time tspan[1] + dt
        u .= uprevprev

        shift!(uprev, uprevprev)
    elseif cache.step == 2
        # increase step count
        cache.step += 1

        # evaluate production matrix at tspan[1] + dt
        evaluate_pds!(P, d, f, uprev, p, t)
        integrator.stats.nf += 1

        # u at time tspan[1] + 2*dt (this was computed in step 1)
        u .= uprev3

        shift!(uprev, uprevprev, uprev3)
    elseif cache.step == 3
        # increase step count
        cache.step += 1

        # evaluate production matrix at tspan[1] + 2*dt
        evaluate_pds!(P, d, f, uprev, p, t)
        integrator.stats.nf += 1

        # u at time tspan[1] + 3*dt (this was computed in step 1)
        u .= uprev4

        shift!(uprev, uprevprev, uprev3, uprev4)
    elseif cache.step == 4
        # increase step count
        cache.step += 1

        # evaluate production matrix at tspan[1] + 3*dt
        evaluate_pds!(P, d, f, uprev, p, t)
        integrator.stats.nf += 1

        # u at time tspan[1] + 4*dt (this was computed in step 1)
        u .= uprev5

        shift!(uprev, uprevprev, uprev3, uprev4, uprev5)
    elseif cache.step == 5
        # increase step count
        cache.step += 1

        # evaluate production matrix at tspan[1] + 4*dt
        evaluate_pds!(P, d, f, uprev, p, t)
        integrator.stats.nf += 1

        # u at time tspan[1] + 5*dt (this was computed in step 1)
        u .= uprev6

        shift!(uprev, uprevprev, uprev3, uprev4, uprev5, uprev6)
    elseif cache.step == 6
        # increase step count
        cache.step += 1

        # evaluate production matrix at tspan[1] + 5*dt
        evaluate_pds!(P, d, f, uprev, p, t)
        integrator.stats.nf += 1

        # u at time tspan[1] + 6*dt (this was computed in step 1)
        u .= v

        shift!(uprev, uprevprev, uprev3, uprev4, uprev5, uprev6, uprev7)
    else
        # increase step count
        cache.step += 1

        # evaluate production matrix
        evaluate_pds!(P, d, f, uprev, p, t)
        integrator.stats.nf += 1

        P_tup = (P, P2, P3, P4, P5, P6, P7)
        d_tup = (d, d2, d3, d4, d5, d6, d7)
        u_tup = (uprev, uprevprev, uprev3, uprev4, uprev5, uprev6, uprev7)

        perform_step_MPLM75!(u, P_tup, d_tup, d_tmp, dt, u_tup, σ, linsolve, αβ,
                             small_constant)
        integrator.stats.nsolve += 5

        shift!(uprev, uprevprev, uprev3, uprev4, uprev5, uprev6, uprev7)
    end

    shift!(P, P2, P3, P4, P5, P6, P7)
    shift!(d, d2, d3, d4, d5, d6, d7)
end

#### MPLM106 ############################################################################
@muladd function start_MPLM106(P, d, t, dt, uprev, f, p, small_constant, linsolve;
                               substep_exp = 2)
    αβ75 = get_constant_parameters(MPLM75(), eltype(uprev))

    substeps = 2^substep_exp
    dts = dt / substeps

    ### 1.5 macro steps ###############################################################
    v, t, nf, ns = start_MPLM75(P, d, t, dts, uprev, f, p, small_constant, linsolve;
                                substep_exp)

    uprev6 = uprev
    P6 = P
    d6 = d

    uprev5 = v[1]
    P5, d5 = evaluate_pds(f, uprev5, p, t - 5 * dts)
    nf += 1

    uprev4 = v[2]
    P4, d4 = evaluate_pds(f, uprev4, p, t - 4 * dts)
    nf += 1

    uprev3 = v[3]
    P3, d3 = evaluate_pds(f, uprev3, p, t - 3 * dts)
    nf += 1

    uprevprev = v[4]
    P2, d2 = evaluate_pds(f, uprevprev, p, t - 2 * dts)
    nf += 1

    uprev = v[5]
    P, d = evaluate_pds(f, uprev, p, t - dts)
    nf += 1

    u = v[6]

    #v1 = v[4]
    local v1
    if substeps == 4
        v1 = v[4]
    end

    ### second macro step ############################################################
    # substeps 7 - 2*substeps
    for i in 1:(2 * substeps - 6)
        uprev7, uprev6, uprev5, uprev4, uprev3, uprevprev, uprev = (uprev6, uprev5, uprev4,
                                                                    uprev3, uprevprev,
                                                                    uprev, u)
        P7, P6, P5, P4, P3, P2 = (P6, P5, P4, P3, P2, P)
        d7, d6, d5, d4, d3, d2 = (d6, d5, d4, d3, d2, d)

        P, d = evaluate_pds(f, uprev, p, t)
        nf += 1

        P_tup = (P, P2, P3, P4, P5, P6, P7)
        d_tup = (d, d2, d3, d4, d5, d6, d7)
        u_tup = (uprev, uprevprev, uprev3, uprev4, uprev5, uprev6, uprev7)

        u = perform_step_MPLM75(P_tup, d_tup, dts, u_tup, linsolve, αβ75,
                                small_constant)
        t += dts
        ns += 5

        if substeps > 4 && i == (substeps - 6)
            v1 = u
        end
    end

    v2 = u

    ### third macro step ############################################################
    for _ in 1:substeps
        uprev7, uprev6, uprev5, uprev4, uprev3, uprevprev, uprev = (uprev6, uprev5, uprev4,
                                                                    uprev3, uprevprev,
                                                                    uprev, u)
        P7, P6, P5, P4, P3, P2 = (P6, P5, P4, P3, P2, P)
        d7, d6, d5, d4, d3, d2 = (d6, d5, d4, d3, d2, d)

        P, d = evaluate_pds(f, uprev, p, t)
        nf += 1

        P_tup = (P, P2, P3, P4, P5, P6, P7)
        d_tup = (d, d2, d3, d4, d5, d6, d7)
        u_tup = (uprev, uprevprev, uprev3, uprev4, uprev5, uprev6, uprev7)

        u = perform_step_MPLM75(P_tup, d_tup, dts, u_tup, linsolve, αβ75,
                                small_constant)
        t += dts
        ns += 5
    end

    v3 = u

    ### fourth macro step ############################################################
    for _ in 1:substeps
        uprev7, uprev6, uprev5, uprev4, uprev3, uprevprev, uprev = (uprev6, uprev5, uprev4,
                                                                    uprev3, uprevprev,
                                                                    uprev, u)
        P7, P6, P5, P4, P3, P2 = (P6, P5, P4, P3, P2, P)
        d7, d6, d5, d4, d3, d2 = (d6, d5, d4, d3, d2, d)

        P, d = evaluate_pds(f, uprev, p, t)
        nf += 1

        P_tup = (P, P2, P3, P4, P5, P6, P7)
        d_tup = (d, d2, d3, d4, d5, d6, d7)
        u_tup = (uprev, uprevprev, uprev3, uprev4, uprev5, uprev6, uprev7)

        u = perform_step_MPLM75(P_tup, d_tup, dts, u_tup, linsolve, αβ75,
                                small_constant)
        t += dts
        ns += 5
    end

    v4 = u

    ### fifth macro step ############################################################
    for _ in 1:substeps
        uprev7, uprev6, uprev5, uprev4, uprev3, uprevprev, uprev = (uprev6, uprev5, uprev4,
                                                                    uprev3, uprevprev,
                                                                    uprev, u)
        P7, P6, P5, P4, P3, P2 = (P6, P5, P4, P3, P2, P)
        d7, d6, d5, d4, d3, d2 = (d6, d5, d4, d3, d2, d)

        P, d = evaluate_pds(f, uprev, p, t)
        nf += 1

        P_tup = (P, P2, P3, P4, P5, P6, P7)
        d_tup = (d, d2, d3, d4, d5, d6, d7)
        u_tup = (uprev, uprevprev, uprev3, uprev4, uprev5, uprev6, uprev7)

        u = perform_step_MPLM75(P_tup, d_tup, dts, u_tup, linsolve, αβ75,
                                small_constant)
        t += dts
        ns += 5
    end

    v5 = u

    ### sixth macro step ############################################################
    for _ in 1:substeps
        uprev7, uprev6, uprev5, uprev4, uprev3, uprevprev, uprev = (uprev6, uprev5, uprev4,
                                                                    uprev3, uprevprev,
                                                                    uprev, u)
        P7, P6, P5, P4, P3, P2 = (P6, P5, P4, P3, P2, P)
        d7, d6, d5, d4, d3, d2 = (d6, d5, d4, d3, d2, d)

        P, d = evaluate_pds(f, uprev, p, t)
        nf += 1

        P_tup = (P, P2, P3, P4, P5, P6, P7)
        d_tup = (d, d2, d3, d4, d5, d6, d7)
        u_tup = (uprev, uprevprev, uprev3, uprev4, uprev5, uprev6, uprev7)

        u = perform_step_MPLM75(P_tup, d_tup, dts, u_tup, linsolve, αβ75,
                                small_constant)
        t += dts
        ns += 5
    end

    v6 = u

    ### seventh macro step ############################################################
    for _ in 1:substeps
        uprev7, uprev6, uprev5, uprev4, uprev3, uprevprev, uprev = (uprev6, uprev5, uprev4,
                                                                    uprev3, uprevprev,
                                                                    uprev, u)
        P7, P6, P5, P4, P3, P2 = (P6, P5, P4, P3, P2, P)
        d7, d6, d5, d4, d3, d2 = (d6, d5, d4, d3, d2, d)

        P, d = evaluate_pds(f, uprev, p, t)
        nf += 1

        P_tup = (P, P2, P3, P4, P5, P6, P7)
        d_tup = (d, d2, d3, d4, d5, d6, d7)
        u_tup = (uprev, uprevprev, uprev3, uprev4, uprev5, uprev6, uprev7)

        u = perform_step_MPLM75(P_tup, d_tup, dts, u_tup, linsolve, αβ75,
                                small_constant)
        t += dts
        ns += 5
    end

    v7 = u

    ### eighth macro step ############################################################
    for _ in 1:substeps
        uprev7, uprev6, uprev5, uprev4, uprev3, uprevprev, uprev = (uprev6, uprev5, uprev4,
                                                                    uprev3, uprevprev,
                                                                    uprev, u)
        P7, P6, P5, P4, P3, P2 = (P6, P5, P4, P3, P2, P)
        d7, d6, d5, d4, d3, d2 = (d6, d5, d4, d3, d2, d)

        P, d = evaluate_pds(f, uprev, p, t)
        nf += 1

        P_tup = (P, P2, P3, P4, P5, P6, P7)
        d_tup = (d, d2, d3, d4, d5, d6, d7)
        u_tup = (uprev, uprevprev, uprev3, uprev4, uprev5, uprev6, uprev7)

        u = perform_step_MPLM75(P_tup, d_tup, dts, u_tup, linsolve, αβ75,
                                small_constant)
        t += dts
        ns += 5
    end

    v8 = u

    ### ninth macro step ############################################################
    for _ in 1:substeps
        uprev7, uprev6, uprev5, uprev4, uprev3, uprevprev, uprev = (uprev6, uprev5, uprev4,
                                                                    uprev3, uprevprev,
                                                                    uprev, u)
        P7, P6, P5, P4, P3, P2 = (P6, P5, P4, P3, P2, P)
        d7, d6, d5, d4, d3, d2 = (d6, d5, d4, d3, d2, d)

        P, d = evaluate_pds(f, uprev, p, t)
        nf += 1

        P_tup = (P, P2, P3, P4, P5, P6, P7)
        d_tup = (d, d2, d3, d4, d5, d6, d7)
        u_tup = (uprev, uprevprev, uprev3, uprev4, uprev5, uprev6, uprev7)

        u = perform_step_MPLM75(P_tup, d_tup, dts, u_tup, linsolve, αβ75,
                                small_constant)
        t += dts
        ns += 5
    end

    return (v1, v2, v3, v4, v5, v6, v7, v8, u), t, nf, ns
end

@muladd function start_MPLM106!(v, vprev, vprev2, vprev3, vprev4, vprev5, vprev6, vprev7,
                                tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8,
                                P, P2, P3, P4, P5, P6, P7, d, d2, d3, d4, d5, d6, d7, d_tmp,
                                t, dt, σ, f, p, small_constant, linsolve; substep_exp = 2)
    αβ75 = get_constant_parameters(MPLM75(), eltype(vprev))

    tmps = (tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8)
    P_tup = (P, P2, P3, P4, P5, P6, P7)
    d_tup = (d, d2, d3, d4, d5, d6, d7)
    v_tup = (vprev, vprev2, vprev3, vprev4, vprev5, vprev6, vprev7)

    substeps = 2^substep_exp
    dts = dt / substeps

    # save current P and d
    P6 .= P
    !isnothing(d) && (d6 .= d)

    ### 1.5 macro steps ###############################################################
    t, nf, ns = start_MPLM75!(v, vprev, vprev2, vprev3, vprev4, vprev5,
                              tmp, tmp2, tmp3, tmp4, tmp5,
                              P, P2, P3, P4, P5, d, d2, d3, d4, d5, d_tmp,
                              t, dts, σ, f, p, small_constant, linsolve; substep_exp)

    # initialize MPLM75                           
    # vprev6 must be initialized as uprev
    vprev5 .= tmp
    evaluate_pds!(P5, d5, f, vprev5, p, t - 5 * dts)
    nf += 1

    vprev4 .= tmp2
    evaluate_pds!(P4, d4, f, vprev4, p, t - 4 * dts)
    nf += 1

    vprev3 .= tmp3
    evaluate_pds!(P3, d3, f, vprev3, p, t - 3 * dts)
    nf += 1

    vprev2 .= tmp4
    evaluate_pds!(P2, d2, f, vprev2, p, t - 2 * dts)
    nf += 1

    vprev .= tmp5
    evaluate_pds!(P, d, f, vprev, p, t - dts)
    nf += 1

    if substeps == 4
        tmp .= tmp4
    end

    sub_steps = (2 * substeps - 6, substeps, substeps, substeps, substeps, substeps,
                 substeps, substeps)

    for (step_idx, n_iter) in enumerate(sub_steps)
        # step_idx 1 corresponds to "second macro step"
        # we need to use tmps[step_idx + 1] because tmp was already handled 

        for i in 1:n_iter
            shift!(v, v_tup...)
            shift!(P_tup...)
            shift!(d_tup...)

            evaluate_pds!(P, d, f, vprev, p, t)
            nf += 1

            perform_step_MPLM75!(v, P_tup, d_tup, d_tmp, dts, v_tup, σ,
                                 linsolve, αβ75,
                                 small_constant)
            t += dts
            ns += 5

            if step_idx == 1 && substeps > 4 && i == (substeps - 6)
                tmp .= v
            end
        end

        # save initial data
        if step_idx < length(sub_steps)
            tmps[step_idx + 1] .= v
        end
    end

    return nf, ns
end

@muladd function perform_step_MPLM106(P_tup, d_tup, dt, u_tup, linsolve, αβ,
                                      small_constant)
    P, P2, P3, P4, P5, P6, P7, P8, P9, P10 = P_tup
    d, d2, d3, d4, d5, d6, d7, d8, d9, d10 = d_tup
    uprev, uprevprev, uprev3, uprev4, uprev5, uprev6, uprev7, uprev8, uprev9, uprev10 = u_tup
    α1, α2, α3, α4, α5, α6, α7, α8, α9, α10, β1, β2, β3, β4, β5, β6, β7, β8, β9, β10 = αβ

    # σ approximations 
    σ = sigma_approx_1(uprev, P, d, dt, linsolve, small_constant)
    σ = sigma_approx_2(σ, uprevprev, P, d, dt, linsolve, small_constant)
    σ = sigma_approx_3(σ, uprev, uprev3, P_tup, d_tup, dt, linsolve, small_constant)
    σ = sigma_approx_4(σ, uprev5, P_tup, d_tup, dt, linsolve, small_constant)
    σ = sigma_approx_5(σ, uprev7, P_tup, d_tup, dt, linsolve, small_constant)

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

@muladd function perform_step_MPLM106!(u, P_tup, d_tup, d_tmp, dt, u_tup, σ, linsolve, αβ,
                                       small_constant)
    P, P2, P3, P4, P5, P6, P7, P8, P9, P10 = P_tup
    d, d2, d3, d4, d5, d6, d7, d8, d9, d10 = d_tup
    uprev, uprevprev, uprev3, uprev4, uprev5, uprev6, uprev7, uprev8, uprev9, uprev10 = u_tup
    α1, α2, α3, α4, α5, α6, α7, α8, α9, α10, β1, β2, β3, β4, β5, β6, β7, β8, β9, β10 = αβ

    # σ approximations
    sigma_approx_1!(σ, uprev, P, d, dt, linsolve, small_constant)
    sigma_approx_2!(σ, uprevprev, P, d, dt, linsolve, small_constant)
    sigma_approx_3!(σ, uprev, uprev3, P_tup, d_tup, d_tmp, dt, linsolve, small_constant)
    sigma_approx_4!(σ, uprev5, P_tup, d_tup, d_tmp, dt, linsolve, small_constant)
    sigma_approx_5!(σ, uprev7, P_tup, d_tup, d_tmp, dt, linsolve, small_constant)

    # Main step 
    @.. broadcast=false σ=σ + small_constant
    lincomb!(P10, β1, P, β2, P2, β3, P3, β4, P4, β5, P5, β6, P6, β7, P7, β8, P8, β9, P9,
             β10, P10)
    lincomb!(d10, β1, d, β2, d2, β3, d3, β4, d4, β5, d5, β6, d6, β7, d7, β8, d8, β9, d9,
             β10, d10)
    @.. broadcast=false linsolve.b=α1 * uprev + α2 * uprevprev + α3 * uprev3 + α4 * uprev4 +
                                   α5 * uprev5 + α6 * uprev6 + α7 * uprev7 + α8 * uprev8 +
                                   α9 * uprev9 + α10 * uprev10
    basic_patankar_step!(u, linsolve.b, P10, d10, linsolve.A, σ, dt, linsolve)

    # statistics: 6 nsolves

    return nothing
end

@muladd function perform_step!(integrator, cache::MPLM106oopCache, repeat_step = false)
    (; alg, t, dt, uprev, uprev2, f, p) = integrator
    (; uprevprev, uprev3, uprev4, uprev5, uprev6, uprev7, uprev8, uprev9, uprev10,
    P2, P3, P4, P5, P6, P7, P8, P9, P10, d2, d3, d4, d5, d6, d7, d8, d9, d10,
    αβ, small_constant) = cache

    # TODO: Check that this actually works
    if integrator.derivative_discontinuity
        cache.step = 1
    end

    if cache.step == 1
        # increase step count
        cache.step += 1

        # evaluate production matrix at tspan[1]
        P, d = evaluate_pds(f, uprev, p, t)
        integrator.stats.nf += 1

        # compute initial values
        v, _, nf, ns = start_MPLM106(P, d, t, dt, uprev, f, p, small_constant,
                                     alg.linsolve; substep_exp = alg.substep_level + 1)
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

        cache.uprev3, cache.uprevprev = (uprevprev, uprev)
    elseif cache.step == 3
        # increase step count
        cache.step += 1

        # evaluate production matrix at tspan[1] + 2*dt
        P, d = evaluate_pds(f, uprev, p, t)
        integrator.stats.nf += 1

        # u at time tspan[1] + 3*dt (this was computed in step 1)
        u = cache.uprev4

        cache.uprev4, cache.uprev3, cache.uprevprev = (uprev3, uprevprev, uprev)
    elseif cache.step == 4
        # increase step count
        cache.step += 1

        # evaluate production matrix at tspan[1] + 3*dt
        P, d = evaluate_pds(f, uprev, p, t)
        integrator.stats.nf += 1

        # u at time tspan[1] + 4*dt (this was computed in step 1)
        u = cache.uprev5

        cache.uprev5, cache.uprev4, cache.uprev3, cache.uprevprev = (uprev4, uprev3,
                                                                     uprevprev, uprev)
    elseif cache.step == 5
        # increase step count
        cache.step += 1

        # evaluate production matrix at tspan[1] + 4*dt
        P, d = evaluate_pds(f, uprev, p, t)
        integrator.stats.nf += 1

        # u at time tspan[1] + 5*dt (this was computed in step 1)
        u = cache.uprev6

        cache.uprev6, cache.uprev5, cache.uprev4, cache.uprev3, cache.uprevprev = (uprev5,
                                                                                   uprev4,
                                                                                   uprev3,
                                                                                   uprevprev,
                                                                                   uprev)
    elseif cache.step == 6
        # increase step count
        cache.step += 1

        # evaluate production matrix at tspan[1] + 5*dt
        P, d = evaluate_pds(f, uprev, p, t)
        integrator.stats.nf += 1

        # u at time tspan[1] + 6*dt (this was computed in step 1)
        u = cache.uprev7

        cache.uprev7, cache.uprev6, cache.uprev5, cache.uprev4, cache.uprev3,
        cache.uprevprev = (uprev6,
                           uprev5,
                           uprev4,
                           uprev3,
                           uprevprev,
                           uprev)
    elseif cache.step == 7
        # increase step count
        cache.step += 1

        # evaluate production matrix at tspan[1] + 6*dt
        P, d = evaluate_pds(f, uprev, p, t)
        integrator.stats.nf += 1

        # u at time tspan[1] + 7*dt (this was computed in step 1)
        u = cache.uprev8

        cache.uprev8, cache.uprev7, cache.uprev6, cache.uprev5, cache.uprev4, cache.uprev3,
        cache.uprevprev = (uprev7,
                           uprev6,
                           uprev5,
                           uprev4,
                           uprev3,
                           uprevprev,
                           uprev)
    elseif cache.step == 8
        # increase step count
        cache.step += 1

        # evaluate production matrix at tspan[1] + 7*dt
        P, d = evaluate_pds(f, uprev, p, t)
        integrator.stats.nf += 1

        # u at time tspan[1] + 8*dt (this was computed in step 1)
        u = cache.uprev9

        cache.uprev9, cache.uprev8, cache.uprev7, cache.uprev6, cache.uprev5, cache.uprev4,
        cache.uprev3, cache.uprevprev = (uprev8,
                                         uprev7,
                                         uprev6,
                                         uprev5,
                                         uprev4,
                                         uprev3,
                                         uprevprev,
                                         uprev)
    elseif cache.step == 9
        # increase step count
        cache.step += 1

        # evaluate production matrix at tspan[1] + 8*dt
        P, d = evaluate_pds(f, uprev, p, t)
        integrator.stats.nf += 1

        # u at time tspan[1] + 9*dt (this was computed in step 1)
        u = cache.uprev10

        cache.uprev10, cache.uprev9, cache.uprev8, cache.uprev7, cache.uprev6, cache.uprev5,
        cache.uprev4, cache.uprev3, cache.uprevprev = (uprev9,
                                                       uprev8,
                                                       uprev7,
                                                       uprev6,
                                                       uprev5,
                                                       uprev4,
                                                       uprev3,
                                                       uprevprev,
                                                       uprev)
    else

        # evaluate production matrix
        P, d = evaluate_pds(f, uprev, p, t)
        integrator.stats.nf += 1

        P_tup = (P, P2, P3, P4, P5, P6, P7, P8, P9, P10)
        d_tup = (d, d2, d3, d4, d5, d6, d7, d8, d9, d10)
        u_tup = (uprev, uprevprev, uprev3, uprev4, uprev5, uprev6, uprev7, uprev8, uprev9,
                 uprev10)

        u = perform_step_MPLM106(P_tup, d_tup, dt, u_tup, alg.linsolve, αβ,
                                 small_constant)
        integrator.stats.nsolve += 6

        cache.uprev10, cache.uprev9, cache.uprev8, cache.uprev7, cache.uprev6, cache.uprev5,
        cache.uprev4, cache.uprev3, cache.uprevprev = (uprev9, uprev8, uprev7, uprev6,
                                                       uprev5, uprev4, uprev3, uprevprev,
                                                       uprev)
    end

    integrator.u = u

    cache.P10, cache.P9, cache.P8, cache.P7, cache.P6, cache.P5, cache.P4, cache.P3,
    cache.P2 = (P9,
                P8,
                P7,
                P6,
                P5,
                P4,
                P3,
                P2,
                P)
    cache.d10, cache.d9, cache.d8, cache.d7, cache.d6, cache.d5, cache.d4, cache.d3,
    cache.d2 = (d9,
                d8,
                d7,
                d6,
                d5,
                d4,
                d3,
                d2,
                d)
end

@muladd function perform_step!(integrator, cache::MPLM106Cache, repeat_step = false)
    (; alg, t, dt, uprev, uprev2, u, f, p) = integrator
    (; uprevprev, uprev3, uprev4, uprev5, uprev6, uprev7, uprev8, uprev9, uprev10,
    v, vprev, vprev2, vprev3, vprev4, vprev5, vprev6, vprev7, P, P2, P3, P4, P5, P6, P7, P8, P9, P10,
    d, d2, d3, d4, d5, d6, d7, d8, d9, d10, d_tmp, σ, αβ, small_constant, linsolve) = cache
    #TODO Check if number of v-vectors can be reduced. 
    # vprevX are only used in the initialization phase.

    # TODO: Check that this actually works
    if integrator.derivative_discontinuity
        cache.step = 1
    end

    if cache.step == 1
        # increase step count
        cache.step += 1

        # initilialze v vectors 
        vprev .= uprev
        vprev2 .= uprev
        vprev3 .= uprev
        vprev4 .= uprev
        vprev5 .= uprev
        vprev6 .= uprev
        vprev7 .= uprev

        # evaluate production matrix at tspan[1]
        evaluate_pds!(P, d, f, uprev, p, t)
        integrator.stats.nf += 1

        # save current P and d
        P10 .= P
        !isnothing(d) && (d10 .= d)

        # compute initial values
        # we use uprevprev as temporary storage for the value of u needed in step 1.
        # we use uprev3 as temporary storage for the value of u needed in step 2.
        # we use uprev4 as temporary storage for the value of u needed in step 3.
        # we use uprev5 as temporary storage for the value of u needed in step 4.
        # we use uprev6 as temporary storage for the value of u needed in step 5.
        # we use uprev7 as temporary storage for the value of u needed in step 6.
        # we use uprev8 as temporary storage for the value of u needed in step 7.
        # we use uprev9 as temporary storage for the value of u needed in step 8.
        # we use v as temporary storage for the value of u needed in step 9.  
        nf, ns = start_MPLM106!(v, vprev, vprev2, vprev3, vprev4, vprev5, vprev6, vprev7,
                                uprevprev, uprev3, uprev4, uprev5, uprev6, uprev7, uprev8,
                                uprev9,
                                P, P2, P3, P4, P5, P6, P7, d, d2, d3, d4, d5, d6, d7, d_tmp,
                                t, dt, σ, f, p, small_constant, linsolve;
                                substep_exp = alg.substep_level + 1)
        integrator.stats.nf += nf
        integrator.stats.nsolve += ns

        # reset P
        P .= P10
        !isnothing(d10) && (d .= d10)

        # u at time tspan[1] + dt
        u .= uprevprev

        shift!(uprev, uprevprev)
    elseif cache.step == 2
        # increase step count
        cache.step += 1

        # evaluate production matrix at tspan[1] + dt
        evaluate_pds!(P, d, f, uprev, p, t)
        integrator.stats.nf += 1

        # u at time tspan[1] + 2*dt (this was computed in step 1)
        u .= uprev3

        shift!(uprev, uprevprev, uprev3)
    elseif cache.step == 3
        # increase step count
        cache.step += 1

        # evaluate production matrix at tspan[1] + 2*dt
        evaluate_pds!(P, d, f, uprev, p, t)
        integrator.stats.nf += 1

        # u at time tspan[1] + 3*dt (this was computed in step 1)
        u .= uprev4

        shift!(uprev, uprevprev, uprev3, uprev4)
    elseif cache.step == 4
        # increase step count
        cache.step += 1

        # evaluate production matrix at tspan[1] + 3*dt
        evaluate_pds!(P, d, f, uprev, p, t)
        integrator.stats.nf += 1

        # u at time tspan[1] + 4*dt (this was computed in step 1)
        u .= uprev5

        shift!(uprev, uprevprev, uprev3, uprev4, uprev5)
    elseif cache.step == 5
        # increase step count
        cache.step += 1

        # evaluate production matrix at tspan[1] + 4*dt
        evaluate_pds!(P, d, f, uprev, p, t)
        integrator.stats.nf += 1

        # u at time tspan[1] + 5*dt (this was computed in step 1)
        u .= uprev6

        shift!(uprev, uprevprev, uprev3, uprev4, uprev5, uprev6)
    elseif cache.step == 6
        # increase step count
        cache.step += 1

        # evaluate production matrix at tspan[1] + 5*dt
        evaluate_pds!(P, d, f, uprev, p, t)
        integrator.stats.nf += 1

        # u at time tspan[1] + 6*dt (this was computed in step 1)
        u .= uprev7

        shift!(uprev, uprevprev, uprev3, uprev4, uprev5, uprev6, uprev7)
    elseif cache.step == 7
        # increase step count
        cache.step += 1

        # evaluate production matrix at tspan[1] + 6*dt
        evaluate_pds!(P, d, f, uprev, p, t)
        integrator.stats.nf += 1

        # u at time tspan[1] + 7*dt (this was computed in step 1)
        u .= uprev8

        shift!(uprev, uprevprev, uprev3, uprev4, uprev5, uprev6, uprev7, uprev8)
    elseif cache.step == 8
        # increase step count
        cache.step += 1

        # evaluate production matrix at tspan[1] + 7*dt
        evaluate_pds!(P, d, f, uprev, p, t)
        integrator.stats.nf += 1

        # u at time tspan[1] + 8*dt (this was computed in step 1)
        u .= uprev9

        shift!(uprev, uprevprev, uprev3, uprev4, uprev5, uprev6, uprev7, uprev8, uprev9)
    elseif cache.step == 9
        # increase step count
        cache.step += 1

        # evaluate production matrix at tspan[1] + 8*dt
        evaluate_pds!(P, d, f, uprev, p, t)
        integrator.stats.nf += 1

        # u at time tspan[1] + 9*dt (this was computed in step 1)
        u .= v

        shift!(uprev, uprevprev, uprev3, uprev4, uprev5, uprev6, uprev7, uprev8, uprev9,
               uprev10)
    else
        # increase step count
        cache.step += 1

        # evaluate production matrix
        evaluate_pds!(P, d, f, uprev, p, t)
        integrator.stats.nf += 1

        P_tup = (P, P2, P3, P4, P5, P6, P7, P8, P9, P10)
        d_tup = (d, d2, d3, d4, d5, d6, d7, d8, d9, d10)
        u_tup = (uprev, uprevprev, uprev3, uprev4, uprev5, uprev6, uprev7, uprev8, uprev9,
                 uprev10)

        perform_step_MPLM106!(u, P_tup, d_tup, d_tmp, dt, u_tup, σ, linsolve, αβ,
                              small_constant)
        integrator.stats.nsolve += 6

        shift!(u_tup...)
    end

    shift!(P, P2, P3, P4, P5, P6, P7, P8, P9, P10)
    shift!(d, d2, d3, d4, d5, d6, d7, d8, d9, d10)
end
