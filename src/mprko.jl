### MPRKO22 #####################################################################################
"""
    MPRKO22(α, β; [linsolve = ..., small_constant = ...])

A family of second-order modified Patankar-Runge-Kutta-Oliver algorithms for
production-destruction systems. Each member of this family is an adaptive, one-step, two-stage method which is
second-order accurate, unconditionally positivity-preserving, and linearly
implicit. The stage-values are conservative as well.
The parameters `α` and `β` are described by Ávila et al. (2021).

For an autonomous PDS, the scheme MPRKO22(α, β) coincides with MPRK22(α).

This method supports adaptive time stepping, using the Patankar-weight denominators
``σ_i``, see Kopecz and Meister (2018), as first order approximations to estimate the error.

The scheme was introduced by Ávila et al. for conservative production-destruction systems.
For nonconservative production–destruction systems we use a straight forward extension
analogous to [`MPE`](@ref).

This modified Patankar-Runge-Kutta method requires the special structure of a
[`PDSProblem`](@ref) or a [`ConservativePDSProblem`](@ref).

You can optionally choose the linear solver to be used by passing an
algorithm from [LinearSolve.jl](https://github.com/SciML/LinearSolve.jl)
as keyword argument `linsolve`.
You can also choose the parameter `small_constant` which is added to all Patankar-weight denominators
to avoid divisions by zero. You can pass a value explicitly, otherwise `small_constant` is set to
`floatmin` of the floating point type used.

## References

- Andrés I. Ávila, Galo Javier González, Stefan Kopecz, and Andreas Meister.
  "Extension of modified Patankar–Runge–Kutta schemes to nonautonomous production–destruction systems based on Oliver’s approach."
  Journal of Computational and Applied Mathematics 389 (2021): 113350.
  [DOI: 10.1016/j.cam.2020.113350](https://doi.org/10.1016/j.cam.2020.113350)
"""
struct MPRKO22{T, F, T2} <: OrdinaryDiffEqAdaptiveAlgorithm
    alpha::T
    beta::T
    linsolve::F
    small_constant_function::T2
end

function MPRKO22(alpha, beta; linsolve = LUFactorization(),
                 small_constant = nothing)
    if isnothing(small_constant)
        small_constant_function = floatmin
    elseif small_constant isa Number
        small_constant_function = Returns(small_constant)
    else # assume small_constant isa Function
        small_constant_function = small_constant
    end
    MPRKO22{typeof(alpha), typeof(linsolve), typeof(small_constant_function)}(alpha, beta,
                                                                              linsolve,
                                                                              small_constant_function)
end

alg_order(::MPRKO22) = 2
isfsal(::MPRKO22) = false

function get_constant_parameters(alg::MPRKO22)
    if !(1 / 2 ≤ alg.alpha)
        throw(ArgumentError("MPRKO22 requires α ≥ 1/2."))
    end
    lb = (alg.alpha - 1) / (2 * alg.alpha - 1)
    ub = (alg.alpha) / (2 * alg.alpha - 1)
    if (1 / 2 ≤ alg.alpha ≤ 1) && !(0 ≤ alg.beta ≤ 1)
        throw(ArgumentError("For α=$(alg.alpha) MPRKO22 requires 0 ≤ β ≤ 1."))
    elseif (1 < alg.alpha) && !(lb ≤ alg.beta ≤ ub)
        throw(ArgumentError("For α=$(alg.alpha) MPRKO22 requires $lb ≤ β ≤ $ub."))
    end

    a21 = alg.alpha
    b2 = 1 / (2 * a21)
    b1 = 1 - b2
    c1 = alg.beta
    c2 = alg.alpha - 2 * alg.alpha * alg.beta + alg.beta

    # This should never happen
    if !all((a21, b1, b2, c1, c2) .≥ 0)
        throw(ArgumentError("MPRKO22 requires nonnegative RK coefficients."))
    end
    return a21, b1, b2, c1, c2
end

struct MPRKO22ConstantCache{T} <: OrdinaryDiffEqConstantCache
    a21::T
    b1::T
    b2::T
    c1::T
    c2::T
    small_constant::T
end

# Out-of-place
function alg_cache(alg::MPRKO22, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
                   uprev, uprev2, f, t, dt, reltol, p, calck,
                   ::Val{false},
                   verbose) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    if !(f isa PDSFunction || f isa ConservativePDSFunction)
        throw(ArgumentError("MPRKO22 can only be applied to production-destruction systems"))
    end

    a21, b1, b2, c1, c2 = get_constant_parameters(alg)
    MPRKO22ConstantCache(a21, b1, b2, c1, c2, alg.small_constant_function(uEltypeNoUnits))
end

function initialize!(integrator, cache::MPRKO22ConstantCache)
end

@muladd function perform_step!(integrator, cache::MPRKO22ConstantCache, repeat_step = false)
    (; alg, t, dt, uprev, f, p) = integrator
    (; a21, b1, b2, c1, c2, small_constant) = cache

    # evaluate production matrix
    P, d = evaluate_pds(f, uprev, p, t + c1 * dt)
    integrator.stats.nf += 1

    Ptmp, dtmp = lincomb(a21, P, d)

    # avoid division by zero due to zero Patankar weights
    σ = add_small_constant(uprev, small_constant)

    u = basic_patankar_step(uprev, Ptmp, σ, dt, alg.linsolve, dtmp)
    integrator.stats.nsolve += 1

    # compute Patankar weight denominator
    if isone(a21)
        σ = u
    else
        # σ = σ .* (u ./ σ) .^ (1 / a21) # generated Infs when solving brusselator
        σ = σ .^ (1 - 1 / a21) .* u .^ (1 / a21)
    end
    # avoid division by zero due to zero Patankar weights
    σ = add_small_constant(σ, small_constant)

    P2, d2 = evaluate_pds(f, u, p, t + c2 * dt)
    integrator.stats.nf += 1

    Ptmp, dtmp = lincomb(b1, P, d, b2, P2, d2)

    u = basic_patankar_step(uprev, Ptmp, σ, dt, alg.linsolve, dtmp)
    integrator.stats.nsolve += 1

    # If a21 = 1 and c1 = 0, then σ is the MPE approximation, i.e. suited for stiff problems.
    # Otherwise, this is not clear.
    tmp = u - σ
    atmp = calculate_residuals(tmp, uprev, u, integrator.opts.abstol,
                               integrator.opts.reltol, integrator.opts.internalnorm, t)
    set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))

    integrator.u = u
end

struct MPRKO22Cache{uType, PType, dType, tabType, F} <: MPRKMutableCache
    tmp::uType
    P::PType
    P2::PType
    D::dType
    D2::dType
    σ::uType
    tab::tabType
    linsolve::F
end

get_tmp_cache(integrator, ::MPRKO22, cache::OrdinaryDiffEqMutableCache) = (cache.σ,)

# In-place
function alg_cache(alg::MPRKO22, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
                   uprev, uprev2, f, t, dt, reltol, p, calck,
                   ::Val{true},
                   verbose) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    a21, b1, b2, c1, c2 = get_constant_parameters(alg)
    tab = MPRKO22ConstantCache(a21, b1, b2, c1, c2,
                               alg.small_constant_function(uEltypeNoUnits))
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

        MPRKO22Cache(tmp, P, P2, nothing, nothing, σ,
                     tab, #MPRKO22ConstantCache
                     linsolve)
    elseif f isa PDSFunction
        linprob = LinearProblem(P2, _vec(tmp))
        linsolve = init(linprob, alg.linsolve,
                        alias = LinearSolve.LinearAliasSpecifier(; alias_A = true,
                                                                 alias_b = true),
                        assumptions = LinearSolve.OperatorAssumptions(true))

        MPRKO22Cache(tmp, P, P2,
                     similar(u), # D
                     similar(u), # D2
                     σ,
                     tab, #MPRKO22ConstantCache
                     linsolve)
    else
        throw(ArgumentError("MPRKO22 can only be applied to production-destruction systems"))
    end
end

function initialize!(integrator, cache::MPRKO22Cache)
end

@muladd function perform_step!(integrator, cache::MPRKO22Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; tmp, P, P2, D, D2, σ, linsolve) = cache
    (; a21, b1, b2, c1, c2, small_constant) = cache.tab

    # We use P2 to store the last evaluation of the PDS
    # as well as to store the system matrix of the linear system

    # evaluate production/destruction terms
    evaluate_pds!(P, D, f, uprev, p, t + c1 * dt)
    integrator.stats.nf += 1

    lincomb!(P2, a21, P)
    lincomb!(D2, a21, D)

    # avoid division by zero due to zero Patankar weights
    @.. broadcast=false σ=uprev + small_constant

    tmp .= uprev
    basic_patankar_step!(u, tmp, P2, D2, σ, dt, linsolve)
    integrator.stats.nsolve += 1

    if isone(a21)
        @.. broadcast=false σ=u + small_constant
    else
        @.. broadcast=false σ=σ^(1 - 1 / a21) * u^(1 / a21) + small_constant
    end

    # evaluate production/destruction terms
    evaluate_pds!(P2, D2, f, u, p, t + c2 * dt)
    integrator.stats.nf += 1

    lincomb!(P2, b1, P, b2, P2)
    lincomb!(D2, b1, D, b2, D2)

    tmp .= uprev
    basic_patankar_step!(u, tmp, P2, D2, σ, dt, linsolve)
    integrator.stats.nsolve += 1

    # Now σ stores the error estimate
    # If a21 = 1 and c1 = 0, then σ is the MPE approximation, i.e. suited for stiff problems.
    # Otherwise, this is not clear.
    @.. broadcast=false σ=u - σ

    # Now tmp stores error residuals
    calculate_residuals!(tmp, σ, uprev, u, integrator.opts.abstol,
                         integrator.opts.reltol, integrator.opts.internalnorm, t,
                         Serial())
    set_EEst!(integrator, integrator.opts.internalnorm(tmp, t))
end
