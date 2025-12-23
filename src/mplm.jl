###################################################################################################################
### GENERAL OBSERVATIONS ##########################################################################################
###################################################################################################################
### 1. Use of (; t, dt, uprev, u, f, p) = integrator instead of @unpack alg, t, dt, uprev, f, p = integrator
### 2. In mprk.jl and sspmprk.jl f = integrator.f is used although f is already unpacked from integrator
### 3. We should have a function for a single step of MPE. This function should be the building block for all MPRK
###    and SSPMPRK schemes. It can also be used in MPLM, see
###    perform_step!(integrator, cache::MPLM22ConstantCache, repeat_step = false)
###################################################################################################################

struct MPLM22{F, T} <: OrdinaryDiffEqAlgorithm
    linsolve::F
    small_constant_function::T
end

alg_order(alg::MPLM22) = 2
alg_extrapolates(alg::MPLM22) = true

###################################################################################################################
### The following started from a copy of AB3 
###################################################################################################################
@cache mutable struct MPLM22ConstantCache{uType,rateType,T} <: OrdinaryDiffEqConstantCache
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
    MPLM22ConstantCache(u,k2, 1, alg.small_constant_function(uEltypeNoUnits))
end

function initialize!(integrator,
        cache::Union{MPLM22ConstantCache})
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
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


        #TODO: The following part is copied from perform_step!(integrator, cache::MPEConstantCache, repeat_step = false)
        # We should have a function _perform_step_MPE! for this.
        # This _perform_step_MPE! should than also be used within the other MPRK and SSPMPRK schemes.  

        # evaluate production matrix
        P = f.p(uprev, p, t)
        integrator.stats.nf += 1

        # avoid division by zero due to zero Patankar weights
        σ = add_small_constant(uprev, small_constant)

        # build linear system matrix and rhs
        if f isa PDSFunction
            d = f.d(uprev, p, t)  # evaluate nonconservative destruction terms
            rhs = uprev + dt * diag(P)
            M = build_mprk_matrix(P, σ, dt, d)
            linprob = LinearProblem(M, rhs)
        else
            # f isa ConservativePDSFunction
            M = build_mprk_matrix(P, σ, dt)
            linprob = LinearProblem(M, uprev)
        end

        # solve linear system
        sol = solve(linprob, alg.linsolve)
        u = sol.u
        integrator.stats.nsolve += 1
    else
        #TODO: Again, we need a function _perform_step_MPE! for this. This will make things much clearer.

        # evaluate production matrix
        P = f.p(uprev, p, t)
        integrator.stats.nf += 1

        # avoid division by zero due to zero Patankar weights
        σ = add_small_constant(uprev, small_constant)

        # build linear system matrix and rhs
        if f isa PDSFunction
            d = f.d(uprev, p, t)  # evaluate nonconservative destruction terms
            rhs = uprev + dt * diag(P)
            M = build_mprk_matrix(P, σ, dt, d)
            linprob = LinearProblem(M, rhs)
        else
            # f isa ConservativePDSFunction
            M = build_mprk_matrix(P, σ, dt)
            linprob = LinearProblem(M, uprev)
        end

        # solve linear system
        sol = solve(linprob, alg.linsolve)
        σ = sol.u
        integrator.stats.nsolve += 1

        # avoid division by zero due to zero Patankar weights
        σ = add_small_constant(σ, small_constant)

        # build linear system matrix and rhs
        if f isa PDSFunction
            d = f.d(uprev, p, t)  # evaluate nonconservative destruction terms
            rhs = uprevprev + dt * diag(P)
            M = build_mprk_matrix(P, σ, 2*dt, d)
            linprob = LinearProblem(M, rhs)
        else
            # f isa ConservativePDSFunction
            M = build_mprk_matrix(P, σ, 2*dt)
            linprob = LinearProblem(M, uprevprev)
        end

        # solve linear system
        sol = solve(linprob, alg.linsolve)
        u = sol.u
        integrator.stats.nsolve += 1

    end

    #TODO: Should be possible to use uprev2. But uprev2 is currently not updated.
    cache.uprevprev = uprev

    #TODO: Don't need k1 and k2 and initialize! must be changed accordingly.
    k2 = k1
    cache.k2 = k2
    integrator.fsallast = f(u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
end