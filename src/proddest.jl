# New problem type PDSProblem
abstract type AbstractPDSProblem end

"""
    PDSProblem(P, D, u0, tspan, p = NullParameters();
               p_prototype = nothing,
               analytic = nothing,
               std_rhs = nothing,
               linear_invariants = nothing)

A structure describing a system of ordinary differential equations in form of a production-destruction system (PDS).
`P` denotes the function defining the production matrix ``P``.
The diagonal of ``P`` contains production terms without destruction counterparts.
`D` is the function defining the vector of destruction terms ``D`` without production counterparts.
`u0` is the vector of initial conditions and `tspan` the time span
`(t_initial, t_final)` of the problem. The optional argument `p` can be used
to pass additional parameters to the functions `P` and `D`.

The functions `P` and `D` can be used either in the out-of-place form with signature
`production_terms = P(u, p, t)` or the in-place form `P(production_terms, u, p, t)`.

### Keyword arguments: ###

- `p_prototype`: If `P` is given in in-place form, `p_prototype` or copies thereof are used to store evaluations of `P`.
  If `p_prototype` is not specified explicitly and `P` is in-place, then `p_prototype` will be internally
  set to `zeros(eltype(u0), (length(u0), length(u0)))`.
- `analytic`: The analytic solution of a PDS must be given in the form `f(u0,p,t)`.
  Specifying the analytic solution can be useful for plotting and convergence tests.
- `std_rhs`: The standard ODE right-hand side evaluation function callable
  as `du = std_rhs(u, p, t)` for the out-of-place form and
  as `std_rhs(du, u, p, t)` for the in-place form. Solvers that do not rely on
  the production-destruction representation of the ODE, will use this function
  instead to compute the solution. 
  
  `std_rhs` can also be passed as a [`SciMLBase.ODEFunction`](https://docs.sciml.ai/DiffEqDocs/stable/types/ode_types/#SciMLBase.ODEFunction). 
  This allows attaching additional information—such as a Jacobian prototype (`jac_prototype`) or sparsity pattern (`sparsity`)—for improving or accelerating the usage of `std_rhs`. 
  See the [`SciMLBase.ODEFunction` documentation](https://docs.sciml.ai/DiffEqDocs/stable/types/ode_types/#SciMLBase.ODEFunction) for details.

  If not specified, a default implementation calling `P` and `D` is used.
- `linear_invariants`: The rows of this matrix contain the linear invariants of the ODE. 
  Certain solvers or callbacks might require this matrix.
  Note that this feature is experimental and its API may change in future releases.

## References

- Hans Burchard, Eric Deleersnijder, and Andreas Meister.
  "A high-order conservative Patankar-type discretisation for stiff systems of
  production-destruction equations."
  Applied Numerical Mathematics 47.1 (2003): 1-30.
  [DOI: 10.1016/S0168-9274(03)00101-6](https://doi.org/10.1016/S0168-9274(03)00101-6)
"""
struct PDSProblem{iip} <: AbstractPDSProblem end

# New ODE function PDSFunction
struct PDSFunction{iip, specialize, P, D, PrototypeP, PrototypeD, StdRHS, Ta, LI} <:
       AbstractODEFunction{iip}
    p::P
    d::D
    p_prototype::PrototypeP
    d_prototype::PrototypeD
    std_rhs::StdRHS
    analytic::Ta
    linear_invariants::LI
end

# define behavior of PDSFunctions for non-existing fields
function Base.getproperty(obj::PDSFunction, sym::Symbol)
    if sym === :mass_matrix
        return I
    elseif sym in fieldnames(typeof(obj))
        return getfield(obj, sym)
    elseif obj.std_rhs isa SciMLBase.AbstractODEFunction
        val = getproperty(obj.std_rhs, sym)
        if sym === :sparsity && val === nothing
            return getproperty(obj.std_rhs, :jac_prototype)
        end
        return val
    elseif sym === :jac_prototype || sym === :sparsity || sym === :sys
        return nothing
    else
        return getfield(obj, sym)
    end
end

# Most general constructor for PDSProblems
function PDSProblem(P, D, u0, tspan, p = NullParameters(); kwargs...)
    Piip = isinplace(P, 4)
    Diip = isinplace(D, 4)
    if Piip == Diip
        iip = Piip
    else
        error("Conflict due to the joint use of in-place and out-of-place functions.")
    end
    return PDSProblem{iip}(P, D, u0, tspan, p; kwargs...)
end

# Specialized constructor for PDSProblems setting `iip` manually
# (arbitrary functions)
function PDSProblem{iip}(P,
                         D,
                         u0,
                         tspan,
                         p = NullParameters();
                         p_prototype = nothing,
                         analytic = nothing,
                         std_rhs = nothing,
                         linear_invariants = nothing,
                         kwargs...,) where {iip}

    # p_prototype is used to store evaluations of P, if P is in-place.
    if isnothing(p_prototype) && iip
        p_prototype = zeros(eltype(u0), (length(u0), length(u0))) / oneunit(first(tspan))
    end
    # If a PDSFunction is to be evaluated and D is in-place, then d_prototype is used to store
    # evaluations of D.
    d_prototype = similar(u0 ./ oneunit(first(tspan)))

    PD = PDSFunction{iip}(P, D; p_prototype, d_prototype, analytic, std_rhs,
                          linear_invariants)
    PDSProblem{iip}(PD, u0, tspan, p; kwargs...)
end

# Specialized constructor for PDSProblems setting `iip` manually
# (PDSFunction)
function PDSProblem{iip}(PD::PDSFunction, u0, tspan, p = NullParameters();
                         kwargs...) where {iip}
    ODEProblem(PD, u0, tspan, p, PDSProblem{iip}(); kwargs...)
end

# Specialized constructor for PDSFunction setting `iip` manually
function PDSFunction{iip}(P, D; kwargs...) where {iip}
    PDSFunction{iip, FullSpecialize}(P, D; kwargs...)
end

# Most specific constructor for PDSFunction
function PDSFunction{iip, FullSpecialize}(P,
                                          D;
                                          p_prototype = nothing,
                                          d_prototype = nothing,
                                          analytic = nothing,
                                          std_rhs = nothing,
                                          linear_invariants = nothing,) where {iip}
    if std_rhs === nothing
        std_rhs = PDSStdRHS(P, D, p_prototype, d_prototype)
    end
    PDSFunction{iip,
                FullSpecialize,
                typeof(P),
                typeof(D),
                typeof(p_prototype),
                typeof(d_prototype),
                typeof(std_rhs),
                typeof(analytic),
                typeof(linear_invariants)}(P, D, p_prototype, d_prototype, std_rhs,
                                           analytic, linear_invariants)
end

# Evaluation of a PDSFunction
function (PD::PDSFunction)(u, p, t)
    return PD.std_rhs(u, p, t)
end

function (PD::PDSFunction)(du, u, p, t)
    return PD.std_rhs(du, u, p, t)
end

# Default implementation of the standard right-hand side evaluation function
struct PDSStdRHS{P, D, CacheP, CacheD, TMP, TMP2} <: Function
    p::P
    d::D
    p_cache::CacheP
    d_cache::CacheD
    tmp::TMP
    tmp2::TMP2
end

function PDSStdRHS(P, D, p_prototype, d_prototype)
    p_cache = isnothing(p_prototype) ? nothing : DiffCache(p_prototype)
    d_cache = isnothing(d_prototype) ? nothing : DiffCache(d_prototype)

    if p_prototype isa AbstractSparseMatrix
        tmp_vec = zeros(eltype(p_prototype), size(p_prototype, 1))
        tmp_cache = DiffCache(tmp_vec)
        tmp2_cache = DiffCache(tmp_vec / oneunit(first(tmp_vec)))
    else
        tmp_cache = nothing
        tmp2_cache = nothing
    end

    return PDSStdRHS(P, D, p_cache, d_cache, tmp_cache,
                     tmp2_cache)
end

# Evaluation of a PDSStdRHS (out-of-place)
function (PD::PDSStdRHS)(u, p, t)
    P = PD.p(u, p, t)
    D = PD.d(u, p, t)
    diag(P) + vec(sum(P; dims = 2)) - vec(sum(P; dims = 1)) - vec(D)
end

# Evaluation of a PDSStdRHS (in-place)
function (PD::PDSStdRHS)(du, u, p, t)
    P_matrix = get_tmp(PD.p_cache, du)
    D_vector = get_tmp(PD.d_cache, du)
    tmp = PD.tmp === nothing ? nothing : get_tmp(PD.tmp, du)
    tmp2 = PD.tmp2 === nothing ? nothing :
           get_tmp(PD.tmp2, zero(eltype(u)) / oneunit(eltype(u)))

    PD.p(P_matrix, u, p, t)

    if P_matrix isa AbstractSparseMatrix
        # row sum coded as matrix-vector product 
        fill!(tmp2, one(eltype(tmp2)))
        mul!(vec(du), P_matrix, tmp2)

        for i in 1:length(u)
            du[i] += P_matrix[i, i]
        end
        sum!(tmp', P_matrix)
        vec(du) .-= tmp
        PD.d(D_vector, u, p, t)
        vec(du) .-= D_vector
    else
        PD.d(D_vector, u, p, t)
        for i in 1:length(u)
            du[i] = P_matrix[i, i] - D_vector[i]
            for j in 1:length(u)
                du[i] += P_matrix[i, j] - P_matrix[j, i]
            end
        end
    end
    return nothing
end

# New problem type ConservativePDSProblem
"""
    ConservativePDSProblem(P, u0, tspan, p = NullParameters();
                           p_prototype = nothing,
                           analytic = nothing,
                           std_rhs = nothing,
                           linear_invariants = nothing)

A structure describing a conservative system of ordinary differential equation in form of a production-destruction system (PDS).
`P` denotes the function defining the production matrix ``P``.
The diagonal of ``P`` contains production terms without destruction counterparts.
`u0` is the vector of initial conditions and `tspan` the time span
`(t_initial, t_final)` of the problem. The optional argument `p` can be used
to pass additional parameters to the function `P`.

The function `P` can be given either in the out-of-place form with signature
`production_terms = P(u, p, t)` or the in-place form `P(production_terms, u, p, t)`.

### Keyword arguments: ###

- `p_prototype`: If `P` is given in in-place form, `p_prototype` or copies thereof are used to store evaluations of `P`.
  If `p_prototype` is not specified explicitly and `P` is in-place, then `p_prototype` will be internally
  set to `zeros(eltype(u0), (length(u0), length(u0)))`.
- `analytic`: The analytic solution of a PDS must be given in the form `f(u0,p,t)`.
  Specifying the analytic solution can be useful for plotting and convergence tests.
- `std_rhs`: The standard ODE right-hand side evaluation function callable
  as `du = std_rhs(u, p, t)` for the out-of-place form and
  as `std_rhs(du, u, p, t)` for the in-place form. Solvers that do not rely on
  the production-destruction representation of the ODE, will use this function
  instead to compute the solution. 

  `std_rhs` can also be passed as a [`SciMLBase.ODEFunction`](https://docs.sciml.ai/DiffEqDocs/stable/types/ode_types/#SciMLBase.ODEFunction). 
  This allows attaching additional information—such as a Jacobian prototype (`jac_prototype`) or sparsity pattern (`sparsity`)—for improving or accelerating the usage of `std_rhs`. 
  See the [`SciMLBase.ODEFunction` documentation](https://docs.sciml.ai/DiffEqDocs/stable/types/ode_types/#SciMLBase.ODEFunction) for details.

  If not specified, a default implementation calling `P` is used
- `linear_invariants`: The rows of this matrix contain the linear invariants of the ODE. 
  Certain solvers or callbacks might require this matrix.
  Note that this feature is experimental and its API may change in future releases.

## References

- Hans Burchard, Eric Deleersnijder, and Andreas Meister.
  "A high-order conservative Patankar-type discretisation for stiff systems of
  production-destruction equations."
  Applied Numerical Mathematics 47.1 (2003): 1-30.
  [DOI: 10.1016/S0168-9274(03)00101-6](https://doi.org/10.1016/S0168-9274(03)00101-6)
"""
struct ConservativePDSProblem{iip} <: AbstractPDSProblem end

# New ODE function ConservativePDSFunction
struct ConservativePDSFunction{iip, specialize, P, PrototypeP, StdRHS, Ta, LI} <:
       AbstractODEFunction{iip}
    p::P # production terms
    p_prototype::PrototypeP # prototype for production terms
    std_rhs::StdRHS # standard right-hand side evaluation function
    analytic::Ta # analytic solution (or nothing)
    linear_invariants::LI
end

# define behavior of ConservativePDSFunction for non-existing fields
function Base.getproperty(obj::ConservativePDSFunction, sym::Symbol)
    if sym === :mass_matrix
        return I
    elseif sym in fieldnames(typeof(obj))
        return getfield(obj, sym)
    elseif obj.std_rhs isa SciMLBase.AbstractODEFunction
        val = getproperty(obj.std_rhs, sym)
        if sym === :sparsity && val === nothing
            return getproperty(obj.std_rhs, :jac_prototype)
        end
        return val
    elseif sym === :jac_prototype || sym === :sparsity || sym === :sys
        return nothing
    else
        return getfield(obj, sym)
    end
end

# `remake(::ODEProblem)` calls `remake(prob.f; f, initialization_data)`.
# The PDS function types carry no initialization data and cannot be rebuilt
# from the keyword arguments of an `ODEFunction`, so only a replacement by
# another PDS function of the same type and with the same `iip` is supported.
function SciMLBase.remake(PD::Union{PDSFunction, ConservativePDSFunction};
                          f = missing, initialization_data = nothing)
    name = nameof(typeof(PD))
    if initialization_data !== nothing
        throw(ArgumentError("$name does not support initialization data."))
    end
    f === missing && return PD
    if !(f isa SciMLBase.parameterless_type(PD)) || isinplace(f) != isinplace(PD)
        throw(ArgumentError("The function of a problem with a $name can only be replaced by a $name with the same `iip` parameter."))
    end
    return f
end

# The default prototypes of a PDS function and the caches of its default
# `std_rhs` are allocated with the element type of `u0 / oneunit(first(tspan))`.
# If `remake` changes this element type without replacing `f`, the PDS function
# is rebuilt with prototypes of the new element type. Passing `f = prob.f` does
# not count as a replacement; `solve` does this when it updates `u0`, `p`, or
# `tspan`.
function SciMLBase.remake(prob::ODEProblem{<:Any, <:Any, iip, <:Any, <:Any, <:Any,
                                           <:AbstractPDSProblem};
                          f = missing, kwargs...) where {iip}
    new_prob = invoke(SciMLBase.remake, Tuple{ODEProblem}, prob; f, kwargs...)
    f === missing || f === prob.f || return new_prob
    new_f = remake_pds_function(prob.f, prob.u0, prob.tspan, new_prob.u0, new_prob.tspan)
    new_f === new_prob.f && return new_prob
    return ODEProblem{iip}(new_f, new_prob.u0, new_prob.tspan, new_prob.p,
                           new_prob.problem_type; new_prob.kwargs...)
end

# Return `PD` if it can be used with `u0` and `tspan`, otherwise a copy of `PD`
# with prototypes of the new element type.
function remake_pds_function(PD, u0_old, tspan_old, u0, tspan)
    if size(u0) != size(u0_old) && has_size_dependent_data(PD)
        throw(ArgumentError("`remake` cannot change the size of `u0` for a problem whose $(nameof(typeof(PD))) stores prototypes or linear invariants; construct a new problem instead."))
    end
    T_old = typeof(zero(eltype(u0_old)) / oneunit(first(tspan_old)))
    T = typeof(zero(eltype(u0)) / oneunit(first(tspan)))
    T == T_old && return PD
    return remake_prototypes(PD, T)
end

function has_size_dependent_data(PD::PDSFunction)
    return PD.p_prototype !== nothing || PD.d_prototype !== nothing ||
           PD.linear_invariants !== nothing
end

function has_size_dependent_data(PD::ConservativePDSFunction)
    return PD.p_prototype !== nothing || PD.linear_invariants !== nothing
end

function remake_prototypes(PD::PDSFunction{iip}, T) where {iip}
    p_prototype = similar_prototype(PD.p_prototype, T)
    d_prototype = similar_prototype(PD.d_prototype, T)
    # A user-supplied `std_rhs` is kept; the default one is rebuilt from the
    # new prototypes by the constructor.
    std_rhs = PD.std_rhs isa PDSStdRHS ? nothing : PD.std_rhs
    return PDSFunction{iip, FullSpecialize}(PD.p, PD.d; p_prototype, d_prototype,
                                            PD.analytic, std_rhs, PD.linear_invariants)
end

function remake_prototypes(PD::ConservativePDSFunction{iip}, T) where {iip}
    p_prototype = similar_prototype(PD.p_prototype, T)
    std_rhs = PD.std_rhs isa ConservativePDSStdRHS ? nothing : PD.std_rhs
    return ConservativePDSFunction{iip, FullSpecialize}(PD.p; p_prototype, PD.analytic,
                                                        std_rhs, PD.linear_invariants)
end

similar_prototype(::Nothing, T) = nothing
# `fill!` keeps the sparsity pattern of sparse prototypes.
similar_prototype(prototype, T) = fill!(similar(prototype, T), zero(T))

# Most general constructor for ConservativePDSProblems
function ConservativePDSProblem(P, u0, tspan, p = NullParameters(); kwargs...)
    iip = isinplace(P, 4)
    return ConservativePDSProblem{iip}(P, u0, tspan, p; kwargs...)
end

# Specialized constructor for ConservativePDSProblems setting `iip` manually
# (arbitrary function)  
function ConservativePDSProblem{iip}(P,
                                     u0,
                                     tspan,
                                     p = NullParameters();
                                     p_prototype = nothing,
                                     analytic = nothing,
                                     std_rhs = nothing,
                                     linear_invariants = nothing,
                                     kwargs...,) where {iip}

    # p_prototype is used to store evaluations of P, if P is in-place.
    if isnothing(p_prototype) && iip
        p_prototype = zeros(eltype(u0), (length(u0), length(u0))) / oneunit(first(tspan))
    end

    PD = ConservativePDSFunction{iip}(P; p_prototype, analytic, std_rhs, linear_invariants)
    ConservativePDSProblem{iip}(PD, u0, tspan, p; kwargs...)
end

# Specialized constructor for ConservativePDSProblems setting `iip` manually
# (ConservativePDSFunction)
function ConservativePDSProblem{iip}(PD::ConservativePDSFunction, u0, tspan,
                                     p = NullParameters(); kwargs...) where {iip}
    ODEProblem(PD, u0, tspan, p, ConservativePDSProblem{iip}(); kwargs...)
end

# Specialized constructor for ConservativePDSFunction setting `iip` manually
function ConservativePDSFunction{iip}(P; kwargs...) where {iip}
    ConservativePDSFunction{iip, FullSpecialize}(P; kwargs...)
end

# Most specific constructor for ConservativePDSFunction
function ConservativePDSFunction{iip, FullSpecialize}(P; p_prototype = nothing,
                                                      analytic = nothing, std_rhs = nothing,
                                                      linear_invariants = nothing) where {iip}
    if std_rhs === nothing
        std_rhs = ConservativePDSStdRHS(P, p_prototype)
    end
    ConservativePDSFunction{iip,
                            FullSpecialize,
                            typeof(P),
                            typeof(p_prototype),
                            typeof(std_rhs),
                            typeof(analytic),
                            typeof(linear_invariants)}(P, p_prototype, std_rhs, analytic,
                                                       linear_invariants)
end

# Evaluation of a ConservativePDSFunction
function (PD::ConservativePDSFunction)(u, p, t)
    return PD.std_rhs(u, p, t)
end

function (PD::ConservativePDSFunction)(du, u, p, t)
    return PD.std_rhs(du, u, p, t)
end

# Default implementation of the standard right-hand side evaluation function
struct ConservativePDSStdRHS{P, CacheP, TMP, TMP2} <: Function
    p::P
    p_cache::CacheP
    tmp::TMP
    tmp2::TMP2
end

function ConservativePDSStdRHS(P, p_prototype)
    p_cache = isnothing(p_prototype) ? nothing : DiffCache(p_prototype)

    if p_prototype isa AbstractSparseMatrix
        tmp_vec = zeros(eltype(p_prototype), size(p_prototype, 1))
        tmp_cache = DiffCache(tmp_vec)
        tmp2_cache = DiffCache(tmp_vec / oneunit(first(tmp_vec))) # drop units
    else
        tmp_cache = nothing
        tmp2_cache = nothing
    end

    return ConservativePDSStdRHS(P, p_cache, tmp_cache, tmp2_cache)
end

# Evaluation of a ConservativePDSStdRHS (out-of-place)
function (PD::ConservativePDSStdRHS)(u, p, t)
    #vec(sum(PD.p(u, p, t), dims = 2)) - vec(sum(PD.p(u, p, t), dims = 1))
    P = PD.p(u, p, t)

    f = zeros(eltype(P), size(u))
    @fastmath @inbounds @simd for I in CartesianIndices(P)
        if !iszero(P[I])
            f[I[1]] += P[I]
            f[I[2]] -= P[I]
        end
    end
    return f
end

function (PD::ConservativePDSStdRHS)(u::SVector, p, t)
    P = PD.p(u, p, t)

    f = similar(P[:, 1]) #constructs MVector
    zeroT = zero(eltype(P))
    for i in eachindex(f)
        f[i] = zeroT
    end

    @fastmath @inbounds @simd for I in CartesianIndices(P)
        if !iszero(P[I])
            f[I[1]] += P[I]
            f[I[2]] -= P[I]
        end
    end

    return SVector(f)
end

# Evaluation of a ConservativePDSStdRHS (in-place)
function (PD::ConservativePDSStdRHS)(du, u, p, t)
    P_matrix = get_tmp(PD.p_cache, du)
    tmp = PD.tmp === nothing ? nothing : get_tmp(PD.tmp, du)

    tmp2 = PD.tmp2 === nothing ? nothing :
           get_tmp(PD.tmp2, zero(eltype(u)) / oneunit(eltype(u)))

    # Populate the matrix and sum up terms
    PD.p(P_matrix, u, p, t)
    sum_terms!(du, tmp, tmp2, P_matrix)
    return nothing
end

# Generic fallback (for dense arrays)
# This implementation does not need any auxiliary vectors
@inline function sum_terms!(du, tmp, tmp2, P)
    for i in eachindex(du)
        du[i] = zero(eltype(du))
        for j in eachindex(du)
            du[i] += P[i, j] - P[j, i]
        end
    end
    return nothing
end

# Same result but more efficient - at least currently for SparseMatrixCSC
@inline function sum_terms!(du, tmp, tmp2, P::AbstractSparseMatrix)
    # row sum coded as matrix vector product
    fill!(tmp2, one(eltype(tmp2)))
    mul!(vec(du), P, tmp2)
    #sum!(vec(du), P)

    #column sum
    sum!(tmp', P)

    vec(du) .-= tmp
    return nothing
end

@inline function sum_terms!(du, tmp, tmp2, P::Tridiagonal)
    Base.require_one_based_indexing(du, P.dl, P.du)
    @assert length(du) == length(P.dl) + 1 == length(P.du) + 1

    let i = 1
        Pij = P.du[i]
        Pji = P.dl[i]
        du[i] = Pij - Pji
    end
    for i in 2:(length(du) - 1)
        Pij = P.dl[i - 1] + P.du[i]
        Pji = P.du[i - 1] + P.dl[i]
        du[i] = Pij - Pji
    end
    let i = lastindex(du)
        Pij = P.dl[i - 1]
        Pji = P.du[i - 1]
        du[i] = Pij - Pji
    end
    return nothing
end
