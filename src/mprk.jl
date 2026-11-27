#####################################################################
### Helper functions
#####################################################################

### add_small_constant ##############################################
# add_small_constant is used in out-of-place implementations to add
# the safety parameter small_constant to avoid divisions by zero.
add_small_constant(v, small_constant) = v .+ small_constant

function add_small_constant(v::SVector{N, T}, small_constant::T) where {N, T}
    v + SVector{N, T}(ntuple(i -> small_constant, N))
end

### p_prototype ######################################################
# p_prototype is used to allocate memory for production matrices
# in in-place implementations.
p_prototype(u, f) = zeros(eltype(u), length(u), length(u))

p_prototype(u, f::ConservativePDSFunction) = _p_prototype(f.p_prototype)
p_prototype(u, f::PDSFunction) = _p_prototype(f.p_prototype)
_p_prototype(prototype) = zero(prototype)

function _p_prototype(prototype::AbstractSparseMatrix)
    # We need to ensure that we store all structural nonzeros that
    # are required for the linear system. In particular, we need to
    # store the diagonal entries.
    p = copy(prototype)
    fill!(nonzeros(p), one(eltype(p)))
    return p + spdiagm(0 => ones(eltype(p), size(p, 1)))
end

### evaluate_pds ######################################################
# evaluate_pds is used in out-of-place implementations to evaluate
# the production and destruction terms
@inline function evaluate_pds(f::PDSFunction, u, p, t)
    P = f.p(u, p, t)
    d = f.d(u, p, t)
    return P, d
end

@inline function evaluate_pds(f::ConservativePDSFunction, u, p, t)
    P = f.p(u, p, t)
    return P, nothing
end

@inline function evaluate_pds!(P, d, f::PDSFunction, u, p, t)
    f.p(P, u, p, t)
    f.d(d, u, p, t)
    return nothing
end

@inline function evaluate_pds!(P, d::Nothing, f::ConservativePDSFunction, u, p, t)
    f.p(P, u, p, t)
    return nothing
end

### lincomb and lincomb! ###############################################
# These functions are used to compute linear combinations of production
# matrices and/or destruction vectors.

# out-of-place versions
# --- Base cases (End of recursion) ---
# For PDSFunctions (with destruction vectors)
@inline lincomb(c1::Number, P1, d1::AbstractArray) = (c1 .* P1, c1 .* d1)

# For ConservativePDSFunctions (without destruction vectors)
@inline lincomb(c1::Number, P1, d1::Nothing) = (c1 .* P1, nothing)

# --- Recursive steps ---
# For PDSFunctions: Processes triplets (coeff, P, d)
@muladd @inline function lincomb(c1::Number, P1, d1::AbstractArray, c2, P2, d2, args...)
    P_tail, d_tail = lincomb(c2, P2, d2, args...)
    return (c1 .* P1 .+ P_tail, c1 .* d1 .+ d_tail)
end

# For ConservativePDSFunctions: Processes triplets (coeff, P, nothing)
@muladd @inline function lincomb(c1::Number, P1, d1::Nothing, c2, P2, d2, args...)
    P_tail, _ = lincomb(c2, P2, d2, args...)
    return (c1 .* P1 .+ P_tail, nothing)
end

# in-place version for general arrays
@muladd @generated function lincomb!(dest::AbstractArray, pairs...)
    n = length(pairs) ÷ 2
    expr = :(pairs[1] * pairs[2])
    for i in 2:n
        expr = :($expr + pairs[$(2 * i - 1)] * pairs[$(2 * i)])
    end
    return quote
        @.. broadcast=false dest=$expr
        return dest
    end
end

# in-place version for destruction terms of conservative PDS
@inline function lincomb!(dest::Nothing, pairs...)
    return nothing
end

# in-place version for sparse matrices
@muladd @generated function lincomb!(dest::SparseMatrixCSC, pairs...)
    n = length(pairs) ÷ 2
    nz_names = [Symbol("nz_", i) for i in 1:n]
    setup_block = Expr(:block)
    # We need to keep the structural nonzeros of the production terms.
    # However, this is not guaranteed by broadcasting, see
    # https://github.com/JuliaSparse/SparseArrays.jl/issues/190
    for i in 1:n
        push!(setup_block.args, :($(nz_names[i]) = nonzeros(pairs[$(2 * i)])))
    end
    expr = :(pairs[1] * $(nz_names[1]))
    for i in 2:n
        expr = :($expr + pairs[$(2 * i - 1)] * $(nz_names[i]))
    end
    return quote
        $setup_block
        nz_dest = nonzeros(dest)
        @.. broadcast=false nz_dest=$expr
        return dest
    end
end

### Patankar linear systems ###################################################################
# A Patankar step solves a linear system for the new values x. Its classical form is
#   (I + B diag(1 ./ σ)) x = b,
# where σ contains the Patankar weight denominators and B holds dt times the production
# and destruction terms (multiplied by the coefficients of the scheme), i.e.,
# B[i, j] = -dt * P[i, j] for i ≠ j and B[j, j] = dt * (d[j] + sum_{i ≠ j} P[i, j]) for
# an MPRK stage. All schemes solve the equivalent system for the Patankar weights
# w = x ./ σ instead,
#   M w = b,  M = diag(σ) + B,  x = σ .* w.
# M is assembled without divisions by σ, so tiny denominators do not create huge or
# infinite entries. For nonnegative production and destruction terms and σ > 0, M has a
# positive diagonal, nonpositive off-diagonal entries, and column sums ≥ σ[j] > 0. Thus,
# M is a nonsingular M-matrix and w ≥ 0 for b ≥ 0. The small constant added to the
# denominators by the schemes keeps all columns of M nonzero.
#
# A denominator σ[j] = Inf (e.g., after an overflow) makes column j of the classical
# matrix the unit vector e_j. Column j of M is then set to oneunit(σ[j]) * e_j, and
# x[j] = w[j] (see `patankar_weight`).
#
# The entries of single columns of M can be many orders of magnitude smaller than the
# others, e.g., σ[j] ≈ floatmin for a component that receives mass but has no outflow.
# Then the Patankar weight w[j] = x[j] / σ[j] can overflow although x[j] is moderate,
# and sparse LU factorizations (e.g., UMFPACK and KLU), which scale the rows and choose
# pivots by threshold criteria, can fail or return inaccurate solutions. Thus, columns of
# M are multiplied by powers of two s[j] before the solve, and x = (σ .* s) .* y with the
# solution y of the scaled system:
# - Sparse matrices: all columns, with 1 ≤ s[j] * max_i |M[i, j]| < 2. Sparse LU can
#   return inaccurate solutions with a successful return code also if no diagonal entry
#   is extreme but the column norms differ by many orders of magnitude, so scaling only
#   the columns with extreme diagonal entries is not sufficient.
# - Dense, `Tridiagonal`, and static matrices: only columns whose diagonal entry, which
#   is the maximal entry of the column of an M-matrix, lies outside [2^-e, 2^e] with
#   e = exponent(floatmax) ÷ 2; they are scaled such that 1 ≤ s[j] * M[j, j] < 2.
#   LU factorizations with partial pivoting choose the same pivots for the scaled and
#   the unscaled matrix, and the scaling changes their results only where the unscaled
#   solve overflows or underflows; the other components of x are bitwise the same.
#   The in-place builders check the diagonal entries while computing them.
# The scaling and the computation of x avoid divisions.
#
# A solve is successful if the linear solver reports success and x is finite. The
# finiteness of x is checked in the loop that computes x.

# Linear solvers that do not set a return code report `ReturnCode.Default`.
function linsolve_successful(linres)
    return SciMLBase.successful_retcode(linres.retcode) ||
           linres.retcode == SciMLBase.ReturnCode.Default
end

# LinearSolve's default solver for `Tridiagonal` matrices throws an exception if the LU
# factorization fails: on Julia 1.10, its QR fallback writes entries outside the band, and
# on newer Julia versions, it uses `ldiv!`, which throws a `SingularException`.
# `LUFactorization` reports the failure with a return code instead.
linsolve_algorithm(linsolve, M) = linsolve
linsolve_algorithm(::Nothing, M::Tridiagonal) = LUFactorization()

# The factor of w[j] in the solution x[j] of a Patankar system, see above
@inline patankar_weight(σ) = ifelse(isinf(σ), oneunit(σ), σ)

# The power of two s with 1 ≤ s * c < 2 for the maximum c ≥ 0 of the absolute values of
# the entries in a column, and one if c is zero or not finite. The exponent is bounded
# such that s is finite for subnormal c.
function column_scale(c::AbstractFloat)
    if iszero(c) || !isfinite(c)
        return one(c)
    end
    return ldexp(one(c), min(-exponent(c), exponent(floatmax(c))))
end
# For IEEE floating point numbers c in the binades [2^(1 - b), 2^b), with the exponent
# bias b, s is obtained directly from the exponent bits of c.
function column_scale(c::T) where {T <: Base.IEEEFloat}
    U = Base.uinttype(T)
    nbits = Base.significand_bits(T)
    e_biased = Int((reinterpret(U, c) & Base.exponent_mask(T)) >> nbits)
    e_max = Int(Base.exponent_mask(T) >> nbits) - 1
    if 1 <= e_biased < e_max
        return reinterpret(T, U(e_max - e_biased) << nbits)
    end
    # zero, subnormal, and not finite c, and c in the largest binade
    if iszero(c) || !isfinite(c)
        return one(c)
    end
    return ldexp(one(c), min(-exponent(c), exponent(floatmax(c))))
end
# number types without `exponent`, e.g., dual numbers
column_scale(c) = one(c)

# The bounds 2^-e and 2^e with e = exponent(floatmax) ÷ 2 of the diagonal entries of
# dense, `Tridiagonal`, and static matrices whose columns are not scaled; `nothing` for
# number types without `exponent`
@inline function extreme_bounds(::Type{T}) where {T <: AbstractFloat}
    e = exponent(floatmax(T)) ÷ 2
    return ldexp(one(T), -e), ldexp(one(T), e)
end
extreme_bounds(::Type) = nothing

# Whether the absolute value of the diagonal entry c lies outside the range of
# `extreme_bounds`. A NaN entry gives a NaN solution, which is rejected anyway.
@inline function is_extreme(c, (lo, hi))
    a = abs(c)
    return (a < lo) | (a > hi)
end
is_extreme(c, ::Nothing) = false

# The power of two of `column_scale` for a diagonal entry c outside the range of
# `extreme_bounds`, and one otherwise
function extreme_column_scale(c::AbstractFloat)
    return is_extreme(c, extreme_bounds(typeof(c))) ? column_scale(abs(c)) : one(c)
end
extreme_column_scale(c) = one(c)

# Whether a diagonal entry of the dense or `Tridiagonal` matrix M is extreme (see
# `is_extreme`). The in-place builders `build_mprk_matrix!` for these matrices return
# the same value.
function has_extreme_diagonal(M)
    bounds = extreme_bounds(eltype(M))
    extreme = false
    for j in axes(M, 2)
        extreme |= is_extreme(M[j, j], bounds)
    end
    return extreme
end
function has_extreme_diagonal(M::Tridiagonal)
    bounds = extreme_bounds(eltype(M))
    extreme = false
    for c in M.d
        extreme |= is_extreme(c, bounds)
    end
    return extreme
end

# Multiply the columns of the dense matrix M with diagonal entries outside the range of
# `extreme_bounds` by powers of two s[j] and store patankar_weight(σ[j]) * s[j] in x[j].
# x may alias σ.
function scale_columns!(x, M::AbstractMatrix, σ)
    Base.require_one_based_indexing(x, M, σ)
    length(x) == length(σ) == size(M, 1) == size(M, 2) ||
        throw(DimensionMismatch("$(size(M)) matrix, but $(length(σ)) weights and $(length(x)) outputs"))
    for j in axes(M, 2)
        s = extreme_column_scale(abs(M[j, j]))
        if !isone(s)
            for i in axes(M, 1)
                M[i, j] *= s
            end
        end
        x[j] = patankar_weight(σ[j]) * s
    end
    return nothing
end

function scale_columns!(x, M::Tridiagonal, σ)
    (; dl, d, du) = M
    Base.require_one_based_indexing(x, d, σ)
    length(x) == length(σ) == length(d) ||
        throw(DimensionMismatch("$(length(d)) columns, but $(length(σ)) weights and $(length(x)) outputs"))
    for j in eachindex(d)
        s = extreme_column_scale(abs(d[j]))
        if !isone(s)
            # column j contains du[j - 1], d[j], and dl[j]
            d[j] *= s
            if j > firstindex(d)
                du[j - 1] *= s
            end
            if j < lastindex(d)
                dl[j] *= s
            end
        end
        x[j] = patankar_weight(σ[j]) * s
    end
    return nothing
end

# Multiply column j of the sparse matrix M by the power of two s of `column_scale(c)`,
# where c is the maximal absolute entry of the column, and return
# patankar_weight(σj) * s.
@inline function scale_column!(M::AbstractSparseMatrixCSC, j, c, σj)
    s = column_scale(c)
    if !isone(s)
        vals = nonzeros(M)
        for idx in nzrange(M, j)
            vals[idx] *= s
        end
    end
    return patankar_weight(σj) * s
end

# Multiply the columns of the sparse matrix M by the powers of two s[j] of `column_scale`
# of their maximal absolute entries and store patankar_weight(σ[j]) * s[j] in x[j].
# x may alias σ.
function scale_columns!(x, M::AbstractSparseMatrixCSC, σ)
    Base.require_one_based_indexing(x, σ)
    length(x) == length(σ) == size(M, 2) ||
        throw(DimensionMismatch("$(size(M, 2)) columns, but $(length(σ)) weights and $(length(x)) outputs"))
    vals = nonzeros(M)
    for j in axes(M, 2)
        c = zero(abs(zero(eltype(M))))
        for idx in nzrange(M, j)
            c = max(c, abs(vals[idx]))
        end
        x[j] = scale_column!(M, j, c, σ[j])
    end
    return nothing
end

# Out-of-place versions returning the scaled matrix and the factors
# patankar_weight.(σ) .* s, or `nothing` if no column is scaled; M can be modified.
function scale_columns(M, σ)
    has_extreme_diagonal(M) || return M, nothing
    factors = similar(σ)
    scale_columns!(factors, M, σ)
    return M, factors
end

function scale_columns(M::AbstractSparseMatrixCSC, σ)
    factors = similar(σ)
    scale_columns!(factors, M, σ)
    return M, factors
end

@inline function scale_columns(M::StaticMatrix{N, N}, σ) where {N}
    bounds = extreme_bounds(eltype(M))
    d = diag(M)
    if any(is_extreme.(d, Ref(bounds)))
        s = extreme_column_scale.(d)
        return M .* s', patankar_weight.(σ) .* s
    end
    return M, nothing
end

# Solve the Patankar system M w = rhs and return x = patankar_weight.(σ) .* w and whether
# the solve was successful. M can be modified.
@inline function patankar_solve(M, rhs, σ, linsolve)
    M_scaled, factors = scale_columns(M, σ)
    sol = solve(LinearProblem(M_scaled, rhs), linsolve_algorithm(linsolve, M_scaled))
    if factors === nothing
        x = patankar_weight.(σ) .* sol.u
    else
        x = factors .* sol.u
    end
    return x, linsolve_successful(sol) && all(isfinite, x)
end

# Store x[i] = f(a[i]) * y[i] and return whether x is finite. x may alias a.
# Long vectors use a vectorized loop, in which the sum of x[i] - x[i] (zero for finite
# x[i] and NaN otherwise) replaces `isfinite`; its setup cost dominates for short
# vectors.
@inline function map_multiply_finite!(f, x, a, y)
    if length(x) < 64
        finite = true
        for i in eachindex(x, a, y)
            xi = f(a[i]) * y[i]
            x[i] = xi
            finite &= isfinite(xi)
        end
        return finite
    end
    acc = zero(eltype(x))
    @simd for i in eachindex(x, a, y)
        xi = f(a[i]) * y[i]
        x[i] = xi
        acc += xi - xi
    end
    return !isnan(acc)
end

# x .= factors .* y and x .= patankar_weight.(σ) .* y, returning whether x is finite
multiply_finite!(x, factors, y) = map_multiply_finite!(identity, x, factors, y)
weights_finite!(x, σ, y) = map_multiply_finite!(patankar_weight, x, σ, y)

# In-place version: solve the Patankar system M w = linsolve.b, where M is the matrix
# aliased by `linsolve`, store x = patankar_weight.(σ) .* w in x, and return whether the
# solve was successful. M is overwritten. x may alias σ but not `linsolve.b`. For dense
# and `Tridiagonal` matrices, `extreme` is the value of `has_extreme_diagonal(M)`.
@inline function patankar_solve!(x, M, σ, linsolve, extreme = has_extreme_diagonal(M))
    if extreme
        # x stores the factors patankar_weight(σ[j]) * s[j] until w is known
        scale_columns!(x, M, σ)
        return solve_scaled!(x, M, linsolve)
    end
    linsolve.A = M
    sol = solve!(linsolve)
    finite = weights_finite!(x, σ, sol.u)
    return linsolve_successful(sol) && finite
end

function patankar_solve!(x, M::AbstractSparseMatrixCSC, σ, linsolve)
    # x stores the factors patankar_weight(σ[j]) * s[j] until w is known
    scale_columns!(x, M, σ)
    return solve_scaled!(x, M, linsolve)
end

# Solve the column-scaled system with the matrix M and store x .* y in x, where x
# contains the factors patankar_weight(σ[j]) * s[j] and y is the solution. Return
# whether the solve was successful.
@inline function solve_scaled!(x, M, linsolve)
    linsolve.A = M
    sol = solve!(linsolve)
    finite = multiply_finite!(x, x, sol.u)
    return linsolve_successful(sol) && finite
end

### basic_patankar_step, basic_patankar_step! ##################################################
# These functions implement the core Patankar step by building and solving the linear system.
# They also return whether the linear solve succeeded; if it did not, the calling
# `perform_step!` rejects the step with `reject_step!`.

# out-of-place implementations
# Version for PDSFunctions (non-conservative PDS)
@muladd function basic_patankar_step(v, P, σ, dt, linsolve, d::AbstractVector, P2 = P)
    rhs = v + dt * diag(P2)
    M = build_mprk_matrix(P, σ, dt, d)
    return patankar_solve(M, rhs, σ, linsolve)
end

# Version for ConservativePDSFunctions (conservative PDS)
function basic_patankar_step(v, P, σ, dt, linsolve, d::Nothing, P2 = P)
    M = build_mprk_matrix(P, σ, dt)
    return patankar_solve(M, v, σ, linsolve)
end

# in-place implementations
# non-conservative PDS
@muladd @inline function add_diagonal!(b, P, dt)
    @inbounds for i in eachindex(b)
        b[i] = b[i] + dt * P[i, i]
    end
    return nothing
end

# in-place version for PDSProblems; u must not alias `linsolve.b`
@inline function basic_patankar_step!(u, v, P, d::AbstractArray, σ, dt, linsolve)
    b = linsolve.b
    b .= v

    add_diagonal!(b, P, dt)

    return build_and_solve!(u, P, σ, dt, d, linsolve)
end

# in-place version for ConservativePDSProblems; u must not alias `linsolve.b`
@inline function basic_patankar_step!(u, v, P, d::Nothing, σ, dt, linsolve)
    b = linsolve.b
    b .= v

    return build_and_solve!(u, P, σ, dt, d, linsolve)
end

# Build the Patankar system in place of P and solve it with `patankar_solve!`
@inline function build_and_solve!(u, P, σ, dt, d, linsolve)
    extreme = build_mprk_matrix!(P, P, σ, dt, d)
    return patankar_solve!(u, P, σ, linsolve, extreme)
end

# Sparse matrices are built and column-scaled in a single pass; u stores the factors
# patankar_weight(σ[j]) * s[j] until w is known.
function build_and_solve!(u, P::AbstractSparseMatrixCSC, σ, dt, d, linsolve)
    build_mprk_matrix!(P, P, σ, dt, d; factors = u)
    return solve_scaled!(u, P, linsolve)
end

### build_mprk_matrix ############################################################
# These functions build the system matrix M = diag(σ) + B of the weight form of the
# Patankar step, see above.

# out-of-place for dense matrices
function build_mprk_matrix(P, sigma, dt, d = nothing)
    # re-use the in-place version implemented below
    M = similar(P)

    build_mprk_matrix!(M, P, sigma, dt, d)

    return M
end

# out-of-place for sparse matrices
function build_mprk_matrix(P::AbstractSparseMatrix, sigma, dt, d = nothing)
    # M must store all diagonal entries, but P need not store them. For example,
    # broadcasting drops stored zeros of the production terms, see
    # https://github.com/JuliaSparse/SparseArrays.jl/issues/190
    M = _p_prototype(P)

    build_mprk_matrix!(M, P, sigma, dt, d)

    return M
end

# out-of-place for StaticMatrix
@inline function sum_offdiagonal_col(P::StaticMatrix{N, N, T}, j) where {N, T}
    res = zero(T)
    for i in 1:N
        if i != j
            res += P[i, j]
        end
    end
    return res
end
@muladd @inline function build_mprk_matrix(P::StaticMatrix{N, N, T}, sigma, dt,
                                           d = nothing) where {N, T}
    M = SMatrix{N, N, T}((i == j) ?
                         # diagonal: sum over P[k, i] where k != i
                         sigma[i] +
                         dt *
                         (sum_offdiagonal_col(P, i) + (d === nothing ? zero(T) : d[i])) :
                         # off-diagonal
                         -(dt * P[i, j])
                         for i in 1:N, j in 1:N)
    any(isinf, sigma) || return M
    return SMatrix{N, N, T}(isinf(sigma[j]) ? (i == j ? oneunit(T) : zero(T)) : M[i, j]
                            for i in 1:N, j in 1:N)
end

# in-place for dense arrays; returns whether a diagonal entry is extreme (`is_extreme`)
@muladd function build_mprk_matrix!(M, P, sigma, dt, d = nothing)
    # M[i,i] = sigma[i] + dt * (d[i] + sum_j≠i P[j,i])
    # M[i,j] = -dt * P[i,j]
    # TODO: the performance of this can likely be improved
    Base.require_one_based_indexing(M, P, sigma)
    @assert size(M, 1) == size(M, 2) == size(P, 1) == size(P, 2) == length(sigma)
    if !isnothing(d)
        Base.require_one_based_indexing(d)
        @assert length(sigma) == length(d)
    end

    zeroM = zero(eltype(P))

    # The diagonal accumulates dt * (d[i] + sum_j≠i P[j,i]) first
    if isnothing(d)
        @inbounds for i in eachindex(sigma)
            M[i, i] = zeroM
        end
    else
        # nonconservative destruction terms (PDSFunctions only!)
        @inbounds for i in eachindex(d)
            M[i, i] = dt * d[i]
        end
    end

    # Run through P and fill M accordingly.
    # If P[i,j] ≠ 0 set M[i,j] = -dt*P[i,j] and add dt*P[i,j] to M[j,j].
    @inbounds @simd for I in CartesianIndices(P)
        if I[1] != I[2]
            if !iszero(P[I])
                dtP = dt * P[I]
                M[I] = -dtP
                M[I[2], I[2]] += dtP
            else
                M[I] = zeroM
            end
        end
    end

    # Add the Patankar weight denominators
    bounds = extreme_bounds(eltype(M))
    extreme = false
    @inbounds @simd for i in eachindex(sigma)
        M[i, i] = sigma[i] + M[i, i]
        extreme |= is_extreme(M[i, i], bounds)
    end

    set_unit_columns!(M, sigma)

    return extreme
end

# in-place optimized version for Tridiagonal matrices; returns whether a diagonal entry
# is extreme (`is_extreme`)
@muladd function build_mprk_matrix!(M::Tridiagonal, P::Tridiagonal, σ, dt, d = nothing)
    # M[i,i] = σ[i] + dt * (d[i] + sum_j≠i P[j,i])
    # M[i,j] = -dt * P[i,j]
    Base.require_one_based_indexing(M.dl, M.d, M.du,
                                    P.dl, P.d, P.du, σ)
    @assert length(M.dl) + 1 == length(M.d) == length(M.du) + 1 ==
            length(P.dl) + 1 == length(P.d) == length(P.du) + 1 == length(σ)
    if !isnothing(d)
        Base.require_one_based_indexing(d)
        @assert length(σ) == length(d)
    end

    if isnothing(d)
        fill!(M.d, zero(eltype(M.d)))
    else
        for i in eachindex(M.d, d)
            M.d[i] = dt * d[i]
        end
    end

    for i in eachindex(M.dl, P.dl)
        dtP = dt * P.dl[i]
        M.dl[i] = -dtP
        M.d[i] += dtP
    end

    for i in eachindex(M.du, P.du)
        dtP = dt * P.du[i]
        M.du[i] = -dtP
        M.d[i + 1] += dtP
    end

    bounds = extreme_bounds(eltype(M))
    extreme = false
    for i in eachindex(M.d, σ)
        M.d[i] = σ[i] + M.d[i]
        extreme |= is_extreme(M.d[i], bounds)
    end

    set_unit_columns!(M, σ)

    return extreme
end

# in-place optimized version for sparse matrices in compressed sparse column format;
# other sparse matrices use the generic version above
# With a vector `factors`, the columns of M are also scaled as in `scale_columns!`, and
# `factors` stores patankar_weight(σ[j]) * s[j]; `factors` may alias σ.
@muladd function build_mprk_matrix!(M::AbstractSparseMatrixCSC, P::AbstractSparseMatrixCSC,
                                    σ, dt, d = nothing; factors = nothing)
    # M[i,i] = σ[i] + dt * (d[i] + sum_j≠i P[j,i])
    # M[i,j] = -dt * P[i,j]
    # M must store all diagonal entries. P need not store them unless M === P,
    # as in the in-place Patankar steps.
    Base.require_one_based_indexing(M, P, σ)
    @assert size(M, 1) == size(M, 2) == size(P, 1) == size(P, 2) == length(σ)
    if !isnothing(d)
        Base.require_one_based_indexing(d)
        @assert length(σ) == length(d)
    end

    if M === P
        # M overwrites P in place. Column j of M depends only on column j of P,
        # so a single pass over the stored entries suffices.
        vals = nonzeros(M)
        rows = rowvals(M)
        if factors !== nothing
            Base.require_one_based_indexing(factors)
            length(factors) == length(σ) ||
                throw(DimensionMismatch("$(length(σ)) weights, but $(length(factors)) factors"))
        end
        for j in axes(M, 2)
            val = zero(eltype(M))
            c = zero(abs(zero(eltype(M))))
            idx_diag = 0
            for idx in nzrange(M, j)
                if rows[idx] == j
                    idx_diag = idx
                else
                    val += vals[idx]
                    vals[idx] = -dt * vals[idx]
                    if factors !== nothing
                        c = max(c, abs(vals[idx]))
                    end
                end
            end
            if idx_diag == 0
                throw(ArgumentError("The MPRK matrix does not store the diagonal entry in column $j."))
            end
            if isnothing(d)
                vals[idx_diag] = σ[j] + dt * val
            else
                vals[idx_diag] = σ[j] + dt * (d[j] + val)
            end
            # Columns are completed one at a time, so that `factors` may alias σ.
            if isinf(σ[j])
                set_unit_column!(M, j, σ[j])
                c = oneunit(c)
            else
                c = max(c, abs(vals[idx_diag]))
            end
            if factors !== nothing
                factors[j] = scale_column!(M, j, c, σ[j])
            end
        end
        return nothing
    end

    M_rows = rowvals(M)
    M_vals = nonzeros(M)
    P_rows = rowvals(P)
    P_vals = nonzeros(P)
    n = size(M, 2)

    # Diagonal entries
    for j in 1:n
        for idx_M in nzrange(M, j)
            # First, we look for the index corresponding to the
            # diagonal entry of M
            i = M_rows[idx_M]
            if i != j
                continue
            end

            # Next, we sum over all production terms in the
            # i-th column of P except for the diagonal entry
            val = zero(eltype(P))
            for idx_P in nzrange(P, i)
                if P_rows[idx_P] != i
                    val += P_vals[idx_P]
                end
            end

            # Finally, we set the diagonal entry of M
            if isnothing(d)
                M_vals[idx_M] = σ[i] + dt * val
            else
                M_vals[idx_M] = σ[i] + dt * (d[i] + val)
            end
        end
    end

    # Fill remaining entries
    for j in 1:n
        iter_M = iterate(nzrange(M, j))
        iter_P = iterate(nzrange(P, j))

        while iter_M !== nothing
            idx_M = first(iter_M)
            i = M_rows[idx_M]

            if iter_P === nothing
                # If there are no more nonzeros in P, we set the
                # entry in M to zero unless we are on the diagonal,
                # where we do nothing since we have already set the
                # diagonal entries.
                if i != j
                    M_vals[idx_M] = zero(eltype(P))
                end
            else
                idx_P = first(iter_P)
                i_P = P_rows[idx_P]
                if i_P < i
                    throw(ArgumentError("The production matrix has more non-zero entries than the MPRK matrix."))
                elseif i_P == i
                    # P has a nonzero entry at this position
                    if i != j
                        M_vals[idx_M] = -dt * P_vals[idx_P]
                    end
                    iter_P = iterate(nzrange(P, j), last(iter_P))
                else
                    # P has a zero entry at this position
                    if i != j
                        M_vals[idx_M] = zero(eltype(P))
                    end
                end
            end

            iter_M = iterate(nzrange(M, j), last(iter_M))

            if (iter_M === nothing) && (iter_P !== nothing)
                throw(ArgumentError("The production matrix has more non-zero entries than the MPRK matrix."))
            end
        end
    end

    set_unit_columns!(M, σ)
    if factors !== nothing
        scale_columns!(factors, M, σ)
    end

    return nothing
end

# Set the columns j of M with σ[j] = Inf to oneunit(σ[j]) * e_j, see
# "Patankar linear systems" above.
function set_unit_columns!(M, σ)
    any(isinf, σ) || return nothing
    for j in eachindex(σ)
        if isinf(σ[j])
            for i in axes(M, 1)
                M[i, j] = i == j ? oneunit(σ[j]) : zero(eltype(M))
            end
        end
    end
    return nothing
end

function set_unit_columns!(M::Tridiagonal, σ)
    any(isinf, σ) || return nothing
    # column j contains du[j - 1], d[j], and dl[j]
    for j in eachindex(σ)
        if isinf(σ[j])
            M.d[j] = oneunit(σ[j])
            if j > firstindex(σ)
                M.du[j - 1] = zero(eltype(M.du))
            end
            if j < lastindex(σ)
                M.dl[j] = zero(eltype(M.dl))
            end
        end
    end
    return nothing
end

function set_unit_columns!(M::AbstractSparseMatrixCSC, σ)
    any(isinf, σ) || return nothing
    for j in eachindex(σ)
        if isinf(σ[j])
            set_unit_column!(M, j, σ[j])
        end
    end
    return nothing
end

# Set column j of the sparse matrix M to oneunit(σj) * e_j
@inline function set_unit_column!(M::AbstractSparseMatrixCSC, j, σj)
    rows = rowvals(M)
    vals = nonzeros(M)
    for idx in nzrange(M, j)
        vals[idx] = rows[idx] == j ? oneunit(σj) : zero(eltype(vals))
    end
    return nothing
end

# We use MPRKMutableCache as supertype for all MPRK scheme caches
abstract type MPRKMutableCache <: OrdinaryDiffEqMutableCache end
get_fsalfirstlast(cache::MPRKMutableCache, rate_prototype) = (nothing, nothing)

### Rejection of step attempts ##################################################################
# OrdinaryDiffEqCore retries a step rejected by `force_stepfail` with a smaller
# time step if `adaptive = true` and stops with `ReturnCode.ConvergenceFailure`
# otherwise. It does not count such a rejection in `stats.nreject`.
# After a rejection, `integrator.u` equals the last accepted value `uprev`, which is
# saved as the final value if the solve stops. In-place methods store stages in `u`,
# and out-of-place methods may have assigned `u` in an earlier attempt that the error
# test rejected.
# With `adaptive = true`, the solve stops with `ReturnCode.Unstable` after
# `MAX_INVALID_ATTEMPTS` consecutive rejected attempts. Each rejection divides the time
# step by `failfactor` (default 2).
const MAX_INVALID_ATTEMPTS = 10

# The counter of consecutive rejected step attempts of adaptive schemes; `nothing` for
# schemes without error estimate, whose attempts are not counted
function invalid_attempts(cache::OrdinaryDiffEqConstantCache)
    return hasfield(typeof(cache), :invalid_attempts) ? cache.invalid_attempts : nothing
end
invalid_attempts(cache::MPRKMutableCache) = invalid_attempts(cache.tab)

function reject_step!(integrator, cache)
    integrator.force_stepfail = true
    integrator.stats.nreject += 1
    if isinplace(integrator.f)
        integrator.u .= integrator.uprev
    else
        integrator.u = integrator.uprev
    end

    count = invalid_attempts(cache)
    if integrator.opts.adaptive && count !== nothing
        count[] += 1
        if count[] >= MAX_INVALID_ATTEMPTS
            count[] = 0
            report_instability(integrator)
            SciMLBase.terminate!(integrator, SciMLBase.ReturnCode.Unstable)
        end
    end
    return nothing
end

# Emit the message of SciMLBase's `check_error` for `ReturnCode.Unstable`, which
# respects the `verbose` setting of the integrator (SciMLBase v3.49 and newer).
function report_instability(integrator)
    if isdefined(SciMLBase, :report_integrator_failure)
        SciMLBase.report_integrator_failure(integrator, Val(:instability))
    end
    return nothing
end

# Store the error estimate `EEst = internalnorm(atmp, t)` of the current step attempt,
# computed from the scaled error components `atmp`, and reset the counter of consecutive
# rejected attempts. If a component of `atmp` is not finite, `EEst` is set to infinity:
# the default `internalnorm` uses fast-math operations, which allow the compiler to
# assume finite values. The step size controller rejects an attempt with `EEst = Inf`
# like any attempt with `EEst > 1` and reduces the time step by the maximal factor
# allowed by `qmin`; such attempts do not count towards `MAX_INVALID_ATTEMPTS`.
function set_error_estimate!(integrator, cache, atmp, t)
    EEst = integrator.opts.internalnorm(atmp, t)
    if !all(isfinite, atmp)
        EEst = oftype(EEst, Inf)
    end
    set_EEst!(integrator, EEst)
    invalid_attempts(cache)[] = 0
    return nothing
end

### MPE #####################################################################################
"""
    MPE([linsolve = ..., small_constant = ...])

The first-order modified Patankar-Euler algorithm for production-destruction systems. This one-step, one-stage method is
first-order accurate, unconditionally positivity-preserving, and
linearly implicit.

The scheme was introduced by Burchard et al. for conservative production-destruction systems.
For nonconservative production–destruction systems we use the straight forward extension

``u_i^{n+1} = u_i^n + Δt \\sum_{j, j≠i} \\biggl(p_{ij}^n \\frac{u_j^{n+1}}{u_j^n}-d_{ij}^n \\frac{u_i^{n+1}}{u_i^n}\\biggr) + {\\Delta}t p_{ii}^n - Δt d_{ii}^n\\frac{u_i^{n+1}}{u_i^n}``,

where ``p_{ij}^n = p_{ij}(t^n,\\mathbf u^n)`` and ``d_{ij}^n = d_{ij}(t^n,\\mathbf u^n)``.

The modified Patankar-Euler method requires the special structure of a [`PDSProblem`](@ref) or a [`ConservativePDSProblem`](@ref).

You can optionally choose the linear solver to be used by passing an
algorithm from [LinearSolve.jl](https://github.com/SciML/LinearSolve.jl)
as keyword argument `linsolve`.
You can also choose the parameter `small_constant` which is added to all Patankar-weight denominators
to avoid divisions by zero. You can pass a value explicitly, otherwise `small_constant` is set to
`floatmin` of the floating point type used.

The current implementation only supports fixed time steps.

## References

- Hans Burchard, Eric Deleersnijder, and Andreas Meister.
  "A high-order conservative Patankar-type discretisation for stiff systems of
  production-destruction equations."
  Applied Numerical Mathematics 47.1 (2003): 1-30.
  [DOI: 10.1016/S0168-9274(03)00101-6](https://doi.org/10.1016/S0168-9274(03)00101-6)
"""
struct MPE{F, T} <: OrdinaryDiffEqAlgorithm
    linsolve::F
    small_constant_function::T
end

function MPE(; linsolve = LUFactorization(), small_constant = nothing)
    if isnothing(small_constant)
        small_constant_function = floatmin
    elseif small_constant isa Number
        small_constant_function = Returns(small_constant)
    else # assume small_constant isa Function
        small_constant_function = small_constant
    end
    MPE(linsolve, small_constant_function)
end

alg_order(::MPE) = 1
isfsal(::MPE) = false

struct MPEConstantCache{T} <: OrdinaryDiffEqConstantCache
    small_constant::T
end

# Out-of-place
function alg_cache(alg::MPE, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
                   uprev, uprev2, f, t, dt, reltol, p, calck,
                   ::Val{false},
                   verbose) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    if !(f isa PDSFunction || f isa ConservativePDSFunction)
        throw(ArgumentError("MPE can only be applied to production-destruction systems"))
    end
    MPEConstantCache(alg.small_constant_function(uEltypeNoUnits))
end

function initialize!(integrator, cache::MPEConstantCache)
end

@muladd function perform_step!(integrator, cache::MPEConstantCache, repeat_step = false)
    (; alg, t, dt, uprev, f, p) = integrator
    (; small_constant) = cache

    # evaluate production matrix and destruction vector
    P, d = evaluate_pds(f, uprev, p, t)
    integrator.stats.nf += 1

    # avoid division by zero due to zero Patankar weights
    σ = add_small_constant(uprev, small_constant)

    u, solved = basic_patankar_step(uprev, P, σ, dt, alg.linsolve, d)
    integrator.stats.nsolve += 1
    solved || return reject_step!(integrator, cache)

    integrator.u = u
end

struct MPECache{PType, dType, uType, tabType, F} <: MPRKMutableCache
    P::PType
    D::dType
    σ::uType
    tab::tabType
    linsolve::F
end

get_tmp_cache(integrator, ::MPE, cache::OrdinaryDiffEqMutableCache) = (cache.σ,)

# In-place
function alg_cache(alg::MPE, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
                   uprev, uprev2, f, t, dt, reltol, p, calck,
                   ::Val{true},
                   verbose) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    P = p_prototype(u, f)
    σ = zero(u)
    tab = MPEConstantCache(alg.small_constant_function(uEltypeNoUnits))

    if f isa ConservativePDSFunction
        # We use P to store the evaluation of the PDS
        # as well as to store the system matrix of the linear system
        # Right hand side of linear system is always uprev
        linprob = LinearProblem(P, _vec(uprev))
        linsolve = init(linprob, linsolve_algorithm(alg.linsolve, linprob.A),
                        alias = LinearSolve.LinearAliasSpecifier(; alias_A = true,
                                                                 alias_b = true),
                        assumptions = LinearSolve.OperatorAssumptions(true))

        MPECache(P, nothing, σ, tab, linsolve)
    elseif f isa PDSFunction
        linsolve_rhs = zero(u)
        # We use P to store the evaluation of the PDS
        # as well as to store the system matrix of the linear system
        linprob = LinearProblem(P, _vec(linsolve_rhs))
        linsolve = init(linprob, linsolve_algorithm(alg.linsolve, linprob.A),
                        alias = LinearSolve.LinearAliasSpecifier(; alias_A = true,
                                                                 alias_b = true),
                        assumptions = LinearSolve.OperatorAssumptions(true))

        MPECache(P, similar(u), σ, tab, linsolve)
    else
        throw(ArgumentError("MPE can only be applied to production-destruction systems"))
    end
end

function initialize!(integrator, cache::MPECache)
end

@muladd function perform_step!(integrator, cache::MPECache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; P, D, σ, linsolve) = cache
    (; small_constant) = cache.tab

    # We use P to store the evaluation of the PDS
    # as well as to store the system matrix of the linear system

    # We require the users to set unused entries to zero!

    # evaluate production/destruction terms
    evaluate_pds!(P, D, f, uprev, p, t)
    integrator.stats.nf += 1

    # avoid division by zero due to zero Patankar weights
    @.. broadcast=false σ=uprev + small_constant

    solved = basic_patankar_step!(u, uprev, P, D, σ, dt, linsolve)
    integrator.stats.nsolve += 1
    solved || return reject_step!(integrator, cache)
end

### MPRK22 #####################################################################################
"""
    MPRK22(α; [linsolve = ..., small_constant = ...])

A family of second-order modified Patankar-Runge-Kutta algorithms for
production-destruction systems. Each member of this family is an adaptive, one-step, two-stage method which is
second-order accurate, unconditionally positivity-preserving, and linearly
implicit. In this implementation the stage-values are conservative as well.
The parameter `α` is described by Kopecz and Meister (2018) and
studied by Izgin, Kopecz and Meister (2022) as well as
Torlo, Öffner and Ranocha (2022).

This method supports adaptive time stepping, using the Patankar-weight denominators
``σ_i``, see Kopecz and Meister (2018), as first order approximations to estimate the error.

The scheme was introduced by Kopecz and Meister for conservative production-destruction systems.
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

- Hans Burchard, Eric Deleersnijder, and Andreas Meister.
  "A high-order conservative Patankar-type discretisation for stiff systems of
  production-destruction equations."
  Applied Numerical Mathematics 47.1 (2003): 1-30.
  [DOI: 10.1016/S0168-9274(03)00101-6](https://doi.org/10.1016/S0168-9274(03)00101-6)
- Stefan Kopecz and Andreas Meister.
  "On order conditions for modified Patankar-Runge-Kutta schemes."
  Applied Numerical Mathematics 123 (2018): 159-179.
  [DOI: 10.1016/j.apnum.2017.09.004](https://doi.org/10.1016/j.apnum.2017.09.004)
- Thomas Izgin, Stefan Kopecz, and Andreas Meister.
  "On Lyapunov stability of positive and conservative time integrators and application
  to second order modified Patankar-Runge-Kutta schemes."
  ESAIM: Mathematical Modelling and Numerical Analysis 56.3 (2022): 1053-1080.
  [DOI: 10.1051/m2an/2022031](https://doi.org/10.1051/m2an/2022031)
- Davide Torlo, Philipp Öffner, and Hendrik Ranocha.
  "Issues with positivity-preserving Patankar-type schemes."
  Applied Numerical Mathematics 182 (2022): 117-147.
  [DOI: 10.1016/j.apnum.2022.07.014](https://doi.org/10.1016/j.apnum.2022.07.014)
"""
struct MPRK22{T, F, T2} <: OrdinaryDiffEqAdaptiveAlgorithm
    alpha::T
    linsolve::F
    small_constant_function::T2
end

function MPRK22(alpha; linsolve = LUFactorization(),
                small_constant = nothing)
    if isnothing(small_constant)
        small_constant_function = floatmin
    elseif small_constant isa Number
        small_constant_function = Returns(small_constant)
    else # assume small_constant isa Function
        small_constant_function = small_constant
    end
    alg = MPRK22{typeof(alpha), typeof(linsolve), typeof(small_constant_function)}(alpha,
                                                                                   linsolve,
                                                                                   small_constant_function)
    check_parameters(alg)
    return alg
end

alg_order(::MPRK22) = 2
isfsal(::MPRK22) = false

# The parameters of the methods are checked in their own type, independent of the
# type of the problem. Limits such as 1/3 are computed in this type as well, so that
# the values a user writes for them, e.g., `1 / 3`, are on the boundary.
function check_parameters(alg::MPRK22)
    (; alpha) = alg
    if !(alpha ≥ one(alpha) / 2)
        throw(ArgumentError("MPRK22 requires α ≥ 1/2."))
    end
    return nothing
end

function get_constant_parameters(alg::MPRK22, type)
    check_parameters(alg)

    a21 = convert(type, alg.alpha)
    b2 = 1 / (2 * a21)
    b1 = 1 - b2

    # This should never happen
    if !all((a21, b1, b2) .≥ 0)
        throw(ArgumentError("MPRK22 requires nonnegative RK coefficients."))
    end
    return a21, b1, b2
end

struct MPRK22ConstantCache{T} <: OrdinaryDiffEqConstantCache
    a21::T
    b1::T
    b2::T
    small_constant::T
    # number of consecutive rejected step attempts
    invalid_attempts::Base.RefValue{Int}
end

# Out-of-place
function alg_cache(alg::MPRK22, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
                   uprev, uprev2, f, t, dt, reltol, p, calck,
                   ::Val{false},
                   verbose) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    if !(f isa PDSFunction || f isa ConservativePDSFunction)
        throw(ArgumentError("MPRK22 can only be applied to production-destruction systems"))
    end

    a21, b1, b2 = get_constant_parameters(alg, uEltypeNoUnits)
    MPRK22ConstantCache(a21, b1, b2, alg.small_constant_function(uEltypeNoUnits), Ref(0))
end

function initialize!(integrator, cache::MPRK22ConstantCache)
    invalid_attempts(cache)[] = 0
end

@muladd function perform_step!(integrator, cache::MPRK22ConstantCache, repeat_step = false)
    (; alg, t, dt, uprev, f, p) = integrator
    (; a21, b1, b2, small_constant) = cache

    # evaluate production matrix
    P, d = evaluate_pds(f, uprev, p, t)
    integrator.stats.nf += 1

    Ptmp, dtmp = lincomb(a21, P, d)

    # avoid division by zero due to zero Patankar weights
    σ = add_small_constant(uprev, small_constant)

    u, solved = basic_patankar_step(uprev, Ptmp, σ, dt, alg.linsolve, dtmp)
    integrator.stats.nsolve += 1
    solved || return reject_step!(integrator, cache)

    # compute Patankar weight denominator
    if isone(a21)
        σ = u
    else
        # σ = σ .* (u ./ σ) .^ (1 / a21) # generated Infs when solving brusselator
        σ = σ .^ (1 - 1 / a21) .* u .^ (1 / a21)
    end
    # avoid division by zero due to zero Patankar weights
    σ = add_small_constant(σ, small_constant)

    P2, d2 = evaluate_pds(f, u, p, t + a21 * dt)
    integrator.stats.nf += 1

    Ptmp, dtmp = lincomb(b1, P, d, b2, P2, d2)

    u, solved = basic_patankar_step(uprev, Ptmp, σ, dt, alg.linsolve, dtmp)
    integrator.stats.nsolve += 1
    solved || return reject_step!(integrator, cache)

    # If a21 = 1, then σ is the MPE approximation, i.e. suited for stiff problems.
    # Otherwise, this is not clear.
    tmp = u - σ
    atmp = calculate_residuals(tmp, uprev, u, integrator.opts.abstol,
                               integrator.opts.reltol, integrator.opts.internalnorm, t)
    set_error_estimate!(integrator, cache, atmp, t)

    integrator.u = u
end

struct MPRK22Cache{uType, PType, dType, tabType, F} <: MPRKMutableCache
    tmp::uType
    P::PType
    P2::PType
    D::dType
    D2::dType
    σ::uType
    tab::tabType
    linsolve::F
end

get_tmp_cache(integrator, ::MPRK22, cache::OrdinaryDiffEqMutableCache) = (cache.σ,)

# In-place
function alg_cache(alg::MPRK22, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
                   uprev, uprev2, f, t, dt, reltol, p, calck,
                   ::Val{true},
                   verbose) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    a21, b1, b2 = get_constant_parameters(alg, uEltypeNoUnits)
    tab = MPRK22ConstantCache(a21, b1, b2, alg.small_constant_function(uEltypeNoUnits),
                              Ref(0))
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
        linsolve = init(linprob, linsolve_algorithm(alg.linsolve, linprob.A),
                        alias = LinearSolve.LinearAliasSpecifier(; alias_A = true,
                                                                 alias_b = true),
                        assumptions = LinearSolve.OperatorAssumptions(true))

        MPRK22Cache(tmp, P, P2, nothing, nothing, σ,
                    tab, #MPRK22ConstantCache
                    linsolve)
    elseif f isa PDSFunction
        linprob = LinearProblem(P2, _vec(tmp))
        linsolve = init(linprob, linsolve_algorithm(alg.linsolve, linprob.A),
                        alias = LinearSolve.LinearAliasSpecifier(; alias_A = true,
                                                                 alias_b = true),
                        assumptions = LinearSolve.OperatorAssumptions(true))

        MPRK22Cache(tmp, P, P2,
                    similar(u), # D
                    similar(u), # D2
                    σ,
                    tab, #MPRK22ConstantCache
                    linsolve)
    else
        throw(ArgumentError("MPRK22 can only be applied to production-destruction systems"))
    end
end

function initialize!(integrator, cache::MPRK22Cache)
    invalid_attempts(cache)[] = 0
end

@muladd function perform_step!(integrator, cache::MPRK22Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; tmp, P, P2, D, D2, σ, linsolve) = cache
    (; a21, b1, b2, small_constant) = cache.tab

    # We use P2 to store the last evaluation of the PDS
    # as well as to store the system matrix of the linear system

    # evaluate production/destruction terms
    evaluate_pds!(P, D, f, uprev, p, t)
    integrator.stats.nf += 1

    lincomb!(P2, a21, P)
    lincomb!(D2, a21, D)

    # avoid division by zero due to zero Patankar weights
    @.. broadcast=false σ=uprev + small_constant

    tmp .= uprev
    solved = basic_patankar_step!(u, tmp, P2, D2, σ, dt, linsolve)
    integrator.stats.nsolve += 1
    solved || return reject_step!(integrator, cache)

    if isone(a21)
        @.. broadcast=false σ=u + small_constant
    else
        @.. broadcast=false σ=σ^(1 - 1 / a21) * u^(1 / a21) + small_constant
    end

    # evaluate production/destruction terms
    evaluate_pds!(P2, D2, f, u, p, t + a21 * dt)
    integrator.stats.nf += 1

    lincomb!(P2, b1, P, b2, P2)
    lincomb!(D2, b1, D, b2, D2)

    tmp .= uprev
    solved = basic_patankar_step!(u, tmp, P2, D2, σ, dt, linsolve)
    integrator.stats.nsolve += 1
    solved || return reject_step!(integrator, cache)

    # Now σ stores the error estimate
    # If a21 = 1, then σ is the MPE approximation, i.e. suited for stiff problems.
    # Otherwise, this is not clear.
    @.. broadcast=false σ=u - σ

    # Now tmp stores error residuals
    calculate_residuals!(tmp, σ, uprev, u, integrator.opts.abstol,
                         integrator.opts.reltol, integrator.opts.internalnorm, t,
                         Serial())
    set_error_estimate!(integrator, cache, tmp, t)
    return nothing
end

### MPRK43 #####################################################################################
"""
    MPRK43I(α, β; [linsolve = ..., small_constant = ...])

A family of third-order modified Patankar-Runge-Kutta schemes for
production-destruction systems. This family of methods is based on the
two-parameter family of third order explicit Runge--Kutta schemes.
Each member of this family is an adaptive, one-step method with four-stages which is
third-order accurate, unconditionally positivity-preserving, conservative and linearly
implicit. In this implementation the stage-values are conservative as well.
The parameters `α` and `β` must be chosen such that the Runge--Kutta coefficients are nonnegative,
see Kopecz and Meister (2018) for details.

These methods support adaptive time stepping, using the Patankar-weight denominators
``σ_i``, see Kopecz and Meister (2018), as second order approximations to estimate the error.

The scheme was introduced by Kopecz and Meister for conservative production-destruction systems.
For nonconservative production–destruction systems we use the straight forward extension
analogous to [`MPE`](@ref).

These modified Patankar-Runge-Kutta methods require the special structure of a
[`PDSProblem`](@ref) or a [`ConservativePDSProblem`](@ref).

You can optionally choose the linear solver to be used by passing an
algorithm from [LinearSolve.jl](https://github.com/SciML/LinearSolve.jl)
as keyword argument `linsolve`.
You can also choose the parameter `small_constant` which is added to all Patankar-weight denominators
to avoid divisions by zero. You can pass a value explicitly, otherwise `small_constant` is set to
`floatmin` of the floating point type used.

## References

- Stefan Kopecz and Andreas Meister.
  "Unconditionally positive and conservative third order modified Patankar–Runge–Kutta
   discretizations of production–destruction systems."
   BIT Numerical Mathematics 58 (2018): 691–728.
  [DOI: 10.1007/s10543-018-0705-1](https://doi.org/10.1007/s10543-018-0705-1)
"""
struct MPRK43I{T, F, T2} <: OrdinaryDiffEqAdaptiveAlgorithm
    alpha::T
    beta::T
    linsolve::F
    small_constant_function::T2
end

function MPRK43I(alpha, beta; linsolve = LUFactorization(),
                 small_constant = nothing)
    if isnothing(small_constant)
        small_constant_function = floatmin
    elseif small_constant isa Number
        small_constant_function = Returns(small_constant)
    else # assume small_constant isa Function
        small_constant_function = small_constant
    end
    alpha, beta = promote(alpha, beta)
    alg = MPRK43I{typeof(alpha), typeof(linsolve), typeof(small_constant_function)}(alpha,
                                                                                    beta,
                                                                                    linsolve,
                                                                                    small_constant_function)
    check_parameters(alg)
    return alg
end

alg_order(::MPRK43I) = 3
isfsal(::MPRK43I) = false

function check_parameters(alg::MPRK43I)
    (; alpha, beta) = alg
    # `3 * alpha == 2` holds for 2/3 rounded to the type of `alpha`.
    if !(alpha ≥ one(alpha) / 3 && 3 * alpha ≠ 2)
        throw(ArgumentError("MPRK43I requires α ≥ 1/3 and α ≠ 2/3."))
    end
    two_thirds = 2 * one(alpha) / 3
    if alpha < two_thirds
        if !(two_thirds ≤ beta ≤ 3 * alpha * (1 - alpha))
            throw(ArgumentError("For this choice of α MPRK43I requires 2/3 ≤ β ≤ 3α(1-α)."))
        end
    else
        # For α > 2/3, the lower bound of β is the larger one of 3α(1-α) and
        # (3α-2)/(6α-3). They coincide at the only real root α0 ≈ 0.8926 of
        # 18α³ - 27α² + 12α - 2, and 3α(1-α) is the larger one for α < α0.
        lb_1 = 3 * alpha * (1 - alpha)
        lb_2 = (3 * alpha - 2) / (6 * alpha - 3)
        if lb_1 ≥ lb_2
            if !(lb_1 ≤ beta ≤ two_thirds)
                throw(ArgumentError("For this choice of α MPRK43I requires 3α(1-α) ≤ β ≤ 2/3."))
            end
        else
            if !(lb_2 ≤ beta ≤ two_thirds)
                throw(ArgumentError("For this choice of α MPRK43I requires (3α-2)/(6α-3) ≤ β ≤ 2/3."))
            end
        end
    end
    return nothing
end

function get_constant_parameters(alg::MPRK43I, type)
    check_parameters(alg)
    alpha = convert(type, alg.alpha)
    beta = convert(type, alg.beta)

    a21 = alpha
    a31 = (3 * alpha * beta * (1 - alpha) - beta^2) /
          (alpha * (2 - 3 * alpha))
    a32 = (beta * (beta - alpha)) / (alpha * (2 - 3 * alpha))
    b1 = 1 + (2 - 3 * (alpha + beta)) / (6 * alpha * beta)
    b2 = (3 * beta - 2) / (6 * alpha * (beta - alpha))
    b3 = (2 - 3 * alpha) / (6 * beta * (beta - alpha))
    c2 = alpha
    c3 = beta

    beta2 = 1 / (2 * a21)
    beta1 = 1 - beta2

    q1 = 1 / (3 * a21 * (a31 + a32) * b3)
    q2 = 1 / a21

    #This should never happen
    if !all((a21, a31, a32, b1, b2, b3, c2, c3, beta1, beta2) .≥ 0)
        throw(ArgumentError("MPRK43I requires nonnegative RK coefficients."))
    end
    return a21, a31, a32, b1, b2, b3, c2, c3, beta1, beta2, q1, q2
end

"""
    MPRK43II(γ; [linsolve = ..., small_constant = ...])

A family of third-order modified Patankar-Runge-Kutta schemes for
production-destruction systems. This family of methods is based on the
one-parameter family of third order explicit Runge--Kutta schemes with
non-negative Runge--Kutta coefficients.
Each member of this family is an adaptive, one-step method with four stages which is
third-order accurate, unconditionally positivity-preserving, conservative and linearly
implicit. In this implementation the stage-values are conservative as well. The parameter `γ` must satisfy
`3/8 ≤ γ ≤ 3/4`.
Further details are given in Kopecz and Meister (2018).

This method supports adaptive time stepping, using the Patankar-weight denominators
``σ_i``, see Kopecz and Meister (2018), as second order approximations to estimate the error.

The scheme was introduced by Kopecz and Meister for conservative production-destruction systems.
For nonconservative production–destruction systems we use the straight forward extension
analogous to [`MPE`](@ref).

These modified Patankar-Runge-Kutta methods require the special structure of a
[`PDSProblem`](@ref) or a [`ConservativePDSProblem`](@ref).

You can optionally choose the linear solver to be used by passing an
algorithm from [LinearSolve.jl](https://github.com/SciML/LinearSolve.jl)
as keyword argument `linsolve`.
You can also choose the parameter `small_constant` which is added to all Patankar-weight denominators
to avoid divisions by zero. To display the default value for data type `type` evaluate
`MPRK43II(gamma).small_constant_function(type)`, where `type` can be, e.g.,
`Float64`.

## References

- Stefan Kopecz and Andreas Meister.
  "Unconditionally positive and conservative third order modified Patankar–Runge–Kutta
   discretizations of production–destruction systems."
   BIT Numerical Mathematics 58 (2018): 691–728.
  [DOI: 10.1007/s10543-018-0705-1](https://doi.org/10.1007/s10543-018-0705-1)
"""
struct MPRK43II{T, F, T2} <: OrdinaryDiffEqAdaptiveAlgorithm
    gamma::T
    linsolve::F
    small_constant_function::T2
end

function small_constant_function_MPRK43II(type)
    if type == Float64
        small_constant = 1e-50
    else
        small_constant = floatmin(type)
    end
    return small_constant
end

function MPRK43II(gamma; linsolve = LUFactorization(),
                  small_constant = small_constant_function_MPRK43II)
    if small_constant isa Number
        small_constant_function = Returns(small_constant)
    else # assume small_constant isa Function
        small_constant_function = small_constant
    end
    alg = MPRK43II{typeof(gamma), typeof(linsolve), typeof(small_constant_function)}(gamma,
                                                                                     linsolve,
                                                                                     small_constant_function)
    check_parameters(alg)
    return alg
end

alg_order(::MPRK43II) = 3
isfsal(::MPRK43II) = false

function check_parameters(alg::MPRK43II)
    (; gamma) = alg
    if !(3 * one(gamma) / 8 ≤ gamma ≤ 3 * one(gamma) / 4)
        throw(ArgumentError("MPRK43II requires 3/8 ≤ γ ≤ 3/4."))
    end
    return nothing
end

function get_constant_parameters(alg::MPRK43II, type)
    check_parameters(alg)
    gamma = convert(type, alg.gamma)

    a21 = 2 * one(gamma) / 3
    a31 = a21 - 1 / (4 * gamma)
    a32 = 1 / (4 * gamma)
    b1 = one(gamma) / 4
    b2 = 3 * one(gamma) / 4 - gamma
    b3 = gamma
    c2 = a21
    c3 = a21

    beta2 = 1 / (2 * a21)
    beta1 = 1 - beta2

    q1 = 1 / (3 * a21 * (a31 + a32) * b3)
    q2 = 1 / a21

    #This should never happen
    if !all((a21, a31, a32, b1, b2, b3, c2, c3, beta1, beta2) .≥ 0)
        throw(ArgumentError("MPRK43II requires nonnegative RK coefficients."))
    end
    return a21, a31, a32, b1, b2, b3, c2, c3, beta1, beta2, q1, q2
end

struct MPRK43ConstantCache{T} <: OrdinaryDiffEqConstantCache
    a21::T
    a31::T
    a32::T
    b1::T
    b2::T
    b3::T
    c2::T
    c3::T
    beta1::T
    beta2::T
    q1::T
    q2::T
    small_constant::T
    # number of consecutive rejected step attempts
    invalid_attempts::Base.RefValue{Int}
end

# Out-of-place
function alg_cache(alg::Union{MPRK43I, MPRK43II}, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
                   uprev, uprev2, f, t, dt, reltol, p, calck,
                   ::Val{false},
                   verbose) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    if !(f isa PDSFunction || f isa ConservativePDSFunction)
        throw(ArgumentError("MPRK43 can only be applied to production-destruction systems"))
    end
    a21, a31, a32, b1, b2, b3, c2, c3, beta1, beta2, q1, q2 = get_constant_parameters(alg,
                                                                                      uEltypeNoUnits)
    MPRK43ConstantCache(a21, a31, a32, b1, b2, b3, c2, c3,
                        beta1, beta2, q1, q2, alg.small_constant_function(uEltypeNoUnits),
                        Ref(0))
end

function initialize!(integrator, cache::MPRK43ConstantCache)
    invalid_attempts(cache)[] = 0
end

@muladd function perform_step!(integrator, cache::MPRK43ConstantCache, repeat_step = false)
    (; alg, t, dt, uprev, f, p) = integrator
    (; a21, a31, a32, b1, b2, b3, c2, c3, beta1, beta2, q1, q2, small_constant) = cache

    # evaluate production and destruction terms
    P, d = evaluate_pds(f, uprev, p, t)
    integrator.stats.nf += 1

    Ptmp, dtmp = lincomb(a21, P, d)

    # avoid division by zero due to zero Patankar weights
    σ = add_small_constant(uprev, small_constant)
    if !(q1 ≈ q2)
        σ0 = σ
    end

    u2, solved = basic_patankar_step(uprev, Ptmp, σ, dt, alg.linsolve, dtmp)
    u = u2
    integrator.stats.nsolve += 1
    solved || return reject_step!(integrator, cache)

    # compute Patankar weight denominator
    σ = σ .^ (1 - q1) .* u .^ q1
    # avoid division by zero due to zero Patankar weights
    σ = add_small_constant(σ, small_constant)

    P2, d2 = evaluate_pds(f, u, p, t + c2 * dt)
    integrator.stats.nf += 1

    Ptmp, dtmp = lincomb(a31, P, d, a32, P2, d2)

    u, solved = basic_patankar_step(uprev, Ptmp, σ, dt, alg.linsolve, dtmp)
    integrator.stats.nsolve += 1
    solved || return reject_step!(integrator, cache)

    # compute Patankar weight denominator
    if !(q1 ≈ q2)
        σ = σ0 .^ (1 - q2) .* u2 .^ q2
        # avoid division by zero due to zero Patankar weights
        σ = add_small_constant(σ, small_constant)
    end

    Ptmp, dtmp = lincomb(beta1, P, d, beta2, P2, d2)

    σ, solved = basic_patankar_step(uprev, Ptmp, σ, dt, alg.linsolve, dtmp)
    integrator.stats.nsolve += 1
    solved || return reject_step!(integrator, cache)
    # avoid division by zero due to zero Patankar weights
    σ = add_small_constant(σ, small_constant)

    P3, d3 = evaluate_pds(f, u, p, t + c3 * dt)
    integrator.stats.nf += 1

    Ptmp, dtmp = lincomb(b1, P, d, b2, P2, d2, b3, P3, d3)

    u, solved = basic_patankar_step(uprev, Ptmp, σ, dt, alg.linsolve, dtmp)
    integrator.stats.nsolve += 1
    solved || return reject_step!(integrator, cache)

    tmp = u - σ
    atmp = calculate_residuals(tmp, uprev, u, integrator.opts.abstol,
                               integrator.opts.reltol, integrator.opts.internalnorm, t)
    set_error_estimate!(integrator, cache, atmp, t)

    integrator.u = u
end

struct MPRK43Cache{uType, PType, dType, tabType, F} <: MPRKMutableCache
    tmp::uType
    tmp2::uType
    P::PType
    P2::PType
    P3::PType
    D::dType
    D2::dType
    D3::dType
    σ::uType
    tab::tabType
    linsolve::F
end

function get_tmp_cache(integrator, ::Union{MPRK43I, MPRK43II},
                       cache::OrdinaryDiffEqMutableCache)
    (cache.σ,)
end

# In-place
function alg_cache(alg::Union{MPRK43I, MPRK43II}, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
                   uprev, uprev2, f, t, dt, reltol, p, calck,
                   ::Val{true},
                   verbose) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    a21, a31, a32, b1, b2, b3, c2, c3, beta1, beta2, q1, q2 = get_constant_parameters(alg,
                                                                                      uEltypeNoUnits)
    tab = MPRK43ConstantCache(a21, a31, a32, b1, b2, b3, c2, c3,
                              beta1, beta2, q1, q2,
                              alg.small_constant_function(uEltypeNoUnits), Ref(0))
    tmp = zero(u)
    tmp2 = zero(u)
    P = p_prototype(u, f)
    P2 = p_prototype(u, f)
    # We use P3 to store the last evaluation of the PDS
    # as well as to store the system matrix of the linear system
    P3 = p_prototype(u, f)
    σ = zero(u)

    if f isa ConservativePDSFunction
        # The right hand side of the linear system is always uprev. But using
        # tmp instead of uprev for the rhs we allow `alias_b=true`. uprev must
        # not be altered, since it is needed to compute the adaptive time step
        # size.
        linprob = LinearProblem(P3, _vec(tmp))
        linsolve = init(linprob, linsolve_algorithm(alg.linsolve, linprob.A),
                        alias = LinearSolve.LinearAliasSpecifier(; alias_A = true,
                                                                 alias_b = true),
                        assumptions = LinearSolve.OperatorAssumptions(true))
        MPRK43Cache(tmp, tmp2, P, P2, P3, nothing, nothing, nothing, σ, tab, linsolve)
    elseif f isa PDSFunction
        D = similar(u)
        D2 = similar(u)
        D3 = similar(u)

        linprob = LinearProblem(P3, _vec(tmp))
        linsolve = init(linprob, linsolve_algorithm(alg.linsolve, linprob.A),
                        alias = LinearSolve.LinearAliasSpecifier(; alias_A = true,
                                                                 alias_b = true),
                        assumptions = LinearSolve.OperatorAssumptions(true))

        MPRK43Cache(tmp, tmp2, P, P2, P3, D, D2, D3, σ, tab, linsolve)
    else
        throw(ArgumentError("MPRK43 can only be applied to production-destruction systems"))
    end
end

function initialize!(integrator, cache::MPRK43Cache)
    invalid_attempts(cache)[] = 0
end

@muladd function perform_step!(integrator, cache::MPRK43Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; tmp, tmp2, P, P2, P3, D, D2, D3, σ, linsolve) = cache
    (; a21, a31, a32, b1, b2, b3, c2, c3, beta1, beta2, q1, q2, small_constant) = cache.tab

    # We use P3 to store the last evaluation of the PDS
    # as well as to store the system matrix of the linear system

    # evaluate production/destruction terms
    evaluate_pds!(P, D, f, uprev, p, t)
    integrator.stats.nf += 1

    lincomb!(P3, a21, P)
    lincomb!(D3, a21, D)

    # avoid division by zero due to zero Patankar weights
    @.. broadcast=false σ=uprev + small_constant

    tmp .= uprev
    solved = basic_patankar_step!(u, tmp, P3, D3, σ, dt, linsolve)
    integrator.stats.nsolve += 1
    solved || return reject_step!(integrator, cache)

    if !(q1 ≈ q2)
        tmp2 .= u #u2 in out-of-place version
    end

    @.. broadcast=false σ=σ^(1 - q1) * u^q1 + small_constant

    # evaluate production/destruction terms
    evaluate_pds!(P2, D2, f, u, p, t + c2 * dt)
    integrator.stats.nf += 1

    lincomb!(P3, a31, P, a32, P2)
    lincomb!(D3, a31, D, a32, D2)

    tmp .= uprev
    solved = basic_patankar_step!(u, tmp, P3, D3, σ, dt, linsolve)
    integrator.stats.nsolve += 1
    solved || return reject_step!(integrator, cache)

    if !(q1 ≈ q2)
        @.. broadcast=false σ=(uprev + small_constant)^(1 - q2) * tmp2^q2 + small_constant
    end

    lincomb!(P3, beta1, P, beta2, P2)
    lincomb!(D3, beta1, D, beta2, D2)

    tmp .= uprev
    solved = basic_patankar_step!(σ, tmp, P3, D3, σ, dt, linsolve)
    integrator.stats.nsolve += 1
    solved || return reject_step!(integrator, cache)

    # avoid division by zero due to zero Patankar weights
    @.. broadcast=false σ=σ + small_constant

    # evaluate production/destruction terms
    evaluate_pds!(P3, D3, f, u, p, t + c3 * dt)
    integrator.stats.nf += 1

    lincomb!(P3, b1, P, b2, P2, b3, P3)
    lincomb!(D3, b1, D, b2, D2, b3, D3)

    tmp .= uprev
    solved = basic_patankar_step!(u, tmp, P3, D3, σ, dt, linsolve)
    integrator.stats.nsolve += 1
    solved || return reject_step!(integrator, cache)

    # Now tmp stores the error estimate
    @.. broadcast=false tmp=u - σ

    # Now tmp2 stores error residuals
    calculate_residuals!(tmp2, tmp, uprev, u, integrator.opts.abstol,
                         integrator.opts.reltol, integrator.opts.internalnorm, t,
                         Serial())
    set_error_estimate!(integrator, cache, tmp2, t)
    return nothing
end

########################################################################################
