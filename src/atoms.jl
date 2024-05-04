### DCP atom rules

add_dcprule(+, ℝ, AnySign, Affine, Increasing)

# function dcprule(::typeof(-), x, y)
#     if y >= 0
#         if x >= 0 && x >= y 
#             return makerule((HalfLine(), HalfLine()), Positive, Affine, (Increasing, Decreasing))
#         elseif x >= 0 && x < y
#             return makerule((HalfLine(), HalfLine()), Negative, Affine, (Increasing, Decreasing))
#         elseif x < 0+
#             return makerule((NegativeHalfLine(), HalfLine()), Negative, Affine, (Increasing, Decreasing))
#         end
#     elseif y < 0
#         if x >= 0
#             return makerule((HalfLine(), NegativeHalfLine()), Positive, Affine, (Increasing, Decreasing))
#         else
#             if x >= y
#                 return makerule((NegativeHalfLine(), NegativeHalfLine()), Positive, Affine, (Increasing, Decreasing))
#             else
#                 return makerule((NegativeHalfLine(), NegativeHalfLine()), Negative, Affine, (Increasing, Decreasing))
#             end
#         end
#     end
# end

add_dcprule(Base.Ref, ℝ, AnySign, Affine, Increasing)

add_dcprule(dot, (array_domain(ℝ), array_domain(ℝ)), AnySign, Affine, Increasing)

"""
    dotsort(x, y)

Sorts `x` and `y` and returns the dot product of the sorted vectors.

# Arguments
    - `x::AbstractVector`: A vector.
    - `y::AbstractVector`: A vector.
"""
function dotsort(x::AbstractVector, y::AbstractVector)
    if length(x) != length(y)
        throw(DimensionMismatch("AbstractVectors must have same length"))
    end
    dot(sort.(x, y))
end
Symbolics.@register_symbolic dotsort(x::AbstractVector, y::AbstractVector)
add_dcprule(dotsort, (array_domain(ℝ,1), array_domain(ℝ,1)), AnySign, Convex, (AnyMono, increasing_if_positive ∘ minimum))

add_dcprule(StatsBase.geomean, array_domain(HalfLine{Number, :open}(),1), Positive, Concave, Increasing)
add_dcprule(StatsBase.harmmean, array_domain(HalfLine{Number, :open}(),1), Positive, Concave, Increasing)

"""
    invprod(x::AbstractVector)

Returns the inverse of the product of the elements of `x`.

# Arguments
    - `x::AbstractVector`: A vector.
"""
function invprod(x::AbstractVector)
    if any(iszero(x))
        throw(DivideError())
    end
    inv(prod(x))
end
Symbolics.@register_symbolic invprod(x::AbstractVector)

add_dcprule(invprod, array_domain(HalfLine{Number, :open}()), Positive, Convex, Decreasing)

add_dcprule(eigmax, symmetric_domain(), AnySign, Convex, AnyMono)

add_dcprule(eigmin, symmetric_domain(), AnySign, Concave, AnyMono)

"""
    eigsummax(m::Symmetric, k)

Returns the sum of the `k` largest eigenvalues of `m`.

# Arguments
    - `m::Symmetric`: A symmetric matrix.
    - `k::Int`: The number of largest eigenvalues to sum.
"""
function eigsummax(m::Symmetric, k::Int)
    if k < 1 || k > size(m, 1)
        throw(DomainError(k, "k must be between 1 and size(m, 1)"))
    end
    sum(eigvals(m)[end-k+1:end])
end
Symbolics.@register_symbolic eigsummax(m::Symmetric, k::Int)
add_dcprule(eigsummax, (array_domain(ℝ, 2), ℝ), AnySign, Convex, AnyMono)

"""
    eigsummin(m::Symmetric, k)

Returns the sum of the `k` smallest eigenvalues of `m`.

# Arguments
    - `m::Symmetric`: A symmetric matrix.
    - `k::Int`: The number of smallest eigenvalues to sum.
"""
function eigsummin(m::Symmetric, k::Int)
    if k < 1 || k > size(m, 1)
        throw(DomainError(k, "k must be between 1 and size(m, 1)"))
    end
    sum(eigvals(m)[1:k])
end
Symbolics.@register_symbolic eigsummin(m::Symmetric, k::Int)
add_dcprule(eigsummin, (array_domain(ℝ, 2), ℝ), AnySign, Concave, AnyMono)

add_dcprule(logdet, semidefinite_domain(), AnySign, Concave, AnyMono)

add_dcprule(LogExpFunctions.logsumexp, array_domain(ℝ,2), AnySign, Convex, Increasing)

"""
    matrix_frac(x::AbstractVector, P::AbstractMatrix)

Returns the quadratic form `x' * P^{-1} * x`.

# Arguments
    - `x::AbstractVector`: A vector.
    - `P::AbstractMatrix`: A matrix.
"""
function matrix_frac(x::AbstractVector, P::AbstractMatrix)
    if length(x) != size(P, 1)
        throw(DimensionMismatch("x and P must have same length"))
    end
    return x' * inv(P) * x
end
Symbolics.@register_symbolic AbstractMatrix_frac(x::AbstractVector, P::AbstractMatrix)
add_dcprule(matrix_frac, (array_domain(ℝ,1), definite_domain()), AnySign, Convex, AnyMono)

add_dcprule(maximum, array_domain(ℝ), AnySign, Convex, Increasing)

add_dcprule(minimum, array_domain(ℝ), AnySign, Concave, Increasing)

#incorrect for p<1
add_dcprule(norm, (array_domain(ℝ), Interval{:closed, :open}(1, Inf)), Positive, Convex, increasing_if_positive)
add_dcprule(norm, (array_domain(ℝ), Interval{:closed, :open}(0, 1)), Positive, Convex, increasing_if_positive)

"""
    perspective(f::Function, x, s::Number)

Returns the perspective function `s * f(x / s)`.

# Arguments
    - `f::Function`: A function.
    - `x`: A number.
    - `s::Number`: A positive number.
"""
function perspective(f::Function, x, s::Number)
    if s < 0
        throw(DomainError(s, "s must be positive"))
    end
    if s == 0
        return zero(typeof(f(x)))
    end
    s * f(x / s)
end
Symbolics.@register_symbolic perspective(f::Function, x, s::Number)
add_dcprule(perspective, (function_domain(), ℝ, Positive), getsign, getcurvature, AnyMono)

"""
    quad_form(x::AbstractVector, P::AbstractMatrix)

Returns the quadratic form `x' * P * x`.

# Arguments
    - `x::AbstractVector`: A vector.
    - `P::AbstractMatrix`: A matrix.
"""
function quad_form(x::AbstractVector, P::AbstractMatrix)
    if length(x) != size(P, 1)
        throw(DimensionMismatch("x and P must have same length"))
    end
    return x' * P * x
end
Symbolics.@register_symbolic quad_form(x::AbstractVector, P::AbstractMatrix)
add_dcprule(quad_form, (array_domain(ℝ,1), semidefinite_domain()), Positive, Convex, (increasing_if_positive, Increasing))

function quad_over_lin(x::AbstractArray, y::Number)
    if y < 0
        throw(DomainError(y, "y must be positive"))
    end
    return sum(x.^2) / y
end

Symbolics.@register_symbolic quad_over_lin(x::Symbolics.Arr, y::Num)

"""
    quad_over_lin(x, y::Number)

Returns the quadratic over linear form `x^2 / y`.

# Arguments
    - `x`: A number or a vector.
    - `y::Number`: A positive number.
"""
function quad_over_lin(x::Real, y::Real)
    if getsign(y) == Negative
        throw(DomainError(y, "y must be positive"))
    end
    return x^2 / y
end

Symbolics.@register_symbolic quad_over_lin(x::Number, y::Number)

add_dcprule(quad_over_lin, (array_domain(ℝ), HalfLine{Number, :open}()), Positive, Convex, (increasing_if_positive, Decreasing))

add_dcprule(quad_over_lin, (ℝ, HalfLine{Number, :open}()), Positive, Convex, (increasing_if_positive, Decreasing))

add_dcprule(sum, array_domain(ℝ, 2), AnySign, Affine, Increasing)

"""
    sum_largest(x::AbstractMatrix, k)

Returns the sum of the `k` largest elements of `x`.

# Arguments
    - `x::AbstractMatrix`: A matrix.
    - `k::Int`: The number of largest elements to sum.
"""
function sum_largest(x::AbstractMatrix, k::Integer)
    return sum(sort(vec(x))[end-k:end])
end
Symbolics.@register_symbolic sum_largest(x::AbstractMatrix, k::Integer)
add_dcprule(sum_largest, (array_domain(ℝ,2), ℤ), AnySign, Convex, Increasing)

"""
    sum_smallest(x::AbstractMatrix, k)

Returns the sum of the `k` smallest elements of `x`.

# Arguments
    - `x::AbstractMatrix`: A matrix.
    - `k::Int`: The number of smallest elements to sum.
"""
function sum_smallest(x::AbstractMatrix, k::Integer)
    return sum(sort(vec(x))[1:k])
end

Symbolics.@register_symbolic sum_smallest(x::AbstractArray, k::Integer)
add_dcprule(sum_smallest, (array_domain(ℝ,2), ℤ), AnySign, Concave, Increasing)

add_dcprule(tr, array_domain(ℝ, 2), AnySign, Affine, Increasing)

"""
    trinv(x::AbstractMatrix)

Returns the trace of the inverse of `x`.

# Arguments
    - `x::AbstractMatrix`: A matrix.
"""
function trinv(x::AbstractMatrix)
    return tr(inv(x))
end
Symbolics.@register_symbolic trinv(x::AbstractMatrix)
add_dcprule(trinv, definite_domain(), Positive, Convex, AnyMono)

"""
    tv(x::AbstractVector{<:Number})

Returns the total variation of `x`, defined as `\sum_i |x_{i+1} - x_i|`.

# Arguments
    - `x::AbstractVector`: A vector.
"""
function tv(x::AbstractVector{<:Number})
    return sum(abs.(x[2:end] - x[1:end-1]))
end
Symbolics.@register_symbolic tv(x::AbstractVector{<:Number})
add_dcprule(tv, array_domain(ℝ,1), Positive, Convex, AnyMono)

"""
    tv(x::AbstractVector{<:AbstractMatrix})

Returns the total variation of `x`, defined as `\sum_{i,j} |x_{k+1}[i,j] - x_k[i,j]|`.

# Arguments
    - `x::AbstractVector`: A vector of matrices.
"""
function tv(x::AbstractVector{<:AbstractMatrix})
    return sum(
        map(1:size(x, 1)-1) do i
            map(1:size(x, 2)-1) do j
                norm([x[k][i+1, j] - x[k][i, j] for k in eachindex(x)])
            end
        end
        )
end
Symbolics.@register_symbolic tv(x::AbstractVector{<:AbstractMatrix})
add_dcprule(tv, array_domain(array_domain(ℝ,2), 1), Positive, Convex, AnyMono)

add_dcprule(abs, ℂ, Positive, Convex, increasing_if_positive)

add_dcprule(conj, ℂ, AnySign, Affine, AnyMono)

add_dcprule(exp, ℝ, Positive, Convex, Increasing)

Symbolics.@register_symbolic LogExpFunctions.xlogx(x::Number)
add_dcprule(xlogx, ℝ, AnySign, Convex, AnyMono)

"""
    huber(x, M=1)

Returns the Huber loss function of `x` with threshold `M`.

# Arguments
    - `x::Number`: A number.
    - `M::Number`: The threshold.
"""
function huber(x::Number, M::Number = 1)
    if M < 0
        throw(DomainError(M, "M must be positive"))
    end
    
    if abs(x) <= M
        return x^2
    else
        return 2* M * abs(x) - M^2
    end
end
Symbolics.@register_symbolic huber(x::Number, M::Number)
add_dcprule(huber, (ℝ, HalfLine()), Positive, Convex, increasing_if_positive)

add_dcprule(imag, ℂ, AnySign, Affine, AnyMono)

add_dcprule(inv, HalfLine{Number, :open}(), Positive, Convex, Decreasing)
add_dcprule(log, HalfLine{Number, :open}(), AnySign, Concave, Increasing)

@register_symbolic Base.log(A::Symbolics.Arr)

add_dcprule(log, array_domain(ℝ, 2), Positive, Concave, Increasing)

@register_symbolic LinearAlgebra.inv(A::Symbolics.Arr)

add_dcprule(inv, semidefinite_domain(), AnySign, Convex, Decreasing)

@register_symbolic LinearAlgebra.sqrt(A::Symbolics.Arr)

add_dcprule(sqrt, semidefinite_domain(), Positive, Concave, Increasing)

add_dcprule(kldivergence, (array_domain(HalfLine{Number, :open},1), array_domain(HalfLine{Number, :open},1)), Positive, Convex, AnyMono)

"""
    lognormcdf(x::Number)

Returns the log of the normal cumulative distribution function of `x`.

# Arguments
    - `x::Number`: A number.
"""
function lognormcdf(x::Number)
    return logcdf(Normal, x)
end
Symbolics.@register_symbolic lognormcdf(x::Number)
add_dcprule(lognormcdf, ℝ, Negative, Concave, Increasing)

add_dcprule(log1p, Interval{:open, :open}(-1, Inf), Negative, Concave, Increasing)

add_dcprule(logistic, ℝ, Positive, Convex, Increasing)

add_dcprule(max, (ℝ, ℝ), AnySign, Convex, Increasing)
add_dcprule(min, (ℝ, ℝ), AnySign, Concave, Increasing)

# special cases which depend on arguments:
function dcprule(::typeof(^), x::Symbolic, i)
    args = (x, i)
    if isone(i)
        return makerule(ℝ, AnySign, Affine, Increasing), args
    elseif isinteger(i) && iseven(i)
        return makerule(ℝ, Positive, Convex, increasing_if_positive), args
    elseif isinteger(i) && isodd(i)
        return makerule(HalfLine(), Positive, Convex, Increasing), args
    elseif i >= 1
        return makerule(HalfLine(), Positive, Convex, Increasing), args
    elseif i > 0 && i < 1
        return makerule(HalfLine(), Positive, Concave, Increasing), args
    elseif i < 0
        return makerule(HalfLine{Float64, :closed}(), Positive, Convex, Increasing), args
    end
end
dcprule(::typeof(Base.literal_pow), f, x...) = dcprule(^, x...)

hasdcprule(::typeof(^)) = true

add_dcprule(real, ℂ, AnySign, Affine, Increasing)

function rel_entr(x::Number, y::Number)
    if x < 0 || y < 0
        throw(DomainError((x, y), "x and y must be positive"))
    end
    if x == 0
        return 0
    end
    x * log(x / y)
end
Symbolics.@register_symbolic rel_entr(x::Number, y::Number)
add_dcprule(rel_entr, (HalfLine{Number, :open}(), HalfLine{Number, :open}()), AnySign, Convex, (AnyMono, Decreasing))

add_dcprule(sqrt, HalfLine(), Positive, Concave, Increasing)

add_dcprule(xexpx, HalfLine, Positive, Convex, Increasing)

add_dcprule(conv, (array_domain(ℝ,1), array_domain(ℝ,1)), AnySign, Affine, AnyMono)

add_dcprule(cumsum, array_domain(ℝ), AnySign, Affine, Increasing)

add_dcprule(diagm, array_domain(ℝ,1), AnySign, Affine, Increasing)

add_dcprule(diag, array_domain(ℝ,2), AnySign, Affine, Increasing)

add_dcprule(diff, array_domain(ℝ), AnySign, Affine, Increasing)

add_dcprule(hcat, array_domain(array_domain(ℝ,1), 1), AnySign, Affine, Increasing)

add_dcprule(kron, (array_domain(ℝ,2), array_domain(ℝ,2)), AnySign, Affine, Increasing)

add_dcprule(reshape, array_domain(ℝ, 2), AnySign, Affine, Increasing)

add_dcprule(triu, array_domain(ℝ, 2), AnySign, Affine, Increasing)

add_dcprule(vec, array_domain(ℝ, 2), AnySign, Affine, Increasing)

add_dcprule(vcat, array_domain(array_domain(ℝ,1), 1), AnySign, Affine, Increasing)

function dcprule(::typeof(broadcast), f, x...)
    return dcprule(f, x...)
end
hasdcprule(::typeof(broadcast)) = true

# add_dcprule(broadcast, (function_domain, array_domain(ℝ)), AnySign, Affine, (AnyMono, AnyMono))