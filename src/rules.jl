### DCP atom rules

struct SDP end

struct ExpSDP end

# add_dcprule(+, ℝ, AnySign, Affine, Increasing)

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

add_dcprule(Base.Ref, ℝ, AnySign, Affine, Increasing, MOI.Reals)

add_dcprule(dot, (array_domain(ℝ), array_domain(ℝ)), AnySign, Affine, Increasing, MOI.Reals)

function dotsort(x::AbstractVector, y::AbstractVector)
    if length(x) != length(y)
        throw(DimensionMismatch("AbstractVectors must have same length"))
    end
    dot(sort.(x, y))
end
Symbolics.@register_symbolic dotsort(x::AbstractVector, y::AbstractVector)
add_dcprule(dotsort, (array_domain(ℝ,1), array_domain(ℝ,1)), AnySign, Vex, (AnyMono, increasing_if_positive ∘ minimum), MOI.Reals)

add_dcprule(StatsBase.geomean, array_domain(HalfLine{Number, :open}(),1), Positive, Cave, Increasing, MOI.GeometricMeanCone)
add_dcprule(StatsBase.harmmean, array_domain(HalfLine{Number, :open}(),1), Positive, Cave, Increasing, MOI.SecondOrderCone) ##maybe wrong cone

# function invprod(x::AbstractVector)
#     if any(iszero(x))
#         throw(DivideError())
#     end
#     inv(prod(x))
# end
# Symbolics.@register_symbolic invprod(x::AbstractVector)

# add_dcprule(invprod, array_domain(HalfLine{Number, :open}()), Positive, Vex, Decreasing, )

add_dcprule(eigmax, symmetric_domain(), AnySign, Vex, AnyMono, SDP)

add_dcprule(eigmin, symmetric_domain(), AnySign, Cave, AnyMono, SDP)

function eigsummax(m::Symmetric, k::Int)
    if k < 1 || k > size(m, 1)
        throw(DomainError(k, "k must be between 1 and size(m, 1)"))
    end
    sum(eigvals(m)[end-k+1:end])
end
Symbolics.@register_symbolic eigsummax(m::Symmetric, k::Int)
add_dcprule(eigsummax, (array_domain(ℝ, 2), ℝ), AnySign, Vex, AnyMono, SDP)

function eigsummin(m::Symmetric, k::Int)
    if k < 1 || k > size(m, 1)
        throw(DomainError(k, "k must be between 1 and size(m, 1)"))
    end
    sum(eigvals(m)[1:k])
end
Symbolics.@register_symbolic eigsummin(m::Symmetric, k::Int)
add_dcprule(eigsummin, (array_domain(ℝ, 2), ℝ), AnySign, Cave, AnyMono, SDP)

add_dcprule(logdet, semidefinite_domain(), AnySign, Cave, AnyMono, ExpSDP)

add_dcprule(LogExpFunctions.logsumexp, array_domain(ℝ,2), AnySign, Vex, Increasing, MOI.ExponentialCone)

function matrix_frac(x::AbstractVector, P::AbstractMatrix)
    if length(x) != size(P, 1)
        throw(DimensionMismatch("x and P must have same length"))
    end
    return x' * inv(P) * x
end
Symbolics.@register_symbolic matrix_frac(x::AbstractVector, P::AbstractMatrix)
add_dcprule(matrix_frac, (array_domain(ℝ,1), definite_domain()), AnySign, Vex, AnyMono, SDP)

add_dcprule(maximum, array_domain(ℝ), AnySign, Vex, Increasing, MOI.Reals)

add_dcprule(minimum, array_domain(ℝ), AnySign, Cave, Increasing, MOI.Reals)

#incorrect for p<1
add_dcprule(norm, (array_domain(ℝ), Interval{:closed, :open}(1, Inf)), Positive, Vex, increasing_if_positive, MOI.NormCone)
add_dcprule(norm, (array_domain(ℝ), Interval{:closed, :open}(0, 1)), Positive, Vex, increasing_if_positive, MOI.NormCone)

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
function getcone() end
add_dcprule(perspective, (function_domain(), ℝ, Positive), getsign, getcurvature, AnyMono, getcone)

function quad_form(x::AbstractVector, P::AbstractMatrix)
    if length(x) != size(P, 1)
        throw(DimensionMismatch("x and P must have same length"))
    end
    return x' * P * x
end
Symbolics.@register_symbolic quad_form(x::AbstractVector, P::AbstractMatrix)
add_dcprule(quad_form, (array_domain(ℝ,1), semidefinite_domain()), Positive, Vex, increasing_if_positive, MOI.SecondOrderCone)
add_dcprule(quad_form, (array_domain(ℝ,1), negsemidefinite_domain()), Negative, Cave, increasing_if_positive ∘ -, MOI.SecondOrderCone)

function quad_over_lin(x::Union{Number, AbstractArray}, y::Number)
    return sum(x.^2) / y
end

Symbolics.@register_symbolic quad_over_lin(x::Union{Number, AbstractArray}, y::Number)
add_dcprule(quad_over_lin, (array_domain(ℝ), HalfLine{Number, :open}()), Positive, Vex, (increasing_if_positive, Decreasing), MOI.SecondOrderCone)


add_dcprule(sum, array_domain(ℝ, 2), AnySign, Affine, Increasing, MOI.Reals)

function sum_largest(x::AbstractMatrix, k::Integer)
    return sum(sort(vec(x))[end-k:end])
end
Symbolics.@register_symbolic sum_largest(x::AbstractMatrix, k::Integer)
add_dcprule(sum_largest, (array_domain(ℝ,2), ℤ), AnySign, Vex, Increasing, MOI.Reals)

function sum_smallest(x::AbstractMatrix, k::Integer)
    return sum(sort(vec(x))[1:k])
end

Symbolics.@register_symbolic sum_smallest(x::AbstractArray, k::Integer)
add_dcprule(sum_smallest, (array_domain(ℝ,2), ℤ), AnySign, Cave, Increasing, MOI.Reals)

add_dcprule(tr, array_domain(ℝ, 2), AnySign, Affine, Increasing)

function trinv(x::AbstractMatrix)
    return tr(inv(x))
end
Symbolics.@register_symbolic trinv(x::AbstractMatrix)
add_dcprule(trinv, definite_domain(), Positive, Vex, AnyMono)

function tv(x::AbstractVector{<:Number})
    return sum(abs.(x[2:end] - x[1:end-1]))
end
Symbolics.@register_symbolic tv(x::AbstractVector{<:Number})
add_dcprule(tv, array_domain(ℝ,1), Positive, Vex, AnyMono)

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
add_dcprule(tv, array_domain(array_domain(ℝ,2), 1), Positive, Vex, AnyMono)

add_dcprule(abs, ℂ, Positive, Vex, increasing_if_positive)

add_dcprule(conj, ℂ, AnySign, Affine, AnyMono)

add_dcprule(exp, ℝ, Positive, Vex, Increasing)

Symbolics.@register_symbolic LogExpFunctions.xlogx(x::Number)
add_dcprule(xlogx, ℝ, AnySign, Vex, AnyMono)

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
add_dcprule(huber, (ℝ, HalfLine()), Positive, Vex, increasing_if_positive)

add_dcprule(imag, ℂ, AnySign, Affine, AnyMono)

add_dcprule(inv, HalfLine{Number, :open}(), Positive, Vex, Decreasing)
add_dcprule(log, HalfLine{Number, :open}(), AnySign, Cave, Increasing)

add_dcprule(kldivergence, (array_domain(HalfLine{Number, :open},1), array_domain(HalfLine{Number, :open},1)), Positive, Vex, AnyMono)

function lognormcdf(x::Number)
    return logcdf(Normal, x)
end
Symbolics.@register_symbolic lognormcdf(x::Number)
add_dcprule(lognormcdf, ℝ, Negative, Cave, Increasing)

add_dcprule(log1p, Interval{:open, :open}(-1, Inf), Negative, Cave, Increasing)

add_dcprule(logistic, ℝ, Positive, Vex, Increasing)

add_dcprule(max, (ℝ, ℝ), AnySign, Vex, Increasing)
add_dcprule(min, (ℝ, ℝ), AnySign, Cave, Increasing)

# special cases which depend on arguments:
function dcprule(::typeof(^), x::Symbolic, i)
    args = (x, i)
    if isone(i)
        return makerule(ℝ, AnySign, Affine, Increasing), args
    elseif isinteger(i) && iseven(i)
        return makerule(ℝ, Positive, Vex, increasing_if_positive), args
    elseif isinteger(i) && isodd(i)
        return makerule(HalfLine(), Positive, Vex, Increasing), args
    elseif i >= 1
        return makerule(HalfLine(), Positive, Vex, Increasing), args
    elseif i > 0 && i < 1
        return makerule(HalfLine(), Positive, Cave, Increasing), args
    elseif i < 0
        return makerule(HalfLine{Float64, :closed}(), Positive, Vex, Increasing), args
    end
end
dcprule(::typeof(Base.literal_pow), f, x...) = dcprule(^, x...)

hasdcprule(::typeof(^)) = true

add_dcprule(sqrt, HalfLine(), Positive, Cave, Increasing)

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
add_dcprule(rel_entr, (HalfLine{Number, :open}(), HalfLine{Number, :open}()), AnySign, Vex, (AnyMono, Decreasing))

add_dcprule(sqrt, HalfLine(), Positive, Cave, Increasing)

add_dcprule(xexpx, HalfLine, Positive, Vex, Increasing)

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