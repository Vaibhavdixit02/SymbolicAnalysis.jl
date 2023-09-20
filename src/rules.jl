
add_dcprule(abs, ℝ, Positive, Vex, increasing_if_positive)
add_dcprule(exp, ℝ, Positive, Vex, Increasing)
Symbolics.@register_symbolic LogExpFunctions.xlogx(x::Number)
add_dcprule(xlogx, ℝ, AnySign, Vex, AnyMono)

add_dcprule(dot, (array_domain(ℝ), array_domain(ℝ)), AnySign, Affine, Increasing)

function dotsort(x::Vector, y::Vector)
    if length(x) != length(y)
        throw(DimensionMismatch("Vectors must have same length"))
    end
    dot(sort.(x, y))
end
Symbolics.@register_symbolic dotsort(x::Vector, y::Vector)
add_dcprule(dotsort, (array_domain(ℝ,1), array_domain(ℝ,1)), AnySign, Vex, (AnyMono, increasing_if_positive ∘ minimum))

add_dcprule(StatsBase.geomean, array_domain(HalfLine{Float64, :open}(),1), Positive, Cave, Increasing)
add_dcprule(StatsBase.harmonicmean, array_domain(HalfLine{Float64, :open}(),1), Positive, Cave, Increasing)

function invprod(x::Vector)
    if any(iszero(x))
        throw(DivideError())
    end
    inv(prod(x))
end
Symbolics.@register_symbolic invprod(x::Vector)

add_dcprule(invprod, array_domain(HalfLine{Float64, :open}()), Positive, Vex, Decreasing)

add_dcprule(eigmax, symmetric_domain(), AnySign, Vex, AnyMono)

add_dcprule(eigmin, symmetric_domain(), AnySign, Cave, AnyMono)

function eigsummax(m::Symmetric, k::Int)
    if k < 1 || k > size(m, 1)
        throw(DomainError(k, "k must be between 1 and size(m, 1)"))
    end
    sum(eigvals(m)[end-k+1:end])
end
Symbolics.@register_symbolic eigsummax(m::Symmetric, k::Int)
add_dcprule(eigsummax, (array_domain(ℝ, 2), ℝ), AnySign, Vex, AnyMono)

function eigsummin(m::Symmetric, k::Int)
    if k < 1 || k > size(m, 1)
        throw(DomainError(k, "k must be between 1 and size(m, 1)"))
    end
    sum(eigvals(m)[1:k])
end
Symbolics.@register_symbolic eigsummin(m::Symmetric, k::Int)
add_dcprule(eigsummin, (array_domain(ℝ, 2), ℝ), AnySign, Cave, AnyMono)

add_dcprule(LogExpFunctions.logdet, semidefinite_domain(), AnySign, Cave, AnyMono)

add_dcprule(LogExpFunctions.logsumexp, array_domain(ℝ,2), AnySign, Vex, Increasing)

function matrix_frac(x::Vector, P::Matrix)
    if length(x) != size(P, 1)
        throw(DimensionMismatch("x and P must have same length"))
    end
    return x' * P * x
end
Symbolics.@register_symbolic matrix_frac(x::Vector, P::Matrix)
add_dcprule(matrix_frac, (array_domain(ℝ,1), definite_domain()), AnySign, Vex, AnyMono)

add_dcprule(maximum, array_domain(ℝ), AnySign, Vex, Increasing)

add_dcprule(minimum, array_domain(ℝ), AnySign, Cave, Increasing)

#incorrect for p<1
add_dcprule(norm, (array_domain(ℝ), Interval{:closed, :open}(1, Inf)), Positive, Vex, increasing_if_positive)
add_dcprule(norm, (array_domain(ℝ), Interval{:closed, :open}(0, 1)), Positive, Vex, increasing_if_positive)

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

function quad_form(x::Vector, P::Matrix)
    if length(x) != size(P, 1)
        throw(DimensionMismatch("x and P must have same length"))
    end
    return x' * P * x
end
Symbolics.@register_symbolic quad_form(x::Vector, P::Matrix)
add_dcprule(quad_form, (array_domain(ℝ,1), semidefinite_domain()), Positive, Vex, increasing_if_positive)
add_dcprule(quad_form, (array_domain(ℝ,1), negsemidefinite_domain()), Negative, Cave, increasing_if_positive ∘ -)

function quad_over_lin(x::Matrix, y::Number)
    if y < 0
        throw(DomainError(y, "y must be positive"))
    end
    return sum(x.^2) / y
end

Symbolics.@register_symbolic quad_over_lin(x::Matrix, y::Number)
add_dcprule(quad_over_lin, (array_domain(ℝ,2), HalfLine{Float64, :open}()), Positive, Vex, (increasing_if_positive, Decreasing))

add_dcprule(sum, array_domain(ℝ, 2), AnySign, Affine, Increasing)

function sum_largest(x::Matrix, k::Integer)
    return sum(sort(vec(x))[end-k:end])
end
Symbolics.@register_symbolic sum_largest(x::Matrix, k::Integer)
add_dcprule(sum_largest, (array_domain(ℝ,2), ℤ), AnySign, Vex, Increasing)

function sum_smallest(x::Matrix, k::Integer)
    return sum(sort(vec(x))[1:k])
end
Symbolics.@register_symbolic sum_smallest(x::Matrix, k::Integer)
add_dcprule(sum_smallest, (array_domain(ℝ,2), ℤ), AnySign, Cave, Increasing)

add_dcprule(tr, array_domain(ℝ, 2), AnySign, Affine, Increasing)

function trinv(x::Matrix)
    return tr(inv(x))
end
Symbolics.@register_symbolic trinv(x::Matrix)
add_dcprule(trinv, definite_domain(), Positive, Vex, AnyMono)

function tv(x::Vector{<:Number})
    return sum(abs.(x[2:end] - x[1:end-1]))
end
Symbolics.@register_symbolic tv(x::Vector{<:Number})
add_dcprule(tv, array_domain(ℝ,1), Positive, Vex, AnyMono)

function tv(x::Vector{<:Matrix})
    return sum(
        map(1:size(x, 1)-1) do i
            map(1:size(x, 2)-1) do j
                norm([x[k][i+1, j] - x[k][i, j] for k in eachindex(x)])
            end
        end
        )
end
Symbolics.@register_symbolic tv(x::Vector{<:Matrix})
add_dcprule(tv, array_domain(array_domain(ℝ,2), 1), Positive, Vex, AnyMono)


# add_dcprule(huber, ℝ, Positive, Vex, increasing_if_positive)
add_dcprule(inv, HalfLine{Float64, :closed}(), Positive, Vex, increasing_if_positive)
add_dcprule(log, HalfLine{Float64, :closed}(), AnySign, Cave, AnyMono)
#add_dcprule(log_sum_exp, array_domain(ℝ,1), AnySign, Cave, AnyMono)


#add_dcprule(positive, ℝ, Positive, Vex, Increasing)
#add_dcprule(^, ℝ, Positive, Vex, Increasing) # Requires special handling based on 2nd arg
#add_dcprule(quad_over_lin, (ℝ, HalfLine{Float64, :closed}()), Positive, Vex, (increasing_if_positive, Decreasing)) # Requires special handling based on 2nd arg
add_dcprule(sqrt, HalfLine(), Positive, Cave, Increasing)
