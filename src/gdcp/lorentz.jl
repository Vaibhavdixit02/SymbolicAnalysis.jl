# DGCP Atoms for Lorentz Model (Hyperbolic Space)
#
# This file implements geodesically convex atoms for the Lorentz model,
# a Cartan-Hadamard manifold of constant negative curvature.
# Based on the results from Ferreira, Németh, and Zhu (2022, 2023).

using Manifolds
using LinearAlgebra
using Symbolics: Symbolic, @register_symbolic, unwrap, variables

@register_symbolic Manifolds.distance(
    M::Manifolds.Lorentz,
    p::AbstractVector,
    q::Union{Symbolics.Arr,AbstractVector},
) false
add_gdcprule(Manifolds.distance, Manifolds.Lorentz, Positive, GConvex, GAnyMono)

"""
    lorentz_log_barrier(a, p)

Computes the log-barrier function for the Lorentz model: `-log(-1 - <a, p>_L)`.

# Arguments
    - `a`: The vector (0, ..., 0, 1) in R^(d+1).
    - `p`: A point on the Lorentz manifold.
"""
function lorentz_log_barrier(p::AbstractVector)
    # Lorentzian inner product: a⋅p_L = a1*p1 + ... + a_d*p_d - a_{d+1}*p_{d+1}
    inner_prod = a[end] * p[end]
    return -log(-1 + inner_prod)
end

@register_symbolic lorentz_log_barrier(p::Union{Symbolics.Arr,AbstractVector})
add_gdcprule(lorentz_log_barrier, Manifolds.Lorentz, Positive, GConvex, GIncreasing)

"""
    lorentz_homogeneous_quadratic(A::AbstractMatrix, p::AbstractVector)

Computes the homogeneous quadratic function `p'Ap` for the Lorentz model.
For geodesic convexity, A should be positive semidefinite.

# Arguments
    - `A::AbstractMatrix`: A (d+1)×(d+1) positive semidefinite matrix.
    - `p::AbstractVector`: A point on the Lorentz manifold.
"""
function lorentz_homogeneous_quadratic(A::AbstractMatrix, p::AbstractVector)
    return p' * A * p
end

@register_symbolic lorentz_homogeneous_quadratic(
    A::AbstractMatrix,
    p::Union{Symbolics.Arr,AbstractVector},
)
add_gdcprule(lorentz_homogeneous_quadratic, Manifolds.Lorentz, Positive, GConvex, GAnyMono)

"""
    lorentz_homogeneous_diagonal(a::AbstractVector, p::AbstractVector)

Computes the homogeneous diagonal quadratic function `∑(a_i * p_i^2)`.
For geodesic convexity, min(a_1,...,a_d) + a_{d+1} ≥ 0.

# Arguments
    - `a::AbstractVector`: A (d+1)-vector where min(a_1,...,a_d) + a_{d+1} ≥ 0.
    - `p::AbstractVector`: A point on the Lorentz manifold.
"""
function lorentz_homogeneous_diagonal(a::AbstractVector, p::AbstractVector)
    if length(a) != length(p)
        throw(DimensionMismatch("Vectors must have same length"))
    end

    if minimum(a[1:end-1]) + a[end] < 0
        throw(
            ArgumentError(
                "For geodesic convexity, min(a[1:end-1]) + a[end] ≥ 0 is required",
            ),
        )
    end

    return sum(a .* p .^ 2)
end

@register_symbolic lorentz_homogeneous_diagonal(
    a::AbstractVector,
    p::Union{Symbolics.Arr,AbstractVector},
)
add_gdcprule(lorentz_homogeneous_diagonal, Manifolds.Lorentz, Positive, GConvex, GAnyMono)

"""
    lorentz_least_squares(X::AbstractMatrix, y::AbstractVector, p::AbstractVector)

Computes the least squares function `‖y - Xp‖²_2 = y'y - 2y'Xp + p'X'Xp` for the Lorentz model.
For geodesic convexity, the following conditions must be satisfied:
1. ∑^d_i=1(X'y)^2_i ≤ (X'y)^2_{d+1}
2. (X'y)_{d+1} ≤ 0

# Arguments
    - `X::AbstractMatrix`: A matrix in R^(n×(d+1)).
    - `y::AbstractVector`: A vector in R^n.
    - `p::AbstractVector`: A point on the Lorentz manifold.
"""
function lorentz_least_squares(X::AbstractMatrix, y::AbstractVector, p::AbstractVector)
    Xty = X' * y
    
    # Check conditions for geodesic convexity
    condition1 = sum(Xty[1:end-1].^2) <= Xty[end]^2
    condition2 = Xty[end] <= 0
    
    if !(condition1 && condition2)
        throw(ArgumentError("Conditions for geodesic convexity not satisfied: " * 
                           "∑^d_i=1(X'y)^2_i ≤ (X'y)^2_{d+1} and (X'y)_{d+1} ≤ 0"))
    end

    return norm(X * p - y)^2
end

@register_symbolic lorentz_least_squares(
    A::AbstractMatrix,
    b::AbstractVector,
    p::Union{Symbolics.Arr,AbstractVector},
)
add_gdcprule(lorentz_least_squares, Manifolds.Lorentz, Positive, GConvex, GAnyMono)

"""
    lorentz_transform(O::AbstractMatrix, p::AbstractVector)

Applies a Lorentz transform to a point on the Lorentz manifold.
The matrix O must be an element of the orthochronous Lorentz group O⁺(1,d).

# Arguments
    - `O::AbstractMatrix`: An element of the orthochronous Lorentz group.
    - `p::AbstractVector`: A point on the Lorentz manifold.
"""
function lorentz_transform(O::AbstractMatrix, p::AbstractVector)
    d = length(p) - 1
    J = Diagonal([ones(d)..., -1])

    # Check if O is in the Lorentz group
    if !isapprox(O' * J * O, J, rtol = 1e-10)
        throw(ArgumentError("Matrix is not in the Lorentz group"))
    end

    # Check if O preserves the positive time direction (orthochronous)
    if (O*[zeros(d)..., 1])[end] <= 0
        throw(ArgumentError("Matrix does not preserve the positive time direction"))
    end

    return O * p
end

@register_symbolic lorentz_transform(
    O::AbstractMatrix,
    p::Union{Symbolics.Arr,AbstractVector},
)
# Not adding a rule since this preserves geodesic convexity but doesn't have a specific curvature

# Export functions
export lorentz_log_barrier, lorentz_homogeneous_quadratic
export lorentz_homogeneous_diagonal, lorentz_least_squares, lorentz_transform
