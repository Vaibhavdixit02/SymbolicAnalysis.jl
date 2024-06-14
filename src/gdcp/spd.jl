@register_symbolic LinearAlgebra.logdet(X::Union{Symbolics.Arr, Matrix{Num}})
add_gdcprule(LinearAlgebra.logdet, SymmetricPositiveDefinite, Positive, GLinear, GIncreasing)

"""
    conjugation(X, B)

Conjugation of a matrix `X` by a matrix `B` is defined as `B'X*B`.

# Arguments
    - `X::Matrix`: A symmetric positive definite matrix.
    - `B::Matrix`: A matrix.
"""
function conjugation(X, B)
    return B'*X*B
end

@register_array_symbolic conjugation(X::Matrix, B::Union{Symbolics.Arr, Matrix{Num}}) begin
    size=(size(B,2), size(B,2))
end

add_gdcprule(conjugation, SymmetricPositiveDefinite, Positive, GConvex, GIncreasing)

@register_symbolic LinearAlgebra.tr(X::Union{Symbolics.Arr, Matrix{Num}})
add_gdcprule(LinearAlgebra.tr, SymmetricPositiveDefinite, Positive, GConvex, GIncreasing)

add_gdcprule(sum, SymmetricPositiveDefinite, Positive, GConvex, GIncreasing)

add_gdcprule(adjoint, SymmetricPositiveDefinite, Positive, GLinear, GIncreasing)

"""
    scalar_mat(X, k=size(X, 1))

Scalar matrix of a symmetric positive definite matrix `X` is defined as `tr(X)*I(k)`.

# Arguments
    - `X::Matrix`: A symmetric positive definite matrix.
    - `k::Int`: The size of the identity matrix.
"""
function scalar_mat(X, k = size(X, 1))
    return tr(X)*I(k)
end

@register_symbolic scalar_mat(X::Union{Symbolics.Arr, Matrix{Num}}, k::Int)

add_gdcprule(scalar_mat, SymmetricPositiveDefinite, Positive, GConvex, GIncreasing)

add_gdcprule(LinearAlgebra.diag, SymmetricPositiveDefinite, Positive, GConvex, GIncreasing)

# """
#     pinching(X, Ps)

# Pinching of a symmetric positive definite matrix `X` by a set of symmetric positive definite matrices `Ps` is defined as `sum(Ps[i]*X*Ps[i])`.

# # Arguments
#     - `X::Matrix`: A symmetric positive definite matrix.
#     - `Ps::Vector`: A vector of symmetric positive definite matrices.
# """
# function pinching(X, Ps)
#     return sum(Ps[i]*X*Ps[i] for i in eachindex(Ps); dims = 1)
# end

@register_symbolic pinching(X::Matrix{Num}, Ps::Vector{Union{Symbolics.Arr, Matrix{Num}}})

add_gdcprule(pinching, SymmetricPositiveDefinite, Positive, GConvex, GIncreasing)

"""
    sdivergence(X, Y)

Symmetric divergence of two symmetric positive definite matrices `X` and `Y` is defined as `logdet((X+Y)/2) - 1/2*logdet(X*Y)`.

# Arguments
    - `X::Matrix`: A symmetric positive definite matrix.
    - `Y::Matrix`: A symmetric positive definite matrix.
"""
function sdivergence(X, Y)
    return logdet((X+Y)/2) - 1/2*logdet(X*Y)
end

@register_symbolic sdivergence(X::Matrix{Num}, Y::Matrix)
add_gdcprule(sdivergence, SymmetricPositiveDefinite, Positive, GConvex, GIncreasing)

@register_symbolic Manifolds.distance(M::Manifolds.SymmetricPositiveDefinite, X::AbstractMatrix, Y::Union{Symbolics.Arr, Matrix{Num}})
add_gdcprule(Manifolds.distance, SymmetricPositiveDefinite, Positive, GConvex, GIncreasing)

# @register_symbolic LinearAlgebra.exp(X::Union{Symbolics.Arr, Matrix{Num}})
# add_gdcprule(exp, SymmetricPositiveDefinite, Positive, GConvex, GIncreasing)

# add_gdcprule(sqrt, SymmetricPositiveDefinite, Positive, GConvex, GIncreasing)

add_gdcprule(SymbolicAnalysis.quad_form, SymmetricPositiveDefinite, Positive, GConvex, GIncreasing)

add_gdcprule(LinearAlgebra.eigmax, SymmetricPositiveDefinite, Positive, GConvex, GIncreasing)

"""
    log_quad_form(y, X)

Log of the quadratic form of a symmetric positive definite matrix `X` and a vector `y` is defined as `log(y'*X*y)`.

# Arguments
    - `y::Vector`: A vector.
    - `X::Matrix`: A symmetric positive definite matrix.
"""
function log_quad_form(y, X)
    return log(y'*X*y)
end

@register_symbolic log_quad_form(y::Vector, X::Union{Symbolics.Arr, Matrix{Num}})
add_gdcprule(log_quad_form, SymmetricPositiveDefinite, Positive, GConvex, GIncreasing)

add_gdcprule(inv, SymmetricPositiveDefinite, Positive, GConvex, GDecreasing)

# add_gdcprule(diag)