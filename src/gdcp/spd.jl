@register_symbolic LinearAlgebra.logdet(X::Symbolics.Arr)

add_gdcprule(logdet, SymmetricPositiveDefinite, GPositive, GLinear, GIncreasing)

function conjugation(X, B)
    return B'*X*B
end

@register_symbolic conjugation(X::Matrix, B::Symbolics.Arr)

add_gdcprule(conjugation, SymmetricPositiveDefinite, GPositive, GVex, GIncreasing)

add_gdcprule(tr, SymmetricPositiveDefinite, GPositive, GVex, GIncreasing)

add_gdcprule(sum, SymmetricPositiveDefinite, GPositive, GVex, GIncreasing)

add_gdcprule(adjoint, SymmetricPositiveDefinite, GPositive, GLinear, GIncreasing)

function scalar_mat(X, k = size(X, 1))
    return tr(X)*I(k)
end

@register_symbolic scalar_mat(X::Symbolics.Arr, k::Int)

add_gdcprule(scalar_mat, SymmetricPositiveDefinite, GPositive, GVex, GIncreasing)

add_gdcprule(diag, SymmetricPositiveDefinite, GPositive, GVex, GIncreasing)

function pinching(X, Ps)
    return sum(Ps[i]*X*Ps[i] for i in eachindex(Ps); dims = 1)
end

@register_symbolic pinching(X::Symbolics.Arr, Ps::Vector{Symbolics.Arr})

add_gdcprule(pinching, SymmetricPositiveDefinite, GPositive, GVex, GIncreasing)

function sdivergence(X, Y)
    return logdet((X+Y)/2) - 1/2*logdet(X) - 1/2*logdet(Y)
end

@register_symbolic sdivergence(X::Symbolics.Arr, Y::Matrix)
add_gdcprule(sdivergence, SymmetricPositiveDefinite, GPositive, GVex, GIncreasing)


@register_symbolic Manifolds.distance(M::Manifolds.SymmetricPositiveDefinite, X::AbstractMatrix, Y::Symbolics.Arr)
add_gdcprule(Manifolds.distance, SymmetricPositiveDefinite, GPositive, GVex, GIncreasing)
