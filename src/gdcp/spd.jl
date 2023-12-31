@register_symbolic LinearAlgebra.logdet(X::Symbolics.Arr)

add_gdcprule(logdet, SymmetricPositiveDefinite, GPositive, GLinear, GIncreasing)

function conjugation(X, B)
    return B'*X*B
end

@register_symbolic conjugation(X::Symbolics.Arr, B::Symbolics.Arr)

add_gdcprule(conjugation, SymmetricPositiveDefinite, GPositive, GLinear, GIncreasing)

add_gdcprule(tr, SymmetricPositiveDefinite, GPositive, GLinear, GIncreasing)

add_gdcprule(sum, SymmetricPositiveDefinite, GPositive, GLinear, GIncreasing)

function scalar_mat(X, k = size(X, 1))
    return tr(X)*I(k)
end

@register_symbolic scalar_mat(X::Symbolics.Arr, k::Int)

add_gdcprule(scalar_mat, SymmetricPositiveDefinite, GPositive, GLinear, GIncreasing)

add_gdcprule(diag, SymmetricPositiveDefinite, GPositive, GLinear, GIncreasing)

function pinching(X, Ps)
    return sum(Ps[i]*X*Ps[i] for i in eachindex(Ps); dims = 1)
end

@register_symbolic pinching(X::Symbolics.Arr, Ps::Vector{Symbolics.Arr})

add_gdcprule(pinching, SymmetricPositiveDefinite, GPositive, GLinear, GIncreasing)
