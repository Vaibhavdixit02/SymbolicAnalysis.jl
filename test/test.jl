using SymbolicAnalysis
using Symbolics, LogExpFunctions
using Symbolics: unwrap
using LinearAlgebra

@syms x y

ex1 = exp(x^2) - log(x) |> unwrap
ex1 = propagate_curvature(propagate_sign(ex1))
getcurvature(ex1)
getsign(ex1)

ex2 = -sqrt(x^2)
ex2 = propagate_curvature(propagate_sign(ex2))
getcurvature(ex2)
getsign(ex2)

# ex2 = -2*norm([1,x]) - x*(x-3) + y
ex2 = propagate_curvature(propagate_sign(ex2))
getcurvature(ex2)
getsign(ex2)


ex = -1*xlogx(x)
ex = propagate_curvature(propagate_sign(ex))
getcurvature(ex)
getsign(ex)

ex = 2*abs(x) -1
ex = propagate_curvature(propagate_sign(ex))
getcurvature(ex)
getsign(ex)

ex = abs(x)^2
ex = propagate_curvature(propagate_sign(ex))
getcurvature(ex)
getsign(ex)

ex = abs(x)^2 + abs(x)^3
ex = propagate_curvature(propagate_sign(ex))
getcurvature(ex)
getsign(ex)

@variables x[1:3] y
ex = x .- y
ex = propagate_curvature(propagate_sign(ex))
getcurvature(ex)
getsign(ex)

ex = exp.(x) |> unwrap
ex = propagate_curvature(propagate_sign(ex))
getcurvature(ex)
getsign(ex)

##vector * scalar gets simplified

@syms x y z
obj = x^2 + y^2 + z^2

ex = propagate_curvature(propagate_sign(obj))
getcurvature(ex)


cons = [
    x + y + z ~ 10
    log1p(x)^2 - log1p(z) ≲ 0
]

ex = propagate_curvature(propagate_sign(cons[1].lhs))
getcurvature(ex)

ex = propagate_curvature(propagate_sign(cons[2].lhs))
getcurvature(ex)

using Manifolds, Symbolics, SymbolicAnalysis

@variables x[1:5, 1:5]
using LinearAlgebra, PDMats
M = Manifolds.SymmetricPositiveDefinite(5)

A = rand(5,5)
A = A*A'
X = x*x'
# A == A'
# Manifolds.check_point(M, A)

using Symbolics: @register_symbolic, unwrap

@register_symbolic Manifolds.log(M::SymmetricPositiveDefinite, p::Symbolics.Arr, q::Symbolics.Arr)

Manifolds.log(M, A, X)

Asqrtinv = sqrt(inv(A))


ex = norm(log(Asqrtinv*X*Asqrtinv), 2)

ex = propagate_curvature(propagate_sign(ex))
getcurvature(ex)

ex = logdet(SymbolicAnalysis.conjugation(A', X)) - logdet(X)
ex = SymbolicAnalysis.propagate_gcurvature(ex)

SymbolicAnalysis.getgcurvature(ex)

ex = logdet(A'*X*A) - logdet(X)
ex = SymbolicAnalysis.propagate_curvature(ex)
SymbolicAnalysis.getcurvature(ex)

using Convex

X = Convex.Variable(5, 5)
ex = logdet(A'*X*A) - logdet(X)
vexity(ex)

