using Manifolds, Symbolics, SymbolicAnalysis
using LinearAlgebra, PDMats
using Symbolics: unwrap

@variables X[1:5, 1:5]

M = Manifolds.SymmetricPositiveDefinite(5)

A = rand(5,5)
A = A*A'

##Brascamplieb Problem
ex = logdet(SymbolicAnalysis.conjugation(A', X)) - logdet(X)
ex = SymbolicAnalysis.propagate_gcurvature(ex)

SymbolicAnalysis.getgcurvature(ex)

# using Convex

# X = Convex.Variable(5, 5)
# Y = Convex.Variable(5, 5)
# ex = exp(X'*Y)
# vexity(ex)

## Karcher Mean
As = [rand(5,5) for i in 1:5]
As = [As[i]*As[i]' for i in 1:5]

ex = SymbolicAnalysis.sdivergence(X, As[1])
ex = SymbolicAnalysis.propagate_gcurvature(ex)

SymbolicAnalysis.getgcurvature(ex)

ex = sum(SymbolicAnalysis.sdivergence(X, As[i]) for i in 1:5) |> unwrap
ex = SymbolicAnalysis.propagate_sign(ex)
ex = SymbolicAnalysis.propagate_gcurvature(ex)

SymbolicAnalysis.getgcurvature(ex)

ex = Manifolds.distance(M, As[1], X)^2 |> unwrap
ex = SymbolicAnalysis.propagate_sign(ex)
ex = SymbolicAnalysis.propagate_gcurvature(ex)

SymbolicAnalysis.getgcurvature(ex)

ex = sum(Manifolds.distance(M, As[i], X)^2 for i in 1:5) |> unwrap
ex = SymbolicAnalysis.propagate_sign(ex)
ex = SymbolicAnalysis.propagate_gcurvature(ex)

SymbolicAnalysis.getgcurvature(ex)

@variables Y[1:5, 1:5]
ex = sqrt(X*Y) |> unwrap
ex = SymbolicAnalysis.propagate_sign(ex)
ex = SymbolicAnalysis.propagate_gcurvature(ex)

SymbolicAnalysis.getgcurvature(ex)

ex = exp(X*Y)
ex = SymbolicAnalysis.propagate_sign(ex)
ex = SymbolicAnalysis.propagate_gcurvature(ex)

SymbolicAnalysis.getgcurvature(ex)