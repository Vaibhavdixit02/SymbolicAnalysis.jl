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

using Convex

X = Convex.Variable(5, 5)
ex = logdet(A'*X*A) - logdet(X)
vexity(ex)

## Karcher Mean
As = [rand(5,5) for i in 1:5]
As = [As[i]*As[i]' for i in 1:5]

ex = SymbolicAnalysis.sdivergence(X, As[1])
ex = SymbolicAnalysis.propagate_gcurvature(ex)

SymbolicAnalysis.getgcurvature(ex)

ex = sum(SymbolicAnalysis.sdivergence(X, As[i]) for i in 1:5) |> unwrap
ex = SymbolicAnalysis.propagate_gcurvature(ex)

SymbolicAnalysis.getgcurvature(ex)