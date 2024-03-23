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
# ex = sqrt(X*Y)
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

using Manopt, Manifolds, Random, LinearAlgebra, ManifoldDiff
using ManifoldDiff: grad_distance, prox_distance
Random.seed!(42);

m = 100
σ = 0.005
q = Matrix{Float64}(I, 5, 5) .+ 2.0
data2 = [exp(M, q, σ * rand(M; vector_at=q)) for i in 1:m];

f(M, x) = sum(distance(M, x, data2[i])^2 for i in 1:m)
f(x) = sum(distance(M, x, data2[i])^2 for i in 1:m)

using FiniteDifferences

r_backend = ManifoldDiff.RiemannianProjectionBackend(
    ManifoldDiff.FiniteDifferencesBackend()
)
gradf1_FD(M, p) = ManifoldDiff.gradient(M, f, p, r_backend)

m1 = gradient_descent(M, f, gradf1_FD, data2[1]; maxiter=1000)

################################
using Optimization, ModelingToolkit, OptimizationManopt, Manifolds, Random, LinearAlgebra

M = SymmetricPositiveDefinite(5)
m = 100
σ = 0.005
q = Matrix{Float64}(I, 5, 5) .+ 2.0
data2 = [exp(M, q, σ * rand(M; vector_at=q)) for i in 1:m];
# f(M, x, p = nothing) = sum(SymbolicAnalysis.distance(M, data2[i], x)^2 for i in 1:m)
# f(x, p = nothing) = sum(SymbolicAnalysis.distance(M, data2[i], x)^2 for i in 1:m)

# optf = OptimizationFunction(f, Optimization.AutoModelingToolkit())
# optprob = Optimization.instantiate_function(optf, data2[1], Optimization.AutoModelingToolkit(), nothing)
# opt = OptimizationManopt.GradientDescentOptimizer(M)
# prob = OptimizationProblem(, data2[1])
# sol = solve(prob, opt)

# @variables X[1:5, 1:5]
# obj = sum(Manifolds.distance(M, data2[i], X)^2 for i in 1:5) |> unwrap

# optsys = complete(OptimizationSystem(obj, X, [], name = :opt1))
# prob = OptimizationProblem(optsys, data2[1])
# opt = OptimizationManopt.NelderMeadOptimizer(M)
# sol = solve(prob, opt, maxiters = 1000)