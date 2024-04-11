using Manifolds, Symbolics, SymbolicAnalysis
using LinearAlgebra, PDMats
using Symbolics: unwrap
using Test

@variables X[1:5, 1:5]

M = Manifolds.SymmetricPositiveDefinite(5)

A = rand(5,5)
A = A*A'

##Brascamplieb Problem
ex = logdet(SymbolicAnalysis.conjugation(A', X)) - logdet(X)
ex = SymbolicAnalysis.propagate_gcurvature(ex)

@test SymbolicAnalysis.getgcurvature(ex) == SymbolicAnalysis.GLinear

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

@test SymbolicAnalysis.getgcurvature(ex) == SymbolicAnalysis.GVex

ex = sum(SymbolicAnalysis.sdivergence(X, As[i]) for i in 1:5) |> Symbolics.unwrap
ex = SymbolicAnalysis.propagate_sign(ex)
ex = SymbolicAnalysis.propagate_gcurvature(ex)

@test SymbolicAnalysis.getgcurvature(ex) == SymbolicAnalysis.GVex

ex = SymbolicAnalysis.distance(M, As[1], X)^2 |> Symbolics.unwrap
ex = SymbolicAnalysis.propagate_sign(ex)
ex = SymbolicAnalysis.propagate_gcurvature(ex)

@test SymbolicAnalysis.getgcurvature(ex) == SymbolicAnalysis.GVex

ex = sum(Manifolds.distance(M, As[i], X)^2 for i in 1:5) |> Symbolics.unwrap
ex = SymbolicAnalysis.propagate_sign(ex)
ex = SymbolicAnalysis.propagate_gcurvature(ex)

@test SymbolicAnalysis.getgcurvature(ex) == SymbolicAnalysis.GVex

@variables Y[1:5, 1:5]
ex = sqrt(X*Y) |> unwrap
ex = SymbolicAnalysis.propagate_sign(ex)
ex = SymbolicAnalysis.propagate_gcurvature(ex)

@test SymbolicAnalysis.getgcurvature(ex) == SymbolicAnalysis.GCave

ex = exp(X*Y)
ex = SymbolicAnalysis.propagate_sign(ex)
ex = SymbolicAnalysis.propagate_gcurvature(ex)

@test SymbolicAnalysis.getgcurvature(ex) == SymbolicAnalysis.GVex

# using Manopt, Manifolds, Random, LinearAlgebra, ManifoldDiff
# using ManifoldDiff: grad_distance, prox_distance
# Random.seed!(42);

# m = 100
# σ = 0.005
# q = Matrix{Float64}(I, 5, 5) .+ 2.0
# data2 = [exp(M, q, σ * rand(M; vector_at=q)) for i in 1:m];

# f(M, x) = sum(distance(M, x, data2[i])^2 for i in 1:m)
# f(x) = sum(distance(M, x, data2[i])^2 for i in 1:m)

# using FiniteDifferences

# r_backend = ManifoldDiff.RiemannianProjectionBackend(
#     ManifoldDiff.FiniteDifferencesBackend()
# )
# gradf1_FD(M, p) = ManifoldDiff.gradient(M, f, p, r_backend)

# m1 = gradient_descent(M, f, gradf1_FD, data2[1]; maxiter=1000)

# ################################
using Optimization, ModelingToolkit, OptimizationManopt, Manifolds, Random, LinearAlgebra, SymbolicAnalysis

# M = SymmetricPositiveDefinite(5)
m = 100
σ = 0.005
q = Matrix{Float64}(I, 5, 5) .+ 2.0

# f(M, x, p = nothing) = sum(SymbolicAnalysis.distance(M, data2[i], x)^2 for i in 1:m)


# optf = OptimizationFunction(f, Optimization.AutoForwardDiff())
# opt = OptimizationManopt.GradientDescentOptimizer()
# prob = OptimizationProblem(optf, data2[2]; manifold = M)
# sol = solve(prob, opt)

M = SymmetricPositiveDefinite(5)
@variables X[1:5, 1:5]
data2 = [exp(M, q, σ * rand(M; vector_at=q)) for i in 1:m];

obj = sum(SymbolicAnalysis.distance(M, data2[i], X)^2 for i in 1:5)
optsys = complete(OptimizationSystem(obj, X, [], name = :opt1))
prob = OptimizationProblem(optsys, data2[1])

f(x, p = nothing) = sum(SymbolicAnalysis.distance(M, data2[i], x)^2 for i in 1:5)
optf = OptimizationFunction(f, Optimization.AutoForwardDiff(); expr = prob.f.expr, sys = optsys)
prob = OptimizationProblem(optf, data2[1]; manifold = M)

opt = OptimizationManopt.GradientDescentOptimizer()
sol = solve(prob, opt, maxiters = 1000)
@test sol.objective < 1e-2

@variables Sigma[1:5, 1:5]
xs = [rand(5) for i in 1:3]
Siginv = inv(Sigma)
ex = sum(SymbolicAnalysis.quad_form(x, Siginv) for x in xs) + 1/5*logdet(Sigma) |> Symbolics.unwrap
ex = SymbolicAnalysis.propagate_sign(ex)
ex = SymbolicAnalysis.propagate_gcurvature(ex)

@test SymbolicAnalysis.getgcurvature(ex) == SymbolicAnalysis.GVex