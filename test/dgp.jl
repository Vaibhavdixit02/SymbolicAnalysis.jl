using Manifolds, Symbolics, SymbolicAnalysis, LinearAlgebra
using LinearAlgebra, PDMats
using Symbolics: unwrap
using Test, Zygote, ForwardDiff

@variables X[1:5, 1:5]

M = Manifolds.SymmetricPositiveDefinite(5)

A = rand(5,5)
A = A*A'

@variables Sigma[1:5, 1:5]
xs = [rand(5) for i in 1:2]
Siginv = inv(Sigma)
ex = sum(SymbolicAnalysis.log_quad_form(x, Siginv) for x in xs) + 1/5*logdet(Sigma) |> Symbolics.unwrap
ex = SymbolicAnalysis.propagate_sign(ex)
ex = SymbolicAnalysis.propagate_gcurvature(ex)

@test SymbolicAnalysis.getgcurvature(ex) == SymbolicAnalysis.GConvex

##Brascamplieb Problem
ex = logdet(SymbolicAnalysis.conjugation(A, X)) - logdet(X) |> unwrap
ex = SymbolicAnalysis.propagate_sign(ex)
ex = SymbolicAnalysis.propagate_gcurvature(ex)

@test SymbolicAnalysis.getgcurvature(ex) == SymbolicAnalysis.GLinear

ex = SymbolicAnalysis.tr(SymbolicAnalysis.conjugation(A, X))

# using Convex

# X = Convex.Variable(5, 5)
# Y = Convex.Variable(5, 5)
# ex = sqrt(X*Y)
# vexity(ex)

## Karcher Mean
As = [rand(5,5) for i in 1:5]
As = [As[i]*As[i]' for i in 1:5]

ex = SymbolicAnalysis.sdivergence(X, As[1]) |> unwrap
ex = SymbolicAnalysis.propagate_sign(ex)
ex = SymbolicAnalysis.propagate_gcurvature(ex)

@test SymbolicAnalysis.getgcurvature(ex) == SymbolicAnalysis.GConvex

ex = sum(SymbolicAnalysis.sdivergence(X, As[i]) for i in 1:5) |> Symbolics.unwrap
ex = SymbolicAnalysis.propagate_sign(ex)
ex = SymbolicAnalysis.propagate_gcurvature(ex)

@test SymbolicAnalysis.getgcurvature(ex) == SymbolicAnalysis.GConvex

ex = SymbolicAnalysis.distance(M, As[1], X)^2 |> Symbolics.unwrap
ex = SymbolicAnalysis.propagate_sign(ex)
ex = SymbolicAnalysis.propagate_gcurvature(ex)

@test SymbolicAnalysis.getgcurvature(ex) == SymbolicAnalysis.GConvex

ex = sum(Manifolds.distance(M, As[i], X)^2 for i in 1:5) |> Symbolics.unwrap
ex = SymbolicAnalysis.propagate_sign(ex)
ex = SymbolicAnalysis.propagate_gcurvature(ex)

@test SymbolicAnalysis.getgcurvature(ex) == SymbolicAnalysis.GConvex

@variables Y[1:5, 1:5]
ex = sqrt(X*Y) |> unwrap
ex = SymbolicAnalysis.propagate_sign(ex)
@test_throws SymbolicUtils.RuleRewriteError SymbolicAnalysis.propagate_gcurvature(ex)

# ex = exp(X*Y) |> unwrap
# ex = SymbolicAnalysis.propagate_sign(ex)
# @test_throws SymbolicUtils.RuleRewriteError SymbolicAnalysis.propagate_gcurvature(ex)


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
using Optimization, OptimizationManopt, Symbolics, Manifolds, Random, LinearAlgebra, SymbolicAnalysis

M = SymmetricPositiveDefinite(5)
m = 100
σ = 0.005
q = Matrix{Float64}(LinearAlgebra.I(5)) .+ 2.0

# f(M, x, p = nothing) = sum(SymbolicAnalysis.distance(M, data2[i], x)^2 for i in 1:m)


# optf = OptimizationFunction(f, Optimization.AutoForwardDiff())
# opt = OptimizationManopt.GradientDescentOptimizer()
# prob = OptimizationProblem(optf, data2[2]; manifold = M)
# sol = solve(prob, opt)

M = SymmetricPositiveDefinite(5)
@variables X[1:5, 1:5]
data2 = [exp(M, q, σ * rand(M; vector_at=q)) for i in 1:m];

f(x, p = nothing) = sum(SymbolicAnalysis.distance(M, data2[i], x)^2 for i in 1:5)
optf = OptimizationFunction(f, Optimization.AutoZygote(); expr = prob.f.expr, sys = prob.f.sys)
prob = OptimizationProblem(optf, data2[1]; manifold = M)

opt = OptimizationManopt.GradientDescentOptimizer()
@time sol = solve(prob, opt, maxiters = 1000)
@test sol.objective < 1e-2

M = SymmetricPositiveDefinite(5)
xs = [rand(5) for i in 1:5]

f(S, p = nothing) = 1/length(xs) * sum(SymbolicAnalysis.log_quad_form(x, S) for x in xs) + 1/5*logdet(inv(S))

optf = OptimizationFunction(f, Optimization.AutoZygote(); expr = prob.f.expr, sys = prob.f.sys)
prob = OptimizationProblem(optf, Array{Float64}(LinearAlgebra.I(5)); manifold = M)

opt = OptimizationManopt.GradientDescentOptimizer()
sol = solve(prob, opt, maxiters = 1000)


A = randn(5,5) #initialize random matrix
A = A*A' #make it a SPD matrix

function matsqrt(X, p = nothing)
    return SymbolicAnalysis.sdivergence(X, A) + SymbolicAnalysis.sdivergence(X, Matrix{Float64}(LinearAlgebra.I(5))) 
end

optf = OptimizationFunction(matsqrt, Optimization.AutoForwardDiff())
prob = OptimizationProblem(optf, A/2, manifold = M)
sol = solve(prob, GradientDescentOptimizer(), maxiters = 1000)

sqrt(A) ≈ sol.minimizer

ex = matsqrt(X) |> unwrap
ex = SymbolicAnalysis.propagate_sign(ex)
ex = SymbolicAnalysis.propagate_gcurvature(ex)

@test SymbolicAnalysis.getgcurvature(ex) == SymbolicAnalysis.GConvex

##Diagonal loading
@variables X[1:5, 1:5]

ex = tr(inv(X)) + logdet(X) |> unwrap
ex = SymbolicAnalysis.propagate_sign(ex)
ex = SymbolicAnalysis.propagate_gcurvature(ex)
@test SymbolicAnalysis.getgcurvature(ex) == SymbolicAnalysis.GConvex

γ = 1/2
ex = (tr(X + γ*I(5)))^(2) |> unwrap
ex = SymbolicAnalysis.propagate_sign(ex)
ex = SymbolicAnalysis.propagate_gcurvature(ex)
@test SymbolicAnalysis.getgcurvature(ex) == SymbolicAnalysis.GConvex

d = 10
n = 50
@variables X[1:d, 1:d]

