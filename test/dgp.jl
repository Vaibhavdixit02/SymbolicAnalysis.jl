using Manifolds, Symbolics, SymbolicAnalysis, LinearAlgebra
using LinearAlgebra, PDMats
using Symbolics: unwrap
using Test, Zygote, ForwardDiff

@variables X[1:5, 1:5]

M = Manifolds.SymmetricPositiveDefinite(5)

A = rand(5, 5)
A = A * A'

ex = SymbolicAnalysis.logdet(SymbolicAnalysis.conjugation(inv(X), A)) |> unwrap
ex = propagate_sign(ex)
ex = propagate_curvature(ex)
ex = propagate_gcurvature(ex, M)
SymbolicAnalysis.getcurvature(ex)
@test SymbolicAnalysis.getgcurvature(ex) == SymbolicAnalysis.GConvex

ex = SymbolicAnalysis.logdet(tr(inv(X))) |> unwrap
ex = propagate_sign(ex)
ex = propagate_curvature(ex)
ex = propagate_gcurvature(ex, M)
@test SymbolicAnalysis.getgcurvature(ex) == SymbolicAnalysis.GConvex
SymbolicAnalysis.getcurvature(ex)

@variables Sigma[1:5, 1:5]
xs = [rand(5) for i = 1:2]
ex =
    sum(SymbolicAnalysis.log_quad_form(x, inv(Sigma)) for x in xs) +
    1 / 5 * logdet(Sigma) |> Symbolics.unwrap
analyze_res = SymbolicAnalysis.analyze(ex, M)
@test analyze_res.gcurvature == SymbolicAnalysis.GConvex

##Brascamplieb Problem
M = SymmetricPositiveDefinite(5)
objective_expr = logdet(SymbolicAnalysis.conjugation(X, A)) - logdet(X) |> unwrap
objective_expr = SymbolicAnalysis.propagate_sign(objective_expr)
analyze_res = analyze(objective_expr, M)
@test analyze_res.gcurvature == SymbolicAnalysis.GConvex

objective_expr = SymbolicAnalysis.propagate_gcurvature(objective_expr, M)
@test SymbolicAnalysis.getgcurvature(objective_expr) == SymbolicAnalysis.GConvex

ex = SymbolicAnalysis.tr(SymbolicAnalysis.conjugation(X, A)) |> unwrap
ex = propagate_sign(ex)
ex = propagate_curvature(ex)
ex = propagate_gcurvature(ex, M)

@test analyze(ex, M).gcurvature == SymbolicAnalysis.GConvex

# using Convex

# X = Convex.Variable(5, 5)
# Y = Convex.Variable(5, 5)
# ex = sqrt(X*Y)
# vexity(ex)

## Karcher Mean
As = [rand(5, 5) for i = 1:5]
As = [As[i] * As[i]' for i = 1:5]

ex = SymbolicAnalysis.sdivergence(X, As[1]) |> unwrap
ex = SymbolicAnalysis.propagate_sign(ex)
ex = SymbolicAnalysis.propagate_gcurvature(ex, M)

@test SymbolicAnalysis.getgcurvature(ex) == SymbolicAnalysis.GConvex

ex = sum(SymbolicAnalysis.sdivergence(X, As[i]) for i = 1:5) |> Symbolics.unwrap
ex = SymbolicAnalysis.propagate_sign(ex)
ex = SymbolicAnalysis.propagate_gcurvature(ex, M)

@test SymbolicAnalysis.getgcurvature(ex) == SymbolicAnalysis.GConvex

ex = SymbolicAnalysis.distance(M, As[1], X)^2 |> Symbolics.unwrap
ex = SymbolicAnalysis.propagate_sign(ex)
ex = SymbolicAnalysis.propagate_gcurvature(ex, M)

@test SymbolicAnalysis.getgcurvature(ex) == SymbolicAnalysis.GConvex

M = SymmetricPositiveDefinite(5)
objective_expr = sum(Manifolds.distance(M, As[i], X)^2 for i = 1:5) |> Symbolics.unwrap
analyze_res = analyze(objective_expr, M)
@test analyze_res.gcurvature == SymbolicAnalysis.GConvex

@variables Y[1:5, 1:5]
ex = sqrt(X * Y) |> unwrap
ex = SymbolicAnalysis.propagate_sign(ex)
@test_throws SymbolicUtils.RuleRewriteError SymbolicAnalysis.propagate_gcurvature(ex, M)

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
using Optimization,
    OptimizationManopt, Symbolics, Manifolds, Random, LinearAlgebra, SymbolicAnalysis

M = SymmetricPositiveDefinite(5)
m = 100
σ = 0.005
q = Matrix{Float64}(LinearAlgebra.I(5)) .+ 2.0

data2 = [exp(M, q, σ * rand(M; vector_at = q)) for i = 1:m];

f(x, p = nothing) = sum(SymbolicAnalysis.distance(M, data2[i], x)^2 for i = 1:5)
optf = OptimizationFunction(f, Optimization.AutoZygote())
prob = OptimizationProblem(optf, data2[1]; manifold = M)

opt = OptimizationManopt.GradientDescentOptimizer()
@time sol = solve(prob, opt, maxiters = 10)
@test sol.objective < 1e-2

M = SymmetricPositiveDefinite(5)
xs = [rand(5) for i = 1:5]

f(S, p = nothing) =
    1 / length(xs) * sum(SymbolicAnalysis.log_quad_form(x, S) for x in xs) +
    1 / 5 * logdet(inv(S))

optf =
    OptimizationFunction(f, Optimization.AutoZygote(); expr = prob.f.expr, sys = prob.f.sys)
prob = OptimizationProblem(optf, Array{Float64}(LinearAlgebra.I(5)); manifold = M)

opt = OptimizationManopt.GradientDescentOptimizer()
sol = solve(prob, opt, maxiters = 10)


A = randn(5, 5) #initialize random matrix
A = A * A' #make it a SPD matrix

function matsqrt(X, p = nothing) #setup objective function
    return SymbolicAnalysis.sdivergence(X, A) +
           SymbolicAnalysis.sdivergence(X, Matrix{Float64}(LinearAlgebra.I(5)))
end

optf = OptimizationFunction(matsqrt, Optimization.AutoZygote()) #setup oracles
prob = OptimizationProblem(optf, A / 2, manifold = M) #setup problem with manifold and initial point

sol = solve(prob, GradientDescentOptimizer()) #solve the problem
@test sqrt(A) ≈ sol.minimizer

ex = matsqrt(X) |> unwrap
ex = SymbolicAnalysis.propagate_sign(ex)
ex = SymbolicAnalysis.propagate_gcurvature(ex, M)

@test SymbolicAnalysis.getgcurvature(ex) == SymbolicAnalysis.GConvex

##Diagonal loading
@variables X[1:5, 1:5]

ex = tr(inv(X)) + logdet(X) |> unwrap
@test analyze(ex, M).gcurvature == SymbolicAnalysis.GConvex

γ = 1 / 2
ex = (tr(X + γ * I(5)))^(2) |> unwrap
@test analyze(ex, M).gcurvature == SymbolicAnalysis.GConvex

d = 10
n = 50
@variables X[1:d, 1:d]

@variables x[1:5] X[1:5, 1:5]
ex = SymbolicAnalysis.log_quad_form(x, inv(X)) |> unwrap
ex = SymbolicAnalysis.propagate_sign(ex)
ex = SymbolicAnalysis.propagate_gcurvature(ex, M)
@test SymbolicAnalysis.getgcurvature(ex) == SymbolicAnalysis.GConvex
