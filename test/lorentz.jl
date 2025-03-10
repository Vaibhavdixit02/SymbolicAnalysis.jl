using Manifolds, Symbolics, SymbolicAnalysis, LinearAlgebra
using Test
using Symbolics: unwrap
using SymbolicAnalysis: propagate_sign, propagate_curvature, propagate_gcurvature

@testset "Lorentz Manifold" begin
    # Create a Lorentz manifold of dimension 2 (3D ambient space)
    M = Lorentz(2)
    
    # Define symbolic variables
    @variables p[1:3]
    
    # Test lorentz_distance
    q = [0.0, 0.0, 1.0]  # A point on the Lorentz model
    ex = Manifolds.distance(M, q, p) |> unwrap
    ex = propagate_sign(ex)
    ex = propagate_gcurvature(ex, M)
    @test SymbolicAnalysis.getgcurvature(ex) == SymbolicAnalysis.GConvex
    
    # Test analyze function
    analyze_res = analyze(ex, M)
    @test analyze_res.gcurvature == SymbolicAnalysis.GConvex
    
    # Test lorentz_log_barrier
    ex = SymbolicAnalysis.lorentz_log_barrier(p) |> unwrap
    ex = propagate_sign(ex)
    ex = propagate_gcurvature(ex, M)
    @test SymbolicAnalysis.getgcurvature(ex) == SymbolicAnalysis.GConvex
    
    # Test lorentz_homogeneous_quadratic
    A = [2.0 0.0 0.0; 0.0 2.0 0.0; 0.0 0.0 1.0]  # Positive definite matrix
    ex = SymbolicAnalysis.lorentz_homogeneous_quadratic(A, p) |> unwrap
    ex = propagate_sign(ex)
    ex = propagate_gcurvature(ex, M)
    @test SymbolicAnalysis.getgcurvature(ex) == SymbolicAnalysis.GConvex
    
    # Test lorentz_homogeneous_diagonal
    a = [2.0, 2.0, 1.0]  # min(a[1:2]) + a[3] = 2 + 1 = 3 ≥ 0
    ex = SymbolicAnalysis.lorentz_homogeneous_diagonal(a, p) |> unwrap
    ex = propagate_sign(ex)
    ex = propagate_gcurvature(ex, M)
    @test SymbolicAnalysis.getgcurvature(ex) == SymbolicAnalysis.GConvex
    
    # Test lorentz_least_squares
    # Create a valid A and b that satisfy the geodesic convexity condition
    A = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    b = [0.0, 0.0, -1.0]  # (b'A)_{d+1} = -1, and ‖b'A‖₂ = 1 ≤ sqrt(1/2) * 1
    ex = SymbolicAnalysis.lorentz_least_squares(A, b, p) |> unwrap
    ex = propagate_sign(ex)
    ex = propagate_gcurvature(ex, M)
    @test SymbolicAnalysis.getgcurvature(ex) == SymbolicAnalysis.GConvex
    
    # Test composition of functions
    ex = 2.0 * Manifolds.distance(M, q, p) + 
         SymbolicAnalysis.lorentz_log_barrier(p) |> unwrap
    ex = propagate_sign(ex)
    ex = propagate_gcurvature(ex, M)
    @test SymbolicAnalysis.getgcurvature(ex) == SymbolicAnalysis.GConvex
    
    # Test lorentz_transform (should preserve geodesic convexity)
    # Create a Lorentz boost in the x-direction
    cosh_phi = 1.2
    sinh_phi = sqrt(cosh_phi^2 - 1)
    O = [1.0 0.0 0.0; 0.0 cosh_phi -sinh_phi; 0.0 -sinh_phi cosh_phi]
    
    # Create a compound expression with the transform
    q_transformed = O * q
    ex = Manifolds.distance(M, q_transformed, p) |> unwrap
    ex = propagate_sign(ex)
    ex = propagate_gcurvature(ex, M)
    @test SymbolicAnalysis.getgcurvature(ex) == SymbolicAnalysis.GConvex
end