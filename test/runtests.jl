using SymbolicAnalysis: propagate_curvature, propagate_sign, propagate_gcurvature
using SafeTestsets, Test

@testset "DCP" begin
    include("test.jl")
end

@testset "DGCP" begin
    include("dgp.jl")
end
