using SymbolicAnalysis, Test
using SafeTestsets

@testset "DCP" begin
    include("test.jl")
end

@testset "DGCP" begin
    include("dgp.jl")
end
