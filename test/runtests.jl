using SymbolicAnalysis, Test
using SafeTestsets, Pkg

@testset "SymbolicAnalysis.jl" begin
    include("test.jl")
    include("dgp.jl")
end
