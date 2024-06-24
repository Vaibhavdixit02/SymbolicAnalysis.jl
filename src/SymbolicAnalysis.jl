module SymbolicAnalysis

using DomainSets
using LinearAlgebra
using LogExpFunctions
using StatsBase
using Distributions
using DSP, DataStructures

using Symbolics
import Symbolics: Symbolic, issym, Term
using SymbolicUtils: iscall
using SymbolicUtils.Rewriters
SymbolicUtils.inspect_metadata[] = true

struct VarDomain end

include("rules.jl")
include("atoms.jl")
include("gdcp/gdcp_rules.jl")
include("gdcp/spd.jl")
include("canon.jl")

struct AnalysisResult
    curvature::SymbolicAnalysis.Curvature
    sign::SymbolicAnalysis.Sign
    gcurvature::Union{SymbolicAnalysis.GCurvature,Nothing}
end

function analyze(ex, M::Union{AbstractManifold,Nothing} = nothing)
    ex = canonize(ex)
    ex = propagate_sign(ex)
    ex = propagate_curvature(ex)
    if isnothing(M)
        return AnalysisResult(getcurvature(ex), getsign(ex), nothing)
    else
        ex = propagate_gcurvature(ex, M)
        return AnalysisResult(getcurvature(ex), getsign(ex), getgcurvature(ex))
    end
end

export propagate_curvature,
    propagate_sign, getcurvature, getsign, propagate_gcurvature, analyze
end
