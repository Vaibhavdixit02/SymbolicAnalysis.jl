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

export propagate_curvature, propagate_sign, getcurvature, getsign
end