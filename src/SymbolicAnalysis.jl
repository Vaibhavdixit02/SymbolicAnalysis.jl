module SymbolicAnalysis

using Symbolics
using DomainSets
using LinearAlgebra
using LogExpFunctions
using StatsBase
using Distributions
using DSP, DataStructures

import Symbolics: Symbolic, issym, istree, Term
using Symbolics.Rewriters
using SymbolicUtils
SymbolicUtils.inspect_metadata[] = true

struct VarDomain end

include("rules.jl")
include("atoms.jl")
include("gdcp/gdcp_rules.jl")
include("gdcp/spd.jl")

export propagate_curvature, propagate_sign, getcurvature, getsign
end