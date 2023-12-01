include("../src/SymbolicsDCP.jl")
using .SymbolicsDCP
using Symbolics, LogExpFunctions
using Symbolics: unwrap

@syms x y

ex = -1*xlogx(x)
ex = propagate_curvature(propagate_sign(ex))
getcurvature(ex)
getsign(ex)

ex = 2*abs(x) -1
ex = propagate_curvature(propagate_sign(ex))
getcurvature(ex)
getsign(ex)

ex = abs(x)^2
ex = propagate_curvature(propagate_sign(ex))
getcurvature(ex)
getsign(ex)

ex = abs(x)^2 + abs(x)^3
ex = propagate_curvature(propagate_sign(ex))
getcurvature(ex)
getsign(ex)

@variables x[1:3] y
ex = SymbolicsDCP.quad_over_lin(x .- y, 1 - y) |> unwrap
ex = propagate_curvature(propagate_sign(ex))
getcurvature(ex)
getsign(ex)

ex = exp.(x) |> unwrap
ex = propagate_curvature(propagate_sign(ex))
getcurvature(ex)
getsign(ex)

##vector * scalar gets simplified