include("../src/convex.jl")

@variables x #only works as @variables because xlogx(::Number)
ex = -1*xlogx(x)
ex = propagate_curvature(propagate_sign(ex))
getcurvature(ex)
getsign(ex) #doesn't work with @variables

@syms x
ex = 2*abs(x) -1
ex = propagate_curvature(propagate_sign(ex))
getcurvature(ex)
getsign(ex)

ex = abs(x)^2
ex = propagate_curvature(propagate_sign(ex))
getcurvature(ex)
getsign(ex)

@variables x y
ex = StatsBase.kldivergence(x, y)