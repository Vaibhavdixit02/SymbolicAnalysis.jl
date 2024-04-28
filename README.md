# SymbolicAnalysis.jl

Symbolics-based function property propagation for optimization

SymbolicAnalysis is a package for implementating the Disciplined Programming approach to optimization,
As demonstrated by the [DCP framework](https://dcp.stanford.edu/), and further followups to it for further classes of
functions https://www.cvxpy.org/tutorial/index.html such as DGP, DQP etc, symbolic representation of problems can be leveraged
to idenstify and facilitate building Convex (or similar function properties) expressions.

This package aims to utlize expression graph rewriting and metadata propagation supported by Symbolics.jl, to support 
propagation of several of these properties - limited right now to Euclidean Convexity and Geodesic Convexity on the Symmetric 
Positive Definite manifold. This package provides an easier to expand implementation of functional properties than the previous 
implementations [CVXPY](https://www.cvxpy.org/index.html) and [Convex.jl](https://github.com/jump-dev/Convex.jl) as well as a 
more performant implementation of the function property propagation. 

Manifold optimization is supported through the Manopt.jl wrapper of Optimization.jl https://github.com/SciML/Optimization.jl that
allows utilizing the Geodesic Convexity propagation to be utilized as a certificate of global optimality with riemannian optimization.
