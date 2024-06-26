# SymbolicAnalysis.jl

Symbolics-based function property propagation for optimization

SymbolicAnalysis is a package for implementing the Disciplined Programming approach to optimization. As demonstrated by the [DCP framework](https://dcp.stanford.edu/), and further followups to it for further classes of functions such as DGP, DQP etc. symbolic representation of problems can be leveraged to certify and construct convex (or similar function properties) expressions.

This package aims to utilize expression graph rewriting and metadata propagation supported by Symbolics.jl, to support propagation of several of these properties - limited right now to Euclidean Convexity and Geodesic Convexity on the Symmetric Positive Definite manifold. This package provides an easy to expand implementation of "atoms" that are functions that have known properties. This allows users to add atoms to the library more easily than the previous implementations [CVXPY](https://www.cvxpy.org/index.html) and [Convex.jl](https://github.com/jump-dev/Convex.jl).

## Installation

To install this package, run the following in the Julia REPL:

```julia
using Pkg
Pkg.add("SymbolicAnalysis")
```

## Usage

The main interface to this package is the `analyze` function.

```@autodocs
Modules=[SymbolicAnalysis]
Pages=["SymbolicAnalysis.jl"]
```
