# SymbolicAnalysis.jl
Symbolics-based function property propagation for optimization

This project was started as a proof of concept for building a replacement for Convex.jl using Symbolics.jl using its rewriting capabilities 
and metadata propagation for registered functions. I found doing the reimplementation of the conic representation very tedious and haven't 
finished it up. Still, this project will be useful in optimization for analyzing programs' convexity automatically in its current form.
I aim to expand it to be more generally useful in determining various program structures for automated solver selection and 
reformulation of problems.
