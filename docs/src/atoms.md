# Atoms

## DCP atoms

| Atom | Curvature | Sign | Monotonicity | Domain |
| :---: | :---: | :---: | :---: | :---: |
| `+` | `Affine` | `AnySign` | `Increasing` | `ℝ` |
| `Base.Ref` | `Affine` | `AnySign` | `Increasing` | `ℝ` |
| `dot` | `Affine` | `AnySign` | `Increasing` | `(array_domain(ℝ), array_domain(ℝ))` |
| `dotsort` | `Convex` | `AnySign` | `(AnyMono, increasing_if_positive ∘ minimum)` | `(array_domain(ℝ,1), array_domain(ℝ,1))` |
| `StatsBase.geomean` | `Concave` | `Positive` | `Increasing` | `array_domain(HalfLine{Number, :open}(),1)` |
| `StatsBase.harmmean` | `Concave` | `Positive` | `Increasing` | `array_domain(HalfLine{Number, :open}(),1)` |
| `invprod` | `Convex` | `Positive` | `Decreasing` | `array_domain(HalfLine{Number, :open}())` |
| `eigmax` | `Convex` | `AnySign` | `AnyMono` | `symmetric_domain()` |
| `eigmin` | `Concave` | `AnySign` | `AnyMono` | `symmetric_domain()` |
| `eigsummax` | `Convex` | `AnySign` | `AnyMono` | `(array_domain(ℝ, 2), ℝ)` |
| `eigsummin` | `Concave` | `AnySign` | `AnyMono` | `(array_domain(ℝ, 2), ℝ)` |
| `logdet` | `Concave` | `AnySign` | `AnyMono` | `semidefinite_domain()` |
| `LogExpFunctions.logsumexp` | `Convex` | `AnySign` | `Increasing` | `array_domain(ℝ,2)` |
| `matrix_frac` | `Convex` | `AnySign` | `AnyMono` | `(array_domain(ℝ,1), definite_domain())` |
| `maximum` | `Convex` | `AnySign` | `Increasing` | `array_domain(ℝ)` |
| `minimum` | `Concave` | `AnySign` | `Increasing` | `array_domain(ℝ)` |
| `norm` | `Convex` | `Positive` | `increasing_if_positive` | `(array_domain(ℝ), Interval{:closed, :open}(1, Inf))` |
| `norm` | `Convex` | `Positive` | `increasing_if_positive` | `(array_domain(ℝ), Interval{:closed, :open}(0, 1))` |
| `perspective` | `getcurvature` | `getsign` | `AnyMono` | `(function_domain(), ℝ, Positive)` |
| `quad_form` | `Convex` | `Positive` | `(increasing_if_positive, Increasing)` | `(array_domain(ℝ,1), semidefinite_domain())` |
| `quad_over_lin` | `Convex` | `Positive` | `(increasing_if_positive, Decreasing)` | `(array_domain(ℝ), HalfLine{Number, :open}())` |
| `quad_over_lin` | `Convex` | `Positive` | `(increasing_if_positive, Decreasing)` | `(ℝ, HalfLine{Number, :open}())` |
| `sum` | `Affine` | `AnySign` | `Increasing` | `array_domain(ℝ, 2)` |
| `sum_largest` | `Convex` | `AnySign` | `Increasing` | `(array_domain(ℝ,2), ℤ)` |
| `sum_smallest` | `Concave` | `AnySign` | `Increasing` | `(array_domain(ℝ,2), ℤ)` |
| `tr` | `Affine` | `AnySign` | `Increasing` | `array_domain(ℝ, 2)` |
| `trinv` | `Convex` | `Positive` | `AnyMono` | `definite_domain()` |
| `tv` | `Convex` | `Positive` | `AnyMono` | `array_domain(ℝ,1)` |
| `tv` | `Convex` | `Positive` | `AnyMono` | `array_domain(array_domain(ℝ,2), 1)` |
| `abs` | `Convex` | `Positive` | `increasing_if_positive` | `ℂ` |
| `conj` | `Affine` | `AnySign` | `AnyMono` | `ℂ` |
| `exp` | `Convex` | `Positive` | `Increasing` | `ℝ` |
| `xlogx` | `Convex` | `AnySign` | `AnyMono` | `ℝ` |
| `huber` | `Convex` | `Positive` | `increasing_if_positive` | `(ℝ, HalfLine())` |
| `imag` | `Affine` | `AnySign` | `AnyMono` | `ℂ` |
| `inv` | `Convex` | `Positive` | `Decreasing` | `HalfLine{Number, :open}()` |
| `log` | `Concave` | `AnySign` | `Increasing` | `HalfLine{Number, :open}()` |
| `log` | `Concave` | `Positive` | `Increasing` | `array_domain(ℝ, 2)` |
| `inv` | `Convex` | `AnySign` | `Decreasing` | `semidefinite_domain()` |
| `sqrt` | `Concave` | `Positive` | `Increasing` | `semidefinite_domain()` |
| `kldivergence` | `Convex` | `Positive` | `AnyMono` | `(array_domain(HalfLine{Number, :open},1), array_domain(HalfLine{Number, :open},1))` |
| `lognormcdf` | `Concave` | `Negative` | `Increasing` | `ℝ` |
| `log1p` | `Concave` | `Negative` | `Increasing` | `Interval{:open, :open}(-1, Inf)` |
| `logistic` | `Convex` | `Positive` | `Increasing` | `ℝ` |
| `max` | `Convex` | `AnySign` | `Increasing` | `(ℝ, ℝ)` |
| `min` | `Concave` | `AnySign` | `Increasing` | `(ℝ, ℝ)` |
| `^` | See special cases in code | See special cases in code | See special cases in code | See special cases in code |
| `real` | `Affine` | `AnySign` | `Increasing` | `ℂ` |
| `rel_entr` | `Convex` | `AnySign` | `(AnyMono, Decreasing)` | `(HalfLine{Number, :open}(), HalfLine{Number, :open}())` |
| `sqrt` | `Concave` | `Positive` | `Increasing` | `HalfLine()`

## DGCP

### Symmetric Positive Definite atoms

| Atom | Geodesic Curvature | Sign | Monotonicity |
| --- | --- | --- | --- |
| `*` | `GLinear` | `Positive` | `GIncreasing` |
| `LinearAlgebra.logdet` | `GLinear` | `Positive` | `GIncreasing` |
| `conjugation` | `GConvex` | `Positive` | `GIncreasing` |
| `LinearAlgebra.tr` | `GConvex` | `Positive` | `GIncreasing` |
| `sum` | `GConvex` | `Positive` | `GIncreasing` |
| `adjoint` | `GLinear` | `Positive` | `GIncreasing` |
| `scalar_mat` | `GConvex` | `Positive` | `GIncreasing` |
| `LinearAlgebra.diag` | `GConvex` | `Positive` | `GIncreasing` |
| `pinching` | `GConvex` | `Positive` | `GIncreasing` |
| `sdivergence` | `GConvex` | `Positive` | `GIncreasing` |
| `Manifolds.distance` | `GConvex` | `Positive` | `GIncreasing` |
| `exp` | `GConvex` | `Positive` | `GIncreasing` |
| `sqrt` | `GConvex` | `Positive` | `GIncreasing` |