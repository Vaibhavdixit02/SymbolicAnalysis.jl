import{_ as t,c as e,o as d,a6 as c}from"./chunks/framework.Te2tIaPL.js";const v=JSON.parse('{"title":"Atoms","description":"","frontmatter":{},"headers":[],"relativePath":"atoms.md","filePath":"atoms.md","lastUpdated":null}'),n={name:"atoms.md"},o=c('<h1 id="Atoms" tabindex="-1">Atoms <a class="header-anchor" href="#Atoms" aria-label="Permalink to &quot;Atoms {#Atoms}&quot;">​</a></h1><h2 id="DCP-atoms" tabindex="-1">DCP atoms <a class="header-anchor" href="#DCP-atoms" aria-label="Permalink to &quot;DCP atoms {#DCP-atoms}&quot;">​</a></h2><table><thead><tr><th style="text-align:center;">Atom</th><th style="text-align:center;">Curvature</th><th style="text-align:center;">Sign</th><th style="text-align:center;">Monotonicity</th><th style="text-align:center;">Domain</th></tr></thead><tbody><tr><td style="text-align:center;"><code>+</code></td><td style="text-align:center;"><code>Affine</code></td><td style="text-align:center;"><code>AnySign</code></td><td style="text-align:center;"><code>Increasing</code></td><td style="text-align:center;"><code>ℝ</code></td></tr><tr><td style="text-align:center;"><code>Base.Ref</code></td><td style="text-align:center;"><code>Affine</code></td><td style="text-align:center;"><code>AnySign</code></td><td style="text-align:center;"><code>Increasing</code></td><td style="text-align:center;"><code>ℝ</code></td></tr><tr><td style="text-align:center;"><code>dot</code></td><td style="text-align:center;"><code>Affine</code></td><td style="text-align:center;"><code>AnySign</code></td><td style="text-align:center;"><code>Increasing</code></td><td style="text-align:center;"><code>(array_domain(ℝ), array_domain(ℝ))</code></td></tr><tr><td style="text-align:center;"><code>dotsort</code></td><td style="text-align:center;"><code>Convex</code></td><td style="text-align:center;"><code>AnySign</code></td><td style="text-align:center;"><code>(AnyMono, increasing_if_positive ∘ minimum)</code></td><td style="text-align:center;"><code>(array_domain(ℝ,1), array_domain(ℝ,1))</code></td></tr><tr><td style="text-align:center;"><code>StatsBase.geomean</code></td><td style="text-align:center;"><code>Concave</code></td><td style="text-align:center;"><code>Positive</code></td><td style="text-align:center;"><code>Increasing</code></td><td style="text-align:center;"><code>array_domain(HalfLine{Number, :open}(),1)</code></td></tr><tr><td style="text-align:center;"><code>StatsBase.harmmean</code></td><td style="text-align:center;"><code>Concave</code></td><td style="text-align:center;"><code>Positive</code></td><td style="text-align:center;"><code>Increasing</code></td><td style="text-align:center;"><code>array_domain(HalfLine{Number, :open}(),1)</code></td></tr><tr><td style="text-align:center;"><code>invprod</code></td><td style="text-align:center;"><code>Convex</code></td><td style="text-align:center;"><code>Positive</code></td><td style="text-align:center;"><code>Decreasing</code></td><td style="text-align:center;"><code>array_domain(HalfLine{Number, :open}())</code></td></tr><tr><td style="text-align:center;"><code>eigmax</code></td><td style="text-align:center;"><code>Convex</code></td><td style="text-align:center;"><code>AnySign</code></td><td style="text-align:center;"><code>AnyMono</code></td><td style="text-align:center;"><code>symmetric_domain()</code></td></tr><tr><td style="text-align:center;"><code>eigmin</code></td><td style="text-align:center;"><code>Concave</code></td><td style="text-align:center;"><code>AnySign</code></td><td style="text-align:center;"><code>AnyMono</code></td><td style="text-align:center;"><code>symmetric_domain()</code></td></tr><tr><td style="text-align:center;"><code>eigsummax</code></td><td style="text-align:center;"><code>Convex</code></td><td style="text-align:center;"><code>AnySign</code></td><td style="text-align:center;"><code>AnyMono</code></td><td style="text-align:center;"><code>(array_domain(ℝ, 2), ℝ)</code></td></tr><tr><td style="text-align:center;"><code>eigsummin</code></td><td style="text-align:center;"><code>Concave</code></td><td style="text-align:center;"><code>AnySign</code></td><td style="text-align:center;"><code>AnyMono</code></td><td style="text-align:center;"><code>(array_domain(ℝ, 2), ℝ)</code></td></tr><tr><td style="text-align:center;"><code>logdet</code></td><td style="text-align:center;"><code>Concave</code></td><td style="text-align:center;"><code>AnySign</code></td><td style="text-align:center;"><code>AnyMono</code></td><td style="text-align:center;"><code>semidefinite_domain()</code></td></tr><tr><td style="text-align:center;"><code>LogExpFunctions.logsumexp</code></td><td style="text-align:center;"><code>Convex</code></td><td style="text-align:center;"><code>AnySign</code></td><td style="text-align:center;"><code>Increasing</code></td><td style="text-align:center;"><code>array_domain(ℝ,2)</code></td></tr><tr><td style="text-align:center;"><code>matrix_frac</code></td><td style="text-align:center;"><code>Convex</code></td><td style="text-align:center;"><code>AnySign</code></td><td style="text-align:center;"><code>AnyMono</code></td><td style="text-align:center;"><code>(array_domain(ℝ,1), definite_domain())</code></td></tr><tr><td style="text-align:center;"><code>maximum</code></td><td style="text-align:center;"><code>Convex</code></td><td style="text-align:center;"><code>AnySign</code></td><td style="text-align:center;"><code>Increasing</code></td><td style="text-align:center;"><code>array_domain(ℝ)</code></td></tr><tr><td style="text-align:center;"><code>minimum</code></td><td style="text-align:center;"><code>Concave</code></td><td style="text-align:center;"><code>AnySign</code></td><td style="text-align:center;"><code>Increasing</code></td><td style="text-align:center;"><code>array_domain(ℝ)</code></td></tr><tr><td style="text-align:center;"><code>norm</code></td><td style="text-align:center;"><code>Convex</code></td><td style="text-align:center;"><code>Positive</code></td><td style="text-align:center;"><code>increasing_if_positive</code></td><td style="text-align:center;"><code>(array_domain(ℝ), Interval{:closed, :open}(1, Inf))</code></td></tr><tr><td style="text-align:center;"><code>norm</code></td><td style="text-align:center;"><code>Convex</code></td><td style="text-align:center;"><code>Positive</code></td><td style="text-align:center;"><code>increasing_if_positive</code></td><td style="text-align:center;"><code>(array_domain(ℝ), Interval{:closed, :open}(0, 1))</code></td></tr><tr><td style="text-align:center;"><code>perspective</code></td><td style="text-align:center;"><code>getcurvature</code></td><td style="text-align:center;"><code>getsign</code></td><td style="text-align:center;"><code>AnyMono</code></td><td style="text-align:center;"><code>(function_domain(), ℝ, Positive)</code></td></tr><tr><td style="text-align:center;"><code>quad_form</code></td><td style="text-align:center;"><code>Convex</code></td><td style="text-align:center;"><code>Positive</code></td><td style="text-align:center;"><code>(increasing_if_positive, Increasing)</code></td><td style="text-align:center;"><code>(array_domain(ℝ,1), semidefinite_domain())</code></td></tr><tr><td style="text-align:center;"><code>quad_over_lin</code></td><td style="text-align:center;"><code>Convex</code></td><td style="text-align:center;"><code>Positive</code></td><td style="text-align:center;"><code>(increasing_if_positive, Decreasing)</code></td><td style="text-align:center;"><code>(array_domain(ℝ), HalfLine{Number, :open}())</code></td></tr><tr><td style="text-align:center;"><code>quad_over_lin</code></td><td style="text-align:center;"><code>Convex</code></td><td style="text-align:center;"><code>Positive</code></td><td style="text-align:center;"><code>(increasing_if_positive, Decreasing)</code></td><td style="text-align:center;"><code>(ℝ, HalfLine{Number, :open}())</code></td></tr><tr><td style="text-align:center;"><code>sum</code></td><td style="text-align:center;"><code>Affine</code></td><td style="text-align:center;"><code>AnySign</code></td><td style="text-align:center;"><code>Increasing</code></td><td style="text-align:center;"><code>array_domain(ℝ, 2)</code></td></tr><tr><td style="text-align:center;"><code>sum_largest</code></td><td style="text-align:center;"><code>Convex</code></td><td style="text-align:center;"><code>AnySign</code></td><td style="text-align:center;"><code>Increasing</code></td><td style="text-align:center;"><code>(array_domain(ℝ,2), ℤ)</code></td></tr><tr><td style="text-align:center;"><code>sum_smallest</code></td><td style="text-align:center;"><code>Concave</code></td><td style="text-align:center;"><code>AnySign</code></td><td style="text-align:center;"><code>Increasing</code></td><td style="text-align:center;"><code>(array_domain(ℝ,2), ℤ)</code></td></tr><tr><td style="text-align:center;"><code>tr</code></td><td style="text-align:center;"><code>Affine</code></td><td style="text-align:center;"><code>AnySign</code></td><td style="text-align:center;"><code>Increasing</code></td><td style="text-align:center;"><code>array_domain(ℝ, 2)</code></td></tr><tr><td style="text-align:center;"><code>trinv</code></td><td style="text-align:center;"><code>Convex</code></td><td style="text-align:center;"><code>Positive</code></td><td style="text-align:center;"><code>AnyMono</code></td><td style="text-align:center;"><code>definite_domain()</code></td></tr><tr><td style="text-align:center;"><code>tv</code></td><td style="text-align:center;"><code>Convex</code></td><td style="text-align:center;"><code>Positive</code></td><td style="text-align:center;"><code>AnyMono</code></td><td style="text-align:center;"><code>array_domain(ℝ,1)</code></td></tr><tr><td style="text-align:center;"><code>tv</code></td><td style="text-align:center;"><code>Convex</code></td><td style="text-align:center;"><code>Positive</code></td><td style="text-align:center;"><code>AnyMono</code></td><td style="text-align:center;"><code>array_domain(array_domain(ℝ,2), 1)</code></td></tr><tr><td style="text-align:center;"><code>abs</code></td><td style="text-align:center;"><code>Convex</code></td><td style="text-align:center;"><code>Positive</code></td><td style="text-align:center;"><code>increasing_if_positive</code></td><td style="text-align:center;"><code>ℂ</code></td></tr><tr><td style="text-align:center;"><code>conj</code></td><td style="text-align:center;"><code>Affine</code></td><td style="text-align:center;"><code>AnySign</code></td><td style="text-align:center;"><code>AnyMono</code></td><td style="text-align:center;"><code>ℂ</code></td></tr><tr><td style="text-align:center;"><code>exp</code></td><td style="text-align:center;"><code>Convex</code></td><td style="text-align:center;"><code>Positive</code></td><td style="text-align:center;"><code>Increasing</code></td><td style="text-align:center;"><code>ℝ</code></td></tr><tr><td style="text-align:center;"><code>xlogx</code></td><td style="text-align:center;"><code>Convex</code></td><td style="text-align:center;"><code>AnySign</code></td><td style="text-align:center;"><code>AnyMono</code></td><td style="text-align:center;"><code>ℝ</code></td></tr><tr><td style="text-align:center;"><code>huber</code></td><td style="text-align:center;"><code>Convex</code></td><td style="text-align:center;"><code>Positive</code></td><td style="text-align:center;"><code>increasing_if_positive</code></td><td style="text-align:center;"><code>(ℝ, HalfLine())</code></td></tr><tr><td style="text-align:center;"><code>imag</code></td><td style="text-align:center;"><code>Affine</code></td><td style="text-align:center;"><code>AnySign</code></td><td style="text-align:center;"><code>AnyMono</code></td><td style="text-align:center;"><code>ℂ</code></td></tr><tr><td style="text-align:center;"><code>inv</code></td><td style="text-align:center;"><code>Convex</code></td><td style="text-align:center;"><code>Positive</code></td><td style="text-align:center;"><code>Decreasing</code></td><td style="text-align:center;"><code>HalfLine{Number, :open}()</code></td></tr><tr><td style="text-align:center;"><code>log</code></td><td style="text-align:center;"><code>Concave</code></td><td style="text-align:center;"><code>AnySign</code></td><td style="text-align:center;"><code>Increasing</code></td><td style="text-align:center;"><code>HalfLine{Number, :open}()</code></td></tr><tr><td style="text-align:center;"><code>log</code></td><td style="text-align:center;"><code>Concave</code></td><td style="text-align:center;"><code>Positive</code></td><td style="text-align:center;"><code>Increasing</code></td><td style="text-align:center;"><code>array_domain(ℝ, 2)</code></td></tr><tr><td style="text-align:center;"><code>inv</code></td><td style="text-align:center;"><code>Convex</code></td><td style="text-align:center;"><code>AnySign</code></td><td style="text-align:center;"><code>Decreasing</code></td><td style="text-align:center;"><code>semidefinite_domain()</code></td></tr><tr><td style="text-align:center;"><code>sqrt</code></td><td style="text-align:center;"><code>Concave</code></td><td style="text-align:center;"><code>Positive</code></td><td style="text-align:center;"><code>Increasing</code></td><td style="text-align:center;"><code>semidefinite_domain()</code></td></tr><tr><td style="text-align:center;"><code>kldivergence</code></td><td style="text-align:center;"><code>Convex</code></td><td style="text-align:center;"><code>Positive</code></td><td style="text-align:center;"><code>AnyMono</code></td><td style="text-align:center;"><code>(array_domain(HalfLine{Number, :open},1), array_domain(HalfLine{Number, :open},1))</code></td></tr><tr><td style="text-align:center;"><code>lognormcdf</code></td><td style="text-align:center;"><code>Concave</code></td><td style="text-align:center;"><code>Negative</code></td><td style="text-align:center;"><code>Increasing</code></td><td style="text-align:center;"><code>ℝ</code></td></tr><tr><td style="text-align:center;"><code>log1p</code></td><td style="text-align:center;"><code>Concave</code></td><td style="text-align:center;"><code>Negative</code></td><td style="text-align:center;"><code>Increasing</code></td><td style="text-align:center;"><code>Interval{:open, :open}(-1, Inf)</code></td></tr><tr><td style="text-align:center;"><code>logistic</code></td><td style="text-align:center;"><code>Convex</code></td><td style="text-align:center;"><code>Positive</code></td><td style="text-align:center;"><code>Increasing</code></td><td style="text-align:center;"><code>ℝ</code></td></tr><tr><td style="text-align:center;"><code>max</code></td><td style="text-align:center;"><code>Convex</code></td><td style="text-align:center;"><code>AnySign</code></td><td style="text-align:center;"><code>Increasing</code></td><td style="text-align:center;"><code>(ℝ, ℝ)</code></td></tr><tr><td style="text-align:center;"><code>min</code></td><td style="text-align:center;"><code>Concave</code></td><td style="text-align:center;"><code>AnySign</code></td><td style="text-align:center;"><code>Increasing</code></td><td style="text-align:center;"><code>(ℝ, ℝ)</code></td></tr><tr><td style="text-align:center;"><code>^</code></td><td style="text-align:center;">See special cases in code</td><td style="text-align:center;">See special cases in code</td><td style="text-align:center;">See special cases in code</td><td style="text-align:center;">See special cases in code</td></tr><tr><td style="text-align:center;"><code>real</code></td><td style="text-align:center;"><code>Affine</code></td><td style="text-align:center;"><code>AnySign</code></td><td style="text-align:center;"><code>Increasing</code></td><td style="text-align:center;"><code>ℂ</code></td></tr><tr><td style="text-align:center;"><code>rel_entr</code></td><td style="text-align:center;"><code>Convex</code></td><td style="text-align:center;"><code>AnySign</code></td><td style="text-align:center;"><code>(AnyMono, Decreasing)</code></td><td style="text-align:center;"><code>(HalfLine{Number, :open}(), HalfLine{Number, :open}())</code></td></tr><tr><td style="text-align:center;"><code>sqrt</code></td><td style="text-align:center;"><code>Concave</code></td><td style="text-align:center;"><code>Positive</code></td><td style="text-align:center;"><code>Increasing</code></td><td style="text-align:center;"><code>HalfLine()</code></td></tr></tbody></table><h2 id="DGCP" tabindex="-1">DGCP <a class="header-anchor" href="#DGCP" aria-label="Permalink to &quot;DGCP {#DGCP}&quot;">​</a></h2><h3 id="Symmetric-Positive-Definite-atoms" tabindex="-1">Symmetric Positive Definite atoms <a class="header-anchor" href="#Symmetric-Positive-Definite-atoms" aria-label="Permalink to &quot;Symmetric Positive Definite atoms {#Symmetric-Positive-Definite-atoms}&quot;">​</a></h3><table><thead><tr><th style="text-align:right;">Atom</th><th style="text-align:right;">Geodesic Curvature</th><th style="text-align:right;">Sign</th><th style="text-align:right;">Monotonicity</th></tr></thead><tbody><tr><td style="text-align:right;"><code>*</code></td><td style="text-align:right;"><code>GLinear</code></td><td style="text-align:right;"><code>Positive</code></td><td style="text-align:right;"><code>GIncreasing</code></td></tr><tr><td style="text-align:right;"><code>LinearAlgebra.logdet</code></td><td style="text-align:right;"><code>GLinear</code></td><td style="text-align:right;"><code>Positive</code></td><td style="text-align:right;"><code>GIncreasing</code></td></tr><tr><td style="text-align:right;"><code>conjugation</code></td><td style="text-align:right;"><code>GConvex</code></td><td style="text-align:right;"><code>Positive</code></td><td style="text-align:right;"><code>GIncreasing</code></td></tr><tr><td style="text-align:right;"><code>LinearAlgebra.tr</code></td><td style="text-align:right;"><code>GConvex</code></td><td style="text-align:right;"><code>Positive</code></td><td style="text-align:right;"><code>GIncreasing</code></td></tr><tr><td style="text-align:right;"><code>sum</code></td><td style="text-align:right;"><code>GConvex</code></td><td style="text-align:right;"><code>Positive</code></td><td style="text-align:right;"><code>GIncreasing</code></td></tr><tr><td style="text-align:right;"><code>adjoint</code></td><td style="text-align:right;"><code>GLinear</code></td><td style="text-align:right;"><code>Positive</code></td><td style="text-align:right;"><code>GIncreasing</code></td></tr><tr><td style="text-align:right;"><code>scalar_mat</code></td><td style="text-align:right;"><code>GConvex</code></td><td style="text-align:right;"><code>Positive</code></td><td style="text-align:right;"><code>GIncreasing</code></td></tr><tr><td style="text-align:right;"><code>LinearAlgebra.diag</code></td><td style="text-align:right;"><code>GConvex</code></td><td style="text-align:right;"><code>Positive</code></td><td style="text-align:right;"><code>GIncreasing</code></td></tr><tr><td style="text-align:right;"><code>pinching</code></td><td style="text-align:right;"><code>GConvex</code></td><td style="text-align:right;"><code>Positive</code></td><td style="text-align:right;"><code>GIncreasing</code></td></tr><tr><td style="text-align:right;"><code>sdivergence</code></td><td style="text-align:right;"><code>GConvex</code></td><td style="text-align:right;"><code>Positive</code></td><td style="text-align:right;"><code>GIncreasing</code></td></tr><tr><td style="text-align:right;"><code>Manifolds.distance</code></td><td style="text-align:right;"><code>GConvex</code></td><td style="text-align:right;"><code>Positive</code></td><td style="text-align:right;"><code>GIncreasing</code></td></tr><tr><td style="text-align:right;"><code>exp</code></td><td style="text-align:right;"><code>GConvex</code></td><td style="text-align:right;"><code>Positive</code></td><td style="text-align:right;"><code>GIncreasing</code></td></tr><tr><td style="text-align:right;"><code>sqrt</code></td><td style="text-align:right;"><code>GConvex</code></td><td style="text-align:right;"><code>Positive</code></td><td style="text-align:right;"><code>GIncreasing</code></td></tr></tbody></table>',6),i=[o];function l(r,a,s,g,y,x){return d(),e("div",null,i)}const h=t(n,[["render",l]]);export{v as __pageData,h as default};
