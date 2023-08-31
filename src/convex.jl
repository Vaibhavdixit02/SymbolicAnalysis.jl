using Symbolics
using DomainSets
using LinearAlgebra

import Symbolics: Symbolic, issym, istree
using Symbolics.Rewriters

using SymbolicUtils
SymbolicUtils.inspect_metadata[] = true

@enum Sign Positive Negative AnySign
@enum Curvature Vex Cave Affine UnknownCurvature
@enum Monotonicity Increasing Decreasing AnyMono

# ispositive(x) --> we'll see

function increasing_if_positive(x)
    sign = getsign(x)
    sign == AnySign ?
    AnyMono : sign == Positive ? Increasing : Decreasing
end

const dcprules_dict = Dict()

function add_dcprule(f, domain, sign, curvature, monotonicity)
    if !(monotonicity isa Tuple)
        monotonicity = (monotonicity,)
    end
    dcprules_dict[f] = makerule(domain, sign, curvature, monotonicity)
end
makerule(domain, sign, curvature, monotonicity) = (domain=domain,
                sign=sign,
                curvature=curvature,
                monotonicity=monotonicity)

hasdcprule(f::Function) = haskey(dcprules_dict, f)
hasdcprule(f) = false
dcprule(f, args...) = dcprules_dict[f]

# special cases which depend on arguments:
function dcprule(::typeof(^), x, i)
    if isone(i)
        return makerule(ℝ, Positive, Vex, Increasing)
    elseif isinteger(i) && iseven(i)
        return makerule(ℝ, Positive, Vex, increasing_if_positive)
    elseif isinteger(i) && isodd(i)
        return makerule(ℝ, Positive, Vex, Increasing)
    elseif i >= 1
        return makerule(HalfLine(), Positive, Vex, Increasing)
    elseif i >= 0 && i < 1
        return makerule(HalfLine(), Positive, Cave, Increasing)
    elseif i < 0
        return makerule(HalfLine{Float64, :open}(), Positive, Cave, Increasing)
    end
end
hasdcprule(::typeof(^)) = true

struct CustomDomain{T} <: Domain{T}
    in::Function
end

Base.in(x, c::CustomDomain) = c.in(x)

function array_domain(element_domain)
    CustomDomain{AbstractArray}() do xs
        all(in(element_domain), xs)
    end
end

function array_domain(element_domain, N)
    CustomDomain{AbstractArray{<:Any, N}}() do xs
        ndims(xs) == N && all(in(element_domain), xs)
    end
end

add_dcprule(abs, ℝ, Positive, Vex, increasing_if_positive)
# add_dcprule(entropy, HalfLine(), AnySign, Cave, AnyMono)
add_dcprule(exp, ℝ, Positive, Vex, Increasing)
#add_dcprule(geomean, array_domain(HalfLine(),1), Positive, Vex, Increasing, vector=true)
# add_dcprule(huber, ℝ, Positive, Vex, increasing_if_positive)
add_dcprule(inv, HalfLine{Float64, :open}(), Positive, Vex, increasing_if_positive)
#add_dcprule(kl_div, (HalfLine{Float64, :open}(),
#                 HalfLine{Float64, :open}()), Positive, Vex, AnyMono)
add_dcprule(log, HalfLine{Float64, :open}(), AnySign, Cave, AnyMono)
#add_dcprule(log_sum_exp, array_domain(ℝ,1), AnySign, Cave, AnyMono)
add_dcprule(maximum, array_domain(ℝ,1), AnySign, Vex, AnyMono)
add_dcprule(minimum, array_domain(ℝ,1), AnySign, Cave, Increasing)
add_dcprule(norm, array_domain(ℝ,1), AnySign, Cave, increasing_if_positive)
#add_dcprule(positive, ℝ, Positive, Vex, Increasing)
#add_dcprule(^, ℝ, Positive, Vex, Increasing) # Requires special handling based on 2nd arg
#add_dcprule(quad_over_lin, (ℝ, HalfLine{Float64, :open}()), Positive, Vex, (increasing_if_positive, Decreasing)) # Requires special handling based on 2nd arg
add_dcprule(sqrt, HalfLine(), Positive, Cave, Increasing)


### Sign ###
#
setsign(ex::Symbolic, sign) = setmetadata(ex, Sign, sign)
setsign(ex, sign) = ex
getsign(ex::Symbolic) = getmetadata(ex, Sign)
getsign(ex) = ex < 0 ? Negative : Positive
hassign(ex::Symbolic) = hasmetadata(ex, Sign)
hassign(ex) = ex isa Real

function add_sign(args)
    signs = getsign.(args)
    if any(==(AnySign), signs)
        AnySign
    elseif all(==(Negative), signs)
        Negative
    elseif all(==(Positive), signs)
        Positive
    else
        AnySign
    end
end

function mul_sign(args)
    signs = getsign.(args)
    if any(==(AnySign), signs)
        AnySign
    elseif isodd(count(==(Negative), signs))
        Negative
    else
        Positive
    end
end

function propagate_sign(ex)
    # Step 1: set the sign of all variables to be AnySign
    r = @rule ~x::issym => hassign(~x) ? ~x : setsign(~x, AnySign)
    ex = Postwalk(PassThrough(r))(ex)

    # Step 2: set the sign of primitve functions
    r = @rule ~x::istree  =>
    setsign(~x, (println(~x); dcprule(operation(~x), arguments(~x)...).sign)) where
        {hasdcprule(operation(~x))}

    ex = Postwalk(PassThrough(r))(ex)

    # Step 3: propagate the sign to top level
    rs = [@rule *(~~x) => setsign(~MATCH, mul_sign(~~x))
          @rule +(~~x) => setsign(~MATCH, add_sign(~~x))]

    Postwalk(Chain(rs))(ex)
end

### Curvature ###
#
setcurvature(ex::Symbolic, curv) = setmetadata(ex, Curvature, curv)
setcurvature(ex, curv) = ex
getcurvature(ex::Symbolic) = getmetadata(ex, Curvature)
getcurvature(ex) = Affine
hascurvature(ex::Symbolic) = hasmetadata(ex, Curvature)
hascurvature(ex) = ex isa Real

function mul_curvature(args)
    # all but one arg is constant
    non_constants = findall(x->issym(x) || istree(x), args)
    constants = findall(x->!issym(x) && !istree(x), args)
    @assert length(non_constants) <= 1
    if !isempty(non_constants)
        expr = args[first(non_constants)]
        curv = find_curvature(expr)
        return if prod(args[constants]) < 0
            # flip
            curv == Vex ? Cave : curv == Cave ? Vex : curv
        else
            curv
        end
    end
    return Affine
end

function add_curvature(args)
    curvs = find_curvature.(args)
    all(==(Affine), curvs) && return Affine
    all(x->x==Vex || x==Affine, curvs) && return Vex
    all(x->x==Cave || x==Affine, curvs) && return Cave
    return UnknownCurvature
end

function propagate_curvature(ex)
    r = [@rule *(~~x) => setcurvature(~MATCH, mul_curvature(~~x))
         @rule +(~~x) => setcurvature(~MATCH, add_curvature(~~x))
         @rule ~x => setcurvature(~x, find_curvature(~x))]
    Postwalk(RestartedChain(r))(ex)
end

function get_arg_property(monotonicity, i, args)
    @label start
    if monotonicity isa Function
        monotonicity(args[i])
    elseif monotonicity isa Tuple
        monotonicity = monotonicity[i]
        @goto start
    else
        monotonicity
    end
end

function find_curvature(ex)
    if hascurvature(ex)
        return getcurvature(ex)
    end

    if istree(ex)
        f, args = operation(ex), arguments(ex)
        rule = dcprule(f, args...)
        f_curvature = rule.curvature
        f_monotonicity = rule.monotonicity

        if f_curvature == Vex || f_curvature == Affine
            if all(enumerate(args)) do (i, arg)
                    arg_curv = find_curvature(arg)
                    m = get_arg_property(f_monotonicity, i, args)
                    if arg_curv == Vex
                        m == Increasing
                    elseif arg_curv == Cave
                        m == Decreasing
                    else
                        arg_curv == Affine
                    end
                end
                return Vex
            end
        elseif f_curvature == Cave || f_curvature == Affine
            if all(enumerate(args)) do (i, arg)
                    arg_curv = find_curvature(arg)
                    m = f_monotonicity[i]
                    if arg_curv == Cave
                        m == Increasing
                    elseif arg_curv == Vex
                        m == Decreasing
                    else
                        arg_curv == Affine
                    end
                end
                return Cave
            end
        elseif f_curvature == Affine
            if all(enumerate(args)) do (i, arg)
                    arg_curv = find_curvature(arg)
                    arg_curv == Affine
                end
                return Affine
            end
        end
        return UnknownCurvature
    else
        return Affine
    end
end
