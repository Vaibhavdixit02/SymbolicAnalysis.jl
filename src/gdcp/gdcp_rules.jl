using Manifolds
using Symbolics: @register_symbolic, unwrap
using LinearAlgebra

# @enum GSign GPositive GNegative GAnySign
@enum GCurvature GVex GCave GLinear GUnknownCurvature
@enum GMonotonicity GIncreasing GDecreasing GAnyMono

const gdcprules_dict = Dict()

function add_gdcprule(f, manifold, sign, curvature, monotonicity)
    if !(monotonicity isa Tuple)
        monotonicity = (monotonicity,)
    end
    gdcprules_dict[f] = makerule(manifold, sign, curvature, monotonicity)
end
makegrule(manifold, sign, curvature, monotonicity) = (manifold=manifold,
                sign=sign,
                gcurvature=curvature,
                gmonotonicity=monotonicity)

hasgdcprule(f::Function) = haskey(gdcprules_dict, f)
hasgdcprule(f) = false
gdcprule(f, args...) = gdcprules_dict[f], args

setgcurvature(ex::Union{Symbolic, Num}, curv) = setmetadata(ex, GCurvature, curv)
setgcurvature(ex, curv) = ex
getgcurvature(ex::Union{Symbolic, Num}) = getmetadata(ex, GCurvature)
getgcurvature(ex) = GLinear
hasgcurvature(ex::Union{Symbolic, Num}) = hasmetadata(ex, GCurvature)
hasgcurvature(ex) = ex isa Real

# function mul_gcurvature(args)
#     # all but one arg is constant
#     non_constants = findall(x->issym(x) || istree(x), args)
#     constants = findall(x->!issym(x) && !istree(x), args)
#     if !isempty(non_constants)
#         expr = args[non_constants]
#         curv = find_gcurvature(expr)
#         return if !isempty(constants) && prod(args[constants]) < 0
#             # flip
#             curv == GVex ? GCave : curv == GCave ? GVex : curv
#         else
#             curv
#         end
#     end
#     return GLinear
# end

function add_gcurvature(args)
    curvs = find_gcurvature.(args)
    all(==(GLinear), curvs) && return GLinear
    all(x->x==GVex || x==GLinear, curvs) && return GVex
    all(x->x==GCave || x==GLinear, curvs) && return GCave
    return GUnknownCurvature
end

function find_gcurvature(ex)
    # @show ex
    # @show hasgcurvature(ex)
    if hasgcurvature(ex)
        return getgcurvature(ex)
    end
    # @show istree(ex)
    if istree(ex)
        f, args = operation(ex), arguments(ex)
        if f in keys(gdcprules_dict)
            rule, args = gdcprule(f, args...)
        else
            rule, args = dcprule(f, args...)
        end
        f_curvature = rule.curvature
        f_monotonicity = rule.monotonicity

        if f_curvature == Vex || f_curvature == Affine
            if all(enumerate(args)) do (i, arg)
                    arg_curv = find_gcurvature(arg)
                    m = get_arg_property(f_monotonicity, i, args)
                    if arg_curv == GVex
                        m == Increasing
                    elseif arg_curv == GCave
                        m == Decreasing
                    else
                        arg_curv == GLinear
                    end
                end
                return GVex
            end
        elseif f_curvature == Cave || f_curvature == Affine
            if all(enumerate(args)) do (i, arg)
                    arg_curv = find_gcurvature(arg)
                    m = f_monotonicity[i]
                    if arg_curv == GCave
                        m == Increasing
                    elseif arg_curv == GVex
                        m == Decreasing
                    else
                        arg_curv == GLinear
                    end
                end
                return GCave
            end
        elseif f_curvature == Affine
            if all(enumerate(args)) do (i, arg)
                    arg_curv = find_gcurvature(arg)
                    arg_curv == GLinear
                end
                return GLinear
            end
        elseif f_curvature isa GCurvature
            return f_curvature
        else
            return GUnknownCurvature
        end
    elseif hasfield(typeof(ex), :val) && operation(ex.val) in keys(gdcprules_dict)
        f, args = operation(ex.val), arguments(ex.val)
        rule, args = gdcprule(f, args...)
        return rule.curvature
    else
        return GLinear
    end
end

function propagate_gcurvature(ex)
    r = [
         @rule +(~~x) => setgcurvature(~MATCH, add_gcurvature(~~x))
         @rule ~x => setgcurvature(~x, find_gcurvature(~x))
        #  @rule *(~~x) => setgcurvature(~MATCH, mul_gcurvature(~~x))
        ]
    return Postwalk(RestartedChain(r))(ex)
end
