using Manifolds
using Symbolics: @register_symbolic, unwrap
using LinearAlgebra

# @enum GSign GPositive GNegative GAnySign
@enum GCurvature GConvex GConcave GLinear GUnknownCurvature
@enum GMonotonicity GIncreasing GDecreasing GAnyMono

const gdcprules_dict = Dict()

function add_gdcprule(f, manifold, sign, curvature, monotonicity)
    if !(monotonicity isa Tuple)
        monotonicity = (monotonicity,)
    end
    gdcprules_dict[f] = makegrule(manifold, sign, curvature, monotonicity)
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

function mul_gcurvature(args)
    # all but one arg is constant
    # non_constants = findall(x->issym(x) || istree(x), args)
    # constants = findall(x->!issym(x) && !istree(x), args)
    # if !isempty(non_constants)
    #     expr = args[non_constants]
    #     curv = find_gcurvature.(expr)
    #     if !isempty(constants) && prod(args[constants]) < 0
    #         # flip
    #         if all(x -> x == GConvex, curv)
    #             return GConcave
    #         elseif all(x -> x == GConvex, curv)
    #             return GConvex
    #         elseif all(x -> x == GLinear, curv)
    #             return GLinear
    #         else
    #             return GUnknownCurvature
    #         end
    #     else
    #         if all(x -> x == GConvex, curv)
    #             return GConvex
    #         elseif all(x -> x == GConcave, curv)
    #             return GConcave
    #         elseif all(x -> x == GLinear, curv)
    #             return GLinear
    #         else
    #             GUnknownCurvature
    #         end
    #     end
    # end
    # return GLinear
    non_constants = findall(x->issym(x) || istree(x), args)
    constants = findall(x->!issym(x) && !istree(x), args)
    try
        @assert length(non_constants) <= 1
    catch
        @warn "DGCP does not support multiple non-constant arguments in multiplication"
        return UnknownGCurvature
    end
    if !isempty(non_constants)
        expr = args[first(non_constants)]
        curv = find_gcurvature(expr)
        return if prod(args[constants]) < 0
            # flip
            curv == GConvex ? GConcave : curv == GConcave ? GConvex : curv
        else
            curv
        end
    end
    return GLinear
end

function add_gcurvature(args)
    curvs = find_gcurvature.(args)
    all(==(GLinear), curvs) && return GLinear
    all(x->x==GConvex || x==GLinear, curvs) && return GConvex
    all(x->x==GConcave || x==GLinear, curvs) && return GConcave
    return GUnknownCurvature
end

function find_gcurvature(ex)
    if hasgcurvature(ex)
        return getgcurvature(ex)
    end
    if istree(ex)
        f, args = operation(ex), arguments(ex)
        @show f
        if f in keys(gdcprules_dict)
            rule, args = gdcprule(f, args...)
            f_curvature = rule.gcurvature
            f_monotonicity = rule.gmonotonicity
        elseif f in keys(dcprules_dict) || Symbol(f) == :^
            rule, args = dcprule(f, args...)
            f_curvature = rule.curvature
            f_monotonicity = rule.monotonicity
        elseif Symbol(f) == :*
            if args[1] isa Number && args[1] > 0
                return find_gcurvature(args[2])
            elseif args[1] isa Number && args[1] < 0
                argscurv = find_gcurvature(args[2])
                if argscurv == GConvex
                    return GConcave
                elseif argscurv == GConcave
                    return GConvex
                else
                    argscurv
                end
            else
                @warn "DCP does not support multiple non-constant arguments in multiplication"
                return UnknownGCurvature
            end
        else
            return GUnknownCurvature
        end
        @show ex
        @show f_curvature
        if f_curvature == Convex || f_curvature == Affine
            if all(enumerate(args)) do (i, arg)
                    arg_curv = find_gcurvature(arg)
                    m = get_arg_property(f_monotonicity, i, args)
                    @show arg
                    if arg_curv == GConvex
                        m == Increasing
                    elseif arg_curv == GConcave
                        m == Decreasing
                    else
                        arg_curv == GLinear
                    end
                end
                return GConvex
            end
        elseif f_curvature == Concave || f_curvature == Affine
            if all(enumerate(args)) do (i, arg)
                    arg_curv = find_gcurvature(arg)
                    m = f_monotonicity[i]
                    if arg_curv == GConcave
                        m == Increasing
                    elseif arg_curv == GConvex
                        m == Decreasing
                    else
                        arg_curv == GLinear
                    end
                end
                return GConcave
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
         @rule *(~~x) => setgcurvature(~MATCH, mul_gcurvature(~~x))
         @rule +(~~x) => setgcurvature(~MATCH, add_gcurvature(~~x))
         @rule ~x => setgcurvature(~x, find_gcurvature(~x))
        ]
    ex= Postwalk(RestartedChain(r))(ex)
    ex = Prewalk(RestartedChain(r))(ex)
    return ex
end
