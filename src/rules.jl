
@enum Sign Positive Negative AnySign
@enum Curvature Vex Cave Affine UnknownCurvature
@enum Monotonicity Increasing Decreasing AnyMono

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

function symmetric_domain()
    CustomDomain{AbstractArray{<:Any, 2}}(issymmetric)
end

function semidefinite_domain()
    CustomDomain{AbstractArray{<:Any, 2}}(isposdef) #not semi so needs to change
end

function negsemidefinite_domain()
    CustomDomain{AbstractArray{<:Any, 2}}(isposdef ∘ -) #not semi so needs to change
end

function definite_domain()
    CustomDomain{AbstractArray{<:Any, 2}}(isposdef)
end

function negdefinite_domain()
    CustomDomain{AbstractArray{<:Any, 2}}(isposdef ∘ -)
end

function function_domain()
    CustomDomain{Function}(x -> typeassert(x, Function))
end

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
dcprule(f, args...) = dcprules_dict[f], args

### Sign ###
setsign(ex::Symbolic, sign) = setmetadata(ex, Sign, sign)
setsign(ex, sign) = ex
function getsign(ex::Symbolic)
    @show ex
    @show istree(ex)
    if issym(ex)
        return getmetadata(ex, Sign)
    elseif istree(ex)
        return getmetadata.(ex, Ref(Sign))
    end
end
getsign(ex::Number) = ex < 0 ? Negative : Positive
hassign(ex::Symbolic) = hasmetadata(ex, Sign)
hassign(ex) = ex isa Real

function add_sign(args)
    @show args
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
    r = @rule ~x::istree  => setsign(~x, (dcprule(operation(~x), arguments(~x)...)[1].sign)) where {hasdcprule(operation(~x))}
    
    ex = Postwalk(PassThrough(r))(ex)

    r = @rule ~x::istree  => setsign(~x, (gdcprule(operation(~x), arguments(~x)...)[1].sign)) where {hasgdcprule(operation(~x))}
    
    ex = Postwalk(PassThrough(r))(ex)

    SymbolicUtils.inspect(ex, metadata=true)
    # Step 3: propagate the sign to top level
    rs = [@rule *(~~x) => setsign(~MATCH, mul_sign(~~x))
          ]
    ex = Prewalk(Chain(rs))(ex)

    rs = [
          @rule +(~~x) => setsign(~MATCH, add_sign(~~x))
          ]
    ex = Postwalk(Chain(rs))(ex)
    # SymbolicUtils.inspect(ex)
    ex
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
        #  @rule broadcast(~f, ~~x) => setcurvature(~MATCH, propagate_curvature(propagate_sign(Symbolics.scalarize((~MATCH)[1]))))
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
    # SymbolicUtils.inspect(ex)
    if istree(ex)
        f, args = operation(ex), arguments(ex)
        rule, args = dcprule(f, args...)
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
