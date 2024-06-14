function canonize(ex)
    rs = [@rule (adjoint(~x) * (~Y * ~x))[1] => quad_form(~x, ~Y)
          @rule (adjoint(~B) * ~X)* ~B => conjugation(~X, ~B)]
    rc = SymbolicUtils.Chain(rs)
    # ex =  SymbolicUtils.Postwalk(rc)(ex)
    ex = SymbolicUtils.Prewalk(rc)(ex)
    SymbolicUtils.inspect(ex, metadata = true)
    return ex
end