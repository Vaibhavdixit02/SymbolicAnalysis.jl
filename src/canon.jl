function canonize(ex)
    rs = [@rule adjoint(~x) * ~Y * ~x => quad_form(x, Y)]
    rc = SymbolicUtils.Chain(rs)
    ex =  SymbolicUtils.Postwalk(rc)(ex)
    # ex = Prewalk(rc)(ex)
    SymbolicUtils.inspect(ex, metadata = true)
    return ex
end