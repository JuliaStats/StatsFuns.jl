# import of Rmath functions

module Rmath

const rmathlib = "libRmath-julia"

### import macro

function _import_rmath(nparams::Int, rname::Symbol, jname::Symbol)
    # C function names
    if rname == :norm
        dfun = Expr(:quote, "dnorm4")
        pfun = Expr(:quote, "pnorm5")
        qfun = Expr(:quote, "qnorm5")
        rfun = Expr(:quote, "rnorm")
    else
        dfun = Expr(:quote, string('d', rname))    # density
        pfun = Expr(:quote, string('p', rname))    # cumulative probability
        qfun = Expr(:quote, string('q', rname))    # quantile
        rfun = Expr(:quote, string('r', rname))    # random sampling
    end

    # Julia function names
    pdf = symbol(string(jname, "pdf"))
    cdf = symbol(string(jname, "cdf"))
    ccdf = symbol(string(jname, "ccdf"))

    logpdf = symbol(string(jname, "logpdf"))
    logcdf = symbol(string(jname, "logcdf"))
    logccdf = symbol(string(jname, "logccdf"))

    invcdf = symbol(string(jname, "invcdf"))
    invccdf = symbol(string(jname, "invccdf"))
    invlogcdf = symbol(string(jname, "invlogcdf"))
    invlogccdf = symbol(string(jname, "invlogccdf"))

    rand = symbol(string(jname, "rand"))

    # arguments & argument types
    _pts = fill(:Cdouble, nparams)

    dtypes = Expr(:tuple, :Cdouble, _pts..., :Cint)
    ptypes = Expr(:tuple, :Cdouble, _pts..., :Cint, :Cint)
    qtypes = Expr(:tuple, :Cdouble, _pts..., :Cint, :Cint)
    rtypes = Expr(:tuple, _pts...)

    pargs = [symbol(string('p', i)) for i in 1:nparams]  # [:p1, :p2, ...]
    pdecls = [Expr(:(::), ps, :Real) for ps in pargs] # [:(p1::Real), :(p2::Real), ...]

    # function implementation
    quote
        $pdf($(pdecls...), x::Real) =
            ccall(($dfun, rmathlib), Float64, $dtypes, x, $(pargs...), 0)

        $logpdf($(pdecls...), x::Real) =
            ccall(($dfun, rmathlib), Float64, $dtypes, x, $(pargs...), 1)

        $cdf($(pdecls...), x::Real) =
            ccall(($pfun, rmathlib), Float64, $ptypes, x, $(pargs...), 1, 0)

        $ccdf($(pdecls...), x::Real) =
            ccall(($pfun, rmathlib), Float64, $ptypes, x, $(pargs...), 0, 0)

        $logcdf($(pdecls...), x::Real) =
            ccall(($pfun, rmathlib), Float64, $ptypes, x, $(pargs...), 1, 1)

        $logccdf($(pdecls...), x::Real) =
            ccall(($pfun, rmathlib), Float64, $ptypes, x, $(pargs...), 0, 1)

        $invcdf($(pdecls...), p::Real) =
            ccall(($qfun, rmathlib), Float64, $qtypes, p, $(pargs...), 1, 0)

        $invccdf($(pdecls...), p::Real) =
            ccall(($qfun, rmathlib), Float64, $qtypes, p, $(pargs...), 0, 0)

        $invlogcdf($(pdecls...), lp::Real) =
            ccall(($qfun, rmathlib), Float64, $qtypes, lp, $(pargs...), 1, 1)

        $invlogccdf($(pdecls...), lp::Real) =
            ccall(($qfun, rmathlib), Float64, $qtypes, lp, $(pargs...), 0, 1)
    end
end

macro import_rmath(nparams, rname, jname)
    esc(_import_rmath(nparams, rname, jname))
end


### Import specific functions

@import_rmath 2 norm norm
@import_rmath 2 beta beta

end
