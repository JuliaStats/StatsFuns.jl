# import of Rmath functions

module RMath

const rmathlib = "libRmath-julia"

### import macro

function _import_rmath(nparams::Int, basename::Symbol)
    # basename: e.g. :norm is the basename for :dnorm, :pnorm, etc
    # nparams: the number of distribution parameters

    # C function names
    if basename == :norm
        dfun = Expr(:quote, "dnorm4")
        pfun = Expr(:quote, "pnorm5")
        qfun = Expr(:quote, "qnorm5")
        rfun = Expr(:quote, "rnorm")
    else
        dfun = Expr(:quote, string('d', basename))    # density
        pfun = Expr(:quote, string('p', basename))    # cumulative probability
        qfun = Expr(:quote, string('q', basename))    # quantile
        rfun = Expr(:quote, string('r', basename))    # random sampling
    end

    # Julia function names
    jbase = basename
    pdf = symbol(string(jbase, "pdf"))
    cdf = symbol(string(jbase, "cdf"))
    ccdf = symbol(string(jbase, "ccdf"))

    logpdf = symbol(string(jbase, "logpdf"))
    logcdf = symbol(string(jbase, "logcdf"))
    logccdf = symbol(string(jbase, "logccdf"))

    invcdf = symbol(string(jbase, "invcdf"))
    invccdf = symbol(string(jbase, "invccdf"))
    invlogcdf = symbol(string(jbase, "invlogcdf"))
    invlogccdf = symbol(string(jbase, "invlogccdf"))

    rand = symbol(string(jbase, "rand"))

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

macro import_rmath(nparams, basename)
    esc(_import_rmath(nparams, basename))
end


### Import specific functions

@import_rmath 2 norm
@import_rmath 2 beta

end
