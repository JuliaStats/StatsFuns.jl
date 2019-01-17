# import of Rmath functions

module RFunctions

import Rmath: libRmath

### import macro

function _import_rmath(rname::Symbol, jname::Symbol, pargs)
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
    pdf = Symbol(jname, "pdf")
    cdf = Symbol(jname, "cdf")
    ccdf = Symbol(jname, "ccdf")

    logpdf = Symbol(jname, "logpdf")
    logcdf = Symbol(jname, "logcdf")
    logccdf = Symbol(jname, "logccdf")

    invcdf = Symbol(jname, "invcdf")
    invccdf = Symbol(jname, "invccdf")
    invlogcdf = Symbol(jname, "invlogcdf")
    invlogccdf = Symbol(jname, "invlogccdf")

    rand = Symbol(jname, "rand")
    has_rand = true
    if rname == :nbeta || rname == :nf || rname == :nt
        has_rand = false
    end

    # arguments & argument types
    _pts = fill(:Cdouble, length(pargs))

    dtypes = Expr(:tuple, :Cdouble, _pts..., :Cint)
    ptypes = Expr(:tuple, :Cdouble, _pts..., :Cint, :Cint)
    qtypes = Expr(:tuple, :Cdouble, _pts..., :Cint, :Cint)
    rtypes = Expr(:tuple, _pts...)

    pdecls = [Expr(:(::), ps, :(Union{Float64,Int})) for ps in pargs] # [:(p1::Union{Float64, Int}), :(p2::Union{...}), ...]

    # function implementation
    quote
        $pdf($(pdecls...), x::Union{Float64,Int}) =
            ccall(($dfun, libRmath), Float64, $dtypes, x, $(pargs...), 0)

        $logpdf($(pdecls...), x::Union{Float64,Int}) =
            ccall(($dfun, libRmath), Float64, $dtypes, x, $(pargs...), 1)

        $cdf($(pdecls...), x::Union{Float64,Int}) =
            ccall(($pfun, libRmath), Float64, $ptypes, x, $(pargs...), 1, 0)

        $ccdf($(pdecls...), x::Union{Float64,Int}) =
            ccall(($pfun, libRmath), Float64, $ptypes, x, $(pargs...), 0, 0)

        $logcdf($(pdecls...), x::Union{Float64,Int}) =
            ccall(($pfun, libRmath), Float64, $ptypes, x, $(pargs...), 1, 1)

        $logccdf($(pdecls...), x::Union{Float64,Int}) =
            ccall(($pfun, libRmath), Float64, $ptypes, x, $(pargs...), 0, 1)

        $invcdf($(pdecls...), q::Union{Float64,Int}) =
            ccall(($qfun, libRmath), Float64, $qtypes, q, $(pargs...), 1, 0)

        $invccdf($(pdecls...), q::Union{Float64,Int}) =
            ccall(($qfun, libRmath), Float64, $qtypes, q, $(pargs...), 0, 0)

        $invlogcdf($(pdecls...), lq::Union{Float64,Int}) =
            ccall(($qfun, libRmath), Float64, $qtypes, lq, $(pargs...), 1, 1)

        $invlogccdf($(pdecls...), lq::Union{Float64,Int}) =
            ccall(($qfun, libRmath), Float64, $qtypes, lq, $(pargs...), 0, 1)

        if $has_rand
            $rand($(pdecls...)) =
                ccall(($rfun, libRmath), Float64, $rtypes, $(pargs...))
        end
    end
end

macro import_rmath(rname, jname, pargs...)
    esc(_import_rmath(rname, jname, pargs))
end


### Import specific functions

@import_rmath beta beta α β
@import_rmath binom binom n p
@import_rmath chisq chisq k
@import_rmath f fdist d1 d2
@import_rmath gamma gamma α β
@import_rmath hyper hyper ms mf n
@import_rmath norm norm μ σ
@import_rmath nbinom nbinom r p
@import_rmath pois pois λ
@import_rmath t tdist k

@import_rmath nbeta nbeta α β λ
@import_rmath nchisq nchisq k λ
@import_rmath nf nfdist k1 k2 λ
@import_rmath nt ntdist k λ

end
