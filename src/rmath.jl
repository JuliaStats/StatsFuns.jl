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

    is_tukey = rname == :tukey

    rand = Symbol(jname, "rand")
    has_rand = true
    if rname == :nbeta || rname == :nf || rname == :nt || is_tukey
        has_rand = false
    end

    # arguments & argument types
    _pts = fill(:Cdouble, length(pargs))

    if is_tukey
        dtypes = Expr(:tuple, :Cdouble, :Cdouble, _pts..., :Cint)
        ptypes = Expr(:tuple, :Cdouble, :Cdouble, _pts..., :Cint, :Cint)
        qtypes = Expr(:tuple, :Cdouble, :Cdouble, _pts..., :Cint, :Cint)
        rtypes = Expr(:tuple, :Cdouble, _pts...)
    else
        dtypes = Expr(:tuple, :Cdouble, _pts..., :Cint)
        ptypes = Expr(:tuple, :Cdouble, _pts..., :Cint, :Cint)
        qtypes = Expr(:tuple, :Cdouble, _pts..., :Cint, :Cint)
        rtypes = Expr(:tuple, _pts...)
    end

    pdecls = [Expr(:(::), ps, :Real) for ps in pargs] # [:(p1::Real), :(p2::Real), ...]

    if is_tukey
        # ptukey and qtukey have an extra literal 1 argument
        pargs = (1, pargs...)
    end

    # Function implementation
    quote
        if $(!is_tukey)
            function $pdf($(pdecls...), x::Real)
                T = float(Base.promote_typeof($(pargs...), x))
                return convert(T, ccall(($dfun, libRmath), Float64, $dtypes, x, $(pargs...), 0))
            end

            function $logpdf($(pdecls...), x::Real)
                T = float(Base.promote_typeof($(pargs...), x))
                return convert(T, ccall(($dfun, libRmath), Float64, $dtypes, x, $(pargs...), 1))
            end
        end

        function $cdf($(pdecls...), x::Real)
            T = float(Base.promote_typeof($(pargs...), x))
            return convert(T, ccall(($pfun, libRmath), Float64, $ptypes, x, $(pargs...), 1, 0))
        end

        function $ccdf($(pdecls...), x::Real)
            T = float(Base.promote_typeof($(pargs...), x))
            return convert(T, ccall(($pfun, libRmath), Float64, $ptypes, x, $(pargs...), 0, 0))
        end

        function $logcdf($(pdecls...), x::Real)
            T = float(Base.promote_typeof($(pargs...), x))
            return convert(T, ccall(($pfun, libRmath), Float64, $ptypes, x, $(pargs...), 1, 1))
        end

        function $logccdf($(pdecls...), x::Real)
            T = float(Base.promote_typeof($(pargs...), x))
            return convert(T, ccall(($pfun, libRmath), Float64, $ptypes, x, $(pargs...), 0, 1))
        end

        function $invcdf($(pdecls...), q::Real)
            T = float(Base.promote_typeof($(pargs...), q))
            return convert(T, ccall(($qfun, libRmath), Float64, $qtypes, q, $(pargs...), 1, 0))
        end

        function $invccdf($(pdecls...), q::Real)
            T = float(Base.promote_typeof($(pargs...), q))
            return convert(T, ccall(($qfun, libRmath), Float64, $qtypes, q, $(pargs...), 0, 0))
        end

        function $invlogcdf($(pdecls...), lq::Real)
            T = float(Base.promote_typeof($(pargs...), lq))
            return convert(T, ccall(($qfun, libRmath), Float64, $qtypes, lq, $(pargs...), 1, 1))
        end

        function $invlogccdf($(pdecls...), lq::Real)
            T = float(Base.promote_typeof($(pargs...), lq))
            return convert(T, ccall(($qfun, libRmath), Float64, $qtypes, lq, $(pargs...), 0, 1))
        end

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

@import_rmath signrank signrank n
@import_rmath beta beta α β
@import_rmath binom binom n p
@import_rmath chisq chisq k
@import_rmath f fdist ν1 ν2
@import_rmath gamma gamma k θ
@import_rmath hyper hyper ms mf n
@import_rmath norm norm μ σ
@import_rmath nbinom nbinom r p
@import_rmath pois pois λ
@import_rmath t tdist k
@import_rmath tukey srdist k ν
@import_rmath wilcox wilcox nx ny

@import_rmath nbeta nbeta α β λ
@import_rmath nchisq nchisq k λ
@import_rmath nf nfdist k1 k2 λ
@import_rmath nt ntdist k λ

end
