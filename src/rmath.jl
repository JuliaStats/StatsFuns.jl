# import of Rmath functions

module RFunctions

using Rmath: Rmath

### import macro

function _import_rmath(rname::Symbol, jname::Symbol, pargs)
    # Rmath function names
    dfun = Symbol(:d, rname)    # density
    pfun = Symbol(:p, rname)    # cumulative probability
    qfun = Symbol(:q, rname)    # quantile
    rfun = Symbol(:r, rname)    # random sampling

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

    # arguments & argument types
    pdecls = [Expr(:(::), ps, :Real) for ps in pargs] # [:(p1::Real), :(p2::Real), ...]

    if rname == :tukey
        # ptukey and qtukey have an extra literal 1 argument
        pargs = (pargs..., 1)
    end

    # Function implementation
    return quote
        if $(isdefined(Rmath, dfun))
            function $pdf($(pdecls...), x::Real)
                T = float(Base.promote_typeof($(pargs...), x))
                return convert(T, Rmath.$dfun(x, $(pargs...), false))
            end

            function $logpdf($(pdecls...), x::Real)
                T = float(Base.promote_typeof($(pargs...), x))
                return convert(T, Rmath.$dfun(x, $(pargs...), true))
            end
        end

        if $(isdefined(Rmath, pfun))
            function $cdf($(pdecls...), x::Real)
                T = float(Base.promote_typeof($(pargs...), x))
                return convert(T, Rmath.$pfun(x, $(pargs...), true, false))
            end

            function $ccdf($(pdecls...), x::Real)
                T = float(Base.promote_typeof($(pargs...), x))
                return convert(T, Rmath.$pfun(x, $(pargs...), false, false))
            end

            function $logcdf($(pdecls...), x::Real)
                T = float(Base.promote_typeof($(pargs...), x))
                return convert(T, Rmath.$pfun(x, $(pargs...), true, true))
            end

            function $logccdf($(pdecls...), x::Real)
                T = float(Base.promote_typeof($(pargs...), x))
                return convert(T, Rmath.$pfun(x, $(pargs...), false, true))
            end
        end

        if $(isdefined(Rmath, qfun))
            function $invcdf($(pdecls...), q::Real)
                T = float(Base.promote_typeof($(pargs...), q))
                return convert(T, Rmath.$qfun(q, $(pargs...), true, false))
            end

            function $invccdf($(pdecls...), q::Real)
                T = float(Base.promote_typeof($(pargs...), q))
                return convert(T, Rmath.$qfun(q, $(pargs...), false, false))
            end

            function $invlogcdf($(pdecls...), lq::Real)
                T = float(Base.promote_typeof($(pargs...), lq))
                return convert(T, Rmath.$qfun(lq, $(pargs...), true, true))
            end

            function $invlogccdf($(pdecls...), lq::Real)
                T = float(Base.promote_typeof($(pargs...), lq))
                return convert(T, Rmath.$qfun(lq, $(pargs...), false, true))
            end
        end
    end
end

macro import_rmath(rname, jname, pargs...)
    return esc(_import_rmath(rname, jname, pargs))
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
