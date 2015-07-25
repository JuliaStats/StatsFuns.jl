using StatsFuns
using Base.Test

function check_rmath(fname, statsfun, rmathfun, params, aname, a, isprob, rtol)
    v = statsfun(params..., a)
    rv = rmathfun(params..., a)
    if isprob
        rd = abs(v / rv - 1.0)
        if rd > rtol
            error("$fname deviates too much from Rmath at " *
                "params = $params, $aname = $a:\n" *
                "  v = $v (rv = $rv)\n  |v/rv - 1| = $rd > $rtol.")
        end
    else
        τ = (1.0 + abs(rv)) * rtol
        ad = abs(v - rv)
        if ad > τ
            error("$fname deviates too much from Rmath at " *
                "params = $params, $aname = $a:\n" *
                "  v = $v (rv = $rv)\n  |v - rv| = $ad > $τ.")
        end
    end
end

get_statsfun(fname) = eval(symbol(fname))
get_rmathfun(fname) = eval(parse(string("Rmath.", fname)))

function rmathcomp(basename, params, X::AbstractArray, rtol=1.0e-12)
    pdf     = string(basename, "pdf")
    logpdf  = string(basename, "logpdf")
    cdf     = string(basename, "cdf")
    ccdf    = string(basename, "ccdf")
    logcdf  = string(basename, "logcdf")
    logccdf = string(basename, "logccdf")

    stats_pdf = get_statsfun(pdf)
    stats_logpdf = get_statsfun(logpdf)

    rmath_pdf = get_statsfun(pdf)
    rmath_logpdf = get_statsfun(logpdf)

    for i = 1:length(X)
        x = X[i]
        check_rmath(pdf, stats_pdf, rmath_pdf,
            params, "x", x, true, rtol)
        check_rmath(logpdf, stats_logpdf, rmath_logpdf,
            params, "x", x, false, rtol)
    end
end


# macro generate_rmath_compfun(basename)
#     # function symbols
#     compfun = symbol(string("rmathcomp_", basename))
#
#     pdf     = symbol(string(basename, "pdf"))
#     logpdf  = symbol(string(basename, "logpdf"))
#     cdf     = symbol(string(basename, "cdf"))
#     ccdf    = symbol(string(basename, "ccdf"))
#     logcdf  = symbol(string(basename, "logcdf"))
#     logccdf = symbol(string(basename, "logccdf"))
#
#     invcdf     = symbol(string(basename, "invcdf"))
#     invccdf    = symbol(string(basename, "invccdf"))
#     invlogcdf  = symbol(string(basename, "invlogcdf"))
#     invlogccdf = symbol(string(basename, "invlogccdf"))
#
#     esc(quote
#         function $(compfun){N}(params::NTuple{N,Real}, X::AbstractArray, rtol=1.0e-12)
#             for i = 1:length(X)
#                 x = X[i]
#                 check_rmath($(string(pdf)), params, "x", x,
#                     $(pdf)(params..., x), Rmath.$(pdf)(params..., x), true, rtol)
#                 check_rmath($(string(logpdf)), params, "x", x,
#                     $(logpdf)(params..., x), Rmath.$(logpdf)(params..., x), false, rtol)
#
#                 check_rmath($(string(cdf)), params, "x", x,
#                     $(cdf)(params..., x), Rmath.$(cdf)(params..., x), true, rtol)
#                 check_rmath($(string(ccdf)), params, "x", x,
#                     $(ccdf)(params..., x), Rmath.$(ccdf)(params..., x), true, rtol)
#                 check_rmath($(string(logcdf)), params, "x", x,
#                     $(logcdf)(params..., x), Rmath.$(logcdf)(params..., x), false, rtol)
#                 check_rmath($(string(logccdf)), params, "x", x,
#                     $(logccdf)(params..., x), Rmath.$(logccdf)(params..., x), true, rtol)
#
#                 p = Rmath.$(cdf)(params..., x)
#                 cp = Rmath.$(ccdf)(params..., x)
#                 lp = log(p)
#                 lcp = log(cp)
#
#                 check_rmath($(string(invcdf)), params, "p", p,
#                     $(invcdf)(params..., p), Rmath.$(invcdf)(params..., p), false, rtol)
#                 check_rmath($(string(invccdf)), params, "p", cp,
#                     $(invccdf)(params..., cp), Rmath.$(invccdf)(params..., cp), false, rtol)
#                 check_rmath($(string(invlogcdf)), params, "lp", lp,
#                     $(invlogcdf)(params..., lp), Rmath.$(invlogcdf)(params..., lp), false, rtol)
#                 check_rmath($(string(invlogccdf)), params, "lp", lcp,
#                     $(invlogccdf)(params..., lcp), Rmath.$(invlogccdf)(params..., lcp), false, rtol)
#             end
#         end
#     end)
# end
#
# @generate_rmath_compfun norm

rmathcomp("norm", (0.0, 1.0), -6.0:0.01:6.0)
rmathcomp("norm", (2.0, 1.0), -3.0:0.01:7.0)
rmathcomp("norm", (0.0, 0.5), -3.0:0.01:3.0)
