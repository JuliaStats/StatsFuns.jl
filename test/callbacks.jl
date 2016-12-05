using StatsFuns
using ForwardDiff: Dual

function check_rmath(fname, statsfun, params, aname, a, isprob, rtol)
    v = statsfun(params..., a)
    rv = statsfun(params..., Dual(a)).value
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

function genericcomp(basename, params, X::AbstractArray, rtol=100eps(float(one(eltype(X)))))
  pdf = string(basename, "pdf")
  logpdf = string(basename, "logpdf")
  stats_pdf = eval(Symbol(pdf))
  stats_logpdf = eval(Symbol(logpdf))
  for i = 1:length(X)
    x = X[i]
    check_rmath(pdf, stats_pdf, params, "x", x, true, rtol)
    check_rmath(logpdf, stats_logpdf, params, "x", x, false, rtol)
  end
end

function genericcomp_tests(basename, configs)
  println("\ttesting $basename ...")
  for (params, data) in configs
      genericcomp(basename, params, data)
  end
end

### Test cases

genericcomp_tests("beta", [
    ((1.0, 1.0), 0.01:0.01:0.99),
    ((2.0, 3.0), 0.01:0.01:0.99),
    ((10.0, 2.0), 0.01:0.01:0.99),
])

genericcomp_tests("chisq", [
    ((1,), 0.0:0.1:8.0),
    ((4,), 0.0:0.1:8.0),
    ((9,), 0.0:0.1:8.0),
])

genericcomp_tests("fdist", [
    ((1, 1), (0.0:0.1:5.0)),
    ((2, 1), (0.0:0.1:5.0)),
    ((5, 2), (0.0:0.1:5.0)),
    ((10, 1), (0.0:0.1:5.0)),
    ((10, 3), (0.0:0.1:5.0)),
])

genericcomp_tests("gamma", [
    ((1.0, 1.0), (0.05:0.05:12.0)),
    ((0.5, 1.0), (0.05:0.05:12.0)),
    ((3.0, 1.0), (0.05:0.05:12.0)),
    ((9.0, 1.0), (0.05:0.05:12.0)),
    ((2.0, 3.0), (0.05:0.05:12.0)),
])

# genericcomp_tests("nbeta", [
#     ((1.0, 1.0, 0.0), 0.01:0.01:0.99),
#     ((2.0, 3.0, 0.0), 0.01:0.01:0.99),
#     ((1.0, 1.0, 2.0), 0.01:0.01:0.99),
#     ((3.0, 4.0, 2.0), 0.01:0.01:0.99),
# ])

# genericcomp_tests("nchisq", [
#     ((2, 1), 0.0:0.2:8.0),
#     ((2, 3), 0.0:0.2:8.0),
#     ((4, 1), 0.0:0.2:8.0),
#     ((4, 3), 0.0:0.2:8.0),
# ])

# genericcomp_tests("nfdist", [
#     ((1.0, 1.0, 0.0), 0.1:0.1:10.0),
#     ((1.0, 1.0, 2.0), 0.1:0.1:10.0),
#     ((2.0, 3.0, 1.0), 0.1:0.1:10.0),
# ])

genericcomp_tests("norm", [
    ((0.0, 1.0), -6.0:0.01:6.0),
    ((2.0, 1.0), -3.0:0.01:7.0),
    ((0.0, 0.5), -3.0:0.01:3.0),
    ((0, 1), -3.0:3.0),
    ((0, 1), -3.0:0.01:3.0)
])

# genericcomp_tests("ntdist", [
#     ((0, 1), -4.0:0.1:10.0),
#     ((0, 4), -4.0:0.1:10.0),
#     ((2, 1), -4.0:0.1:10.0),
#     ((2, 4), -4.0:0.1:10.0),
# ])

genericcomp_tests("tdist", [
    ((1,), -5.0:0.1:5.0),
    ((2,), -5.0:0.1:5.0),
    ((5,), -5.0:0.1:5.0),
])
