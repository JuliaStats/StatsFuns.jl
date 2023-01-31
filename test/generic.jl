using StatsFuns
using StatsFuns: RFunctions
using ForwardDiff: Dual
using Test

include("utils.jl")

function check_rmath(fname, statsfun, rmathfun, params, aname, a, isprob, rtol)
    v = @inferred(rmathfun(params..., a))
    rv = @inferred(statsfun(params..., Dual(a))).value
    @test v isa float(Base.promote_typeof(params..., a))
    @test rv isa float(Base.promote_typeof(params..., a))
    if isprob
        @test v ≈ rv rtol=rtol nans=true
    else
        @test v ≈ rv atol=rtol rtol=rtol nans=true
    end
end

function genericcomp(basename, params, X::AbstractArray, rtol=_default_rtol(params, X))
  pdf = string(basename, "pdf")
  logpdf = string(basename, "logpdf")
  stats_pdf = eval(Symbol(pdf))
  stats_logpdf = eval(Symbol(logpdf))
  rmath_pdf = eval(Meta.parse(string("RFunctions.", pdf)))
  rmath_logpdf = eval(Meta.parse(string("RFunctions.", logpdf)))

  @testset "pdf with params=$params, x=$x" for x in X
    check_rmath(pdf, stats_pdf, rmath_pdf, params, "x", x, true, rtol)
  end

  @testset "logpdf with params=$params, x=$x" for x in X
    check_rmath(logpdf, stats_logpdf, rmath_logpdf, params, "x", x, false, rtol)
  end
end

function genericcomp_tests(basename, configs)
  println("\ttesting $basename ...")
  for (params, data) in configs
      genericcomp(basename, params, data)
  end
end

### Test cases

@testset "Generic" begin
    genericcomp_tests("beta", [
        ((1.0, 1.0), 0.01:0.01:0.99),
        ((2.0, 3.0), 0.01:0.01:0.99),
        ((10.0, 2.0), 0.01:0.01:0.99),
    ])

    genericcomp_tests("binom", [
        ((1, 0.5), 0.0:1.0),
        ((1, 0.7), 0.0:1.0),
        ((8, 0.6), 0.0:8.0),
        ((20, 0.1), 0.0:20.0),
        ((20, 0.9), 0.0:20.0),
        ((20, 0.9), 0:20),
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

    genericcomp_tests("pois", [
        ((1.0,), 0:30),
        ((10.0,), 0:37),
        ((10.0,), 39:42),
    ])
    # Requires slightly larger tolerance: #157
    genericcomp("pois", (10.0,), 38:38, 2.5e-14)

    genericcomp_tests("tdist", [
        ((1,), -5.0:0.1:5.0),
        ((2,), -5.0:0.1:5.0),
        ((5,), -5.0:0.1:5.0),
    ])

    #genericcomp_tests("srdist", [
    #    ((1,2), (0.0:0.1:5.0)),
    #    ((2,2), (0.0:0.1:5.0)),
    #    ((5,3), (0.0:0.1:5.0)),
    #    ((10,2), (0.0:0.1:5.0)),
    #    ((10,5), (0.0:0.1:5.0))
    #])
end
