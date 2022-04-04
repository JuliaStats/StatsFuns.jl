module StatsFuns

using Base: Math.@horner
using IrrationalConstants
using LogExpFunctions
using SpecialFunctions
import ChainRulesCore
import InverseFunctions

export
    # distrs/beta
    betapdf,            # pdf of beta distribution
    betalogpdf,         # logpdf of beta distribution
    betacdf,            # cdf of beta distribution
    betaccdf,           # ccdf of beta distribution
    betalogcdf,         # logcdf of beta distribution
    betalogccdf,        # logccdf of beta distribution
    betainvcdf,         # inverse-cdf of beta distribution
    betainvccdf,        # inverse-ccdf of beta distribution
    betainvlogcdf,      # inverse-logcdf of beta distribution
    betainvlogccdf,     # inverse-logccdf of beta distribution

    # distrs/binom
    binompdf,           # pdf of binomial distribution
    binomlogpdf,        # logpdf of binomial distribution
    binomcdf,           # cdf of binomial distribution
    binomccdf,          # ccdf of binomial distribution
    binomlogcdf,        # logcdf of binomial distribution
    binomlogccdf,       # logccdf of binomial distribution
    binominvcdf,        # inverse-cdf of binomial distribution
    binominvccdf,       # inverse-ccdf of binomial distribution
    binominvlogcdf,     # inverse-logcdf of binomial distribution
    binominvlogccdf,    # inverse-logccdf of binomial distribution

    # distrs/chisq
    chisqpdf,           # pdf of chi-square distribution
    chisqlogpdf,        # logpdf of chi-square distribution
    chisqcdf,           # cdf of chi-square distribution
    chisqccdf,          # ccdf of chi-square distribution
    chisqlogcdf,        # logcdf of chi-square distribution
    chisqlogccdf,       # logccdf of chi-square distribution
    chisqinvcdf,        # inverse-cdf of chi-square distribution
    chisqinvccdf,       # inverse-ccdf of chi-square distribution
    chisqinvlogcdf,     # inverse-logcdf of chi-square distribution
    chisqinvlogccdf,    # inverse-logccdf of chi-square distribution

    # distrs/fdist
    fdistpdf,           # pdf of F distribution
    fdistlogpdf,        # logpdf of F distribution
    fdistcdf,           # cdf of F distribution
    fdistccdf,          # ccdf of F distribution
    fdistlogcdf,        # logcdf of F distribution
    fdistlogccdf,       # logccdf of F distribution
    fdistinvcdf,        # inverse-cdf of F distribution
    fdistinvccdf,       # inverse-ccdf of F distribution
    fdistinvlogcdf,     # inverse-logcdf of F distribution
    fdistinvlogccdf,    # inverse-logccdf of F distribution

    # distrs/gamma
    gammapdf,           # pdf of gamma distribution
    gammalogpdf,        # logpdf of gamma distribution
    gammacdf,           # cdf of gamma distribution
    gammaccdf,          # ccdf of gamma distribution
    gammalogcdf,        # logcdf of gamma distribution
    gammalogccdf,       # logccdf of gamma distribution
    gammainvcdf,        # inverse-cdf of gamma distribution
    gammainvccdf,       # inverse-ccdf of gamma distribution
    gammainvlogcdf,     # inverse-logcdf of gamma distribution
    gammainvlogccdf,    # inverse-logccdf of gamma distribution

    # distrs/hyper
    hyperpdf,           # pdf of hypergeometric distribution
    hyperlogpdf,        # logpdf of hypergeometric distribution
    hypercdf,           # cdf of hypergeometric distribution
    hyperccdf,          # ccdf of hypergeometric distribution
    hyperlogcdf,        # logcdf of hypergeometric distribution
    hyperlogccdf,       # logccdf of hypergeometric distribution
    hyperinvcdf,        # inverse-cdf of hypergeometric distribution
    hyperinvccdf,       # inverse-ccdf of hypergeometric distribution
    hyperinvlogcdf,     # inverse-logcdf of hypergeometric distribution
    hyperinvlogccdf,    # inverse-logccdf of hypergeometric distribution

    # distrs/nbeta
    nbetapdf,           # pdf of noncentral beta distribution
    nbetalogpdf,        # logpdf of noncentral beta distribution
    nbetacdf,           # cdf of noncentral beta distribution
    nbetaccdf,          # ccdf of noncentral beta distribution
    nbetalogcdf,        # logcdf of noncentral beta distribution
    nbetalogccdf,       # logccdf of noncentral beta distribution
    nbetainvcdf,        # inverse-cdf of noncentral beta distribution
    nbetainvccdf,       # inverse-ccdf of noncentral beta distribution
    nbetainvlogcdf,     # inverse-logcdf of noncentral beta distribution
    nbetainvlogccdf,    # inverse-logccdf of noncentral beta distribution

    # distrs/nbinom
    nbinompdf,           # pdf of negative nbinomial distribution
    nbinomlogpdf,        # logpdf of negative nbinomial distribution
    nbinomcdf,           # cdf of negative nbinomial distribution
    nbinomccdf,          # ccdf of negative nbinomial distribution
    nbinomlogcdf,        # logcdf of negative nbinomial distribution
    nbinomlogccdf,       # logccdf of negative nbinomial distribution
    nbinominvcdf,        # inverse-cdf of negative nbinomial distribution
    nbinominvccdf,       # inverse-ccdf of negative nbinomial distribution
    nbinominvlogcdf,     # inverse-logcdf of negative nbinomial distribution
    nbinominvlogccdf,    # inverse-logccdf of negative nbinomial distribution

    # distrs/nchisq
    nchisqpdf,           # pdf of noncentral chi-square distribution
    nchisqlogpdf,        # logpdf of noncentral chi-square distribution
    nchisqcdf,           # cdf of noncentral chi-square distribution
    nchisqccdf,          # ccdf of noncentral chi-square distribution
    nchisqlogcdf,        # logcdf of noncentral chi-square distribution
    nchisqlogccdf,       # logccdf of noncentral chi-square distribution
    nchisqinvcdf,        # inverse-cdf of noncentral chi-square distribution
    nchisqinvccdf,       # inverse-ccdf of noncentral chi-square distribution
    nchisqinvlogcdf,     # inverse-logcdf of noncentral chi-square distribution
    nchisqinvlogccdf,    # inverse-logccdf of noncentral chi-square distribution

    # distrs/nfdist
    nfdistpdf,           # pdf of noncentral F distribution
    nfdistlogpdf,        # logpdf of noncentral F distribution
    nfdistcdf,           # cdf of noncentral F distribution
    nfdistccdf,          # ccdf of noncentral F distribution
    nfdistlogcdf,        # logcdf of noncentral F distribution
    nfdistlogccdf,       # logccdf of noncentral F distribution
    nfdistinvcdf,        # inverse-cdf of noncentral F distribution
    nfdistinvccdf,       # inverse-ccdf of noncentral F distribution
    nfdistinvlogcdf,     # inverse-logcdf of noncentral F distribution
    nfdistinvlogccdf,    # inverse-logccdf of noncentral F distribution

    # distrs/norm
    normpdf,            # pdf of normal distribution
    normlogpdf,         # logpdf of normal distribution
    normcdf,            # cdf of normal distribution
    normccdf,           # ccdf of normal distribution
    normlogcdf,         # logcdf of normal distribution
    normlogccdf,        # logccdf of normal distribution
    norminvcdf,         # inverse-cdf of normal distribution
    norminvccdf,        # inverse-ccdf of normal distribution
    norminvlogcdf,      # inverse-logcdf of normal distribution
    norminvlogccdf,     # inverse-logccdf of normal distribution

    # distrs/ntdist
    ntdistpdf,           # pdf of noncentral t distribution
    ntdistlogpdf,        # logpdf of noncentral t distribution
    ntdistcdf,           # cdf of noncentral t distribution
    ntdistccdf,          # ccdf of noncentral t distribution
    ntdistlogcdf,        # logcdf of noncentral t distribution
    ntdistlogccdf,       # logccdf of noncentral t distribution
    ntdistinvcdf,        # inverse-cdf of noncentral t distribution
    ntdistinvccdf,       # inverse-ccdf of noncentral t distribution
    ntdistinvlogcdf,     # inverse-logcdf of noncentral t distribution
    ntdistinvlogccdf,    # inverse-logccdf of noncentral t distribution

    # distrs/pois
    poispdf,           # pdf of Poisson distribution
    poislogpdf,        # logpdf of Poisson distribution
    poiscdf,           # cdf of Poisson distribution
    poisccdf,          # ccdf of Poisson distribution
    poislogcdf,        # logcdf of Poisson distribution
    poislogccdf,       # logccdf of Poisson distribution
    poisinvcdf,        # inverse-cdf of Poisson distribution
    poisinvccdf,       # inverse-ccdf of Poisson distribution
    poisinvlogcdf,     # inverse-logcdf of Poisson distribution
    poisinvlogccdf,    # inverse-logccdf of Poisson distribution

    # distrs/tdist
    tdistpdf,           # pdf of student's t distribution
    tdistlogpdf,        # logpdf of student's t distribution
    tdistcdf,           # cdf of student's t distribution
    tdistccdf,          # ccdf of student's t distribution
    tdistlogcdf,        # logcdf of student's t distribution
    tdistlogccdf,       # logccdf of student's t distribution
    tdistinvcdf,        # inverse-cdf of student's t distribution
    tdistinvccdf,       # inverse-ccdf of student's t distribution
    tdistinvlogcdf,     # inverse-logcdf of student's t distribution
    tdistinvlogccdf,    # inverse-logccdf of student's t distribution

    # distrs/srdist
    srdistcdf,           # cdf of studentized range distribution
    srdistccdf,          # ccdf of studentized range distribution
    srdistlogcdf,        # logcdf of studentized range distribution
    srdistlogccdf,       # logccdf of studentized range distribution
    srdistinvcdf,        # inverse-cdf of studentized range distribution
    srdistinvccdf,       # inverse-ccdf of studentized range distribution
    srdistinvlogcdf,     # inverse-logcdf of studentized range distribution
    srdistinvlogccdf,    # inverse-logccdf of studentized range distribution

    # misc
    logmvgamma,         # logarithm of multivariate gamma function
    logmvbeta,          # logarithm of multivariate beta function
    lstirling_asym,

    # tvpack
    bvncdf    # bivariate standard normal cdf

## source files
include("misc.jl")
include("rmath.jl")
include("tvpack.jl")

using .RFunctions

include(joinpath("distrs", "beta.jl"))
include(joinpath("distrs", "binom.jl"))
include(joinpath("distrs", "chisq.jl"))
include(joinpath("distrs", "fdist.jl"))
include(joinpath("distrs", "gamma.jl"))
include(joinpath("distrs", "hyper.jl"))
include(joinpath("distrs", "nbeta.jl"))
include(joinpath("distrs", "nbinom.jl"))
include(joinpath("distrs", "nchisq.jl"))
include(joinpath("distrs", "nfdist.jl"))
include(joinpath("distrs", "norm.jl"))
include(joinpath("distrs", "ntdist.jl"))
include(joinpath("distrs", "pois.jl"))
include(joinpath("distrs", "tdist.jl"))
include(joinpath("distrs", "srdist.jl"))

include("chainrules.jl")
include("inverse.jl")

end # module
