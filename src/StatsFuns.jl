module StatsFuns

using Base: Math.@horner
using Reexport: @reexport
using SpecialFunctions: beta_inc, beta_inc_inv, digamma,
    erfc, erfcinv, erfcx, gamma_inc, gamma_inc_inv, logbeta, loggamma

# reexports
@reexport using IrrationalConstants:
    twoπ, # 2π
    fourπ, # 4π
    halfπ, # π / 2
    quartπ, # π / 4
    invπ, # 1 / π
    twoinvπ, # 2 / π
    fourinvπ, # 4 / π
    inv2π, # 1 / (2π)
    inv4π, # 1 / (4π)
    sqrt2, # √2
    sqrt3, # √3
    sqrtπ, # √π
    sqrt2π, # √2π
    sqrt4π, # √4π
    sqrthalfπ, # √(π / 2)
    invsqrt2, # 1 / √2
    invsqrt2π, # 1 / √2π
    loghalf, # log(1 / 2)
    logtwo, # log(2)
    logπ, # log(π)
    log2π, # log(2π)
    log4π       # log(4π)

@reexport using LogExpFunctions:
    xlogx, # x * log(x) for x > 0, or 0 when x == 0
    xlogy, # x * log(y) for x > 0, or 0 when x == 0
    xlog1py, # x * log(1 + y) for x > 0, or 0 when x == 0
    logistic, # 1 / (1 + exp(-x))
    logit, # log(x / (1 - x))
    log1psq, # log(1 + x^2)
    log1pexp, # log(1 + exp(x))
    log1mexp, # log(1 - exp(x))
    log2mexp, # log(2 - exp(x))
    logexpm1, # log(exp(x) - 1)
    softplus, # alias of log1pexp
    invsoftplus, # alias of logexpm1
    log1pmx, # log(1 + x) - x
    logmxp1, # log(x) - x + 1
    logaddexp, # log(exp(x) + exp(y))
    logsubexp, # log(abs(e^x - e^y))
    logsumexp, # log(sum(exp(x)))
    softmax, # exp(x_i) / sum(exp(x)), for i
    softmax!        # inplace softmax

export
    # distrs/beta
    betapdf, # pdf of beta distribution
    betalogpdf, # logpdf of beta distribution
    betalogupdf, # unnormalized logpdf of beta distribution (parameters constant)
    betalogulikelihood, # unnormalized logpdf of beta distribution (data constant)
    betacdf, # cdf of beta distribution
    betaccdf, # ccdf of beta distribution
    betalogcdf, # logcdf of beta distribution
    betalogccdf, # logccdf of beta distribution
    betainvcdf, # inverse-cdf of beta distribution
    betainvccdf, # inverse-ccdf of beta distribution
    betainvlogcdf, # inverse-logcdf of beta distribution
    betainvlogccdf, # inverse-logccdf of beta distribution

    # distrs/binom
    binompdf, # pdf of binomial distribution
    binomlogpdf, # logpdf of binomial distribution
    binomlogupdf, # unnormalized logpdf of binomial distribution (parameters constant)
    binomlogulikelihood, # unnormalized logpdf of binomial distribution (data constant)
    binomcdf, # cdf of binomial distribution
    binomccdf, # ccdf of binomial distribution
    binomlogcdf, # logcdf of binomial distribution
    binomlogccdf, # logccdf of binomial distribution
    binominvcdf, # inverse-cdf of binomial distribution
    binominvccdf, # inverse-ccdf of binomial distribution
    binominvlogcdf, # inverse-logcdf of binomial distribution
    binominvlogccdf, # inverse-logccdf of binomial distribution

    # distrs/chisq
    chisqpdf, # pdf of chi-square distribution
    chisqlogpdf, # logpdf of chi-square distribution
    chisqlogupdf, # unnormalized logpdf of chi-square distribution (parameters constant)
    chisqlogulikelihood, # unnormalized logpdf of chi-square distribution (data constant)
    chisqcdf, # cdf of chi-square distribution
    chisqccdf, # ccdf of chi-square distribution
    chisqlogcdf, # logcdf of chi-square distribution
    chisqlogccdf, # logccdf of chi-square distribution
    chisqinvcdf, # inverse-cdf of chi-square distribution
    chisqinvccdf, # inverse-ccdf of chi-square distribution
    chisqinvlogcdf, # inverse-logcdf of chi-square distribution
    chisqinvlogccdf, # inverse-logccdf of chi-square distribution

    # distrs/fdist
    fdistpdf, # pdf of F distribution
    fdistlogpdf, # logpdf of F distribution
    fdistlogupdf, # unnormalized logpdf of F distribution (parameters constant)
    fdistlogulikelihood, # unnormalized logpdf of F distribution (data constant)
    fdistcdf, # cdf of F distribution
    fdistccdf, # ccdf of F distribution
    fdistlogcdf, # logcdf of F distribution
    fdistlogccdf, # logccdf of F distribution
    fdistinvcdf, # inverse-cdf of F distribution
    fdistinvccdf, # inverse-ccdf of F distribution
    fdistinvlogcdf, # inverse-logcdf of F distribution
    fdistinvlogccdf, # inverse-logccdf of F distribution

    # distrs/gamma
    gammapdf, # pdf of gamma distribution
    gammalogpdf, # logpdf of gamma distribution
    gammalogupdf, # unnormalized logpdf of gamma distribution (parameters constant)
    gammalogulikelihood, # unnormalized logpdf of gamma distribution (data constant)
    gammacdf, # cdf of gamma distribution
    gammaccdf, # ccdf of gamma distribution
    gammalogcdf, # logcdf of gamma distribution
    gammalogccdf, # logccdf of gamma distribution
    gammainvcdf, # inverse-cdf of gamma distribution
    gammainvccdf, # inverse-ccdf of gamma distribution
    gammainvlogcdf, # inverse-logcdf of gamma distribution
    gammainvlogccdf, # inverse-logccdf of gamma distribution

    # distrs/hyper
    hyperpdf, # pdf of hypergeometric distribution
    hyperlogpdf, # logpdf of hypergeometric distribution
    hyperlogupdf, # unnormalized logpdf of hypergeometric distribution (parameters constant)
    hyperlogulikelihood, # unnormalized logpdf of hypergeometric distribution (data constant)
    hypercdf, # cdf of hypergeometric distribution
    hyperccdf, # ccdf of hypergeometric distribution
    hyperlogcdf, # logcdf of hypergeometric distribution
    hyperlogccdf, # logccdf of hypergeometric distribution
    hyperinvcdf, # inverse-cdf of hypergeometric distribution
    hyperinvccdf, # inverse-ccdf of hypergeometric distribution
    hyperinvlogcdf, # inverse-logcdf of hypergeometric distribution
    hyperinvlogccdf, # inverse-logccdf of hypergeometric distribution

    # distrs/nbeta
    nbetapdf, # pdf of noncentral beta distribution
    nbetalogpdf, # logpdf of noncentral beta distribution
    nbetalogupdf, # unnormalized logpdf of noncentral beta distribution (parameters constant)
    nbetalogulikelihood, # unnormalized logpdf of noncentral beta distribution (data constant)
    nbetacdf, # cdf of noncentral beta distribution
    nbetaccdf, # ccdf of noncentral beta distribution
    nbetalogcdf, # logcdf of noncentral beta distribution
    nbetalogccdf, # logccdf of noncentral beta distribution
    nbetainvcdf, # inverse-cdf of noncentral beta distribution
    nbetainvccdf, # inverse-ccdf of noncentral beta distribution
    nbetainvlogcdf, # inverse-logcdf of noncentral beta distribution
    nbetainvlogccdf, # inverse-logccdf of noncentral beta distribution

    # distrs/nbinom
    nbinompdf, # pdf of negative nbinomial distribution
    nbinomlogpdf, # logpdf of negative nbinomial distribution
    nbinomlogupdf, # unnormalized logpdf of negative nbinomial distribution (parameters constant)
    nbinomlogulikelihood, # unnormalized logpdf of negative nbinomial distribution (data constant)
    nbinomcdf, # cdf of negative nbinomial distribution
    nbinomccdf, # ccdf of negative nbinomial distribution
    nbinomlogcdf, # logcdf of negative nbinomial distribution
    nbinomlogccdf, # logccdf of negative nbinomial distribution
    nbinominvcdf, # inverse-cdf of negative nbinomial distribution
    nbinominvccdf, # inverse-ccdf of negative nbinomial distribution
    nbinominvlogcdf, # inverse-logcdf of negative nbinomial distribution
    nbinominvlogccdf, # inverse-logccdf of negative nbinomial distribution

    # distrs/nchisq
    nchisqpdf, # pdf of noncentral chi-square distribution
    nchisqlogpdf, # logpdf of noncentral chi-square distribution
    nchisqlogupdf, # unnormalized logpdf of noncentral chi-square distribution (parameters constant)
    nchisqlogulikelihood, # unnormalized logpdf of noncentral chi-square distribution (data constant)
    nchisqcdf, # cdf of noncentral chi-square distribution
    nchisqccdf, # ccdf of noncentral chi-square distribution
    nchisqlogcdf, # logcdf of noncentral chi-square distribution
    nchisqlogccdf, # logccdf of noncentral chi-square distribution
    nchisqinvcdf, # inverse-cdf of noncentral chi-square distribution
    nchisqinvccdf, # inverse-ccdf of noncentral chi-square distribution
    nchisqinvlogcdf, # inverse-logcdf of noncentral chi-square distribution
    nchisqinvlogccdf, # inverse-logccdf of noncentral chi-square distribution

    # distrs/nfdist
    nfdistpdf, # pdf of noncentral F distribution
    nfdistlogpdf, # logpdf of noncentral F distribution
    nfdistlogupdf, # unnormalized logpdf of noncentral F distribution (parameters constant)
    nfdistlogulikelihood, # unnormalized logpdf of noncentral F distribution (data constant)
    nfdistcdf, # cdf of noncentral F distribution
    nfdistccdf, # ccdf of noncentral F distribution
    nfdistlogcdf, # logcdf of noncentral F distribution
    nfdistlogccdf, # logccdf of noncentral F distribution
    nfdistinvcdf, # inverse-cdf of noncentral F distribution
    nfdistinvccdf, # inverse-ccdf of noncentral F distribution
    nfdistinvlogcdf, # inverse-logcdf of noncentral F distribution
    nfdistinvlogccdf, # inverse-logccdf of noncentral F distribution

    # distrs/norm
    normpdf, # pdf of normal distribution
    normlogpdf, # logpdf of normal distribution
    normlogupdf, # unnormalized logpdf of normal distribution (parameters constant)
    normlogulikelihood, # unnormalized logpdf of normal distribution (data constant)
    normcdf, # cdf of normal distribution
    normccdf, # ccdf of normal distribution
    normlogcdf, # logcdf of normal distribution
    normlogccdf, # logccdf of normal distribution
    norminvcdf, # inverse-cdf of normal distribution
    norminvccdf, # inverse-ccdf of normal distribution
    norminvlogcdf, # inverse-logcdf of normal distribution
    norminvlogccdf, # inverse-logccdf of normal distribution

    # distrs/ntdist
    ntdistpdf, # pdf of noncentral t distribution
    ntdistlogpdf, # logpdf of noncentral t distribution
    ntdistlogupdf, # unnormalized logpdf of noncentral t distribution (parameters constant)
    ntdistlogulikelihood, # unnormalized logpdf of noncentral t distribution (data constant)
    ntdistcdf, # cdf of noncentral t distribution
    ntdistccdf, # ccdf of noncentral t distribution
    ntdistlogcdf, # logcdf of noncentral t distribution
    ntdistlogccdf, # logccdf of noncentral t distribution
    ntdistinvcdf, # inverse-cdf of noncentral t distribution
    ntdistinvccdf, # inverse-ccdf of noncentral t distribution
    ntdistinvlogcdf, # inverse-logcdf of noncentral t distribution
    ntdistinvlogccdf, # inverse-logccdf of noncentral t distribution

    # distrs/pois
    poispdf, # pdf of Poisson distribution
    poislogpdf, # logpdf of Poisson distribution
    poislogupdf, # unnormalized logpdf of Poisson distribution (parameters constant)
    poislogulikelihood, # unnormalized logpdf of Poisson distribution (data constant)
    poiscdf, # cdf of Poisson distribution
    poisccdf, # ccdf of Poisson distribution
    poislogcdf, # logcdf of Poisson distribution
    poislogccdf, # logccdf of Poisson distribution
    poisinvcdf, # inverse-cdf of Poisson distribution
    poisinvccdf, # inverse-ccdf of Poisson distribution
    poisinvlogcdf, # inverse-logcdf of Poisson distribution
    poisinvlogccdf, # inverse-logccdf of Poisson distribution

    # distrs/tdist
    tdistpdf, # pdf of student's t distribution
    tdistlogpdf, # logpdf of student's t distribution
    tdistlogupdf, # unnormalized logpdf of student's t distribution (parameters constant)
    tdistlogulikelihood, # unnormalized logpdf of student's t distribution (data constant)
    tdistcdf, # cdf of student's t distribution
    tdistccdf, # ccdf of student's t distribution
    tdistlogcdf, # logcdf of student's t distribution
    tdistlogccdf, # logccdf of student's t distribution
    tdistinvcdf, # inverse-cdf of student's t distribution
    tdistinvccdf, # inverse-ccdf of student's t distribution
    tdistinvlogcdf, # inverse-logcdf of student's t distribution
    tdistinvlogccdf, # inverse-logccdf of student's t distribution

    # distrs/signrank
    signrankpdf,
    signranklogpdf,
    signranklogupdf,
    signranklogulikelihood,
    signrankcdf,
    signranklogcdf,
    signrankccdf,
    signranklogccdf,
    signrankinvcdf,
    signrankinvccdf,
    signrankinvlogcdf,
    signrankinvlogccdf,

    # distrs/srdist
    srdistcdf, # cdf of studentized range distribution
    srdistccdf, # ccdf of studentized range distribution
    srdistlogcdf, # logcdf of studentized range distribution
    srdistlogccdf, # logccdf of studentized range distribution
    srdistinvcdf, # inverse-cdf of studentized range distribution
    srdistinvccdf, # inverse-ccdf of studentized range distribution
    srdistinvlogcdf, # inverse-logcdf of studentized range distribution
    srdistinvlogccdf, # inverse-logccdf of studentized range distribution

    # distrs/wilcox
    wilcoxpdf,
    wilcoxlogpdf,
    wilcoxlogupdf,
    wilcoxlogulikelihood,
    wilcoxcdf,
    wilcoxlogcdf,
    wilcoxccdf,
    wilcoxlogccdf,
    wilcoxinvcdf,
    wilcoxinvccdf,
    wilcoxinvlogcdf,
    wilcoxinvlogccdf,

    # misc
    logmvgamma, # logarithm of multivariate gamma function
    logmvbeta, # logarithm of multivariate beta function
    lstirling_asym,

    # tvpack
    bvncdf    # bivariate standard normal cdf

## source files
include("misc.jl")
include("rmath.jl")
include("tvpack.jl")

using .RFunctions

include("distrs/beta.jl")
include("distrs/binom.jl")
include("distrs/chisq.jl")
include("distrs/fdist.jl")
include("distrs/gamma.jl")
include("distrs/hyper.jl")
include("distrs/nbeta.jl")
include("distrs/nbinom.jl")
include("distrs/nchisq.jl")
include("distrs/nfdist.jl")
include("distrs/norm.jl")
include("distrs/ntdist.jl")
include("distrs/pois.jl")
include("distrs/tdist.jl")
include("distrs/signrank.jl")
include("distrs/srdist.jl")
include("distrs/wilcox.jl")

end # module
