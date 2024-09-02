# functions related to student's T distribution

tdistpdf(ν::Real, x::Real) = exp(tdistlogpdf(ν, x))

tdistlogpdf(ν::Real, x::Real) = tdistlogpdf(promote(ν, x)...)
function tdistlogpdf(ν::T, x::T) where {T<:Real}
    isinf(ν) && return normlogpdf(x)
    νp12 = (ν + 1) / 2
    return loggamma(νp12) - (logπ + log(ν)) / 2 - loggamma(ν / 2) - νp12 * log1p(x^2 / ν)
end

function tdistcdf(ν::T, x::T) where T<:Real
    if isinf(ν)
        return normcdf(x)
    elseif x < 0
        return fdistccdf(one(ν), ν, x^2)/2
    else
        return 1 - fdistccdf(one(ν), ν, x^2)/2
    end
end
tdistcdf(ν::Real, x::Real) = tdistcdf(map(float, promote(ν, x))...)

tdistccdf(ν::Real, x::Real) = tdistcdf(ν, -x)

function tdistlogcdf(ν::T, x::T) where T<:Real
    if isinf(ν)
        return normlogcdf(x)
    elseif x < 0
        ret = fdistlogccdf(one(ν), ν, x^2)
        return ret - logtwo
    else
        return log1p(-fdistccdf(one(ν), ν, x^2)/2)
    end
end
tdistlogcdf(ν::Real, x::Real) = tdistlogcdf(map(float, promote(ν, x))...)

tdistlogccdf(ν::Real, x::Real) = tdistlogcdf(ν, -x)

function tdistinvcdf(ν::T, p::T) where T<:Real
    if isinf(ν)
        return norminvcdf(p)
    elseif p < 0.5
        return -sqrt(fdistinvccdf(one(ν), ν, 2*p))
    else
        return sqrt(fdistinvccdf(one(ν), ν, 2*(1 - p)))
    end
end
tdistinvcdf(ν::Real, p::Real) = tdistinvcdf(map(float, promote(ν, p))...)

tdistinvccdf(ν::Real, p::Real) = -tdistinvcdf(ν, p)
function tdistinvlogcdf(ν::T, logp::T) where T<:Real
    if isinf(ν)
        return norminvlogcdf(logp)
    else
        logq = logp + logtwo
        if logq < 0
            return -sqrt(fdistinvlogccdf(one(ν), ν, logq))
        else
            return sqrt(fdistinvlogccdf(one(ν), ν, log2mexp(logq)))
        end
    end
end
tdistinvlogcdf(ν::Real, logp::Real) = tdistinvlogcdf(map(float, promote(ν, logp))...)

tdistinvlogccdf(ν::Real, logp::Real) = -tdistinvlogcdf(ν, logp)
