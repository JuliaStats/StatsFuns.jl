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

if VERSION < v"1.7.0-DEV.1172"
    function _expm1(x::Float16)
        if -0.2 < x < 0.1
            return @evalpoly(x, Float16(0), Float16(1), Float16(1/2), Float16(1/6), Float16(1/24), Float16(1/120))
        else
            return exp(x) - 1
        end
    end
end
_expm1(x::Number) = expm1(x)

function tdistinvlogcdf(ν::T, logp::T) where T<:Real
    if isinf(ν)
        return norminvlogcdf(logp)
    elseif logp < -log(2)
        return -sqrt(fdistinvlogccdf(one(ν), ν, logp + logtwo))
    else
        return sqrt(fdistinvccdf(one(ν), ν, -2*_expm1(logp)))
    end
end
tdistinvlogcdf(ν::Real, logp::Real) = tdistinvlogcdf(map(float, promote(ν, logp))...)

tdistinvlogccdf(ν::Real, logp::Real) = -tdistinvlogcdf(ν, logp)
