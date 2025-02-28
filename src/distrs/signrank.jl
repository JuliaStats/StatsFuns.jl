#=
The signrank distribution is equivalent to the problem
of finding the number of subsets of {1,2,...,n} summing to W,
relative to the total number of subsets.
The empty subset is defined to sum to zero.
This can be calculated using the recursion:
either n is in the subset in which case we need to calculate
the number of subsets of {1,2,...,n-1} summing to W-n,
or n is not in the subset in which case we need to calculate
the number of subsets of {1,2,...,n-1} summing to W.
This can be calculated bottom up using dynamic programming.

The i'th element of DP in the j'th outer loop iteration represents:
the number of ways {1,2,...,j} can sum to W-i+1.
 =#

@inline function signrankDP(n, W)
    DP = zeros(Int, W + 1)
    DP[W+1] = 1
    for j in 1:n
        for i in 1:(W+1-j)
            DP[i] += DP[i+j]
        end
    end
    return DP
end

function signrankpdf(n::Int, W::Union{Float64,Int})
    f, _ = modf(W)
    if f != 0
        return 0.0
    end
    W = Int(W)
    max_W = (n * (n + 1)) >> 1
    if W < 0
        return 0.0
    elseif W > max_W >> 1
        return signrankpdf(n, max_W - W)
    end
    DP = signrankDP(n, W)
    return ldexp(float(DP[1]), -n)
end

function signranklogpdf(n::Int, W::Union{Float64,Int})
    return log(signrankpdf(n, W))
end

function signrankcdf(n::Int, W::Union{Float64,Int})
    W = round(Int, W, RoundNearestTiesUp)
    max_W = (n * (n + 1)) >> 1
    if W < 0
        return 0.0
    elseif W > max_W >> 1
        return 1.0 - signrankcdf(n, max_W - W - 1)
    end
    DP = signrankDP(n, W)
    return sum(Base.Fix2(ldexp, -n) ∘ float, DP)
end

function signranklogcdf(n::Int, W::Union{Float64,Int})
    return log(signrankcdf(n, W))
end

function signrankccdf(n::Int, W::Union{Float64,Int})
    max_W = (n * (n + 1)) >> 1
    W = round(Int, W, RoundNearestTiesUp)
    return signrankcdf(n, max_W - W - 1,)
end

function signranklogccdf(n::Int, W::Union{Float64,Int})
    return log(signrankccdf(n, W))
end

function signrankinvcdf(n::Int, p::Float64)
    if p < 0.0 || p > 1.0
        return NaN
    end
    W = 0
    while signrankcdf(n, W) < p # TODO binary search and symmetry
        W += 1
    end
    return float(W)
end

function signrankinvlogcdf(n::Int, logp::Float64)
    if logp > 0.0 || logp == -Inf
        return NaN
    end
    W = 0
    while signranklogcdf(n, W) < logp # TODO binary search and symmetry
        W += 1
    end
    return float(W)
end

function signrankinvccdf(n::Int, p::Float64)
    signrankinvcdf(n, 1 - p)
end

function signrankinvlogccdf(n::Int, logp::Float64)
    if logp == -Inf
        return NaN
    end
    if logp == 0
        return 0.0
    end
    #signrankinvlogcdf(n, log1mexp(logp)) # does not roundtrip well
    signrankinvccdf(n, exp(logp))
end
