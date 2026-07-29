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
the probability that a subset of {1,2,...,j} sums to W-i+1,
i.e. the number of such subsets divided by the 2^j subsets in total.

The counts themselves grow like 2^n and overflow `Int` well before the
distribution stops being of interest, so the recursion is normalised as it goes
rather than divided by 2^n at the end. Halving is exact in binary floating
point, so while every count stays below 2^53 this gives exactly the counts
divided by 2^n; beyond that the additions round, but the entries stay in [0,1]
and can never overflow.
 =#

function signrankDP(n, W)
    DP = zeros(Float64, W + 1)
    DP[W + 1] = 1.0
    for j in 1:n
        for i in 1:(W + 1 - j)
            DP[i] = (DP[i] + DP[i + j]) / 2
        end
        # Sums smaller than j cannot contain j, so their counts are unchanged and
        # the corresponding probabilities only pick up the factor 1/2
        for i in max(1, W + 2 - j):(W + 1)
            DP[i] /= 2
        end
    end
    return DP
end

function signrankpdf(n::Int, W::Float64)
    return isinteger(W) ? signrankpdf(n, Int(W)) : 0.0
end
function signrankpdf(n::Int, W::Int)
    if W < 0
        return 0.0
    end
    max_W = (n * (n + 1)) >> 1
    W2 = max_W - W
    if W2 < W
        return signrankpdf(n, W2)
    end
    DP = signrankDP(n, W)
    return DP[1]
end

function signranklogpdf(n::Int, W::Union{Float64, Int})
    return log(signrankpdf(n, W))
end

function signrankcdf(n::Int, W::Float64)
    return signrankcdf(n, round(Int, W, RoundNearestTiesUp))
end
function signrankcdf(n::Int, W::Int)
    if W < 0
        return 0.0
    end
    max_W = (n * (n + 1)) >> 1
    W2 = max_W - W - 1
    if W2 < W
        return 1.0 - signrankcdf(n, W2)
    end
    DP = signrankDP(n, W)
    return sum(DP)
end

function signranklogcdf(n::Int, W::Union{Float64, Int})
    return log(signrankcdf(n, W))
end

function signrankccdf(n::Int, W::Float64)
    return signrankccdf(n, round(Int, W, RoundNearestTiesUp))
end
function signrankccdf(n::Int, W::Int)
    max_W = (n * (n + 1)) >> 1
    return signrankcdf(n, max_W - W - 1)
end

function signranklogccdf(n::Int, W::Union{Float64, Int})
    return log(signrankccdf(n, W))
end

function signrankinvcdf(n::Int, p::Float64)
    if !(0.0 <= p <= 1.0)
        return NaN
    end
    W = 0
    while signrankcdf(n, W) < p # TODO binary search and symmetry
        W += 1
    end
    return float(W)
end

function signrankinvlogcdf(n::Int, logp::Float64)
    if !(-Inf < logp <= 0.0)
        return NaN
    end
    W = 0
    while signranklogcdf(n, W) < logp # TODO binary search and symmetry
        W += 1
    end
    return float(W)
end

function signrankinvccdf(n::Int, p::Float64)
    if !(0.0 <= p <= 1)
        return NaN
    end
    if p == 0.0
        max_W = (n * (n + 1)) >> 1
        return float(max_W)
    end
    W = 0
    while signrankccdf(n, W) > p # TODO binary search and symmetry
        W += 1
    end
    return float(W)
end

function signrankinvlogccdf(n::Int, logp::Float64)
    if !(-Inf < logp <= 0.0)
        return NaN
    end
    W = 0
    while signranklogccdf(n, W) > logp # TODO binary search and symmetry
        W += 1
    end
    return float(W)
end
