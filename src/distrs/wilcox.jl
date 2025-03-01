#=
The wilcox distribution is equivalent to the problem
of finding the number of subsets with size nx
of {1,2,...,nx,nx+1,...,nx+ny} summing to (U-sum(1:nx))
relative to the total number of subsets with size nx.
The empty subset is defined to sum to zero.
This can be calculated using the recursion:
either nx+ny is in the subset in witch case we need to calculate
the number of subsets with size nx-1 of {1,2,...,nx,nx+1,...,nx+ny-1}
summing to (U-sum(1:nx)-nx-ny),
or nx is not in the subset in which case we need to calculate
the number of subsets with size nx of {1,2,...,nx,nx+1,...,nx+ny-1}
summing to (U-sum(1:nx)),
This can be calculated bottom up using dynamic programming.

The (i,j)'th element of DP in the k'th outer loop iteration represents:
the number of subsets with size k of {1,2,...,k,k+1,...,k+j-1} summing to i.
 =#

@inline function wilcoxDP(nx, ny, U)
    DP = zeros(Int, sum(1:nx)+U, ny + 1)
    for i in 1:(nx+ny)
        for j in i:(ny+1)
            DP[i, j] = 1
        end
    end
    for k in 2:nx
        for j in 1:(ny+1)
            for i in min(U + sum(1:k), sum((j):(k+j-1))):-1:1
                if j > 1
                    DP[i, j] = DP[i, j-1]
                else
                    DP[i, j] = 0
                end
                if i - j - k + 1 >= 1
                    DP[i, j] += DP[i-j-k+1, j]
                end
            end
        end
    end
    return DP
end

function wilcoxpdf(nx::Int, ny::Int, U::Float64)
    return isinteger(U) ? wilcoxpdf(nx, ny, Int(W)) : 0.0
end
function wilcoxpdf(nx::Int, ny::Int, U::Int)
    if U < 0
        return 0.0
    end
    max_U = sum((ny+1):(ny+nx)) - sum(1:nx)
    U2 = max_U - U
    if U2 < U
        return wilcoxpdf(nx, ny, U2)
    end
    DP = wilcoxDP(nx, ny, U)
    DP[sum(1:nx)+U, ny+1] / binomial(nx+ny, nx)
end

function wilcoxlogpdf(nx::Int, ny::Int, U::Union{Float64,Int})
    return log(wilcoxpdf(nx, ny, U))
end

function wilcoxcdf(nx::Int, ny::Int, U::Float64)
    wilcoxcdf(nx, ny, round(Int, U, RoundNearestTiesUp))
end
function wilcoxcdf(nx::Int, ny::Int, U::Int)
    if U < 0
        return 0.0
    end
    max_U = sum((ny+1):(ny+nx)) - sum(1:nx)
    U2 = max_U - U - 1
    if U2 < U
        return 1.0 - wilcoxcdf(nx, ny, U2)
    end
    DP = wilcoxDP(nx, ny, U)
    sum(DP[(sum(1:nx)):(sum(1:nx)+U), ny+1] ./ binomial(nx+ny, nx))
end

function wilcoxlogcdf(nx::Int, ny::Int, U::Union{Float64,Int})
    max_U = sum((ny+1):(ny+nx)) - sum(1:nx)
    U2 = max_U - U - 1
    if U2 < U
        return log1p(-wilcoxcdf(nx, ny, U2))
    else
        return log(wilcoxcdf(nx, ny, U))
    end
end

function wilcoxccdf(nx::Int, ny::Int, U::Float64)
    wilcoxccdf(nx, ny, round(Int, U, RoundNearestTiesUp))
end
function wilcoxccdf(nx::Int, ny::Int, U::Int)
    max_U = sum((ny+1):(nx+ny)) - sum(1:nx)
    wilcoxcdf(nx, ny, max_U - U - 1)
end

function wilcoxlogccdf(nx::Int, ny::Int, U::Union{Float64,Int})
    max_U = sum((ny+1):(ny+nx)) - sum(1:nx)
    U2 = max_U - U - 1
    if U2 < U+2 # +2 is emperical
        return log(wilcoxccdf(nx, ny, U))
    else
        return log1p(-wilcoxccdf(nx, ny, U2))
    end
    return 
end

function wilcoxinvcdf(nx::Int, ny::Int, p::Float64)
    if !(0.0 <= p <= 1.0)
        return NaN
    end
    U = 0
    while wilcoxcdf(nx, ny, U) < p # TODO binary search and symmetry
        U += 1
    end
    return float(U)
end

function wilcoxinvlogcdf(nx::Int, ny::Int, logp::Float64)
    if !(-Inf < logp <= 0.0)
        return NaN
    end
    U = 0
    while wilcoxlogcdf(nx, ny, U) < logp # TODO binary search and symmetry
        U += 1
    end
    return float(U)
end

function wilcoxinvccdf(nx::Int, ny::Int, p::Float64)
    if !(0.0 <= p <= 1)
        return NaN
    end
    if p == 0.0
        max_U = sum((ny+1):(ny+nx)) - sum(1:nx)
        return float(max_U)
    end
    U = 0
    while wilcoxccdf(nx, ny, U) > p # TODO binary search and symmetry
        U += 1
    end
    return float(U)
end

function wilcoxinvlogccdf(nx::Int, ny::Int, logp::Float64)
    if !(-Inf < logp <= 0.0)
        return NaN
    end
    U = 0
    while wilcoxlogccdf(nx, ny, U) > logp # TODO binary search and symmetry
        U += 1
    end
    return float(U)
end
