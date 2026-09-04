#=
Compute the number of sequences of `nx` 0s and `ny` 1s in each of which a 1 precedes a 0 `U` times.

## Details

Due to symmetry, we only have to consider the case `nx ≤ ny`.
Let `m = min(nx, ny)` and `n = max(nx, ny)`, and denote the number of sequences of `m` 0s and `n` 1s
in which a 1 precedes a 0 `U` times by `pₘ,ₙ(U)`.

Mann and Whitney (1947) computed `pₘ,ₙ(U)` by exploiting the recurrence relation

pₘ,ₙ(U) = pₘ₋₁,ₙ(U - n) + pₘ,ₙ₋₁(U)

It can be obtained by considering the sequence with the last element removed:
1. If the last element is a 0, it is preceded by `n` 1s;
   and hence in the sequence with the last element removed, consisting of `m - 1` 0s and `n` 1s, a 1 must precede a 0 `U - n` times.
2. If the last element is a 1, it does not precede any 0;
   and hence in the sequence with the last element removed, consisting of `m` 0s and `n - 1` 1s, a 1 must preced a 0 `U` times.

Löffler (1983) published the recurrence relation

pₘ,ₙ(U) = 1/U \sum_{a = 0}^{U - 1} σₘ,ₙ(U - a) pₘ,ₙ(a)

where

σₘ,ₙ(k) := \sum_{d | k} ϵₘ,ₙ(d) d

with ϵₘ,ₙ(d) = 1 for 1 ≤ d ≤ m, ϵₘ,ₙ(d) = -1 for n < d ≤ m + n, and ϵₘ,ₙ(d) = 0 otherwise.
Implementation of this recurrence relation allows faster computations with fewer memory allocations.

## References

H. B. Mann, D. R. Whitney. "On a Test of Whether one of Two Random Variables is Stochastically Larger than the Other." Ann. Math. Statist. 18 (1) 50 - 60, March, 1947. https://doi.org/10.1214/aoms/1177730491 
A. Löffler: "Über eine Partition der nat. Zahlen und ihre Anwendung beim U-Test." Wissenschaftliche Zeitschrift der Martin-Luther-Universität Halle-Wittenberg; Mathematisch-Naturwissenschaftliche Reihe, XXXII'83 M, Heft 5, 87–89; available as https://upload.wikimedia.org/wikipedia/commons/f/f5/LoefflerWilcoxonMannWhitneyTest.pdf
=#

#=
The recurrence below is linear and homogeneous in `pₘ,ₙ`, so seeding it with
`1 / binomial(nx + ny, nx)` instead of `1` yields the probabilities directly.
This matters because the counts themselves grow like `binomial(nx + ny, nx)` and overflow
`Int` from `nx + ny = 68` onwards, whereas the probabilities are bounded by 1.
=#
function wilcox_probabilities(nx::Int, ny::Int, U::Int)
    # This internal function expects 0 <= U <= nx * ny / 2
    if !(0 <= U <= (nx * ny) / 2)
        throw(ArgumentError("`wilcox_probabilities(nx, ny, U)` is only implemented for 0 <= U <= (nx * ny) / 2"))
    end

    # Due to symmetry, `wilcox_probabilities(nx, ny, U) = wilcox_probabilities(ny, nx, U)`
    # Hence for simplicity we only consider the case `wilcox_probabilities(min(nx, ny), max(nx, ny), U)` here
    m, n = minmax(nx, ny)

    # Compute σ(k) = ∑_{d|k} ϵ(d) d where
    # - ϵ(d) := 1 for 1 ≤ d ≤ m
    # - ϵ(d) := -1 for n < d ≤ m + n
    # - ϵ(d) := 0 otherwise
    sigmas = zeros(Int, U)
    for d in 1:m
        for i in d:d:U
            sigmas[i] += d
        end
    end
    for d in (n + 1):(m + n)
        for i in d:d:U
            sigmas[i] -= d
        end
    end

    # Recursively compute pₘ,ₙ(a) / binomial(nx + ny, nx) for 0 <= a <= U, in units of
    # 2^shift. The logarithm of the coefficient is used since `binomial(nx + ny, nx)` itself
    # overflows `Int` from `nx + ny = 68` onwards.
    #
    # The seed is the probability of the least likely outcome, so it leaves the normal range
    # long before the probabilities of interest do: it is subnormal from `nx + ny = 1030` and
    # zero from `nx + ny = 1090`, and since the recurrence is homogeneous a seed of zero
    # would silently zero every probability. Working in units of 2^shift avoids that. The
    # recurrence is homogeneous so the choice of unit cannot affect the result, `ldexp`
    # removes it from the scalar result exactly, and `shift` is zero whenever the seed is
    # normal.
    logbinom = first(logabsbinomial(nx + ny, nx))
    shift = max(0, ceil(Int, (logbinom - 700) / logtwo))
    probabilities = Vector{Float64}(undef, U + 1)
    probabilities[1] = exp(shift * logtwo - logbinom)
    for a in 1:U
        p = 0.0
        for i in 1:a
            p += probabilities[i] * sigmas[a + 1 - i]
        end
        probabilities[a + 1] = p / a
    end

    return probabilities, shift
end

function wilcoxpdf(nx::Int, ny::Int, U::Float64)
    return isinteger(U) ? wilcoxpdf(nx, ny, Int(U)) : 0.0
end
function wilcoxpdf(nx::Int, ny::Int, U::Int)
    max_U = nx * ny
    if !(0 <= U <= max_U)
        return 0.0
    end
    U = min(U, max_U - U)
    probabilities, shift = wilcox_probabilities(nx, ny, U)
    return ldexp(probabilities[end], -shift)
end

function wilcoxlogpdf(nx::Int, ny::Int, U::Union{Float64, Int})
    return log(wilcoxpdf(nx, ny, U))
end

function wilcoxcdf(nx::Int, ny::Int, U::Float64)
    return wilcoxcdf(nx, ny, round(Int, U, RoundNearestTiesUp))
end
function wilcoxcdf(nx::Int, ny::Int, U::Int)
    max_U = nx * ny
    if U < 0
        return 0.0
    elseif U >= max_U
        return 1.0
    else
        U2 = max_U - U - 1
        probabilities, shift = wilcox_probabilities(nx, ny, min(U, U2))
        p = ldexp(sum(probabilities), -shift)
        return U2 < U ? 1.0 - p : p
    end
end

function wilcoxlogcdf(nx::Int, ny::Int, U::Union{Float64, Int})
    max_U = nx * ny
    U2 = max_U - U - 1
    if U2 < U
        return log1p(-wilcoxcdf(nx, ny, U2))
    else
        return log(wilcoxcdf(nx, ny, U))
    end
end

function wilcoxccdf(nx::Int, ny::Int, U::Float64)
    return wilcoxccdf(nx, ny, round(Int, U, RoundNearestTiesUp))
end
function wilcoxccdf(nx::Int, ny::Int, U::Int)
    max_U = nx * ny
    return wilcoxcdf(nx, ny, max_U - U - 1)
end

function wilcoxlogccdf(nx::Int, ny::Int, U::Union{Float64, Int})
    max_U = nx * ny
    U2 = max_U - U - 1
    if U2 < U + 2 # +2 is empirical
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
        max_U = nx * ny
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
