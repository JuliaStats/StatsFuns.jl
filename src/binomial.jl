# Computation of binomial probability distribution function using Catherine Loader's Saddle point algorithm (http://octave.1599824.n4.nabble.com/attachment/3829107/0/loader2000Fast.pdf)
# binompdf(n, p, x)
# log(binompdf(n, p, x)) = log(binompdf(n, x/n, x)) - Dev(n, p, x)
# Deviance term 'Dev': Dev(n, p, x) = xlog(x/(np)) + (n-x)log((n - x)/(n*(1 - p)))
# Stirling's approximation: n! = √(2*π*n) * (n/e)^n + δ(n), where δ(n) is the error in the stirling approximation
# δ(n) ≈ lstirling_asym(n), the Asymptotic stirling series expansion error, https://dlmf.nist.gov/5.11, lstirling_asym defined in misc.jl
# p(n, x/n, x) = √(n/(2πx(n - x)))*e^(δ(n) - δ(x) - δ(n - x)), log(p(n, x/n, x)) = 0.5*log(n/(2*pi*x*(n - x))) + δ(n) - δ(x) - δ(n - x)
"""
    binompdf(n::Integer, p::Float, x::Integer)

Computes the binomial probability distibution function.
Given probability *`p`* of for success of each trial, returns the probability that *`x`* out of *`n`* trials are successful.

# Examples
```julia-repl
julia> binompdf(13, 0.58, 7)
0.20797396939077062
```
# Arguments
- `n::Integer`: total number of trials.
- `p::Float`: probability of success of each trial.
- `x::Integer`: number of successful trials.
...
"""
function binompdf{ T <: Union{Float16, Float32, Float64} }(n::Integer, p::T, x::Integer)
    n > typemax(Int64) && error("n is too large.")
    n < 1 && throw(ArgumentError("n = $n must be a positive integer"))
    ( x > n || x < 0 ) && throw(ArgumentError("x = $x must ∈ [0, n]"))
    (p < 0.0 || p > 1.0) && throw(ArgumentError("p = $p must ∈ [0, 1]"))
    p == 0.0 && return ( (x == 0) ? T(1.0) : T(0.0) )
    p == 1.0 && return ( (x == n) ? T(1.0) : T(0.0) )
    x == 0 && return exp(n*log1p(-p))
    x == n && return p^n
    p = Float64(p)
    (n, x) = (Int64(n), Int64(x))
    lc = stirling_err(n) - stirling_err(x) - stirling_err(n - x) -D(x, n*p) - D(n - x, n*(1 - p))
    return  T(exp(lc)*sqrt(n/(2π*x*(n - x))))
end

#
const Stirlingerror = Float64[log(factorial(n)) - log(sqrt(2pi * n) * (n/Base.e)^n) for n in big(1):big(20)]
function stirling_err(n)
    n <= 20 && return Stirlingerror[n]
    return lstirling_asym(n)
end
# Deviance term: D(x, np) = x*log(x/np) + np - x
function D(x::Int64, np::Float64)
    if abs(x - np) < 0.1*(x + np)
        s = (x - np)*(x - np)/(x + np)
        v = (x - np)/(x + np)
        ej = 2*x*v
        j = 1.0
        s1 = 0.0
        while true
            ej *= v*v
            s1 = s + ej/(2*j + 1)
            if s1 == s return s1 end
            s = s1
            j += 1
        end
    end
    return x*log(x/np) + np - x
end

# binomlogpdf(n, p, x) = log(binompdf(n, p, x))
# We use the same strategy as above but do not exponentiate the final result
"""
    binomlogpdf(n::Integer, p::Float, x::Integer)

Computes the log of the binomial probability distibution function.

# Examples
```julia-repl
julia> binomlogpdf(13, 0.58, 7)
-1.5703423542721349
```
# Arguments
- `n::Integer`: total number of trials.
- `p::Float`: probability of success of each trial.
- `x::Integer`: number of successful trials.
...
"""
function binomlogpdf{ T <: Union{Float16, Float32, Float64} }(n::Integer, p::T, x::Integer)
    n > typemax(Int64) && error("n is too large.")
    n < 1 && throw(ArgumentError("n = $n must be a positive integer"))
    ( x > n || x < 0 ) && throw(ArgumentError("x = $x must ∈ [0, n]"))
    (p < 0.0 || p > 1.0) && throw(ArgumentError("p = $p must ∈ [0, 1]"))
    p == 0.0 && return ( (x == 0) ? T(0.0) : T(-Inf) )
    p == 1.0 && return ( (x == n) ? T(0.0) : T(-Inf) )
    x == 0 && return n*log1p(-p)
    x == n && return n*log(p)
    p = Float64(p)
    (n, x) = (Int64(n), Int64(x))
    lc = stirling_err(n) - stirling_err(x) - stirling_err(n - x) - D(x, n*p) - D(n - x, n*(1.0 - p))
    return  T(lc + 0.5*log(n/(2π*x*(n - x))))
end
