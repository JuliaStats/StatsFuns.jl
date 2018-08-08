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
    n < 0 && throw(ArgumentError("n = $n must be a non-zero positive integer"))
    ( x > n || x < 0 ) && throw(ArgumentError("x = $x must ∈ [0, n]"))
    (p < 0.0 || p > 1.0) && throw(ArgumentError("p = $p must ∈ [0, 1]"))
    n == 0 && return one(T)
    p == 0.0 && return ( (x == 0) ? one(T): zero(T) )
    p == 1.0 && return ( (x == n) ? zero(T) : zero(T) )
    x == 0 && return exp(n*log1p(-p))
    x == n && return p^n
    p = Float64(p)
    lc = lstirling_asym(n) - lstirling_asym(x) - lstirling_asym(n - x) -D(x, n*p) - D(n - x, n*(1 - p))
    return  T(exp(lc)*sqrt(n/(2π*x*(n - x))))
end

# Deviance term: D(x, np) = x*log(x/np) + np - x
D(x, np) = -x*logmxp1(np/x)

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
    n > typemax(Int64) && error("n is too large")
    n < 0 && throw(ArgumentError("n = $n must be a non-zero positive integer"))
    ( x > n || x < 0 ) && throw(ArgumentError("x = $x must ∈ [0, n]"))
    (p < 0.0 || p > 1.0) && throw(ArgumentError("p = $p must ∈ [0, 1]"))
    n == 0 && return zero(T)
    p == 0.0 && return ( (x == 0) ? zero(T) : T(-Inf) )
    p == 1.0 && return ( (x == n) ? zero(T) : T(-Inf) )
    x == 0 && return n*log1p(-p)
    x == n && return n*log(p)
    p = Float64(p)
    (n, x) = (Int64(n), Int64(x))
    lc = lstirling_asym(n) - lstirling_asym(x) - lstirling_asym(n - x) - D(x, n*p) - D(n - x, n*(1.0 - p))
    return  T(lc + 0.5*log(n/(2π*x*(n - x))))
end
