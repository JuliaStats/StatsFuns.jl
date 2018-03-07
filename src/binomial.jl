#sfe(x) = log(n!*e^n/((n^n)*sqrt(2*pi*n))), for stirlerr
const sfe = Float64[log(factorial(n)) - log(sqrt(2pi * n) * (n/Base.e)^n) for n in big(1):big(20)]
# stirlerr(x) = log(x!) -log(sqrt(2*pi*x)*(n/e)^n)
function stirlerr(n::Float64)
    if (n <= 15)
        @inbounds return(sfe[Int64(n)])
    end
    return lstirling_asym(n)
end

stirlerr(n::Int64) = stirlerr(Float64(n))

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
function binompdf(n::Integer, p::AbstractFloat, x::Integer)
    # Using Saddle Point Algorithm by Catherine Loader which can be found here: http://octave.1599824.n4.nabble.com/attachment/3829107/0/loader2000Fast.pdf
    if (n  < 1 || x > n || x < 0) return NaN end
    if ( p < 0.0 || p > 1.0) return NaN end
    if ( p == 0.0 ) return ( (x == 0) ? 1.0 : 0.0) end
    if ( p == 1.0 ) return ( (x == n) ? 1.0 : 0.0) end
    if ( x == 0 ) return exp(n*log1p(-p)) end
    if ( x == n ) return exp(n*log(p)) end
    if n > typemax(Int64)
        error("n is too large.")
    end
    n = convert(Int64, n)
    p = convert(Float64, p)
    x = convert(Int64, x)
    lc = stirlerr(n) - stirlerr(x) - stirlerr(n-x) -D(x, n*p) - D(n-x, n*(1 - p))
    return  exp(lc)*sqrt(n/(2*pi*x*(n - x)))
end

function D(x::Int64, np::Float64) # Deviance term, x*log(x/np) + np - x
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
function binomlogpdf(n::Integer, p::AbstractFloat, x::Integer)
    if (n  < 1 || x > n || x < 0) return NaN end
    if ( p < 0.0 || p > 1.0) return NaN end
    if ( p == 0.0 ) return ( (x == 0) ? 0.0 : -Inf) end
    if ( p == 1.0 ) return ( (x == n) ? 0.0 : -Inf) end
    if ( x == 0 ) return n*log(1 - p) end
    if ( x == n ) return n*log(p) end
    if n > typemax(Int64)
        error("`n` is too large.")
    end
    n = convert(Int64, n)
    p = convert(Float64, p)
    x = convert(Int64, x)
    lc = stirlerr(n) - stirlerr(x) - stirlerr(n - x) - D(x, n*p) - D(n - x, n*(1.0 - p))
    return  lc + 0.5*log(n/(2*pi*x*(n - x)))
end
