function stirlerr(n::Int64)
    # stirlerr(x) = log(x!) -log(sqrt(2*pi*x)*(n/e)^n)
    const S0 = 0.083333333333333333333        # 1/12
    const S1 = 0.00277777777777777777778      # 1/360
    const S2 = 0.00079365079365079365079365   # 1/1260
    const S3 = 0.000595238095238095238095238  # 1/1680
    const S4 = 0.0008417508417508417508417508 # 1/1188
    # 30 values, sfe(x) = log(n!*e^n/((n^n)*sqrt(2*pi*n)))
    sfe = (0.1534264097200273452913848, 0.0810614667953272582196702, 0.0548141210519176538961390,
    0.0413406959554092940938221, 0.03316287351993628748511048, 0.02767792568499833914878929,
	0.02374616365629749597132920, 0.02079067210376509311152277, 0.01848845053267318523077934,
	0.01664469118982119216319487, 0.01513497322191737887351255, 0.01387612882307074799874573,
	0.01281046524292022692424986, 0.01189670994589177009505572, 0.01110455975820691732662991,
    0.010411265261972096497478567,0.009799416126158803298389475, 0.009255462182712732917728637,
    0.008768700134139385462952823, 0.008330563433362871256469318,0.007934114564314020547248100,
    0.007573675487951840794972024, 0.007244554301320383179543912, 0.006942840107209529865664152,
    0.006665247032707682442354394, 0.006408994188004207068439631,0.006171712263039457647532867,
    0.005951370112758847735624416,0.005746216513010115682023589, 0.005554733551962801371038690)
    if (n <= 15)
        @inbounds return(sfe[n + n])
    end
    nn = n*n
    if (n > 500) return((S0-S1/nn)/n) end
    if (n > 80) return((S0-(S1-S2/nn)/nn)/n) end
    if (n > 35) return((S0-(S1-(S2-S3/nn)/nn)/nn)/n) end
    return((S0-(S1-(S2-(S3-S4/nn)/nn)/nn)/nn)/n)
end

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
    if ( p == 0.0 ) return ( (x == 0) ? Float64(1.0):Float64(0.0)) end
    if ( p == 1.0 ) return ( (x == n) ? Float64(1.0):Float64(0.0)) end
    if ( x == 0 ) return exp(n*log(1-p)) end
    if ( x == n ) return exp(n*log(p)) end
    n = convert(Int64, min(typemax(Int64), n))
    p = convert(Float64, max(eps(Float64), p))
    x = convert(Int64, min(typemax(Int64), x))
    lc = stirlerr(n) - stirlerr(x) - stirlerr(n-x) -D(x, n*p) - D(n-x, n*(1.0-p))
    return  exp(lc)*sqrt(n/(2*pi*x*(n-x)))
end

function D(x::Int64, np::Float64) # Deviance term, x*log(x/np) + np - x
    if( abs(x - np) < 0.1*(x+np))
        s = (x - np)*(x - np)/(x + np)
        v = (x - np)/(x + np)
        ej = 2*x*v
        j = 1
        s1 = 0
        while true
            ej *= v*v
            s1 = s+ej/(2*j+1)
            if( s1 == s) return s1 end
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
    return log(binompdf(n, p, x))
end
