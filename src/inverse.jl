for d in (:norm,)
    @eval begin
        InverseFunctions.inverse(::typeof($(Symbol(d,:cdf)))) = $(Symbol(d, :invcdf))
        InverseFunctions.inverse(::typeof($(Symbol(d,:ccdf)))) = $(Symbol(d, :invccdf))
        InverseFunctions.inverse(::typeof($(Symbol(d,:logcdf)))) = $(Symbol(d, :invlogcdf))
        InverseFunctions.inverse(::typeof($(Symbol(d,:logccdf)))) = $(Symbol(d, :invlogccdf))
        InverseFunctions.inverse(::typeof($(Symbol(d,:invcdf)))) = $(Symbol(d, :cdf))
        InverseFunctions.inverse(::typeof($(Symbol(d,:invccdf)))) = $(Symbol(d, :ccdf))
        InverseFunctions.inverse(::typeof($(Symbol(d,:invlogcdf)))) = $(Symbol(d, :logcdf))
        InverseFunctions.inverse(::typeof($(Symbol(d,:invlogccdf)))) = $(Symbol(d, :logccdf))
    end
end
