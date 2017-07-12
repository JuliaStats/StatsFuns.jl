# functions related to binomial distribution

import .RFunctions:
    binompdf,
    binomlogpdf,
    binomcdf,
    binomccdf,
    binomlogcdf,
    binomlogccdf,
    binominvcdf,
    binominvccdf,
    binominvlogcdf,
    binominvlogccdf

# pdf for numbers with generic types
binompdf(n::Real, p::Real, k::Number) = begin
    # NOTE: fail to do below
    #       prod(map(i -> (n + 1 - i) / i, 1:Int(k.))) * p^k * (1 - p)^(n - k)
    n_choose_k = 1; i = 1;
    while i <= k n_choose_k *= (n + 1 - i) / i; i += 1 end
    n_choose_k * p^k * (1 - p)^(n - k)
end

# logpdf for numbers with generic types
binomlogpdf(n::Real, p::Real, k::Number) = begin
    # NOTE: fail to do below
    #       sum(map(i -> log((n + 1 - i) / i), 1:k)) + k * log(p) + (n - k) * log(1 - p)
    log_n_choose_k = 0; i = 1;
    while i <= k log_n_choose_k += log((n + 1 - i) / i); i += 1 end
    log_n_choose_k + k * log(p) + (n - k) * log(1 - p)
end
