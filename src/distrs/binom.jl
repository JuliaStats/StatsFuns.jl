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
binompdf(n::Real, p::Real, k::Number) = exp(binomlogpdf(n, p, k))

# logpdf for numbers with generic types
binomlogpdf(n::Real, p::Real, k::Number) = begin
    # NOTE: fail to do below
    #       sum(map(i -> log((n + 1 - i) / i), 1:k)) + k * log(p) + (n - k) * log(1 - p)
    log_n_choose_k = 0; i = 1;
    while i <= k log_n_choose_k += log((n + 1 - i) / i); i += 1 end
    log_n_choose_k + k * log(p) + (n - k) * log(1 - p)
end
