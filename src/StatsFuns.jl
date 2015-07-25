module StatsFuns

using Compat
import Compat.@irrational

export
    # basicfuns
    xlogx,       # x * log(x)
    xlogy,       # x * log(y)
    logistic,    # 1 / (1 + exp(-x))
    logit,       # log(x / (1 - x))
    softplus,    # log(1 + exp(x))
    invsoftplus, # log(exp(x) - 1)
    logsumexp,   # log(exp(x) + exp(y)) or log(sum(exp(x)))
    softmax,
    softmax!


## source files
include("constants.jl")
include("basicfuns.jl")

end # module
