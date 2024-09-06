module StatsFunsInverseFunctionsExt

using StatsFuns
import InverseFunctions

InverseFunctions.inverse(::typeof(normcdf)) = norminvcdf
InverseFunctions.inverse(::typeof(norminvcdf)) = normcdf

InverseFunctions.inverse(::typeof(normccdf)) = norminvccdf
InverseFunctions.inverse(::typeof(norminvccdf)) = normccdf

InverseFunctions.inverse(::typeof(normlogcdf)) = norminvlogcdf
InverseFunctions.inverse(::typeof(norminvlogcdf)) = normlogcdf

InverseFunctions.inverse(::typeof(normlogccdf)) = norminvlogccdf
InverseFunctions.inverse(::typeof(norminvlogccdf)) = normlogccdf

InverseFunctions.inverse(::typeof(logit)) = logistic
InverseFunctions.inverse(::typeof(logistic)) = logit

InverseFunctions.inverse(::typeof(log1mexp)) = log1mexp

InverseFunctions.inverse(::typeof(log1pexp)) = logexpm1
InverseFunctions.inverse(::typeof(logexpm1)) = log1pexp

end # module
