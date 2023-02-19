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

end # module
