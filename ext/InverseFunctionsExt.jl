module InverseFunctionsExt

using StatsFuns
import InverseFunctions

InverseFunctions.inverse(::typeof(laplacecdf)) = laplaceinvcdf
InverseFunctions.inverse(::typeof(laplaceinvcdf)) = laplacecdf

InverseFunctions.inverse(::typeof(laplaceccdf)) = laplaceinvccdf
InverseFunctions.inverse(::typeof(laplaceinvccdf)) = laplaceccdf

InverseFunctions.inverse(::typeof(laplacelogcdf)) = laplaceinvlogcdf
InverseFunctions.inverse(::typeof(laplaceinvlogcdf)) = laplacelogcdf

InverseFunctions.inverse(::typeof(laplacelogccdf)) = laplaceinvlogccdf
InverseFunctions.inverse(::typeof(laplaceinvlogccdf)) = laplacelogccdf

InverseFunctions.inverse(::typeof(normcdf)) = norminvcdf
InverseFunctions.inverse(::typeof(norminvcdf)) = normcdf

InverseFunctions.inverse(::typeof(normccdf)) = norminvccdf
InverseFunctions.inverse(::typeof(norminvccdf)) = normccdf

InverseFunctions.inverse(::typeof(normlogcdf)) = norminvlogcdf
InverseFunctions.inverse(::typeof(norminvlogcdf)) = normlogcdf

InverseFunctions.inverse(::typeof(normlogccdf)) = norminvlogccdf
InverseFunctions.inverse(::typeof(norminvlogccdf)) = normlogccdf

end # module
