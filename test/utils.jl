# default relative tolerance for comparisons with Rmath
function _default_rtol(params, X::AbstractArray)
    # has to take into account `params` as well since otherwise e.g. `X::Array{<:Rational}`
    # always uses a tolerance based on `eps(one(Float64))` even when parameters are of type
    # Float32
    return _default_rtol(float(promote_type(Base.promote_typeof(params...), eltype(X))))
end

# We use less sharp tolerances for Float16 and Float32 since the proportion of significant
# digits that are equal when evaluating Rmath and StatsFuns is smaller
# eps^(7//8) means requiring equality of about 7/8 of the significant digits, etc.
_default_rtol(::Type{Float64}) = eps(Float64)^(7//8) # ~2.0e-14
_default_rtol(::Type{Float32}) = eps(Float32)^(3//4) # ~6.4e-6
_default_rtol(::Type{Float16}) = eps(Float16)^(2//3) # ~9.9e-3 
