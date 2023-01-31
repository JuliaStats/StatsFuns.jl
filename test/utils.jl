# default relative tolerance for comparisons with Rmath
function _default_rtol(params, X::AbstractArray)
    # has to take into account `params` as well since otherwise e.g. `X::Array{<:Rational}`
    # always uses a tolerance based on `eps(one(Float64))` even when parameters are of type
    # Float32
    # eps^(7//8) means requiring equality of about 7/8 of the significant digits
    # Corresponds to tolerances of ~2e-14 (Float64), ~9f-7 (Float32) and ~0.002 (Float16)
    return eps(float(one(promote_type(Base.promote_typeof(params...), eltype(X)))))^(7//8)
end
