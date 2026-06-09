# Owen's T function
#
# Ported from Boost owen_t.hpp, which is Copyright Benjamin Sobotta 2012
#      Use, modification and distribution are subject to the
#      Boost Software License, Version 1.0. (http://www.boost.org/LICENSE_1_0.txt)
# Port by Steven G. Johnson, also under Boost Software License, Version 1.0.
#
# Algorithm based on Mike Patefield and David Tandy, "Fast and accurate calculation
# of Owen's T-function," Journal of Statistical Software, 5, pp. 1-25 (2000).

# P(-∞<Z<=x)-0.5 with Z being normally distributed.
owens_t_znorm1(x::Real) = erf(x * invsqrt2) / 2

# P(x<=Z<∞) with Z being normally distributed.
owens_t_znorm2(x::Real) = erfc(x * invsqrt2) / 2

# Auxiliary function, it computes an array key that is used to determine
# the specific computation method for Owen's T and the order thereof
# used in owens_t_dispatch.  Differs from C++ in using 1-based indices
let hrange = (0.02f0, 0.06f0, 0.09f0, 0.125f0, 0.26f0, 0.4f0, 0.6f0, 1.6f0, 1.7f0, 2.33f0, 2.4f0, 3.36f0, 3.4f0, 4.8f0), # length 14
        arange = (0.025f0, 0.09f0, 0.15f0, 0.36f0, 0.5f0, 0.9f0, 0.99999f0), # length 7
        select = (
        1, 1, 2, 13, 13, 13, 13, 13, 13, 13, 13, 16, 16, 16, 9,
        1, 2, 2, 3, 3, 5, 5, 14, 14, 15, 15, 16, 16, 16, 9,
        2, 2, 3, 3, 3, 5, 5, 15, 15, 15, 15, 16, 16, 16, 10,
        2, 2, 3, 5, 5, 5, 5, 7, 7, 16, 16, 16, 16, 16, 10,
        2, 3, 3, 5, 5, 6, 6, 8, 8, 17, 17, 17, 12, 12, 11,
        2, 3, 5, 5, 5, 6, 6, 8, 8, 17, 17, 17, 12, 12, 12,
        2, 3, 4, 4, 6, 6, 8, 8, 17, 17, 17, 17, 17, 12, 12,
        2, 3, 4, 4, 6, 6, 18, 18, 18, 18, 17, 17, 17, 12, 12,
    ) # 1-based indices (in 8x15 "matrix")
    global owens_t_compute_code
    function owens_t_compute_code(h::Real, a::Real)
        ihint = something(findfirst(>=(h), hrange), length(hrange) + 1)
        iaint = something(findfirst(>=(a), arange), length(arange) + 1)
        return select[(iaint - 1) * 15 + ihint]
    end
end

let ord = (2, 3, 4, 5, 7, 10, 12, 18, 10, 20, 30, 0, 4, 7, 8, 20, 0, 0) # 53-bit precision table
    global owens_t_get_order
    function owens_t_get_order(icode::Integer, ::Type{Float64})
        return ord[icode]
    end
end

# compute the value of Owen's T function with method T1 from the reference paper
function owens_t_T1(h::T, a::T, m::Integer) where {T <: Real}
    hs = -h * h / 2
    as = a * a
    aj = a * inv2π
    dj = expm1(hs)
    gj = hs * exp(hs)
    val = atan(a) * inv2π
    j, jj = one(m), T(1)
    while true
        val += dj * aj / jj
        m <= j && break
        j += one(m)
        jj += 2
        aj *= as
        dj = gj - dj
        gj *= hs / j
    end
    return val
end

# compute the value of Owen's T function with method T2 from the reference paper
function owens_t_T2(h::T, a::T, m::Integer, ah::T) where {T <: Real}
    maxii = m + m + one(m)
    hs = h * h
    as = -a * a
    y = inv(hs)
    ii = one(m)
    val = zero(T)
    vi = a * exp(-ah * ah / 2) * invsqrt2π
    z = owens_t_znorm1(ah) / h
    while true
        val += z
        maxii <= ii && break
        z = y * (vi - ii * z)
        vi *= as
        ii += oftype(ii, 2)
    end
    return val * exp(-hs / 2) * invsqrt2π
end

# compute the value of Owen's T function with method T3 from the reference paper
let c2 = (
        0.9999999999999998751, -0.99999999999988796462, 0.99999999998290743652,
        -0.99999999896282500134, 0.99999996660459362918, -0.9999993398627247676,
        0.99999125611136965852, -0.99991777624463387686, 0.99942835555870132569,
        -0.99697311720723000295, 0.98751448037275303682, -0.95915857980572882813,
        0.89246305511006708555, -0.76893425990463999675, 0.5889352846848469325,
        -0.38380345160440256652, 0.20317601701045299653, -0.82813631607004984866e-1,
        0.24167984735759576523e-1, -0.44676566663971825242e-2, 0.39141169402373836468e-3,
    )
    global owens_t_T3
    function owens_t_T3(h::Float64, a::Float64, ah::Float64)
        as = a * a
        hs = h * h
        y = inv(hs)
        ii = one(h)
        vi = a * exp(-ah * ah / 2) * invsqrt2π
        zi = owens_t_znorm1(ah) / h
        val = zero(h)
        i = 1
        while true
            val += zi * c2[i]
            i == lastindex(c2) && break
            i += 1
            zi = y * (ii * zi - vi)
            vi *= as
            ii += 2
        end
        return val * exp(-hs / 2) * invsqrt2π
    end
end

# compute the value of Owen's T function with method T4 from the reference paper
function owens_t_T4(h::T, a::T, m::Integer) where {T <: Real}
    maxii = m + m + one(m)
    hs = h * h
    as = -a * a
    ii = 1
    ai = a * exp(-hs * (1 - as) / 2) * inv2π
    yi = one(T)
    val = zero(T)
    while true
        val += ai * yi
        maxii <= ii && break
        ii += 2
        yi = (1 - hs * yi) / ii
        ai *= as
    end
    return val
end

# compute the value of Owen's T function with method T5 from the reference paper
let pts = (
        0.35082039676451715489e-2, 0.3127904233803075374e-1, 0.8526682628321945109e-1,
        0.16245071730812277011, 0.25851196049125434828, 0.36807553840697533536,
        0.48501092905604697475, 0.60277514152618576821, 0.71477884217753226516,
        0.81475510988760098605, 0.89711029755948965867, 0.95723808085944261843, 0.99178832974629703586,
    ),
        wts = (
        0.18831438115323502887e-1, 0.18567086243977649478e-1, 0.18042093461223385584e-1,
        0.17263829606398753364e-1, 0.1624321997598985673e-1, 0.14994592034116704829e-1,
        0.13535474469662088392e-1, 0.11886351605820165233e-1, 0.10070377242777431897e-1,
        0.81130545742299586629e-2, 0.60419009528470238773e-2, 0.38862217010742057883e-2, 0.16793031084546090448e-2,
    )
    #=
       NOTICE:
       - The pts[] array contains the squares (!) of the abscissas, i.e. the roots of the Legendre
         polynomial P_n(x), instead of the plain roots as required in Gauss-Legendre
         quadrature, because T5(h,a,m) contains only x^2 terms.
       - The wts[] array contains the weights for Gauss-Legendre quadrature scaled with a factor
         of 1/(2*pi) according to T5(h,a,m).
    =#
    global owens_t_T5
    function owens_t_T5(h::Float64, a::Float64)
        as = a * a
        hs = -h * h / 2
        r = 1 .+ as .* pts
        return sum(wts .* exp.(hs .* r) ./ r) * a
    end
end

# compute the value of Owen's T function with method T6 from the reference paper
function owens_t_T6(h::T, a::T) where {T <: Float64}
    normh = owens_t_znorm2(h)
    y = 1 - a
    r = atan(y, 1 + a)
    val = normh * (1 - normh) / 2
    !iszero(r) && (val -= r * exp(-y * h * h / 2r) * inv2π)
    return val
end

# This routine dispatches the call to one of six subroutines, depending on the values of h and a.
# preconditions: h >= 0, 0<=a<=1, ah=a*h.
# Simple main case for 64-bit precision or less, this is as per the Patefield-Tandy paper:
function owens_t_dispatch(h::Float64, a::Float64, ah::Float64)
    # Handle some special cases first, these are from
    # page 1077 of Owen's original paper:
    iszero(h) && return atan(a) * inv2π
    iszero(a) && return zero(h)
    a == 1 && return owens_t_znorm2(-h) * owens_t_znorm2(h) / 2
    @assert a <= 1 # when a>1 we call this routine with 1/a:

    icode = owens_t_compute_code(h, a)
    m = owens_t_get_order(icode, typeof(h))
    meth = (1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 4, 4, 4, 4, 5, 6) # 18 entries
    method = meth[icode]

    # determine the appropriate method, T1 ... T6
    method == 1 && return owens_t_T1(h, a, m)
    method == 2 && return owens_t_T2(h, a, m, ah)
    method == 3 && return owens_t_T3(h, a, ah)
    method == 4 && return owens_t_T4(h, a, m)
    method == 5 && return owens_t_T5(h, a)
    return owens_t_T6(h, a)
end

# compute Owen's T function, T(h,a), for arbitrary values of h and a
function _owens_t(h::Float64, a::Float64)
    # exploit that T(-h,a) == T(h,a)
    h = abs(h)

    # Use equation (2) in the paper to remap the arguments
    # such that h>=0 and 0<=a<=1 for the call of the actual
    # computation routine.
    abs_a = abs(a)
    abs_ah = abs_a * h

    val = if abs_a <= 1
        owens_t_dispatch(h, abs_a, abs_ah)
    elseif h <= oftype(h, 0.67)
        normh = owens_t_znorm1(h)
        normah = owens_t_znorm1(abs_ah)
        1 // 4 - normh * normah - owens_t_dispatch(abs_ah, inv(abs_a), h)
    else
        normh = owens_t_znorm2(h)
        normah = owens_t_znorm2(abs_ah)
        (normh + normah) / 2 - normh * normah - owens_t_dispatch(abs_ah, inv(abs_a), h)
    end
    return a < 0 ? -val : val # exploit that T(h,-a) == -T(h,a)
end

_owens_t(h::Float32, a::Float32) = Float32(Float64(h), Float64(a))
_owens_t(h::Float16, a::Float16) = Float16(Float64(h), Float64(a))

"""
    owens_t(h::Real, a::Real)

Returns Owen's T function

```math
T(h,a) = \\frac{1}{2\\pi} \\int_0^a \\frac{e^{-h^2(1+x^2)/2}{1 + x^2} dx
```

(This is the probability of ``X > h`` and ``0 < Y < aX``, where ``X`` and ``Y`` are
i.i.d. standard normal random variables.)
"""
owens_t(h::Real, a::Real) = _owens_t(promote(float(h), float(a))...)
