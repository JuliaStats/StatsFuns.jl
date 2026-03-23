# functions related to noncentral beta distribution
# Pure Julia implementation based on VBA code by Ian Smith

const _nc_limit = 1_000_000.0
const _cSmall = 5.562684646268003457725581793331e-309

"""
    _beta_nc1(x, y, a, b, nc) -> (cdf, nc_derivative)

Core CDF computation for the noncentral beta distribution via Poisson-weighted
sum of incomplete beta functions. `y = 1 - x` is passed separately for accuracy.
Returns the CDF value and the PDF (nc_derivative) as a tuple.
Based on VBA `beta_nc1` by Ian Smith.
"""
function _beta_nc1(x::Float64, y::Float64, a::Float64, b::Float64, nc::Float64)
    nc_derivative = 0.0

    # Find starting index n via quadratic formula
    bb = (x * nc - 1.0) - a
    if bb < -1.0e150
        n_over_bb = a / bb
        aa = n_over_bb - nc * x * (n_over_bb + b / bb)
        n_temp = bb * (1.0 + sqrt(1.0 - (4.0 * aa / bb)))
        n = floor(2.0 * aa * (bb / n_temp))
    else
        aa = a - nc * x * (a + b)
        if bb < 0.0
            n = bb - sqrt(bb^2 - 4.0 * aa)
            n = floor(2.0 * aa / n)
        else
            n = floor((bb + sqrt(bb^2 - 4.0 * aa)) / 2.0)
        end
    end
    if n < 0.0
        n = 0.0
    end

    aa = n + a
    bb_idx = n
    ptnc = _poisson_term(n, nc, nc - n, 0.0)
    # ptx = b * binomialTerm(aa, b, x, y, b*x - aa*y, 0) which equals
    # (aa+b)*(I(x,aa,b) - I(x,aa+1,b))
    ptx = b * _binomial_term(aa, b, x, y, ab_minus_cd(b, x, aa, y), 0.0)
    aa = aa + 1.0
    bb_idx = bb_idx + 1.0
    p = nc / bb_idx
    ps = p
    nc_derivative = ps
    s = x / aa  # (I(x, aa, b) - I(x, aa+1, b)) / ptx
    w = p
    term = s * w
    result = term

    if ptx > 0
        while ((term > 1.0e-15 * result) && (p > 1.0e-16 * w)) || (ps > 1.0e-16 * nc_derivative)
            s = (aa + b) * s
            aa = aa + 1.0
            bb_idx = bb_idx + 1.0
            p = nc / bb_idx * p
            ps = p * s
            nc_derivative = nc_derivative + ps
            s = x / aa * s
            w = w + p
            term = s * w
            result = result + term
        end
        w = w * ptnc
    else
        w = poisccdf(nc, n - 1.0)
    end

    if x > y
        s = betaccdf(b, a + (bb_idx + 1.0), y)
    else
        s = betacdf(a + (bb_idx + 1.0), b, x)
    end
    cdf_result = result * ptx * ptnc + s * w

    # Downward summation
    ps = 1.0
    nc_dtemp = 0.0
    aa = n + a
    bb_idx = n
    p = 1.0
    s = ptx / (aa + b) # I(x, aa, b) - I(x, aa+1, b)
    if x > y
        w = betaccdf(b, aa, y) # I(x, aa, b)
    else
        w = betacdf(aa, b, x)
    end
    term = p * w
    result = term

    while bb_idx > 0.0 && (((term > 1.0e-15 * result) && (s > 1.0e-16 * w)) || (ps > 1.0e-16 * nc_dtemp))
        s = aa / x * s
        ps = p * s
        nc_dtemp = nc_dtemp + ps
        p = bb_idx / nc * p
        aa = aa - 1.0
        bb_idx = bb_idx - 1.0
        if bb_idx == 0.0
            aa = a
        end
        s = s / (aa + b)
        w = w + s
        term = p * w
        result = result + term
    end

    if n > 0.0
        nc_dtemp = nc_derivative * ptx + nc_dtemp + p * aa / x * s
    elseif b == 0.0
        nc_dtemp = 0.0
    else
        nc_dtemp = _binomial_term(aa, b, x, y, ab_minus_cd(b, x, aa, y), log(b) + log(nc_derivative + aa / (x * (aa + b))))
    end
    nc_dtemp = nc_dtemp / y

    cdf_result = cdf_result + result * ptnc + poiscdf(nc, bb_idx - 1.0) * w

    if nc_dtemp == 0.0
        nc_derivative = 0.0
    else
        nc_derivative = _poisson_term(n, nc, nc - n, log(nc_dtemp))
    end

    return cdf_result, nc_derivative
end

"""
    _comp_beta_nc1(x, y, a, b, nc) -> (ccdf, nc_derivative)

Complementary CDF computation for the noncentral beta distribution.
`y = 1 - x` is passed separately for accuracy.
Returns the complementary CDF value and the PDF (nc_derivative) as a tuple.
Based on VBA `comp_beta_nc1` by Ian Smith.
"""
function _comp_beta_nc1(x::Float64, y::Float64, a::Float64, b::Float64, nc::Float64)
    nc_derivative = 0.0

    # Find starting index n via quadratic formula
    bb = (x * nc - 1.0) - a
    if bb < -1.0e150
        n_over_bb = a / bb
        aa = n_over_bb - nc * x * (n_over_bb + b / bb)
        n_temp = bb * (1.0 + sqrt(1.0 - (4.0 * aa / bb)))
        n = floor(2.0 * aa * (bb / n_temp))
    else
        aa = a - nc * x * (a + b)
        if bb < 0.0
            n = bb - sqrt(bb^2 - 4.0 * aa)
            n = floor(2.0 * aa / n)
        else
            n = floor((bb + sqrt(bb^2 - 4.0 * aa)) / 2.0)
        end
    end
    if n < 0.0
        n = 0.0
    end

    aa = n + a
    bb_idx = n
    ptnc = _poisson_term(n, nc, nc - n, 0.0)
    ptx = b / (aa + b) * _binomial_term(aa, b, x, y, ab_minus_cd(b, x, aa, y), 0.0)

    # Downward sum
    ps = 1.0
    nc_dtemp = 0.0
    p = 1.0
    s = 1.0
    w = p
    term = 1.0
    result = 0.0

    if ptx > 0
        while bb_idx > 0.0 && (((term > 1.0e-15 * result) && (p > 1.0e-16 * w)) || (ps > 1.0e-16 * nc_dtemp))
            s = aa / x * s
            ps = p * s
            nc_dtemp = nc_dtemp + ps
            p = bb_idx / nc * p
            aa = aa - 1.0
            bb_idx = bb_idx - 1.0
            if bb_idx == 0.0
                aa = a
            end
            s = s / (aa + b)
            term = s * w
            result = result + term
            w = w + p
        end
        w = w * ptnc
    else
        w = poiscdf(nc, n)
    end

    if n > 0.0
        nc_dtemp = (nc_dtemp + p * aa / x * s) * ptx
    elseif a == 0.0 || b == 0.0
        nc_dtemp = 0.0
    else
        nc_dtemp = _binomial_term(aa, b, x, y, ab_minus_cd(b, x, aa, y), log(b) + log(aa / (x * (aa + b))))
    end

    if x > y
        s = betacdf(b, aa, y)
    else
        s = betaccdf(aa, b, x)
    end
    ccdf_result = result * ptx * ptnc + s * w

    # Upward sum
    aa = n + a
    bb_idx = n
    p = 1.0
    nc_derivative = 0.0
    s = ptx
    if x > y
        w = betacdf(b, aa, y) # 1 - I(x, aa, b)
    else
        w = betaccdf(aa, b, x)
    end
    term = 0.0
    result = term

    while true
        w = w + s # 1 - I(x, aa, b)
        s = (aa + b) * s
        aa = aa + 1.0
        bb_idx = bb_idx + 1.0
        p = nc / bb_idx * p
        ps = p * s
        nc_derivative = nc_derivative + ps
        s = x / aa * s
        term = p * w
        result = result + term
        if !(((term > 1.0e-15 * result) && (s > 1.0e-16 * w)) || (ps > 1.0e-16 * nc_derivative))
            break
        end
    end

    nc_dtemp = (nc_derivative + nc_dtemp) / y
    ccdf_result = ccdf_result + result * ptnc + poisccdf(nc, bb_idx) * w

    if nc_dtemp == 0.0
        nc_derivative = 0.0
    else
        nc_derivative = _poisson_term(n, nc, nc - n, log(nc_dtemp))
    end

    return ccdf_result, nc_derivative
end

"""
    _inv_beta_nc1(prob, a, b, nc) -> x

Inverse CDF for the noncentral beta distribution via Newton-Raphson.
Based on VBA `inv_beta_nc1` by Ian Smith.
"""
function _inv_beta_nc1(prob::Float64, a::Float64, b::Float64, nc::Float64)
    if prob > 0.5
        return _comp_inv_beta_nc1(1.0 - prob, a, b, nc)
    end

    lop = 0.0
    hip = 1.0
    lox = 0.0
    hix = 1.0
    pr_exp = exp(-nc)

    if pr_exp > prob
        if 2.0 * prob > pr_exp
            x = betainvccdf(a + _cSmall, b, (pr_exp - prob) / pr_exp)
        else
            x = betainvcdf(a + _cSmall, b, prob / pr_exp)
        end
        if x == 0.0
            return 0.0
        else
            oneMinusP = 1.0 - x
            temp = oneMinusP
            y = betainvcdf(a + nc^2 / (a + 2 * nc), b, prob)
            oneMinusP2 = (a + nc) * (1.0 - y) / (a + nc * (1.0 + y))
            if temp > oneMinusP2
                oneMinusP = temp
            else
                x = (a + 2.0 * nc) * y / (a + nc * (1.0 + y))
                oneMinusP = oneMinusP2
            end
        end
    else
        y = betainvcdf(a + nc^2 / (a + 2 * nc), b, prob)
        x = (a + 2.0 * nc) * y / (a + nc * (1.0 + y))
        oneMinusP = (a + nc) * (1.0 - y) / (a + nc * (1.0 + y))
        if oneMinusP < _cSmall
            oneMinusP = _cSmall
            pr, nc_derivative = _beta_nc1(x, oneMinusP, a, b, nc)
            if pr < prob
                return 1.0
            end
        end
    end

    dif = 0.0
    while true
        pr, nc_derivative = _beta_nc1(x, oneMinusP, a, b, nc)
        if pr < 3.0e-308 && nc_derivative == 0.0
            hip = oneMinusP
            lox = x
            dif = dif / 2.0
            x = x - dif
            oneMinusP = oneMinusP + dif
        elseif nc_derivative == 0.0
            lop = oneMinusP
            hix = x
            dif = dif / 2.0
            x = x - dif
            oneMinusP = oneMinusP + dif
        else
            if pr < prob
                hip = oneMinusP
                lox = x
            else
                lop = oneMinusP
                hix = x
            end
            dif = -(pr / nc_derivative) * _logdif(pr, prob)
            if x > oneMinusP
                if oneMinusP - dif < lop
                    dif = (oneMinusP - lop) * 0.9
                elseif oneMinusP - dif > hip
                    dif = (oneMinusP - hip) * 0.9
                end
            elseif x + dif < lox
                dif = (lox - x) * 0.9
            elseif x + dif > hix
                dif = (hix - x) * 0.9
            end
            x = x + dif
            oneMinusP = oneMinusP - dif
        end
        if !((abs(pr - prob) > prob * 1.0e-14) && (abs(dif) > abs(min(x, oneMinusP)) * 1.0e-10))
            break
        end
    end
    return x
end

"""
    _comp_inv_beta_nc1(prob, a, b, nc) -> x

Inverse complementary CDF for the noncentral beta distribution via Newton-Raphson.
Based on VBA `comp_inv_beta_nc1` by Ian Smith.
"""
function _comp_inv_beta_nc1(prob::Float64, a::Float64, b::Float64, nc::Float64)
    if prob > 0.5
        return _inv_beta_nc1(1.0 - prob, a, b, nc)
    end

    lop = 0.0
    hip = 1.0
    lox = 0.0
    hix = 1.0
    pr_exp = exp(-nc)

    if pr_exp > prob
        if 2.0 * prob > pr_exp
            x = betainvcdf(a + _cSmall, b, (pr_exp - prob) / pr_exp)
        else
            x = betainvccdf(a + _cSmall, b, prob / pr_exp)
        end
        oneMinusP = 1.0 - x
        if oneMinusP < _cSmall
            oneMinusP = _cSmall
            pr, nc_derivative = _comp_beta_nc1(x, oneMinusP, a, b, nc)
            if pr > prob
                return 1.0
            end
        else
            temp = oneMinusP
            y = betainvccdf(a + nc^2 / (a + 2 * nc), b, prob)
            oneMinusP2 = (a + nc) * (1.0 - y) / (a + nc * (1.0 + y))
            if temp < oneMinusP2
                oneMinusP = temp
            else
                x = (a + 2.0 * nc) * y / (a + nc * (1.0 + y))
                oneMinusP = oneMinusP2
            end
            if oneMinusP < _cSmall
                oneMinusP = _cSmall
                pr, nc_derivative = _comp_beta_nc1(x, oneMinusP, a, b, nc)
                if pr > prob
                    return 1.0
                end
            elseif x < _cSmall
                x = _cSmall
                pr, nc_derivative = _comp_beta_nc1(x, oneMinusP, a, b, nc)
                if pr < prob
                    return 0.0
                end
            end
        end
    else
        y = betainvccdf(a + nc^2 / (a + 2 * nc), b, prob)
        x = (a + 2.0 * nc) * y / (a + nc * (1.0 + y))
        oneMinusP = (a + nc) * (1.0 - y) / (a + nc * (1.0 + y))
        if oneMinusP < _cSmall
            oneMinusP = _cSmall
            pr, nc_derivative = _comp_beta_nc1(x, oneMinusP, a, b, nc)
            if pr > prob
                return 1.0
            end
        elseif x < _cSmall
            x = _cSmall
            pr, nc_derivative = _comp_beta_nc1(x, oneMinusP, a, b, nc)
            if pr < prob
                return 0.0
            end
        end
    end

    dif = x
    while true
        pr, nc_derivative = _comp_beta_nc1(x, oneMinusP, a, b, nc)
        if pr < 3.0e-308 && nc_derivative == 0.0
            lop = oneMinusP
            hix = x
            dif = dif / 2.0
            x = x - dif
            oneMinusP = oneMinusP + dif
        elseif nc_derivative == 0.0
            hip = oneMinusP
            lox = x
            dif = dif / 2.0
            x = x - dif
            oneMinusP = oneMinusP + dif
        else
            if pr < prob
                lop = oneMinusP
                hix = x
            else
                hip = oneMinusP
                lox = x
            end
            dif = (pr / nc_derivative) * _logdif(pr, prob)
            if x > oneMinusP
                if oneMinusP - dif < lop
                    dif = (oneMinusP - lop) * 0.9
                elseif oneMinusP - dif > hip
                    dif = (oneMinusP - hip) * 0.9
                end
            elseif x + dif < lox
                dif = (lox - x) * 0.9
            elseif x + dif > hix
                dif = (hix - x) * 0.9
            end
            x = x + dif
            oneMinusP = oneMinusP - dif
        end
        if !((abs(pr - prob) > prob * 1.0e-14) && (abs(dif) > abs(min(x, oneMinusP)) * 1.0e-10))
            break
        end
    end
    return x
end

# Public API functions
#
# NOTE: The R/StatsFuns convention uses noncentrality parameter λ, while the VBA
# internal functions use nc = λ/2 (half the noncentrality parameter). All public
# functions convert between these conventions.

function nbetapdf(α::Real, β::Real, λ::Real, x::Real)
    T = float(Base.promote_typeof(α, β, λ, x))
    a = Float64(α)
    b = Float64(β)
    nc = Float64(λ) / 2.0  # VBA convention: nc = λ/2
    xf = Float64(x)

    if a < 0.0 || b < 0.0 || nc < 0.0 || nc > _nc_limit || (a == 0.0 && b == 0.0)
        return convert(T, NaN)
    elseif xf < 0.0 || xf > 1.0
        return convert(T, 0.0)
    elseif xf == 0.0 || nc == 0.0
        return convert(T, exp(-nc) * betapdf(a, b, xf))
    elseif xf == 1.0 && b == 1.0
        return convert(T, a + nc)
    elseif xf == 1.0
        return convert(T, betapdf(a, b, xf))
    else
        if a < 1.0 || xf * b <= (1.0 - xf) * (a + nc)
            _, nc_derivative = _beta_nc1(xf, 1.0 - xf, a, b, nc)
        else
            _, nc_derivative = _comp_beta_nc1(xf, 1.0 - xf, a, b, nc)
        end
        return convert(T, nc_derivative + 0.0) # +0.0 to avoid -0.0
    end
end

function nbetalogpdf(α::Real, β::Real, λ::Real, x::Real)
    return log(nbetapdf(α, β, λ, x))
end

function nbetacdf(α::Real, β::Real, λ::Real, x::Real)
    T = float(Base.promote_typeof(α, β, λ, x))
    a = Float64(α)
    b = Float64(β)
    nc = Float64(λ) / 2.0
    xf = Float64(x)

    if a < 0.0 || b < 0.0 || nc < 0.0 || nc > _nc_limit || (a == 0.0 && b == 0.0)
        return convert(T, NaN)
    elseif xf < 0.0
        return convert(T, 0.0)
    elseif xf >= 1.0
        return convert(T, 1.0)
    elseif xf == 0.0 && a == 0.0
        return convert(T, exp(-nc))
    elseif xf == 0.0
        return convert(T, 0.0)
    elseif nc == 0.0
        return convert(T, betacdf(a, b, xf))
    elseif a < 1.0 || xf * b <= (1.0 - xf) * (a + nc)
        cdf_val, _ = _beta_nc1(xf, 1.0 - xf, a, b, nc)
        return convert(T, cdf_val + 0.0)
    else
        ccdf_val, _ = _comp_beta_nc1(xf, 1.0 - xf, a, b, nc)
        return convert(T, (1.0 - ccdf_val) + 0.0)
    end
end

function nbetaccdf(α::Real, β::Real, λ::Real, x::Real)
    T = float(Base.promote_typeof(α, β, λ, x))
    a = Float64(α)
    b = Float64(β)
    nc = Float64(λ) / 2.0
    xf = Float64(x)

    if a < 0.0 || b < 0.0 || nc < 0.0 || nc > _nc_limit || (a == 0.0 && b == 0.0)
        return convert(T, NaN)
    elseif xf < 0.0
        return convert(T, 1.0)
    elseif xf >= 1.0
        return convert(T, 0.0)
    elseif xf == 0.0 && a == 0.0
        return convert(T, -expm1(-nc))
    elseif xf == 0.0
        return convert(T, 1.0)
    elseif nc == 0.0
        return convert(T, betaccdf(a, b, xf))
    elseif a < 1.0 || xf * b >= (1.0 - xf) * (a + nc)
        ccdf_val, _ = _comp_beta_nc1(xf, 1.0 - xf, a, b, nc)
        return convert(T, ccdf_val + 0.0)
    else
        cdf_val, _ = _beta_nc1(xf, 1.0 - xf, a, b, nc)
        return convert(T, (1.0 - cdf_val) + 0.0)
    end
end

function nbetalogcdf(α::Real, β::Real, λ::Real, x::Real)
    return log(nbetacdf(α, β, λ, x))
end

function nbetalogccdf(α::Real, β::Real, λ::Real, x::Real)
    return log(nbetaccdf(α, β, λ, x))
end

function nbetainvcdf(α::Real, β::Real, λ::Real, q::Real)
    T = float(Base.promote_typeof(α, β, λ, q))
    a = Float64(α)
    b = Float64(β)
    nc = Float64(λ) / 2.0
    prob = Float64(q)

    if a < 0.0 || b <= 0.0 || nc < 0.0 || nc > _nc_limit || prob < 0.0 || prob > 1.0
        return convert(T, NaN)
    elseif prob == 0.0 || (a == 0.0 && prob <= exp(-nc))
        return convert(T, 0.0)
    elseif prob == 1.0
        return convert(T, 1.0)
    elseif nc == 0.0
        return convert(T, betainvcdf(a, b, prob))
    else
        return convert(T, _inv_beta_nc1(prob, a, b, nc) + 0.0)
    end
end

function nbetainvccdf(α::Real, β::Real, λ::Real, q::Real)
    T = float(Base.promote_typeof(α, β, λ, q))
    a = Float64(α)
    b = Float64(β)
    nc = Float64(λ) / 2.0
    prob = Float64(q)

    if a < 0.0 || b <= 0.0 || nc < 0.0 || nc > _nc_limit || prob < 0.0 || prob > 1.0
        return convert(T, NaN)
    elseif prob == 1.0 || (a == 0.0 && prob >= -expm1(-nc))
        return convert(T, 0.0)
    elseif prob == 0.0
        return convert(T, 1.0)
    elseif nc == 0.0
        return convert(T, betainvccdf(a, b, prob))
    else
        return convert(T, _comp_inv_beta_nc1(prob, a, b, nc) + 0.0)
    end
end

function nbetainvlogcdf(α::Real, β::Real, λ::Real, lq::Real)
    T = float(Base.promote_typeof(α, β, λ, lq))
    return convert(T, nbetainvcdf(α, β, λ, exp(lq)))
end

function nbetainvlogccdf(α::Real, β::Real, λ::Real, lq::Real)
    T = float(Base.promote_typeof(α, β, λ, lq))
    return convert(T, nbetainvccdf(α, β, λ, exp(lq)))
end
