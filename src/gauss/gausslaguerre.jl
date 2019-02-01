# Estimate for the number of weights that don't underflow
# This is only an estimate, not an upper bound, so code has to be able to cope
# with the estimate being off
estimate_reduced_n(n, alpha) = round(typeof(n), min(17*sqrt(n), n))

# Our threshold for deciding on underflow
underflow_threshold(x) = underflow_threshold(typeof(x))
underflow_threshold(::Type{T}) where {T <: AbstractFloat} = 10realmin(T)

# We explicitly store the first 11 roots of the Airy function in double precision
const airy_roots = [-2.338107410459767, -4.08794944413097, -5.520559828095551,
    -6.786708090071759, -7.944133587120853, -9.02265085340981, -10.04017434155809,
    -11.00852430373326, -11.93601556323626, -12.828776752865757, -13.69148903521072]


function ensure_allocation(n, k, x, w)
    n_alloc = length(x)
    if k > n_alloc
        n_alloc = min(n, 2*n_alloc)
        x1 = x
        w1 = w
        x = zeros(ELT, n_alloc)
        w = zeros(ELT, n_alloc)
        x[1:k] = x1
        w[1:k] = w1
    end
    x, w
end

"""
Compute the Gauss-Laguerre rule using explicit asymptotic expansions for the nodes
and weights.
Optional parameters are the number of terms in the expansion ('T'), the index where
one transitions from the Bessel region to the bulk (`k_besssel`) and from the bulk
to the Airy region ('k_airy').
"""
function asy_gausslaguerre(n::Integer, alpha, reduced = false,
    T = max(1, ceil(Int, 50/log(n))),           # Heuristic for number of terms
    k_bessel = max(ceil(Int, sqrt(n) ), 7),     # Heuristical indices for Bessel and Airy regions
    k_airy = floor(Int, 0.9*n))

    if alpha^2/n > 1
        warn("A large value of alpha may lead to inaccurate results.")
    end

    n_alloc = reduced ? estimate_reduced_n(n, alpha) : n

    ELT = typeof(float(alpha))
    x = zeros(ELT, n_alloc)
    w = zeros(ELT, n_alloc)

    # The expansions are given in powers of 1/(4n+2α+2)
    d = 1/(4n+2alpha+2)


    # The Bessel region
    # - Compute the roots of the Bessel function of order alpha
    jak_vector = besselroots(alpha, k_bessel)
    # - Iterate over all elements
    for k in 1:k_bessel
        jak = jak_vector[k]
        xk, wk, δ_bessel = asy_gausslaguerre_bessel(n, alpha, jak_vector[k], T, d)

        # Are we producing a compressed representation?
        if reduced
            # We check whether or not the current weight underflows
            if abs(wk) < underflow_threshold(ELT)
                # It does, we can stop here
                x = x[1:k-1]
                w = w[1:k-1]
                return x, w
            end
            # We check whether or not our allocated vector is still large enough
            x, w = ensure_allocation(n, k, x, w)
        end
        x[k] = xk
        w[k] = wk
    end

    # The bulk region
    for k in k_bessel+1:k_airy-1
        xk, wk, δ_bulk = asy_gauslaguerre_bulk(n, alpha, k, d, T)
        if reduced
            if abs(wk) < underflow_threshold(ELT)
                x = x[1:k-1]
                w = w[1:k-1]
                return x, w
            end
            x, w = ensure_allocation(n, k, x, w)
        end
        x[k] = xk
        w[k] = wk
    end

    # The Airy region
    for k in k_airy:n
        xk, wk, δ_airy = asy_gausslaguerre_airy(n, alpha, k, d, T)
        if reduced
            if abs(wk) < underflow_threshold(ELT)
                x = x[1:k-1]
                w = w[1:k-1]
                return x, w
            end
            x, w = ensure_allocation(n, k, x, w)
        end
        x[k] = xk
        w[k] = wk
    end

    # Sanity check
    if ( minimum(x) < 0.0 ) || ( maximum(x) > 4*n + 2*alpha + 2 ) ||  ( minimum(diff(x)) <= 0.0 ) || (minimum(w) < 0.0)
        warn("Unexpected inconsistency in the computation of nodes and weights")
    end

    x, w
end

function asy_gauslaguerre_bulk(n, alpha, k, d, T)
    ELT = eltype(alpha)

    pt = (4n-4k+3)*d
    t = ELT(pi)^2/16*(pt -1)^2
    for it = 1:6
        t = t - (pt*pi +2*sqrt(t-t^2) -acos(2*t-1) )*sqrt(t/(1-t))/2
    end

    x3 = -d/12*(12*alpha^2 + 5*(1-t)^(-2) - 4*(1-t)^(-1) - 4)
    x5 = d^3*(1-t)/t/720*(1600*(1-t)^(-6) - 3815*(1-t)^(-5) + 480*alpha^4 +2814*(1-t)^(-4) - 576*(1-t)^(-3)   - 960*alpha^2 - 48*(15*alpha^4 - 30*alpha^2 + 7)*(1-t)^(-1) -16*(1-t)^(-2) + 224)
    x7 = -d^5/181440*(1-t)^2/t^2*(10797500*(1-t)^(-10) - 43122800*(1-t)^(-9)
        + 66424575*(1-t)^(-8) -48469876*(1-t)^(-7) + 193536*alpha^6 + 16131880*(1-t)^(-6)
        + 80*(315*alpha^4 - 630*alpha^2 -221)*(1-t)^(-4) - 1727136*(1-t)^(-5)
        - 967680*alpha^4 - 320*(63*alpha^4 - 126*alpha^2 +43)*(1-t)^(-3)
        + 384*(945*alpha^6 - 4620*alpha^4 + 6405*alpha^2 - 1346)*(1-t)^(-2)
        + 1354752*alpha^2   - 23040*(21*alpha^6 - 105*alpha^4 + 147*alpha^2 - 31)*(1-t)^(-1) -285696)
    x9 = d^7/10886400*(1-t)^3/t^3*(43222750000*(1-t)^(-14) - 241928673000*(1-t)^(-13)
        + 566519158800*(1-t)^(-12)   -714465642135*(1-t)^(-11) + 518401904799*(1-t)^(-10)
        + 672*(12000*alpha^4 - 24000*alpha^2 +64957561)*(1-t)^(-8) - 212307298152*(1-t)^(-9)
        + 24883200*alpha^8 - 192*(103425*alpha^4 -206850*alpha^2 + 15948182)*(1-t)^(-7)
        + 3360*(4521*alpha^4 - 9042*alpha^2 - 7823)*(1-t)^(-6) -232243200*alpha^6
        - 1792*(3375*alpha^6 - 13905*alpha^4  + 17685*alpha^2 - 1598)*(1-t)^(-5)
        + 16128*(450*alpha^6 - 2155*alpha^4 + 2960*alpha^2 - 641)*(1-t)^(-4)
        + 812851200*alpha^4 -768*(70875*alpha^8 - 631260*alpha^6 + 2163630*alpha^4
        - 2716980*alpha^2 +555239)*(1-t)^(-3)  + 768*(143325*alpha^8 - 1324260*alpha^6
        + 4613070*alpha^4 -5826660*alpha^2 + 1193053)*(1-t)^(-2) - 1028505600*alpha^2
        - 5806080*(15*alpha^8 -140*alpha^6 + 490*alpha^4 - 620*alpha^2 + 127)*(1-t)^(-1) + 210677760)

    w3 = d^2/6*(2*t + 3)/(t-1)^3
    w5 = (1-t)^2/720/t^2*d^4*(8000*(1-t)^(-8) - 24860*(1-t)^(-7) + 27517*(1-t)^(-6)
        - 12408*(1-t)^(-5) + 1712*(1-t)^(-4) +16*(15*alpha^4 - 30*alpha^2 + 7)*(1-t)^(-2) + 32*(1-t)^(-3))
    w7 = -(1-t)^3/90720/t^3*d^6*(43190000*(1-t)^(-12) -204917300*(1-t)^(-11)
        + 393326325*(1-t)^(-10)  - 386872990*(1-t)^(-9) + 201908326*(1-t)^(-8)
        + 80*(315*alpha^4 - 630*alpha^2 + 53752)*(1-t)^(-6)  - 50986344*(1-t)^(-7)
        - 320*(189*alpha^4 -378*alpha^2 - 89)*(1-t)^(-5) + 480*(63*alpha^4 - 126*alpha^2
        + 43)*(1-t)^(-4)  -384*(315*alpha^6 - 1470*alpha^4 + 1995*alpha^2 - 416)*(1-t)^(-3)
        + 2304*(21*alpha^6 -105*alpha^4 + 147*alpha^2 - 31)*(1-t)^(-2) )

    xk = zero(ELT)
    wk = zero(ELT)
    if T >= 9
        xk = x9 + x7 + x5 + x3
        wk = w7 + w5 + w3
        δ = abs(x9)
    elseif T >= 7
        xk = x7 + x5 + x3
        wk = w7 + w5 + w3
        δ = abs(x9)
    elseif T >= 5
        xk = x5 + x3
        wk = w5 + w3
        δ = abs(x7)
    elseif T >= 3
        xk = x3
        wk = w3
        δ = abs(x5)
    else
        δ = abs(x3)
    end
    if T == -1
        if abs(x3) > eps(ELT) && abs(x3) < 1
            xk += x3
            wk += w3
            δ = abs(x5)
            if δ > eps(ELT) && δ < abs(x3)
                xk += x5
                wk += w5
                δ = abs(x7)
                if δ > eps(ELT) && δ < abs(x5)
                    xk += x7
                    wk += w7
                    if abs(x9) > eps(ELT) && abs(x9) < abs(x7)
                        xk += x9
                    end
                end
            end
        end
    end

    xk += + t/d
    wk = xk^alpha * exp(-xk) * 2pi * sqrt(t/(1-t)) * (1+wk)
    xk, wk, δ
end

function asy_gausslaguerre_bessel(n, alpha, jak, T, d)
    ELT = typeof(alpha)

    # First the formulas
    x3 = (jak^2 + 2*alpha^2 - 2)*d^2/3
    x5 = (11*jak^4 +3*jak^2*(11*alpha^2-19) +46*alpha^4 -140*alpha^2 +94)*d^4/45
    x7 = (657*jak^6 +36*jak^4*(73*alpha^2-181) +2*jak^2*(2459*alpha^4 -10750*alpha^2 +14051)
        + 4*(1493*alpha^6 -9303*alpha^4 +19887*alpha^2 - 12077) )*d^6/2835
    x9 = (10644*jak^8 + 60*(887*alpha^2 - 2879)*jak^6 + (125671*alpha^4 -729422*alpha^2 + 1456807)*jak^4
        + 3*(63299*alpha^6 - 507801*alpha^4 + 1678761*alpha^2 - 2201939)*jak^2 + 2*(107959*alpha^8
        - 1146220*alpha^6 + 5095482*alpha^4 -10087180*alpha^2 + 6029959) )*d^8/42525
    x11 = -4/945*(10644*(alpha^2 + 2*alpha + 1)*jak^8 + 60*(3953*alpha^4 + 14038*alpha^3 + 16404*alpha^2
        + 6506*alpha + 187)*jak^6 + (1637671*alpha^6 + 7851662*alpha^5 + 13629209*alpha^4 +
        9709796*alpha^3 + 2158745*alpha^2 + 272654*alpha + 408487)*jak^4 + (4802617*alpha^8 +
        27116434*alpha^7 + 57215614*alpha^6 + 51622634*alpha^5 +14740880*alpha^4 + 1802486*alpha^3
        + 7023106*alpha^2 + 1388686*alpha -1851977)*jak^2 + 2*(3595079*alpha^10 + 21764718*alpha^9
        + 48110139*alpha^8 +    41096360*alpha^7 + 1739502*alpha^6 - 4767756*alpha^5 + 14434702*alpha^4
        + 4980840*alpha^3 - 7744821*alpha^2 - 1142642*alpha + 1796919))

    w3 = (alpha^2 + jak^2 -1)*2/3*d^2
    w5 = (46*alpha^4 + 33*jak^4 +6*jak^2*(11*alpha^2 -19) -140*alpha^2 +94)/45*d^4
    w7 = (1493*alpha^6 + 657*jak^6 + 27*(73*alpha^2 - 181)*jak^4 - 9303*alpha^4
        + (2459*alpha^4 -10750*alpha^2 + 14051)*jak^2 + 19887*alpha^2 - 12077)*4/2835*d^6
    w9 = (215918*alpha^8 + 53220*jak^8 + 240*(887*alpha^2 - 2879)*jak^6 -2292440*alpha^6 +
        3*(125671*alpha^4 - 729422*alpha^2 + 1456807)*jak^4 + 10190964*alpha^4  +
        6*(63299*alpha^6 - 507801*alpha^4 + 1678761*alpha^2 -2201939)*jak^2 -
        20174360*alpha^2 + 12059918)/42525*d^8

    xk = zero(ELT)
    wk = zero(ELT)
    if T >= 11
        xk = x11 + x9 + x7 + x5 + x3
        wk = w9 + w7 + w5 + w3
        δ = abs(x11)
    elseif T >= 9
        xk = x9 + x7 + x5 + x3
        wk = w9 + w7 + w5 + w3
        δ = abs(x11)
    elseif T >= 7
        xk = x7 + x5 + x3
        wk = w7 + w5 + w3
        δ = abs(x9)
    elseif T >= 5
        xk = x5 + x3
        wk = w5 + w3
        δ = abs(x7)
    elseif T >= 3
        xk = x3
        wk = w3
        δ = abs(x5)
    else
        δ = abs(x3)
    end
    if T == -1
        if abs(x3) > eps(ELT) && abs(x3) < 1
            xk += x3
            wk += w3
            δ = abs(x5)
            if δ > eps(ELT) && δ < abs(x3)
                xk += x5
                wk += w5
                δ = abs(x7)
                if δ > eps(ELT) && δ < abs(x5)
                    xk += x7
                    wk += w7
                    δ = abs(x9)
                    if δ > eps(ELT) && δ < abs(x7)
                        xk += x9
                        wk += w9
                        if abs(x11) > eps(ELT) && abs(x11) < abs(x9)
                            xk += x11
                        end
                    end
                end
            end
        end
    end
    xk = jak^2*d*(1 + xk)
    wk = 4*d*xk^alpha*exp(-xk)/besselj(alpha-1, jak)^2*(1+wk)
    xk, wk, δ
end

function asy_gausslaguerre_airy(n, alpha, k, d, T)
    ELT = eltype(alpha)

    # Compute the corresponding root of the Airy function
    index = n - k + 1
    if index <= 11
        ak = airy_roots[index]
    else
        t = 3*pi/2*(index-0.25)
        ak = -t^(2/3)*(1 + 5/48/t^2 - 5/36/t^4 + 77125/82944/t^6 -10856875/6967296/t^8)
    end

    x5 = -(15152/3031875*ak^5+1088/121275*ak^2)*2^(1/3)*d^(7/3)
    x3 = ak^2*(d*16)^(1/3)/5 + (11/35-alpha^2-12/175*ak^3)*d + (16/1575*ak+92/7875*ak^4)*2^(2/3)*d^(5/3)
    x1 = 1/d + ak*(d/4)^(-1/3)

    xk = zero(ELT)
    wk = zero(ELT)
    if T >= 5
        xk = x5 + x3 + x1
        δ = abs(x5)
    elseif T >= 3
        xk = x3 + x1
        δ = abs(x5)
    elseif T >= 1
        xk = x1
        δ = abs(x3)
    end
    if T == -1
        xk = x1
        δ = abs(x3)
        if δ > eps(ELT) && δ < abs(x1)
            xk += x3
            δ = abs(x5)
            if δ > eps(ELT) && δ < abs(x3)
                xk += x5
            end
        end
    end
    wk = 4^(1/3)*xk^(alpha+1/3)*exp(-xk)/(airyaiprime(ak))^2
    xk, wk, δ
end



# Special case: alpha = 0
function asy_gausslaguerre0(n::Integer, reduced = false)

    ELT = Float64
    alpha = zero(ELT)

    T = max(1,ceil(Int, 50/log(n)))          # Heuristic for number of terms
    k_bessel = max(ceil(Int, sqrt(n) ), 7)   # Heuristical indices for Bessel and Airy regions
    k_airy = floor(Int, 0.9*n)

    n_alloc = reduced ? estimate_reduced_n(n, alpha) : n

    x = zeros(ELT, n_alloc)
    w = zeros(ELT, n_alloc)

    # The expansions are given in powers of 1/(4n+2)
    d = 1/(4n+2)


    # The Bessel region
    # - Compute the roots of the Bessel function of order alpha
    jak_vector = besselroots(zero(ELT), k_bessel)
    # - Iterate over all elements
    for k in 1:k_bessel
        xk = zero(ELT)
        wk = zero(ELT)
        jak = jak_vector[k]
        if (T >= 13)
            xk += -8/42525*(410649*jak^10 + 629298*jak^8 + 40737621*jak^6 -364320354*jak^4 + 1871633936*jak^2 - 3655024096)*d^12
        end
        if (T >= 11)
            xk += 1/1403325*( 410649*jak^10 -  9908262*jak^8 + 138902061*jak^6 - 1248722004*jak^4 + 6028914206*jak^2 - 11427291076)*d^10
            wk += 2/1403325*(1231947*jak^10 - 24770655*jak^8 + 277804122*jak^6 - 1873083006*jak^4 + 6028914206*jak^2 - 5713645538)*d^10
        end
        if (T >= 9)
            xk += (10644*jak^8 - 172740*jak^6 + 1456807*jak^4 -  6605817*jak^2 + 12059918)/42525*d^8
            wk += (53220*jak^8 - 690960*jak^6 + 4370421*jak^4 - 13211634*jak^2 + 12059918)/42525*d^8
        end
        if (T >= 7)
            xk += (657*jak^6 - 6516*jak^4 + 28102*jak^2 - 48308)*1/2835*d^6
            wk += (657*jak^6 - 4887*jak^4 + 14051*jak^2 - 12077)*4/2835*d^6
        end
        if (T >= 5)
            xk += (11*jak^4 - 57*jak^2 + 94)/45*d^4
            wk += (33*jak^4 -114*jak^2 + 94)/45*d^4
        end
        if (T >= 3)
            # From here, the terms are also in [Tricomi 1947 pg. 296]
            xk += (jak^2 - 2)*1/3*d^2
            wk += (jak^2 - 1)*2/3*d^2
        end
        xk = jak^2*d*(1 + xk)
        wk = 4*d*exp(-xk)/(besselj(-one(ELT), jak))^2*(1+wk)

        # Are we producing a compressed representation?
        if reduced
            # We check whether or not the current weight underflows
            if abs(wk) < underflow_threshold(ELT)
                # It does, we can stop here
                x = x[1:k-1]
                w = w[1:k-1]
                return x, w
            end
            # We check whether or not our allocated vector is still large enough
            if k > n_alloc
                n_alloc = min(n, 2*n_alloc)
                x1 = x
                w1 = w
                x = zeros(ELT, n_alloc)
                w = zeros(ELT, n_alloc)
                x[1:k] = x1
                w[1:k] = w1
            end
        end

        x[k] = xk
        w[k] = wk
    end

    # The bulk region
    for k in k_bessel+1:k_airy-1
        xk = zero(ELT)
        wk = zero(ELT)
        pt = (4n-4k+3)*d
        t = pi^2/16*(pt -1)^2
        for it = 1:6
            t = t - (pt*pi +2*sqrt(t-t^2) -acos(2*t-1) )*sqrt(t/(1-t))/2
        end
        if (T >= 9)
            xk += d^7/10886400*(1-t)^3/t^3*(43222750000*(1-t)^(-14) - 241928673000*(1-t)^(-13) + 566519158800*(1-t)^(-12)   -714465642135*(1-t)^(-11) + 518401904799*(1-t)^(-10) + 672*64957561*(1-t)^(-8)   - 212307298152*(1-t)^(-9) - 192*15948182*(1-t)^(-7)  - 3360*7823*(1-t)^(-6) + 1792*1598*(1-t)^(-5)+ 16128*(- 641)*(1-t)^(-4)  -768*555239*(1-t)^(-3)  + 768*1193053*(1-t)^(-2)- 5806080*127*(1-t)^(-1) + 210677760)
            # We do not have the corresponding term of the same order for the weights
        end
        if (T >= 7)
            xk += -d^5/181440*(1-t)^2/t^2*(10797500*(1-t)^(-10) - 43122800*(1-t)^(-9) + 66424575*(1-t)^(-8) -48469876*(1-t)^(-7) + 193536*alpha^6 + 16131880*(1-t)^(-6) + 80*(315*alpha^4 - 630*alpha^2 -221)*(1-t)^(-4) - 1727136*(1-t)^(-5) - 967680*alpha^4 - 320*(63*alpha^4 - 126*alpha^2 +43)*(1-t)^(-3)  + 384*(945*alpha^6 - 4620*alpha^4 + 6405*alpha^2 - 1346)*(1-t)^(-2) + 23040*31*(1-t)^(-1) -285696)
            wk += -(1-t)^3/90720/t^3*d^6*(43190000*(1-t)^(-12) -204917300*(1-t)^(-11) + 393326325*(1-t)^(-10)  - 386872990*(1-t)^(-9) + 201908326*(1-t)^(-8) +80*(315*alpha^4 - 630*alpha^2 + 53752)*(1-t)^(-6)  - 50986344*(1-t)^(-7) - 320*(189*alpha^4 -378*alpha^2 - 89)*(1-t)^(-5) + 480*(63*alpha^4 - 126*alpha^2 + 43)*(1-t)^(-4)  -384*(315*alpha^6 - 1470*alpha^4 + 1995*alpha^2 - 416)*(1-t)^(-3) + 2304*(21*alpha^6 -105*alpha^4 + 147*alpha^2 - 31)*(1-t)^(-2) )
        end
        if (T >= 5)
            # These are terms of size O(1/n^5)
            xk += d^3*(1-t)/t/720*(1600*(1-t)^(-6) - 3815*(1-t)^(-5) + 2814*(1-t)^(-4) - 576*(1-t)^(-3) - 48*7*(1-t)^(-1) -16*(1-t)^(-2) + 224)
            wk += (1-t)^2/720/t^2*d^4*(8000*(1-t)^(-8) - 24860*(1-t)^(-7) + 27517*(1-t)^(-6) - 12408*(1-t)^(-5) + 1712*(1-t)^(-4) +16*7*(1-t)^(-2) + 32*(1-t)^(-3))
        end
        if (T >= 3)
            xk += -d/12*(5*(1-t)^(-2) - 4*(1-t)^(-1) - 4)
            wk += d^2/6*(2*t + 3)/(t-1)^3
        end
        xk += + t/d
        wk = exp(-xk) * 2pi * sqrt(t/(1-t)) * (1+wk)

        if reduced
            if abs(wk) < underflow_threshold(ELT)
                x = x[1:k-1]
                w = w[1:k-1]
                return x, w
            end
            if k > n_alloc
                n_alloc = min(n, 2*n_alloc)
                x1 = x
                w1 = w
                x = zeros(ELT, n_alloc)
                w = zeros(ELT, n_alloc)
                x[1:k] = x1
                w[1:k] = w1
            end
        end

        x[k] = xk
        w[k] = wk
    end

    # The Airy region
    for k in k_airy:n
        xk = zero(ELT)
        wk = zero(ELT)

        # Compute the corresponding root of the Airy function
        index = n - k + 1
        if index <= 11
            ak = airy_roots[index]
        else
            t = 3*pi/2*(index-0.25)
            ak = -t^(2/3)*(1 + 5/48/t^2 - 5/36/t^4 + 77125/82944/t^6 -10856875/6967296/t^8)
        end

        if (T >= 5)
            # [Gatteshi 2002 (4.9)], gives O(n^{-4}) relative error
            xk += -(15152/3031875*ak^5+1088/121275*ak^2)*2^(1/3)*d^(7/3)
        end
        if (T >= 3)
            xk = xk + ak^2*(d*16)^(1/3)/5 + (11/35-alpha^2-12/175*ak^3)*d + (16/1575*ak+92/7875*ak^4)*2^(2/3)*d^(5/3)
        end
        xk += 1/d + ak*(d/4)^(-1/3)
        wk = 4^(1/3)*xk^(alpha+1/3)*exp(-xk)/(airyaiprime(ak))^2

        if reduced
            if abs(wk) < underflow_threshold(ELT)
                x = x[1:k-1]
                w = w[1:k-1]
                return x, w
            end
            if k > n_alloc
                n_alloc = min(n, 2*n_alloc)
                x1 = x
                w1 = w
                x = zeros(ELT, n_alloc)
                w = zeros(ELT, n_alloc)
                x[1:k] = x1
                w[1:k] = w1
            end
        end

        x[k] = xk
        w[k] = wk
    end

    # Sanity check
    if ( minimum(x) < 0.0 ) || ( maximum(x) > 4*n + 2 ) ||  ( minimum(diff(x)) <= 0.0 ) || (minimum(w) < 0.0)
        warn("Unexpected inconsistency in the computation of nodes and weights")
    end

    x, w
end




########################## Routines for the forward recurrence ##########################

"Compute Gauss-Laguerre rule based on the recurrence relation, using Newton iterations on an initial guess."
function rec_gausslaguerre(n, alpha, reduced = false)
    T = typeof(float(alpha))

    n_alloc = reduced ? estimate_reduced_n(n, alpha) : n
    w = zeros(T, n_alloc)
    x = zeros(T, n_alloc)

    # We compute up to 7 starting values for the Newton iterations
    n_pre = min(n, 7)

    nu = 4n + 2alpha + 2
    bes = T.(besselroots(alpha, n_pre)).^2 / nu # this is a lower bound by [DLMF 18.16.10]
    x[1:n_pre] = bes

    local pn_deriv
    noUnderflow = true      # this flag turns false once the weights start to underflow
    for k in 1:n
        if k > n_pre
            # Use sextic extrapolation for a new initial guess
            x[k] = 7*x[k-1] -21*x[k-2] +35*x[k-3] -35*x[k-4] +21*x[k-5] -7*x[k-6] +x[k-7]
        end

        step = x[k]
        ov = realmax(T) # Previous/old value
        ox = x[k] # Old x

        l = 0 # Newton-Raphson iteration number
        max_iter = 20
        while (abs(step) > 40eps(T)*x[k]) && (l < max_iter)
            l = l + 1
            pn, pn_deriv = evalLaguereRec(n, alpha, x[k])
            if abs(pn) >= abs(ov)*(1-50eps(T))
                # The function values do not decrease enough any more due to roundoff errors.
                x[k] = ox # Set to the previous value and quit.
                break
            end
            step = pn / pn_deriv
            ox = x[k]
            x[k] = x[k] - step
            ov = pn
        end
        if ( x[k] < 0 ) || ( x[k] > 4*n + 2*alpha + 2 ) || ( l == max_iter ) || ( ( k != 1 ) && ( x[k - 1] >= x[k] ) )
            warn("Newton method may not have converged in laguerreRec($n,$alpha).")
        end

        # Compute the weight
        if noUnderflow
            pn_prev, ~ = evalLaguereRec(n-1, alpha, x[k])
            w[k] = (n^2 +alpha*n)^(-1/2)/pn_prev/pn_deriv
        end
        if noUnderflow && (w[k] < underflow_threshold(T))
            # Frome here on after weights will no longer be computed
            noUnderflow = false
        end

        if reduced
            if (k > 1) && !noUnderflow
                x = x[1:k-1]
                w = w[1:k-1]
                return x, w
            end
            if k == n_alloc
                # We have to allocate a bigger array
                n_alloc *= 2
                x1 = x
                w1 = w
                x = zeros(T, n_alloc)
                w = zeros(T, n_alloc)
                x[1:k] = x1
                w[1:k] = w1
            end
        end
    end
    x, w
end


"""
Evaluate the orthonormal associated Laguerre polynomial with positive leading coefficient,
as well as its derivative, in the point x using the recurrence relation.
"""
function evalLaguereRec(n, alpha, x)
    T = promote_type(typeof(float(alpha)), typeof(float(x)))
    pnprev = zero(T)
    pn = 1/sqrt(gamma(alpha+1))
    pndprev = zero(T)
    pnd = zero(T)
    for k in 1:n
        pnold = pn
        pn = (x -2*k -alpha+1)/sqrt(k*(alpha+k))*pnold-sqrt((k-1+alpha)*(k-1)/k/(k+alpha))*pnprev
        pnprev = pnold
        pndold = pnd
        pnd = (pnold+(x-2*k-alpha+1)*pndold)/sqrt(k*(alpha+k)) -sqrt((k-1+alpha)*(k-1)/k/(alpha+k))*pndprev
        pndprev = pndold
    end
    pn, pnd
end
