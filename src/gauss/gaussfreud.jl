
# Estimate for the number of weights that don't underflow.
# This is only an estimate, not an upper bound, so code has to be able to cope
# with the estimate being off.
# The estimate here is similar to the Laguerre one, but more general with the
# additional argument m.
estimate_freud_reduced_n(n, alpha, m) = round(typeof(n), min(exp(exp(1/m)*1.05)*n^(1-1/2/m), n))


"""
Compute the n-point Gaussian quadrature rule for the Freud-type weight function
`x^α exp(-q_m x^m)`.
The computations are based on asymptotic expansions for the corresponding orthogonal
polynomials.
"""
function asy_gaussfreud(n::Integer, alpha = 0.0, m = 1, qm = 1.0; reduced = false)
    if alpha <= -1
        error("The parameter α <= -1 corresponds to a nonintegrable weight function")
    end
    if n < 0
        error("gaussfreud($n,$alpha,$m,$qm) not defined: n must be positive.")
    end
    if alpha^2/n > 1
        warn("A large alpha may lead to inaccurate results because the weight is low and R(z) is not close to identity.")
    end

    ELT = promote_type(typeof(float(alpha)), typeof(qm))

    T = ceil(Int, 34/log(n) ) # Heuristic for number of terms, should be scaled by the logarithm of eps(ELT) over the machine precision.
    UQ0 = getUQ(alpha, qm, m, T)
    n_alloc = reduced ? estimate_freud_reduced_n(n, alpha, m) : n

    n_pre = max(ceil(Int, 3.6*n^0.188), 7)
    # Heuristics to switch between Bessel, extrapolation and Airy initial guesses.
    igatt = ceil(Int, n_alloc + 1.31*n^0.4 - n)

    A = zeros(ELT, m+1)
    for k in 0:m
        A[k+1] = prod((2*(1:k) .- 1)/2 ./ (1:k))
    end
    softEdge = (n*2/m/qm/A[m+1] )^(1/m)
    # Use finite differences for derivative of polynomial when not x^alpha*exp(-x) and use other initial approximations
    useFinDiff = (m != 1) || (qm != 1.0)
    bes = besselroots(alpha, n_pre).^2 # [Tricomi 1947 pg. 296]
    w = zeros(ELT, n_alloc)

    # Find initial values for x
    if useFinDiff
        x = [bes*(2*m-1)^2/16/m^2/n^2*softEdge ; zeros(ELT, n_alloc-n_pre) ]
    else
        ak = [-13.69148903521072; -12.828776752865757; -11.93601556323626;    -11.00852430373326; -10.04017434155809; -9.02265085340981; -7.944133587120853;    -6.786708090071759; -5.520559828095551; -4.08794944413097; -2.338107410459767]
        t = 3*pi/2*( (igatt:-1:12).-0.25) # [DLMF (9.9.6)]
        ak = [-t.^(2/3).*(1 .+ 5/48 ./ t.^2 - 5/36 ./t.^4 .+ 77125/82944 ./t.^6 .- 10856875/6967296 ./t.^8); ak[max(1,12-igatt):11] ]
        nu = 4*n+2*alpha+2 # [Gatteshi 2002 (4.9)]
        air = (nu .+ ak*(4*nu)^(1/3) .+ ak.^2*(nu/16)^(-1/3)/5 .+ (11/35 .- alpha^2 .- 12/175*ak.^3)/nu + (16/1575*ak .+ 92/7875*ak.^4)*2^(2/3)*nu^(-5/3) .- (15152/3031875*ak.^5 .+ 1088/121275*ak.^2)*2^(1/3)*nu^(-7/3))
        x = [ bes/(4*n + 2*alpha+2) .* (1 .+ (bes .+ 2*(alpha^2 - 1) )/(4*n + 2*alpha+2)^2/3 ) ; zeros(ELT, n_alloc - n_pre -max(igatt,0) ) ; air]
    end

    if !useFinDiff
        UQ1 = getUQ(alpha+1, qm, m, T)
        factor0 = 1-4im*4^alpha*sum((UQ0[1,2,1:(T-1),1,2] + UQ0[1,2,1:(T-1),1,1])./n.^reshape(1:(T-1), (1,1,T-1)) )
        factor1 = 1-4im*4^(alpha+1)*sum((UQ1[1,2,1:(T-1),1,2] + UQ1[1,2,1:(T-1),1,1])./n.^reshape(1:(T-1), (1,1,T-1)) )
        factorx = real(sqrt(factor0/factor1 )/2/(1 - 1/n)^(1+alpha/2))
        factorw = real( -(1 - 1/(n + 1) )^(n + 1+ alpha/2)*(1 - 1/n)^(1 + alpha/2)*exp(1 + 2*log(2) )*4^(1+alpha)*pi*n^alpha*sqrt(factor0*factor1)*(1 + 1/n)^(alpha/2) )
    end

    noUnderflow = true
    for k = 1:n
        if useFinDiff && (k > n_pre)
            # Use linear extrapolation for the initial guesses for robustness in generalised weights.
            x[k] = 2*x[k-1] -x[k-2]
        elseif k > n_pre
            # Use sextic extrapolation for the initial guesses.
            x[k] = 7*x[k-1] -21*x[k-2] +35*x[k-3] -35*x[k-4] +21*x[k-5] -7*x[k-6] +x[k-7]
        end
        step = x[k]
        l = 0 # Newton-Raphson iteration number
        max_iter = 20
        ov = floatmax(ELT) # Previous/old value
        ox = x[k] # Old x
        # Accuracy of the expansions up to machine precision would lower this bound.
        while ( abs(step) > eps(ELT)*40*x[k] ) && ( l < max_iter)
            l = l + 1
            pe = polyAsyRHgen(n, x[k], alpha, T, qm, m, UQ0)
            if abs(pe) >= abs(ov)*(1-35*eps(ELT))
                # The function values do not decrease enough any more due to roundoff errors.
                x[k] = ox # Set to the previous value and quit.
                break
            end
            if useFinDiff
                hh = max(sqrt(eps(ELT))*x[k], sqrt(eps(ELT)) )
                step = pe*hh/(polyAsyRHgen(n, x[k]+hh, alpha, T, qm, m, UQ0) - pe)
            else
                # poly' = (p*exp(-Q/2) )' = exp(-Q/2)*(p' -p/2) with orthonormal p.
                step = pe/(polyAsyRHgen(n-1, x[k], alpha+1, T, qm, m, UQ1)*factorx - pe/2)
            end
            ox = x[k]
            x[k] = x[k] -step
            ov = pe
        end
        if (x[k] < 0) || (l == max_iter) || ( (k != 1) && (x[k-1] >= x[k]) ) || isnan(x[k])
            # Print some debugging information and throw an error.
            print(x[k], " = x[k], k = ", k, ", l = ", l, ", x[k-1] = ", x[k-1], ", x[k-2] = ", x[k-1], ", step = ", step, ", ox=", ox, ", ov = ", ov, ".\n")
            error("Newton method may not have converged.")
        elseif x[k] > softEdge
            warn("Node is outside the support of the measure: inaccuracy is expected.")
        end
        if noUnderflow
            if useFinDiff
                hh = max(sqrt(eps(ELT))*x[k], sqrt(eps(ELT)) )
                w[k] = hh/(polyAsyRHgen(n,x[k]+hh,alpha,T,qm,m,UQ0) -polyAsyRHgen(n,x[k],alpha,T,qm,m,UQ0))/polyAsyRHgen(n-1,x[k],alpha,T,qm,m,UQ0)/exp(qm*x[k]^m) # This leaves out a constant factor, given by a ratio of leading order coefficients and normalising constants
            else noUnderflow
                w[k] = factorw/polyAsyRHgen(n-1, x[k], alpha+1, T, qm, m, UQ1)/polyAsyRHgen(n+1, x[k], alpha, T,qm,m,UQ0)/exp( x[k] )
            end
        end
        if noUnderflow && (w[k] < underflow_threshold(ELT))
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
