# Routines for Laguerre polynomials

# Compute the expansion of the orthonormal polynomial without e^(qm*x^m/2) nor a constant factor based on some heuristics.
function polyAsyRHgen(np, y, alpha, T::Int64, qm, m::Int64, UQ)
    if (qm == 1) && (m == 1)
        z = y/4/np
        mnxi = 2*np*( sqrt(z)*sqrt(1-z) - acos(sqrt(z) ) ) # = -n*xin/i
    else
        A = zeros(m+1)
        for k =0:m
            A[k+1] = prod((2*(1:k)-1)/2/(1:k))
        end
        z = y/(np*2/m/qm/A[m+1] )^(1/m)
        # Also correct but much slower: Hn = 4*m/(2*m-1)*double(hypergeom([1, 1-m], 3/2-m, z))/m
        Hn = 2/A[m+1]*sum(z.^(0:m-1).*A[m-(0:m-1)])/m
        mnxi = np*(sqrt(z+0im).*sqrt(1-z+0im).*Hn/2 -2*acos(sqrt(z+0im)))
    end
    # We could avoid these tests by splitting the loop k=1:mn into three parts with heuristics for the bounding indices.
    if y < sqrt(np)
        # The fixed delta in the Riemann-Hilbert Problem would mean this bound has to be proportional to n, but x(1:k) are O(1/n) so choose the bound in between them to make more use of the (cheap) expansion in the bulk.
        return asyBesselgen(np, z, alpha, T, qm, m, UQ, mnxi + pi*np, true)
    elseif y > 3.7*np
        # Use the expansion in terms of the (expensive) Airy function, although the corresponding weights will start underflowing for n >= 186 for standard associated Laguerre polynomials.
        return asyAirygen(np, z, alpha, T, qm, m, UQ, (mnxi*3im/2)^(2/3), true)
    end
    asyBulkgen(np, z, alpha, T, qm, m, UQ, mnxi)
end

function asyBulkgen(np, z, alpha, T::Int64, qm, m::Int64, UQ, mnxi)
    if T == 1
        return real( 2/(z+0im)^(1/4 + alpha/2)/(1 - z+0im)^(1/4)*cos(acos(2*z - 1+0im)*(1/2 + alpha/2) - mnxi - pi/4) )
    end
    R = [1 0]
    for k = 1:T-1
        for i = 1:ceil(Int64, 3*T/2)
            R = R + reshape((UQ[1,:,k,i,1]/(z-1)^i + UQ[1,:,k,i,2]/z^i)/np^k,(1,2))
        end
    end
    p = real( 2/(z+0im)^(1/4 + alpha/2)*(cos(acos(2*z-1+0im)*(1/2+alpha/2) - mnxi-pi/4)*R[1]-cos(acos(2*z-1+0im)*(-1/2+alpha/2)-mnxi-pi/4)*R[2]*1im*4^alpha)/(1 - z+0im)^(1/4) )
end

function asyBesselgen(np, z, alpha, T::Int64, qm, m::Int64, UQ, npb, useQ::Bool)
    if T == 1
        return real( sqrt(2*pi)*(-1)^np*sqrt(npb)/(z+0im)^(1/4+alpha/2)/(1 - z+0im)^(1/4)*(sin( (alpha + 1)/2*acos(2*z - 1+0im) - pi*alpha/2)*besselj(alpha,npb) + cos( (alpha + 1)/2*acos(2*z - 1+0im) - pi*alpha/2)*(besselj(alpha-1,npb) - alpha/(npb)*besselj(alpha, npb) ) ) )
    end
    R = (1+0im)*[1 0]
    # Use the series expansion of R because it is faster and we use asyBessel only very close to zero to have less calls to besselj.
    if useQ
        for k = 1:T-1
            for i = 1:min(size(UQ,4),9-k)
                R = R + reshape(UQ[1, :, k, i, 4]*z^(i-1)/np^k,(1,2))
            end
        end
    else
        d = z-1+0im
        phi = 2*z-1+2*sqrt(z)*sqrt(d)
        Rko = (1+0im)zeros(T-1,2)
        sL = (1+0im)zeros(2,2,T-1)
		for m = 1:T-1
            for i = 1:ceil(Int64,3*m/2)
                Rko[m,:] += UQ[1, :, m, i, 1]/d^i+UQ[1, :, m, i, 2]/z^i
            end
            sL[:,:,m] = brac(m-1,alpha)/2^(1+2*m)/(npb/2im/np)^m*( [2^(-alpha)  0 ; 0   2^(alpha)]*[sqrt(phi)    1im/sqrt(phi)  ;  -1im/sqrt(phi)   sqrt(phi)]/2/z^(1/2)/d^(1/2)*[(-phi)^(alpha/2)  0 ; 0   (-phi)^(-alpha/2) ]*[((-1)^m)/m*(alpha^2+m/2-1/4)     (m-1/2)*1im  ;  -((-1)^m)*(m-1/2)*1im    (alpha^2+m/2-1/4)/m]*[(-phi)^(-alpha/2)  0 ; 0   (-phi)^(alpha/2) ]*[sqrt(phi)     -1im/sqrt(phi) ;  1im/sqrt(phi)    sqrt(phi)]*[2^(alpha)   0 ; 0 2^(-alpha)] -mod(m+1,2)*(4*alpha^2+2*m-1)/m*eye(2,2) )
        end
        for k = 1:T-1
            R += reshape((Rko[k,:] - sL[1,:,k])/np^k,(1,2))
            for m = 1:k-1
                R -= reshape(reshape(Rko[k-m,:],(1,2))*sL[:,:,m]/np^k,(1,2))
            end
        end
    end
    p = real( sqrt(2*pi)*(-1)^np*sqrt(npb)/(z+0im)^(1/4+alpha/2)/(1 - z+0im)^(1/4)*( (sin( (alpha + 1)/2*acos(2*z - 1+0im) - pi*alpha/2)*R[1] -sin( (alpha - 1)/2*acos(2*z - 1+0im) - pi*alpha/2)*R[2]*1im*4^alpha)*besselj(alpha, npb) + (cos( (alpha + 1)/2*acos(2*z - 1+0im)- pi*alpha/2)*R[1] - cos( (alpha - 1)/2*acos(2*z - 1+0im) - pi*alpha/2)*R[2]*1im*4^alpha)*(besselj(alpha-1, npb) - alpha/npb*besselj(alpha, npb) ) ) )
end

function asyAirygen(np, z, alpha, T::Int64, qm, m::Int64, UQ, fn, useQ::Bool, xin=NaN+NaN*1im)
    d = z - 1.0 +0im
    if T == 1
        return real( 4*sqrt(pi)/(z+0im)^(1/4+alpha/2)/d^(1/4)*(cos( (alpha + 1)/2*acos(2*z - 1+0im) )*fn^(1/4)*airyai(fn) -1im*sin( (alpha + 1)/2*acos(2*z - 1+0im) )*ifelse(angle(z-1) <= 0, -one(z), one(z) )*fn^(-1/4)*airyaiprime(fn) ) )
    end
    R = (1+0im)*[1 0]
    if useQ
        for k = 1:T-1
            for i = 1:min(size(UQ,4),9-k)
                R = R + reshape(UQ[1, :, k, i, 3]*d^(i-1)/np^k,(1,2))
            end
        end
    else
        phi = 2*z-1+2*sqrt(z)*sqrt(d)
        Rko = (1+0im)zeros(T-1,2)
        sR = (1+0im)zeros(2,2,T-1)
	for m = 1:T-1
            for i = 1:ceil(Int64,3*m/2)
                Rko[m,:] += UQ[1, :, m, i, 1]/d^i+UQ[1, :, m, i, 2]/z^i
            end
            sR[:,:,m] = nuk(m)/xin^m*( [2^(-alpha)  0 ; 0   2^(alpha)]*[sqrt(phi)    1im/sqrt(phi)  ;  -1im/sqrt(phi)   sqrt(phi)]/8/z^(1/2)/d^(1/2)*[phi^(alpha/2)  0 ; 0   phi^(-alpha/2) ]*[(-1.0)^m  -m*6im ; 6im*m*(-1)^m    1.0]*[phi^(-alpha/2)  0 ; 0   phi^(alpha/2) ]*[sqrt(phi)     -1im/sqrt(phi) ;  1im/sqrt(phi)    sqrt(phi)]*[2^(alpha)   0 ; 0 2^(-alpha)] - mod(m+1,2)*eye(2,2) )
        end
        for k = 1:T-1
            R += reshape((Rko[k,:] -sR[1,:,k])/np^k,(1,2))
            for m = 1:k-1
                R -= reshape(reshape(Rko[k-m,:],(1,2))*sR[:,:,m]/np^k,(1,2))
            end
        end
    end
    p = real( 4*sqrt(pi)/(z+0im)^(1/4+alpha/2)/d^(1/4)*( (R[1]*cos( (alpha + 1)/2*acos(2*z - 1+0im) ) -cos( (alpha - 1)/2*acos(2*z - 1+0im) )*R[2]*1im*4^alpha)*fn^(1/4)*airyai(fn) + 1im*(-sin( (alpha + 1)/2*acos(2*z - 1+0im) )*R[1] +sin( (alpha - 1)/2*acos(2*z - 1+0im) )*R[2]*1im*4^alpha)*ifelse(angle(z-1) <= 0, -one(z), one(z) )*fn^(-1/4)*airyaiprime(fn) ) )
end

# Additional short functions
# pochhammer does not seem to exist yet in Julia
poch(x,n) = prod(x .+ (0:(n-1)))

function binom(x,n) # binomial only works for integer x
    b = one(x)
    for i = 1:n
        b *= (x-(n-i))/i
    end
    b
end

nuk(n) = -gamma(3*n-1/2)*2^n/27^n/2/n/sqrt(pi)/gamma(n*2)
brac(n,alpha) = prod(4*alpha^2 .- (2*(1:n) .- 1).^2 )/(2^(2*n)*gamma(1.0+n))


# Compute the W or V-matrices to construct the asymptotic expansion of R.
# Input
#   alpha, qm, m - Factors in the weight function w(x) = x^alpha*exp(-qm*x^m)
#   maxOrder     - The maximum order of the error
#   r            - 1 when computing Wright, -1 when computing Wleft
#   isW          - Whether we compute W(iso V)-matrices
# Output
#   WV           - Coefficient matrices for (z + 1/2 \pm 1/2)^m of Delta_k(z) or s_k(z)
function getV(alpha,qm,m::Int64,maxOrder::Int64,r)
    mo = ceil(Int64, 3*maxOrder/2) + 4
    ns = 0:mo
    f = NaN*zeros(mo+1) # Coefficients in the expansion of \bar{phi}_n(z) or \xi_n(z)
    g = NaN*zeros(maxOrder-1,mo+1)

    A = zeros(m+1)
    for k =0:m
        A[k+1] = prod((2*(1:k) .- 1.0)/2/(1:k))
    end
    if (r == 1) # Right disk: near z=1
        f = NaN*zeros(mo+2)
        ns = [ns; ns[mo+1]+1] # Extend by one because f(1) = 0 while not for left
        u = zeros(ns[mo+2]+1,ns[mo+2]+2)
        v = zeros(ns[mo+2]+1,ns[mo+2]+2)
        u[1,1] = 1.0
        v[1,1] = 1.0
        for n = [ns; ns[mo+2]+1]
            u[2,n+1] = binom(1/2,n+1)
            v[2,n+1] = binom(1/2,n+2)
        end
        for kt = 2:ns[mo+2]
            for n = ns
                u[kt+1,n+1] = sum(u[kt,(0:n) .+ 1].*u[2,n .- (0:n) .+ 1])
                v[kt+1,n+1] = sum(v[kt,(0:n) .+ 1].*v[2,n .- (0:n) .+ 1])
            end
        end
        q = zeros(ns[mo+2]+1)
        rr = zeros(ns[mo+2]+1) # Coeffs in the expansion of sqrt(2-2*sqrt(1-w))
        for kt = ns
            for l = 0:kt
                q[kt+1] = q[kt+1] + poch(1/2,kt-l)*u[kt-l+1,l+1]/(-2)^(kt-l)/gamma(1.0+kt-l)/(1+2*(kt-l) )
                rr[kt+1] = rr[kt+1] + binom(1/2,kt-l)*v[kt-l+1,l+1]*2^(kt-l)
            end
        end
        if (m == 1)
            for n = ns
                f[n+1,1] = -2*binom(1/2,n)
                for l = 0:n
                    f[n+1,1] = f[n+1,1] + 2*q[l+1]*rr[n-l+1]
                end
            end
        else
            for j = ns
                f[j+1,1] = 0
                for i=0:min(j,m-1)
                    f[j+1,1] = f[j+1,1] + binom(1/2,j-i)*(-1)^(m-i-1)*gamma(-1/2-i)/gamma(1/2-m)/gamma(m-i)
                end
                f[j+1,1] = -f[j+1,1]/m/A[m+1]
                for l = 0:j
                    f[j+1,1] = f[j+1,1] + 2*q[l+1]*rr[j-l+1]
                end
            end
        end
        if(abs(f[1]) > 10*eps(Float64) )
            error("xi_n should be O( (z-1)^(3/2) ): Expected f[1] to be zero")
        end
        ns = ns[1:mo+1] # Reset ns to its value before computing f's
        g[1,1,1] = -1/f[2,1]
        for n = 1:mo
            g[1,n+1,1] = -sum(reshape(g[1,1:n,1],(n,1)).*reshape(f[(n+2):-1:3,1],(n,1) ) )/f[2,1]
        end
    else # Left disk: near z=0
        if (m == 1)
            for n = ns
                f[n+1,1] = -(binom(1/2,n)*(-1)^n + poch(1/2,n)./(1+2*n)./gamma(1.0+n))
            end
        else
            for n = ns
                f[n+1,1] = 0.0
                for k = 0:min(m-1,n)
                    f[n+1,1] = f[n+1,1] + binom(1/2,n-k)*(-1)^(n-k)*A[m-k]
                end
                f[n+1,1] = -f[n+1,1]/2/m/A[m+1]-poch(1/2,n)./(1+2*n)./gamma(1.0+n)
            end
        end
        g[1,1,1] = 1/f[1,1]
        for n = 1:mo
            g[1,n+1,1] = -sum(reshape(g[1,1:n,1],(n,1)).*reshape(f[(n+1):-1:2,1],(n,1)) )/f[1,1]
        end
    end
    rho = (1+1im)*zeros(2*mo+3,mo+1)
    for n = ns
        rho[2,n+1] = poch(1/2,n)/gamma(1.0+n)/(1+2*n)*(-r)^n
    end
    rho[1,1] = 1
    for i = 2:(maxOrder-1)
        for n = ns
            g[i,n+1] = sum(g[i-1,1:(n+1) ].*g[1,(n+1):-1:1] )
        end
    end
    for i = 2:(mo*2+2)
        for n = ns
            rho[i+1,n+1] = sum(rho[i,1:(n+1) ].*rho[2,(n+1):-1:1] )
        end
    end
    OmOdd = (1+1im)*zeros(mo+1); OmEven = (1+1im)*zeros(mo+1)
    XiOdd = (1+1im)*zeros(mo+1); XiEven = (1+1im)*zeros(mo+1)
    ThOdd = (1+1im)*zeros(mo+1); ThEven = (1+1im)*zeros(mo+1)
    OmO = (1+1im)*zeros(mo+1); OmE = (1+1im)*zeros(mo+1)
    XiO = (1+1im)*zeros(mo+1); XiE = (1+1im)*zeros(mo+1)
    ThO = (1+1im)*zeros(mo+1); ThE = (1+1im)*zeros(mo+1)
    for n = ns
        js = 0:n
        for j = js
            OmOdd[n+1] = OmOdd[n+1] + (-1)^j/gamma(1.0+2.0*j)*(-2*alpha/sqrt(-r+0.0im))^(2*j)*rho[2*j+1,n-j+1]
            XiOdd[n+1] = XiOdd[n+1] + (-1)^j/gamma(1.0+2.0*j)*(-2*(alpha+1)/sqrt(-r+0.0im))^(2*j)*rho[2*j+1,n-j+1]
            ThOdd[n+1] = ThOdd[n+1] + (-1)^j/gamma(1.0+2.0*j)*(-2*(alpha-1)/sqrt(-r+0.0im))^(2*j)*rho[2*j+1,n-j+1]
            OmEven[n+1] = OmEven[n+1] + (-1)^j/gamma(1.0+2*j+1.0)*(-2*alpha/sqrt(-r+0.0im))^(2*j+1)*rho[2*j+2,n-j+1]
            XiEven[n+1] = XiEven[n+1] + (-1)^j/gamma(1.0+2*j+1.0)*(-2*(alpha+1)/sqrt(-r+0.0im))^(2*j+1)*rho[2*j+2,n-j+1]
            ThEven[n+1] = ThEven[n+1] + (-1)^j/gamma(1.0+2*j+1.0)*(-2*(alpha-1)/sqrt(-r+0.0im))^(2*j+1)*rho[2*j+2,n-j+1]
        end
        for j = js
            OmO[n+1] = OmO[n+1] + binom(-1/2,j)*(r)^j*OmOdd[n-j+1]
            XiO[n+1] = XiO[n+1] + binom(-1/2,j)*(r)^j*XiOdd[n-j+1]
            ThO[n+1] = ThO[n+1] + binom(-1/2,j)*(r)^j*ThOdd[n-j+1]
            OmE[n+1] = OmE[n+1] + binom(-1/2,j)*(r)^j*OmEven[n-j+1]
            XiE[n+1] = XiE[n+1] + binom(-1/2,j)*(r)^j*XiEven[n-j+1]
            ThE[n+1] = ThE[n+1] + binom(-1/2,j)*(r)^j*ThEven[n-j+1]
        end
    end
    Ts = zeros(ComplexF64, 2, 2, mo+1) # = G_{k,n}^{odd/even} depending on k, overwritten on each new k
    WV = zeros(ComplexF64, 2, 2, maxOrder-1, mo+1)
    for k = 1:(maxOrder-1)
        Ts[:,:,:] .= 0
        if r == 1
            if mod(k,2) == 1
                for n = 0:mo
                    Ts[:,:,n+1] = nuk(k)*[-2*(2*binom(-1/2,n-1)*(n>0)+binom(-1/2,n))     2im*4^(-alpha)*binom(-1/2,n)    ;    2im*4^(alpha)*binom(-1/2,n)    (2*(2*binom(-1/2,n-1)*(n>0) +binom(-1/2,n)))] -6*k*nuk(k)*[-2*OmO[n+1]   4^(-alpha)*2im*XiO[n+1]  ;   4^(alpha)*2im*ThO[n+1]    2*OmO[n+1]]
                    WV[:,:,k,n+1] = sum(repeat(reshape(g[k,1:(n+1)], (1,1,n+1) ), outer=(2,2,1)) .* Ts[:,:,(n+1):-1:1], dims=3)/8
                end
            else
                for n = 0:mo
                     Ts[:,:,n+1] = nuk(k)*4*(n==0)*I +6*k*nuk(k)*[-2im*OmE[n+1]    -2*4^(-alpha)*XiE[n+1]  ;   -2*4^alpha*ThE[n+1]   2im*OmE[n+1]]
                     WV[:,:,k,n+1] = sum(repeat(reshape(g[k,1:(n+1) ], (1,1,n+1) ), outer=[2,2,1]).*Ts[:,:,(n+1):-1:1], dims=3)/8
                end
            end
        else
            if mod(k,2) == 1
                for n = 0:mo
                    Ts[:,:,n+1] = -(alpha^2+k/2-1/4)/k*[-(-1)^n*(2*binom(-1/2,n-1)*(n>0)+binom(-1/2,n))*2    -1im*4^(-alpha)*2*(-1)^n*binom(-1/2,n)  ;  -1im*4^(alpha)*2*(-1)^n*binom(-1/2,n)     ( (-1)^n*(2*binom(-1/2,n-1)*(n>0) +binom(-1/2,n))*2)] - (k-1/2)*[2*OmO[n+1]   4^(-alpha)*2im*XiO[n+1]  ;   4^(alpha)*2im*ThO[n+1]   -2*OmO[n+1]] # binom(-1/2,-1) should be zero
                    WV[:,:,k,n+1] = -(-1)^(ceil(Int64, k/2)+1)*(1im*sqrt(2))^k*(-2+0im)^(-k/2)/4^(k+1)*brac(k-1,alpha)*sum(repeat(reshape(g[k,1:(n+1) ], (1,1,n+1) ), outer=[2,2,1]).*Ts[:,:,(n+1):-1:1], dims=3)
                end
            else
                for n = 0:mo
                    Ts[:,:,n+1] = (alpha^2+k/2-1/4)/k*4*(n==0)*I  -2*(k-1/2)*[ OmE[n+1]   4^(-alpha)*1im*XiE[n+1]  ;   4^alpha*1im*ThE[n+1]   -OmE[n+1] ]
                    WV[:,:,k,n+1] = -(-1)^(ceil(Int64, k/2)+1)*(1im*sqrt(2))^k*(-2)^(-k/2)/4^(k+1)*brac(k-1,alpha)*sum(repeat(reshape(g[k,1:(n+1) ], (1,1,n+1) ), outer=[2,2,1]).*Ts[:,:,(n+1):-1:1],dims=3)
                end
            end
        end
    end
    WV
end

# Get the U-matrices to construct the asymptotic expansion of R using the procedure with the convolutions with a specified method.
# Input
#   alpha, qm, m - Parts of the weight function
#   maxOrder     - The maximal order of the error
# Output
#   UQ           - Coefficient matrices of R_k(z) for (z-1)^(-m) [Uright], or z^(-m) [Uleft] of R_k^{right}(z) for (z-1)^n [Qright] and of R_k^{left}(z) for z^n [Qleft]
function getUQ(alpha, qm, m::Int64, maxOrder::Int64)
    Vr = getV(alpha, qm, m, maxOrder, 1)
    Vl = getV(alpha, qm, m, maxOrder, -1)
    UQ = (1+1im)*zeros(2,2,maxOrder-1,ceil(Int64,3*maxOrder/2)+2, 4)
    for kt = 0:(maxOrder-2)
        # Uright(:,:,(maxOrder-1)+1,:) will not be used later on because first term in expansions is without U's
        for mt = 0:(ceil(Int64,3*(kt+1)/2)-1)
            UQ[:,:,kt+1,mt+1,1] = Vr[:,:,kt+1,ceil(Int64,3*(kt+1)/2)-mt]
            for j = 0:(kt-1)
                for l = 0:(ceil(Int64,3*(j+1)/2)-mt-1)
                    UQ[:,:,kt+1,mt+1,1] = UQ[:,:,kt+1,mt+1,1] + UQ[:,:,kt-j,l+1,3]*Vr[:,:,j+1,ceil(Int64,3*(j+1)/2)-l-mt]
                end
            end
        end
        for mt = 0:(ceil(Int64,(kt+1)/2)-1)
            UQ[:,:,kt+1,mt+1,2] = Vl[:,:,kt+1,ceil(Int64,(kt+1)/2)-mt]
            for j= 0:(kt-1)
                for l = 0:(ceil(Int64,(j+1)/2)-mt-1)
                    UQ[:,:,kt+1,mt+1,2] = UQ[:,:,kt+1,mt+1,2] + UQ[:,:,kt-j,l+1,4]*Vl[:,:,j+1,ceil(Int64,(j+1)/2)-l-mt]
                end
            end
        end
        for n = 0:(ceil(Int64,3*(maxOrder-kt+1)/2)-1)
            UQ[:,:,kt+1,n+1,3] = -Vr[:,:,kt+1,ceil(Int64,3*(kt+1)/2)+1+n]
            UQ[:,:,kt+1,n+1,4] = -Vl[:,:,kt+1,ceil(Int64,(kt+1)/2)+1+n]
            for i = 0:(ceil(Int64,(kt+1)/2)-1)
                UQ[:,:,kt+1,n+1,3] = UQ[:,:,kt+1,n+1,3] + binom(-i-1,n)*UQ[:,:,kt+1,i+1,2]
            end
            for i = 0:(ceil(Int64,3*(kt+1)/2)-1)
                UQ[:,:,kt+1,n+1,4] = UQ[:,:,kt+1,n+1,4] + binom(-i-1,n)*(-1.0)^(-i-1-n)*UQ[:,:,kt+1,i+1,1]
            end
            for j = 0:(kt-1)
                for l = 0:(ceil(Int64, (j+1)/2)+n)
                    UQ[:,:,kt+1,n+1,4] = UQ[:,:,kt+1,n+1,4] -UQ[:,:,kt-j,l+1,4]*Vl[:,:,j+1,n-l+1+ceil(Int64,(j+1)/2) ]
                end
                for l = 0:(ceil(Int64, 3*(j+1)/2)+n)
                    UQ[:,:,kt+1,n+1,3] = UQ[:,:,kt+1,n+1,3] -UQ[:,:,kt-j,l+1,3]*Vr[:,:,j+1,n-l+1+ceil(Int64, 3*(j+1)/2) ]
                end
            end
        end
    end
    UQ
end
