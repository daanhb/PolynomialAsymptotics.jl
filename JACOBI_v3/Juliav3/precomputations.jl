# precomputations.jl: contains all functions doing the precomputations
# About
#   Author       - Peter Opsomer (Peter.Opsomer@cs.kuleuven.be)
#   History      - Created October 2014, last edit October 2015


# Compute the contour integral of the given function using the trapezoidal rule around the interval.
# Input
#   f   - the integrand
#   rho - parameter for the Bernstein ellipse
#   M   - number of terms, use successive doubling from [F. Bornemann 2011] if NaN]
# Output
#   cont- the contour integral
function trap_rule(f,rho,M=NaN)
if isnan(M)
	M = 8; tol = 1e-15;
	function thi(k,l) 2im*pi*(1:l:k)/k; end
	function bernst(k,l) rho/2*exp(thi(k,l)) + 1/rho/2*exp(-thi(k,l)); end
	s = f(bernst(M,1)).*(rho*1im/2*exp(thi(M,1)) -1im/rho/2*exp(-thi(M,1)) )/1im; cont = mean(s); err1 = NaN;
	while M < 1e6
		M = 2*M;
		s = [s; f(bernst(M,2)).*(rho*1im/2*exp(thi(M,2)) -1im/rho/2*exp(-thi(M,2)) )/1im];
		val = mean(s); kappa = mean(abs(s))/abs(val);
		err0 = abs(val-cont)/abs(val); 
		if (err0^3/err1^2 <= kappa*tol) || (~isfinite(kappa)) break; end
		cont = val; err1 = err0;
	end
elseif true
	ths = linspace(0,2*pi,M+1); # = theta_k, k = 0 to M-1
	ths = ths[1:(end-1)]
	hh = ths[2]-ths[1]; # To be sure to have the right distance
	zs = rho/2*exp(1im*ths) + 1/rho/2*exp(-1im*ths);
	cont = hh*sum(f(zs).*(rho*1im/2*exp(1im*ths) -1im/rho/2*exp(-1im*ths) ) )/(2*pi*1im );
else
	# Below an ellipse, but could also need other parametrisations, for example when evaluating psi(x) for x further from the interval than the point where log(h) is not analytic
    r1 = 1.8;
    r2 = 0.8;
    zs = r1*cos(ths) +r2*1im*sin(ths);
    cont = hh*sum(f(zs).*(-r1*sin(ths) +r2*1im*cos(ths) ) )/(2*pi*1im );
end
if rho < 1 # Bernstein in negative direction but need positive direction
    cont = -cont;
end
return cont
end

# Compute the contour integrals c_n, d_n, D_\infty, psi(z) [, its derivative and its contour integral].
# Input
#   alpha, beta - Parts of the weight w(x) = (1-x).^alpha.*(1+x).^beta.*h(x)
#   h           - Function for analytic function h(x)
#   nrT         - Number of terms, where 1 gives leading order term
#   [rho        - Berstein ellipse parameter: small to not encircle where log(h) is analytic and large to encircle the point x where to evaluate psi(x). Choose it to optimise the conditioning number as in [F. Bornemann, 2011] & [Wang & Huybrechs 2015], so it is a vector with rho(n+2) optimal for log(h)/sqrt(zeta-1)/sqrt(zeta+1)/(z-zeta)^(n+1) with rho+1/rho going to twice the distance of where log h is not analytic to the interval, or it is a scalar which is acceptable for all n.]
#   [M          - Number of points in trapezoidal rules, with M(n+2) optimal for c&d(n+1) and increasing with n, or M a scalar or NaN for successive doubling]
# Output
#   c, d        - The c- and d-coefficients, size(c&d) == ceil(Integer, nrT/2)
#   Dinf        - Limit of the Szego function
#   psi         - Anonymous function for the phase function
#   [dpsi       - Derivative of the phase function]
#   [contpsi	- Contour integral in psi, contpsi(x,p) has (zeta-x)^p]
function contour_integrals(alpha,beta,h,nrT,rho=false,M=NaN)
c = (1.0+1im)*zeros(ceil(Integer, nrT/2),1);
d = (1.0+1im)*zeros(ceil(Integer, nrT/2),1);
ra = 2*rand()-1; mm = 2; cc = 7;
# Checking for known exact solutions by comparing value at one random z = ra
if h(ra) == 1
    print("Using exact result for Legendre weight in contour_integrals",'\n')
    Dinf = 2.0^(-alpha/2-beta/2);
    function psi(x) 1/2.*( (alpha + beta).*acos(x+0.0im)-alpha*pi); end
    function dpsi(x) 1/2.*(alpha+beta).*(-1./(sqrt(1.0+0.0im-x).*sqrt(1.0+0.0im+x) ) ); end
    function contpsi(x,p) 0; end
elseif abs(h(ra) - exp(-cc*ra.^(2*mm)) ) <100*eps()
    print("Using exact result for Gibbs weight in contour_integrals",'\n')
    # binomial is only implemented for integers in Julia version 0.3.1
    function binom(x,n) prod(x:-1:(x-n+1) )/prod(1:n); end
    for n = 0:(2*mm-1)
	for j = 0:floor( (2*mm-n-1)/2)
	    # Actually the second binomial only has integers but seems to be converted to floats
	    c[n+1] = c[n+1] -cc*binom(j-1/2,j)*binom(2*mm-1-2*j,2*mm-n-1-2*j);
	end
	d[n+1] = (-1)^(n+1)*c[n+1];
    end
    Dinf = 2.0^(-alpha/2-beta/2)*exp(-cc/2*binom(mm-1/2,mm));
    degs = (0:(length(c)-1) );
    function psi(x) 1/2.*( (alpha + beta).*acos(x+0.0im)-alpha*pi) +sqrt(1.0+0.0im-x).*sqrt(1.0+0.0im+x)/2.*sum( c.*(x-1.0+0.0im).^(degs) ) ; end
    function dpsi(x) 1/2.*(alpha+beta).*(-1./(sqrt(1.0+0.0im-x).*sqrt(1.0+0.0im+x) ) ) -x./sqrt(1.0+0.0im-x)./sqrt(1.0+0.0im+x)/2.*sum(c.*(x-1.0+0.0im).^degs) +sqrt(1.0+0.0im-x).*sqrt(1.0+0.0im+x)/2.*sum(degs[2:end].*c[2:end].*(x-1.0+0.0im).^(degs[2:end]-1) ); end
   function contpsi(x,p) 2im*pi*sum( (degs[p:end]).^(p-1).*c[p:end].*(0im+x-1).^(degs[p:end]-p+1) ); end
else # Trapezoidal rules
    if abs(h(ra) - 1./sqrt(ra+3)) < 100*eps()
	print("Using trapezoidal rules for Fourier weight in contour_integrals", '\n')
	rho = 4; M = NaN; # M = 80; # But determine M automatically later
    elseif rho == false
	rho = 1.25; # > 1 so positive direction
    end
    lr = length(rho);
    if isnan(M)
	for nr = 1:(ceil(Integer, nrT/2)-1) # nr = n+1 in c_n = c(nr)
		c[nr] = trap_rule(z -> log(h(0.0im+z ) +0.0im)./( sqrt(z -1.0+0.0im).*sqrt(z +1.0+0.0im).*(z -1.0+0.0im).^nr ), rho[min(nr+1,lr)],M);
		d[nr] = trap_rule(z -> log(h(0.0im+z ) +0.0im)./( sqrt(z -1.0+0.0im).*sqrt(z +1.0+0.0im).*(z +1.0+0.0im).^nr ), rho[min(nr+1,lr)],M);
	end
	Dinf = real(2.0^( -(alpha+beta)./2).*exp(trap_rule(z -> log(h(0.0im+z)+0.0im )./sqrt(z -1.0+0.0im)./sqrt(z +1.0+0.0im), rho[1],M)/2) );
	    
	function psi(x) 1/2.*( (alpha + beta).*acos(x+0.0im)-alpha*pi) + sqrt(1.0+0.0im-x).*sqrt(1.0+0.0im+x).*trap_rule(z -> log(h(0.0im+z)+0.0im )./sqrt(z-1.0+0.0im)./sqrt(z+1.0+0.0im)./(0.0im+z-x), rho[min(2,lr)], M)/2; end

	function dpsi(x) 1/2.*(alpha + beta).*(-1)./(sqrt(1.0+0.0im-x).*sqrt(1.0+0.0im +x)) -x/(sqrt(1.0+0.0im-x).*sqrt(1.0+0.0im+x)*2)*trap_rule(z -> log(h(0.0im+z) +0.0im)./sqrt(z-1.0+0.0im)./sqrt(z+1.0+0.0im)./(0.0im+z-x),rho[min(2,lr)],M) + sqrt(1.0+0.0im-x).*sqrt(1.0+0.0im+x)/2*trap_rule(z -> log(h(0.0im+z)+0.0im )./sqrt(z-1.0+0.0im)./sqrt(z+1.0+0.0im)./(z-x+0.0im).^2,rho[min(3,lr)],M); end
	function contpsi(x,p) 2im*pi*trap_rule(z -> log(h(0.0im+z) )./sqrt(0.0im+z-1)./sqrt(0.0im+z+1)./(0.0im+z-x).^p, rho[min(2,lr)], M); end

   else
	lm = length(M);
	for nr = 1:(ceil(Integer, nrT/2)-1) # nr = n+1 in c_n = c(nr)
		c[nr] = trap_rule(z -> log(h(0.0im+z ) +0.0im)./( sqrt(z -1.0+0.0im).*sqrt(z +1.0+0.0im).*(z -1.0+0.0im).^nr ), rho[min(nr+1,lr)],M[min(nr+1, lm)]);
		d[nr] = trap_rule(z -> log(h(0.0im+z ) +0.0im)./( sqrt(z -1.0+0.0im).*sqrt(z +1.0+0.0im).*(z +1.0+0.0im).^nr ), rho[min(nr+1,lr)],M[min(nr+1, lm)]);
	end
	Dinf = real(2.0^( -(alpha+beta)./2).*exp(trap_rule(z -> log(h(0.0im+z)+0.0im )./sqrt(z -1.0+0.0im)./sqrt(z +1.0+0.0im), rho[1],M[1])/2) );
	    
	function psi(x) 1/2.*( (alpha + beta).*acos(x+0.0im)-alpha*pi) + sqrt(1.0+0.0im-x).*sqrt(1.0+0.0im+x).*trap_rule(z -> log(h(0.0im+z)+0.0im )./sqrt(z-1.0+0.0im)./sqrt(z+1.0+0.0im)./(0.0im+z-x), rho[min(2, lm)], M[min(2, lm)])/2; end

	function dpsi(x) 1/2.*(alpha + beta).*(-1)./(sqrt(1.0+0.0im-x).*sqrt(1.0+0.0im +x)) -x/(sqrt(1.0+0.0im-x).*sqrt(1.0+0.0im+x)*2)*trap_rule(z -> log(h(0.0im+z) +0.0im)./sqrt(z-1.0+0.0im)./sqrt(z+1.0+0.0im)./(0.0im+z-x),rho[min(2,lr)],M[min(2, lm)]) + sqrt(1.0+0.0im-x).*sqrt(1.0+0.0im+x)/2*trap_rule(z -> log(h(0.0im+z)+0.0im )./sqrt(z-1.0+0.0im)./sqrt(z+1.0+0.0im)./(z-x+0.0im).^2,rho[min(3,lr)],M[min(3, lm)]); end
	function contpsi(x,p) 2im*pi*trap_rule(z -> log(h(0.0im+z) )./sqrt(0.0im+z-1)./sqrt(0.0im+z+1)./(0.0im+z-x).^p, rho[min(2,lr)], M[min(3, lm)]); end
   end
end
return (c, d, Dinf, psi, dpsi,contpsi)
end

# Compute the W- or V-matrices to construct the asymptotic expansion of R 
# using the procedure with the convolutions as explained in the paper.
# Input
#   q            - alpha (when computing Wright) or beta (when computing Wleft)
#   t            - beta or alpha
#   Dinf         - Limit of the Szego function
#   cd           - c-or d- Coefficients needed for the U-matrices
#   maxOrder     - The maximum order of the error
#   r            - 1 when computing Wright, -1 when computing Wleft
#   isW          - true if the result are the Ws-, false if the V-s
#   [mos         - Maximum orders for each k, default ceil(Integer, (maxOrder-2)/2) for all k]
# Output
#   WVc          - Coefficient matrices for (z \pm 1)^m of s_k(z) or Delta_k(z)
function WV(q,t,Dinf,cd, maxOrder,r,isW, mos=ceil(Integer,  (maxOrder-2)/2)*ones(maxOrder-1) )
if maxOrder <= 1
	return (1+1im)*zeros(2,2,0,0);
end
mosDefault = norm(mos-ceil(Integer,  (maxOrder-2)/2)*ones(maxOrder-1)) == 0
(mo,~) = findmax(mos);
mo = convert(Int64,mo[1]);
ns = 0:mo;
# 'pochhammer' undefined in Julia
function poch(x,n) prod(x+(0:(n-1)) ); end
# binomial is only implemented for integers in version 0.3.1
function binom(x,n) prod(x:-1:(x-n+1) )/prod(1:n); end
# Will use gamma(1.0+n) instead of factorial(n) to avoid overflow
function brac(k) convert(FloatingPoint, (k == 0)) + convert(FloatingPoint, (k != 0))*prod(4*q^2 -(2*(1:k)-1).^2)/(2.0^(2*k)*(gamma(1.0+k) ) ); end

f = zeros(mo+1,1); # Coefficients in the expansion of \varphi(z)
for n = ns
    f[n+1] = poch(1/2,n)./(-r*2).^n./(1+2*n)./gamma(1.0+n);
end
# Or: f = poch(1/2,ns)./(-r*2).^ns./(1+2*ns)./gamma(1.0+ns);
g = zeros(max(mo*2+1,(maxOrder-1)),mo+1);
g[1,1] = 1; # Coefficients in the expansion of \varphi(z)^{-i}
for n = 1:mo
	g[1,n+1] = (-g[1,1:(n+1) ]*f[(n+1):-1:1])[1];
end

st = (1+1im)*zeros(max(mo*2+1,(maxOrder-1)),mo+1); # = \rho_{i,n,alpha+beta}
stp = (1+1im)*zeros(max(mo*2+1,(maxOrder-1)),mo+1); # = \rho_{i,n,alpha+beta+1}
stm = (1+1im)*zeros(max(mo*2+1,(maxOrder-1)),mo+1); # = \rho_{i,n,alpha+beta-1}
for n = ns
    tmp = 0;
    for j = 0:n
       tmp = tmp + binom(1/2,j)/(r*2)^j*cd[n-j+1];
    end
    st[1,n+1] = r*tmp + f[n+1]*(q+t);
    stp[1,n+1] = r*tmp + f[n+1]*(q+t+1);
    stm[1,n+1] = r*tmp + f[n+1]*(q+t-1);
end

for i = 2:max(mo*2+1,(maxOrder-1))
    for n = ns
        g[i,n+1] = sum(g[i-1,1:(n+1) ].*g[1,(n+1):-1:1] );
        st[i,n+1] = sum(st[i-1,1:(n+1) ].*st[1,(n+1):-1:1] );
        stp[i,n+1] = sum(stp[i-1,1:(n+1) ].*stp[1,(n+1):-1:1] );
        stm[i,n+1] = sum(stm[i-1,1:(n+1) ].*stm[1,(n+1):-1:1] );
    end
end

OmOdd = (1+1im)*zeros(mo+1,1); # = H_{n,alpha+beta}^{odd}
OmEven = (1+1im)*zeros(mo+1,1); # = H_{n,alpha+beta}^{even}
XiOdd = (1+1im)*zeros(mo+1,1); # = H_{n,alpha+beta+1}^{odd}
XiEven = (1+1im)*zeros(mo+1,1); # = H_{n,alpha+beta+1}^{even}
ThOdd = (1+1im)*zeros(mo+1,1); # = H_{n,alpha+beta-1}^{odd}
ThEven = (1+1im)*zeros(mo+1,1); # = H_{n,alpha+beta-1}^{even}
for n = ns
    if ( (mod(maxOrder,2)==0) || (n != mo) )
        OmEven[n+1] = OmEven[n+1] + st[1,n+1];
        XiEven[n+1] = XiEven[n+1] + stp[1,n+1];
        ThEven[n+1] = ThEven[n+1] + stm[1,n+1];
        for j = 1:n
            OmOdd[n+1] = OmOdd[n+1] + (r*2)^j*st[2*j,n-j+1]/gamma(1.0+2*j);
            XiOdd[n+1] = XiOdd[n+1] + (r*2)^j*stp[2*j,n-j+1]/gamma(1.0+2*j);
            ThOdd[n+1] = ThOdd[n+1] + (r*2)^j*stm[2*j,n-j+1]/gamma(1.0+2*j);
            OmEven[n+1] = OmEven[n+1] + (r*2)^j*st[2*j+1,n-j+1]/gamma(2.0+2*j);
            XiEven[n+1] = XiEven[n+1] + (r*2)^j*stp[2*j+1,n-j+1]/gamma(2.0+2*j);
            ThEven[n+1] = ThEven[n+1] + (r*2)^j*stm[2*j+1,n-j+1]/gamma(2.0+2*j);
        end
	else # This splitting is needed to avoid wrong OmEven[n+1] etc when n=mo and maxOrder is odd
        for j = 1:n
            OmOdd[n+1] = OmOdd[n+1] + (r*2)^j*st[2*j,n-j+1]/gamma(1.0+2*j);
            XiOdd[n+1] = XiOdd[n+1] + (r*2)^j*stp[2*j,n-j+1]/gamma(1.0+2*j);
            ThOdd[n+1] = ThOdd[n+1] + (r*2)^j*stm[2*j,n-j+1]/gamma(1.0+2*j);
        end
    end
end
OmOdd[1] = 1;
XiOdd[1] = 1; # Overwrite because else wrong
ThOdd[1] = 1;

OmO = (1+1im)*zeros(mo+1,1); # = X_{n,alpha+beta}^{odd}
OmE = (1+1im)*zeros(mo+1,1); # = X_{n,alpha+beta}^{even}
XiO = (1+1im)*zeros(mo+1,1); # = X_{n,alpha+beta+1}^{odd}
XiE = (1+1im)*zeros(mo+1,1); # = X_{n,alpha+beta+1}^{even}
ThO = (1+1im)*zeros(mo+1,1); # = X_{n,alpha+beta-1}^{odd}
ThE = (1+1im)*zeros(mo+1,1); # = X_{n,alpha+beta-1}^{even}
for n = ns
    for j = 0:n
        OmO[n+1] = OmO[n+1] + binom(-1/2,j)/(r*2)^j*OmOdd[n-j+1]/sqrt(r*2);
        XiO[n+1] = XiO[n+1] + binom(-1/2,j)/(r*2)^j*XiOdd[n-j+1]/sqrt(r*2);
        ThO[n+1] = ThO[n+1] + binom(-1/2,j)/(r*2)^j*ThOdd[n-j+1]/sqrt(r*2);
        if ( (mod(maxOrder,2)==0) || (n != mo) )
            OmE[n+1] = OmE[n+1] + binom(-1/2,j)/(r*2)^j*OmEven[n-j+1];
            XiE[n+1] = XiE[n+1] + binom(-1/2,j)/(r*2)^j*XiEven[n-j+1];
            ThE[n+1] = ThE[n+1] + binom(-1/2,j)/(r*2)^j*ThEven[n-j+1];
        end
    end
end

Ts = (1+1im)*zeros(2,2,mo+1);  # = T_{k,n}^{odd} or T_{k,n}^{even} depending on k, overwritten each new k
WVc = (1+1im)*zeros(2,2,maxOrder-1,mo+1);
for k = 1:(maxOrder-1)
    a = (q^2 + k/2-1/4)/k;
    b = -1im*r*(k-1/2);
    if (mosDefault)
        mo = ceil(Integer,  (maxOrder-1-k)/2)-1+ceil(Integer, k/2);
    else
        mo = mos[k];
    end
    Ts[:,:,:] = 0;
    if mod(k,2) == 1
# Actually, for m = (-ceil(Integer, k/2)):(ceil(Integer, (maxOrder-1-k)/2)-1) but it seems to be safer not to use negative indices so shift by ceil(Integer, k/2)+1 and easier to let n go to "mo" always iso ceil(Integer, (maxOrder-1-k)/2)+ceil(Integer, k/2)
        for n = 0:mo
	    # Or a*binom(-1/2,n)*(0im+r)^(n+1/2)/2^(n+1/2)*(2*n+1)/(2*n-1), or a*(-r*binom(-1/2,n)/(0im+r*2)^(n+1/2) -binom(-1/2,n-1)/(0im+r*2)^(n-1/2))			
	    Ts[:,:,n+1] = [ (-a*binom(n-3/2,n)*(-r*2.0)^(-n)*(2*n+1.0)*sqrt(0im+r/2.0) +1im*b*OmO[n+1] ) 	Dinf^2*(1im*a*binom(-1/2,n)*(r*2.0)^(-n)/sqrt(0im+r*2) +r*b*XiO[n+1] )	;	(1im*a*binom(-1/2,n)*(r*2.0)^(-n)/sqrt(0im+r*2.0) +r*b*ThO[n+1] )/Dinf^2	 (a*binom(n-3/2,n)*(-r*2.0)^(-n)*(2*n+1.0)*sqrt(0im+r/2.0) -1im*b*OmO[n+1] )];
            WVc[:,:,k,n+1] = brac(k-1)/(r*2.0)^(3*k/2.0)*sum(repeat(reshape(g[k,1:(n+1) ],1,1,convert(Int64,n+1)),outer=[2,2,1]).*Ts[:,:,(n+1):-1:1],3);
        end
    else
        for n = 0:mo
	    Ts[:,:,n+1] = [(a*convert(FloatingPoint, (n == 0)) +1im*r*b*OmE[n+1])     Dinf^2*b*XiE[n+1] 	; 	b*ThE[n+1]/Dinf^2	(a*convert(FloatingPoint, (n == 0)) -1im*b*r*OmE[n+1] )];
            WVc[:,:,k,n+1] = brac(k-1)/(r*2.0)^(3*k/2.0)*sum(repeat(reshape(g[k,1:(n+1) ],1,1,convert(Int64,n+1)),outer=[2,2,1]).*Ts[:,:,(n+1):-1:1],3);
            if isW # The extra matrix of the theorem:
                WVc[:,:,k,n+1] = WVc[:,:,k,n+1] + brac(k-1)/(r*2.0)^(3*k/2.0)*(-(4*q^2+2*k-1)*g[k,n+1]/2/k*eye(2) );
            end
        end
    end
end
return WVc
end


# Get the U-matrices to construct the asymptotic expansion of R using the procedure with the convolutions as explained in the paper. Optionally, specify method or number of matrices to compute or get the Q-s.
# Input
#   alpha, beta  - Parts of the weight w(x) = (1-x).^alpha.*(1+x).^beta.*h(x)
#   Dinf         - Limit of the Szego function
#   c,d          - Coefficients needed for the U-matrices
#   maxOrder     - The maximum order of the error
#   [method       - 'UW' (default) to get the U-s through the W-s, 'UQW' to get the U- & Q-s through the W-s and 'UQV' to get the U- & Q-s through the V-s]
#   [ns          - If method contains 'Q', specify how many Q-s to compute so compute Q_{k,0:(ns(k)-1)}^{right/left}, also for W/V]
# Output
#   Uright       - Coefficient matrices of R_k(z) for (z-1)^(-m)
#   Uleft        - Coefficient matrices of R_k(z) for (z+1)^(-m)
#   [Qright      - Coefficient matrices of R_k^{right}(z) for (z-1)^n]
#   [Qleft       - Coefficient matrices of R_k^{left}(z) for (z+1)^n]
function UQ(alpha,beta,Dinf,c,d,maxOrder,method = "UW",ns= ceil(Integer,  (maxOrder-1-(1:(maxOrder-1)))/2)-1+ceil(Integer, (1:(maxOrder-1))/2) )

tmp = ceil(Integer,  (maxOrder-1-(1:(maxOrder-1)))/2)-1+ceil(Integer, (1:(maxOrder-1))/2)
existed = norm(ns[:]-tmp) != 0;
if ~existed
    if ismatch(r"Q",method)
	ns = convert(Array{Int64,1},maximum([ (maxOrder-2-(1:(maxOrder-1))+1+ceil(Integer, (1:(maxOrder-1))/2) )     zeros(1,maxOrder-1)[:]],2)[:] );
    else 
        # Could use the following to specify the necessary (exact) number of W- or V-matrices to compute, but not really needed...
#         ns = ceil(Integer,  (maxOrder-1-(1:(maxOrder-1)))/2)-1+ceil(Integer, (1:(maxOrder-1))/2);
    end
end

Uright = (1+1im)*zeros(2,2,maxOrder-1,ceil(Integer,  (maxOrder-1)/2) );
Uleft = (1+1im)*zeros(2,2,maxOrder-1,ceil(Integer,  (maxOrder-1)/2));
# About half of these tensors will not be used
if ismatch(r"Q",method)
    tmp = ceil(Integer,  (maxOrder-1-(1:(maxOrder-1)))/2)-1+ceil(Integer, (1:(maxOrder-1))/2)
    if norm(ns[:]-tmp) != 0
	(mn,ignore) = findmax(ns);
	mn = convert(Int64,mn);
        Qright = (1+1im)*zeros(2,2,maxOrder-1,mn+1 );
        Qleft = (1+1im)*zeros(2,2,maxOrder-1,mn+1 );
    else
        Qright = (1+1im)*zeros(2,2,maxOrder-1,convert(Int64,ceil( (maxOrder-1)/2) ));
        Qleft = (1+1im)*zeros(2,2,maxOrder-1,convert(Int64,ceil( (maxOrder-1)/2) ));
    end
end

function poch(x,n) prod(x+(0:(n-1) ) ); end
# 'pochhammer' undefined in Julia

if ismatch(r"W",method)
    tmp = ceil(Integer,  (maxOrder-1-(1:(maxOrder-1)))/2)-1+ceil(Integer, (1:(maxOrder-1))/2)
    if norm(ns[:]-tmp) != 0
	Wr = WV(alpha,beta,Dinf,c,maxOrder,1,true,ns[:]+2+ceil(Integer,  (1:(maxOrder-1)) /2) ); 
	# Wr = W_{k,m}^{right}, k=1..maxOrder-1, m given by ns; Wl analogous.
	Wl = WV(beta,alpha,Dinf,d,maxOrder,(1im)^2,true,ns[:]+2+ceil(Integer,  (1:(maxOrder-1)) /2));
    else
	Wr = WV(alpha,beta,Dinf,c,maxOrder,1,true);
	Wl = WV(beta,alpha,Dinf,d,maxOrder,(1im)^2,true);
    end
    for k = 1:(maxOrder-1)
        # Actually, for m = (-ceil(Integer, k/2)):1, but better no negative indices so shift by ceil(Integer, k/2)+1
        for m = 1:ceil(Integer, k/2)
            Uright[:,:,k,m] = Wr[:,:,k,ceil(Integer, k/2)+1-m];
            Uleft[:,:,k,m] = Wl[:,:,k,ceil(Integer, k/2)+1-m];
            for j=1:(k-1)
                for l = max(m-ceil(Integer, j/2),1):ceil(Integer, (k-j)/2)
                    Uright[:,:,k,m] = Uright[:,:,k,m] + Uright[:,:,k-j,l]*Wr[:,:,j,ceil(Integer, j/2)+1+l-m];
                    Uleft[:,:,k,m] = Uleft[:,:,k,m] + Uleft[:,:,k-j,l]*Wl[:,:,j,ceil(Integer, j/2)+1+l-m];
                end
            end
            for j=1:(k-1)
                for n = 0:(ceil(Integer, j/2)-m)
                    for i = 1:ceil(Integer, (k-j)/2)
                        Uright[:,:,k,m] = Uright[:,:,k,m] + poch(1-i-n,n)/2.0^(i+n)/gamma(1.0+n)*Uleft[:,:,k-j,i]*Wr[:,:,j,ceil(Integer, j/2)+1-n-m];
                        Uleft[:,:,k,m] = Uleft[:,:,k,m] + poch(1-i-n,n)/(-2.0)^(i+n)/gamma(1.0+n)*Uright[:,:,k-j,i]*Wl[:,:,j,ceil(Integer, j/2)+1-n-m];
                    end
                end
            end
        end
    end
    if method == "UW"
	return (Uright,Uleft);
    end
elseif ismatch(r"V",method)
    tmp = ceil(Integer,  (maxOrder-1-(1:(maxOrder-1)))/2)-1+ceil(Integer, (1:(maxOrder-1))/2)
    if norm(ns[:]-tmp) != 0
		Vr = WV(alpha,beta,Dinf,c,maxOrder,1,false,ns[:]+2+ceil(Integer,  (1:(maxOrder-1)) /2) );
		Vl = WV(beta,alpha,Dinf,d,maxOrder,(1im)^2,false,ns[:]+2+ceil(Integer,  (1:(maxOrder-1)) /2) );
	else
		Vr = WV(alpha,beta,Dinf,c,maxOrder,1,false);
		Vl = WV(beta,alpha,Dinf,d,maxOrder,(1im)^2,false);
	end
else
    error(string("Wrong method format, got ", method));
end

if "UQW" == method
    # Already did the U's, calculate the Q-s from them:
    for k = 1:(maxOrder-1)
        for n = 0:ns[k]
            Qright[:,:,k,n+1] = -Wr[:,:,k,ceil(Integer, k/2)+1+n];
            Qleft[:,:,k,n+1] = -Wl[:,:,k,ceil(Integer, k/2)+1+n];
            for i = 1:ceil(Integer, k/2)
                Qright[:,:,k,n+1] = Qright[:,:,k,n+1] + poch(1-i-n,n)/2.0^(i+n)/gamma(1.0+convert(Int64,n) )*Uleft[:,:,k,i];
                Qleft[:,:,k,n+1] = Qleft[:,:,k,n+1] + poch(1-i-n,n)/(-2.0)^(i+n)/gamma(1.0+convert(Int64,n) )*Uright[:,:,k,i];
            end
            for j = 1:(k-1)
                for i = -ceil(Integer, j/2):n
                    for l = 1:ceil(Integer, (k-j)/2)
                        Qright[:,:,k,n+1] = Qright[:,:,k,n+1] - poch(1-l-n+i,n-i)/2.0^(n-i+l)/gamma(1.0+convert(Int64,n-i) )*Uleft[:,:,k-j,l]*Wr[:,:,j,i+1+ceil(Integer, j/2) ];
                        Qleft[:,:,k,n+1] = Qleft[:,:,k,n+1] - poch(1-l-n+i,n-i)/(-2.0)^(n-i+l)/gamma(1.0+convert(Int64,n-i) )*Uright[:,:,k-j,l]*Wl[:,:,j,i+1+ceil(Integer, j/2) ];
                    end
                end
                for i = (n+1):(n+ceil(Integer, (k-j)/2) )
                    Qright[:,:,k,n+1] = Qright[:,:,k,n+1] -Uright[:,:,k-j,i-n]*Wr[:,:,j,i+1+ceil(Integer, j/2) ];
                    Qleft[:,:,k,n+1] = Qleft[:,:,k,n+1] -Uleft[:,:,k-j,i-n]*Wl[:,:,j,i+1+ceil(Integer, j/2) ];
                end
            end
        end
    end
elseif "UQV"== method
    for k = 1:(maxOrder-1)
        for m = 1:ceil(Integer, k/2)
            Uright[:,:,k,m] = Vr[:,:,k,ceil(Integer, k/2)+1-m];
            Uleft[:,:,k,m] = Vl[:,:,k,ceil(Integer, k/2)+1-m];
            for j=1:(k-1)
                for l = 0:(ceil(Integer, j/2)-m)
                    Uright[:,:,k,m] = Uright[:,:,k,m] + Qright[:,:,k-j,l+1]*Vr[:,:,j,ceil(Integer, j/2)+1-l-m];
                    Uleft[:,:,k,m] = Uleft[:,:,k,m] + Qleft[:,:,k-j,l+1]*Vl[:,:,j,ceil(Integer, j/2)+1-l-m];
                end
            end
        end
        for n = 0:ns[k]
            Qright[:,:,k,n+1] = -Vr[:,:,k,ceil(Integer, k/2)+1+n];
            Qleft[:,:,k,n+1] = -Vl[:,:,k,ceil(Integer, k/2)+1+n];
            for i = 1:ceil(Integer, k/2)
                Qright[:,:,k,n+1] = Qright[:,:,k,n+1] + poch(1-i-n,n)/2.0^(i+n)/gamma(1.0+convert(Int64,n) )*Uleft[:,:,k,i];
                Qleft[:,:,k,n+1] = Qleft[:,:,k,n+1] + poch(1-i-n,n)/(-2.0)^(i+n)/gamma(1.0+convert(Int64,n) )*Uright[:,:,k,i];
            end
            for j = 1:(k-1)
                for l = 0:(ceil(Integer, j/2)+n)
                    Qright[:,:,k,n+1] = Qright[:,:,k,n+1] -Qright[:,:,k-j,l+1]*Vr[:,:,j,n-l+1+ceil(Integer, j/2) ];
                    Qleft[:,:,k,n+1] = Qleft[:,:,k,n+1] -Qleft[:,:,k-j,l+1]*Vl[:,:,j,n-l+1+ceil(Integer, j/2) ];
                end
            end
        end
    end
else
    error(string("Invalid method, got ", method));
end
return (Uright,Uleft,Qright,Qleft)
end
