# asy_expansions.jl: contains all asymptotic expansions, needing results from the precomputations
# About
#   Author       - Peter Opsomer (Peter.Opsomer@cs.kuleuven.be)
#   History      - Created October 2014, last edit October 2015


# Evaluate the asymptotic expansion for the generalised Jacobi polynomial in the lens-shaped region with a certain normalisation. Derivative is only computed when dw is specified, default is monic normalisation without derivatives.
# Input
#   n,z          - Degree and point at which to evaluate
#   alpha,beta,h - Parts of the weight 
#   psi(x)       - Function for analytic phase function psi(x)
#   nrT          - Number of terms. 1 gives leading order term
#   Dinf         - Limit of the Szego function
#   Uright,Uleft - Right and left U-matrices
#   [nor         - Normalisation: 'm' for monic(default), 'o' for orthonormal]
#   [dh,dpsi     - Derivatives of part of weight and phase function]
# Output
#   ipi          - Asymptotic expansion of the polynomial in the lens
#   [dipi        - Derivative of polynomial in the lens]
function asy_lens(n,z,alpha,beta,h,psi,nrT,Dinf,Uright,Uleft,nor = "m",dh = false,dpsi=false)

if ismatch(r"o",nor)
    if (nrT == 1) || (length(Uright) == 0) 	nmz = 1;     else
	# In Julia, 1:(nrT-1) already is a column vector instead of a row
	nmz = sqrt(real(1+2*1im*Dinf^2*sum(reshape(Uright[2,1,1:(nrT-1),1] + Uleft[2,1,1:(nrT-1),1], (nrT-1), 1)./(n+1).^(1:(nrT-1) ) ) ) );
    end
end
function w(x) return (0im+1-x).^alpha.*(0im+1+x).^beta.*h(0im+x); end

RI = (1+0im)*eye(2);
for k = 1:nrT-1
    for m = 1:ceil(Integer, k/2)
	RI = RI + (Uright[:,:,k,m]./(z-1).^m + Uleft[:,:,k,m]./(z+1).^m)/n.^k;
    end
end
if nor == "m"
    ipi = ([2.0.^(1/2-n)./(sqrt(0im+w(z) ).*(1.0+z+0.0im).^(1/4).*(1.0-z+0.0im).^(1/4) )   0;]*RI*[Dinf*cos((n+1/2)*acos(z) +psi(z) -pi/4); -1im/Dinf*cos((n-1/2)*acos(z) +psi(z) -pi/4)])[1,1];
elseif nor == "o"
    ipi = ([(2/pi)^(1/2)./(sqrt(0im+w(z) ).*(1.0+z+0.0im).^(1/4).*(1.0-z+0.0im).^(1/4) )   0;]*RI*[cos((n+1/2)*acos(z) +psi(z) -pi/4)  ;  -1im/Dinf^2*cos((n-1/2)*acos(z) +psi(z) -pi/4)]*nmz)[1,1];
else
    error("Wrong normalisation")
end
if dh == false return ipi end

function dw(x) return -alpha*(1-x).^(alpha-1).*(1+x).^beta.*h(x)  +(1-x).^alpha.*beta*(1+x).^(beta-1).*h(x) +(1-x).^alpha.*(1+x).^beta.*dh(x); end

dRI = (1+1im)*zeros(2,2);
for k = 1:nrT-1
    for m = 1:ceil(Integer, k/2)
	dRI = dRI + (Uright[:,:,k,m]./(z-1).^(m+1)*(-m) + Uleft[:,:,k,m]./(z+1).^(m+1)*(-m))/n.^k;
    end
end

if nor == "m"
    dipi = (([2.0.^(1/2-n).*(dw(z)*(-1/2)./(sqrt(0im+w(z) ).^3.*(1.0+z+0.0im).^(1/4).*(1.0+0.0im-z).^(1/4) ) -1/4*(-2*z)./(sqrt(0im+w(z) ).*(1.0+0.0im+z).^(5/4).*(1.0+0.0im-z).^(5/4) ) )    0;]*RI + [2.0.^(1/2-n)./(sqrt(0im+w(z) ).*(1.0+0.0im+z).^(1/4).*(1.0+0.0im-z).^(1/4) )   0;]*dRI)*[Dinf*cos((n+1/2)*acos(z) +psi(z) -pi/4); -1im/Dinf*cos((n-1/2)*acos(z) +psi(z) -pi/4)]+ [2.0.^(1/2-n)./(sqrt(0im+w(z) ).*(1.0+0.0im+z).^(1/4).*(1.0+0.0im-z).^(1/4) )     0;]*RI*[-Dinf*sin((n+1/2)*acos(z) +psi(z) -pi/4)*(dpsi(z)-(n+1/2)./sqrt(0im+1.0-z)./sqrt(0im+1.0+z) ); 1im/Dinf*sin((n-1/2)*acos(z) +psi(z) -pi/4)*(dpsi(z)-(n-1/2)./sqrt(0im+1.0-z)./sqrt(0im+1.0+z) )])[1,1];
elseif nor == "o"
    dipi = ((([sqrt(0im+2/pi).*(dw(z)*(-1/2)./(sqrt(0im+w(z) ).^3.*(1.0+0.0im+z).^(1/4).*(1.0+0.0im-z).^(1/4) ) -1/4*(-2*z)./(sqrt(0im+w(z) ).*(1.0+0.0im+z).^(5/4).*(1.0+0.0im-z).^(5/4) ) )    0;]*RI+ [sqrt(0im+2/pi)./(sqrt(0im+w(z) ).*(1.0+0.0im+z).^(1/4).*(1.0+0.0im-z).^(1/4) )   0;]*dRI)*[cos((n+1/2)*acos(z) +psi(z) -pi/4); -1im/Dinf^2*cos((n-1/2)*acos(z) +psi(z) -pi/4)] + [sqrt(0im+2/pi)./(sqrt(0im+w(z) ).*(1.0+0.0im+z).^(1/4).*(1.0+0.0im-z).^(1/4) )     0;]*RI* [-sin((n+1/2)*acos(z) +psi(z) -pi/4)*(dpsi(z)-(n+1/2)./sqrt(0im+1.0-z)./sqrt(0im+1.0+z) ); 1im/Dinf^2*sin((n-1/2)*acos(z) +psi(z) -pi/4)*(dpsi(z)-(n-1/2)./sqrt(0im+1.0-z)./sqrt(0im+1.0+z) )])*nmz)[1,1];
end
return (ipi, dipi)
end


# Evaluate the asymptotic expansion for the generalised Jacobi polynomial in the outer region with a certain normalisation. Derivative is only computed when dw is specified, default is monic normalisation without derivatives.
# Input
#   n,z          - Degree and point at which to evaluate
#   alpha,beta,h - Parts of the weight 
#   psi(x)       - Function for analytic phase function psi(x)
#   nrT          - Number of terms. 1 gives leading order term
#   Dinf         - Limit of the Szego function
#   Uright,Uleft - Right and left U-matrices
#   [nor         - Normalisation: 'm' for monic(default), 'o' for orthonormal]
#   [dh,dpsi     - Derivatives of part of weight and phase function]
# Output
#   opi          - Asymptotic expansion of the outer polynomial
#   [dopi        - Derivative of outer polynomial]
function asy_outer(n,z,alpha,beta,h,psi,nrT,Dinf,Uright,Uleft,nor = "m",dh = false,dpsi=false)

if ismatch(r"o",nor)
    if (nrT == 1) || (length(Uright) == 0) nmz = 1; else
	nmz = sqrt(real(1+2*1im*Dinf^2*sum(reshape(Uright[2,1,1:(nrT-1),1] + Uleft[2,1,1:(nrT-1),1], (nrT-1), 1)./(n+1).^(1:(nrT-1) ) ) ) );
    end
end
function w(x) return (1-x).^alpha.*(1+x).^beta.*h(x); end

RI = (1+0im)*eye(2);
for k = 1:nrT-1
    for m = 1:ceil(Integer, k/2)
	RI = RI + (Uright[:,:,k,m]./(z-1).^m + Uleft[:,:,k,m]./(z+1).^m)/n.^k;
    end
end

if nor == "m"
    opi = ([2.0.^(-1/2-n+0.0)./(sqrt(0im+w(z)).*(1.0+0.0im+z).^(1/4).*(1.0+0.0im-z).^(1/4) )   0;]*RI*[Dinf*exp((-1.0+0im)^(angle(z) <= 0)*1im*((n+1/2)*acos(z) +psi(z) -pi/4) );-1im/Dinf*exp((-1.0+0im)^(angle(z) <= 0)*1im*((n-1/2)*acos(z) +psi(z) -pi/4))])[1,1];
elseif nor == "o"
    opi = ([(2*pi).^(-1/2+0.0)./(sqrt(0im+w(z)).*(1.0+0.0im+z).^(1/4).*(1.0+0.0im-z).^(1/4) )    0;]*RI*[exp((-1.0+0im)^(angle(z) <= 0)*1im*((n+1/2)*acos(z) +psi(z) -pi/4) );   -1im/Dinf^2*exp((-1.0+0im)^(angle(z) <= 0)*1im*((n-1/2)*acos(z) +psi(z) -pi/4))]*nmz)[1,1];
else
    error("Wrong normalisation")
end

if dh == false return opi end

function dw(x) return -alpha*(1-x).^(alpha-1).*(1+x).^beta.*h(x)  +(1-x).^alpha.*beta*(1+x).^(beta-1).*h(x) +(1-x).^alpha.*(1+x).^beta.*dh(x); end

dRI = (1+1im)*zeros(2,2);
for k = 1:nrT-1
    for m = 1:ceil(Integer, k/2)
	dRI = dRI + (Uright[:,:,k,m]./(z-1).^(m+1)*(-m) + Uleft[:,:,k,m]./(z+1).^(m+1)*(-m))/n.^k;
    end
end
if nor=="m"
    dopi = (([2.0.^(-1/2-n+0.0).*(dw(z)*(-1/2)./(sqrt(0im+w(z) ).^3.*(1.0+0.0im+z).^(1/4).*(1.0+0.0im-z).^(1/4) ) -1/4*(-2*z)./(sqrt(0im+w(z) ).*(1.0+0.0im+z).^(5/4).*(1.0+0.0im-z).^(5/4) ) )    0;]*RI  + [2.0.^(0.0-1/2-n)./(sqrt(0im+w(z) ).*(1.0+0.0im+z).^(1/4).*(1.0+0.0im-z).^(1/4) )    0;]*dRI)* [Dinf*exp((-1.0+0im)^(angle(z) <= 0)*1im*((n+1/2)*acos(z) +psi(z) -pi/4) );  -1im/Dinf*exp((-1.0+0im)^(angle(z) <= 0)*1im*((n-1/2)*acos(z) +psi(z) -pi/4))]    + [2.0.^(-1/2-n+0.0)./(sqrt(0im+w(z)).*(1.0+0.0im+z).^(1/4).*(1.0+0.0im-z).^(1/4) )    0;]*RI*[Dinf*exp((-1.0+0im)^(angle(z) <= 0)*1im*((n+1/2)*acos(z) +psi(z) -pi/4) )*(-1.0+0im)^(angle(z) <= 0)*1im*(dpsi(0.0im+z)-(n+1/2)./sqrt(0im+1.0-z)./sqrt(0im+1.0+z) );-1im/Dinf*exp((-1.0+0im)^(angle(z) <= 0)*1im*((n-1/2)*acos(z) +psi(z) -pi/4) )*(-1.0+0im)^(angle(z) <= 0)*1im*(dpsi(0.0im+z)-(n-1/2)./sqrt(0im+1.0-z)./sqrt(0im+1.0+z) )])[1,1];
elseif nor=="o"
    dopi = ((([1/sqrt(0im+2*pi).*(dw(z)*(-1/2)./(sqrt(0im+w(z) ).^3.*(1.0+0.0im+z).^(1/4).*(1.0+0.0im-z).^(1/4) ) -1/4*(-2*z)./(sqrt(0im+w(z) ).*(1.0+0.0im+z).^(5/4).*(1.0+0.0im-z).^(5/4) ) )    0;]*RI    + [1/sqrt(0im+2*pi)./(sqrt(0im+w(z) ).*(1.0+0.0im+z).^(1/4).*(1.0+0.0im-z).^(1/4) )    0;]*dRI)*   [exp((-1.0+0im)^(angle(z) <= 0)*1im*((n+1/2)*acos(z) +psi(z) -pi/4) );-1im/Dinf^2*exp((-1.0+0im)^(angle(z) <= 0)*1im*((n-1/2)*acos(z) +psi(z) -pi/4))] + [1/sqrt(0im+2*pi)./(sqrt(0im+w(z)).*(1.0+0.0im+z).^(1/4).*(1.0+0.0im-z).^(1/4) )     0;]*RI* [exp((-1.0+0im)^(angle(z) <= 0)*1im*((n+1/2)*acos(z) +psi(z) -pi/4) )*(-1.0+0im)^(angle(z) <= 0)*1im*(dpsi(0.0im+z)-(n+1/2)./sqrt(0im+1.0-z)./sqrt(0im+1.0+z) );-1im/Dinf^2*exp((-1.0+0im)^(angle(z) <= 0)*1im*((n-1/2)*acos(z) +psi(z) -pi/4) )* (-1.0+0im)^(angle(z) <= 0)*1im*(dpsi(0.0im+z)-(n-1/2)./sqrt(0im+1.0-z)./sqrt(0im+1.0+z) )])*nmz)[1,1];
end
return (opi, dopi)
end


#Evaluate the asymptotic expansion for the generalised Jacobi polynomial in the right boundary region with a certain method. Derivative is only computed when assigned, default is monic normalisation without derivatives or Q-s. When using the Q-s, a series of the cosine transform of the asymptotic expansion is also used and theta=arccos(z).
# Input
#   n,z          - Degree and point at which to evaluate.
#   alpha,beta,h - Parts of the weight 
#   psi(x)       - Function for analytic phase function psi(x)
#   nrT          - Number of terms. 1 gives leading order term
#   Dinf         - Limit of the Szego function
#   Uright,Uleft - Right and left U-matrices
#   [method        - 'm' for monic(default), 'o' for orthonormal, 'mQ' for monic using the Q-s, 'oQ' for orthonormal with Q-s]
#   [dh,dpsi     - Derivatives of part of weight and phase function]
#   [Qright      - Right Q-matrices for when ismatch(r"Q",method)]
#   [theta       - theta=arccos(z) when computing with Q-s to save digits]
#   [contpsi	 - Contour integral in psi, needed for series expansions]
# Output
#   bpiR         - Asymptotic expansion of the right boundary polynomial
#   [dbpiR       - Derivative of right boundary polynomial]
function asy_right(n,z,alpha,beta,h,psi,nrT,Dinf,Uright,Uleft,method = "m",dh=false,dpsi=false,Qright=false,theta=false,contpsi=false)

if ismatch(r"o",method)
    if (nrT == 1) || (length(Uright) == 0) nmz = 1; else
	nmz = sqrt(real(1+2*1im*Dinf^2*sum(reshape(Uright[2,1,1:(nrT-1),1] + Uleft[2,1,1:(nrT-1),1], (nrT-1), 1)./(n+1).^(1:(nrT-1) ) ) ) );
    end
end

if ismatch(r"Q",method)
    if (abs(z-1) >1e-2) || (abs(theta) > 1e-2)
	warn("Theta is not small so using Q-s will probably give an inaccurate result")
    end
    mo = size(Qright,4)-1;
    sbes = 0+0im;	sden = 0+0im;	sthosth = 0+0im;	ssinpl = 0+0im;	ssinmi = 0+0im;	sdbes = 0+0im;
    for nn = 0:mo # Replaced factorial(n) by gamma(1.0+n) to avoid overflow
	sbes = sbes + (-1)^nn*(n/2)^(alpha+2*nn)*theta^(2*nn)/gamma(1+nn+alpha)/gamma(1.0+nn);
	sden = sden + (-1)^nn*theta^(2*nn)/gamma(3+2*nn);
	sthosth = sthosth + (-1)^nn*theta^(2*nn)/gamma(2.0+2*nn);
	ssinpl = ssinpl + (-1)^nn*(psi(z)+alpha*pi/2+theta/2).^(2*nn)/gamma(2.0+2*nn);
	ssinmi = ssinmi + (-1)^nn*(psi(z)+alpha*pi/2-theta/2).^(2*nn)/gamma(2.0+2*nn);
	sdbes = sdbes + (-1)^nn/gamma(1.0+nn)*(n/2)^(2*nn+alpha+1)*(4/n^2*theta^(2*nn)/gamma(nn+alpha)-theta^(2*nn+2)/gamma(2+nn+alpha))/2;
    end
    plcoef = (alpha+beta+1)/2 +sthosth/(4im*pi)*contpsi(z,1);
    micoef = (alpha+beta-1)/2 +sthosth/(4im*pi)*contpsi(z,1);
    sbesowei = sbes*(sden)^(-alpha/2+0.0);
    splsindbes = plcoef*ssinpl*sdbes*sden^(-alpha/2+0.0);
    smisindbes = micoef*ssinmi*sdbes*sden^(-alpha/2+0.0);

    RR = (1+1im)*zeros(2,2);
    for k = (nrT-1):-1:1 # Counting down to avoid roundoff error
	tmp = zeros(2,2);
	for i = 1:size(Qright,4)
	    tmp = tmp + Qright[:,:,k,i]*(-theta^2*sden)^(i-1);
	end
	RR = (RR + tmp)/n;
    end
    RR = RR + eye(2);
else
    function w(x) return (1-x).^alpha.*(1+x).^beta.*h(x); end
    if abs(z-1) < eps()^(1/3)
	warn("z is close to 1 so use Q-s and series expansions in 'method'")
    end
    function brac(k) convert(FloatingPoint, (k == 0)) + convert(FloatingPoint, (k != 0))*prod(4*alpha^2 -(2*(1:k)-1).^2)/(2.0^(2*k)*(gamma(1.0+k) ) ); end
    
    function DeltaR(k,x) brac(k-1)/2.0^k./(log(x + sqrt(0im+x-1).*sqrt(0im+x+1)+0im).^k)./(2*sqrt(0im+x+1.0)*sqrt(0im+x-1.0))*[Dinf 0; 0 Dinf^(-1+0.0)]*		[sqrt(0im+x + sqrt(0im+x-1.0).*sqrt(0im+x+1.0) )    1im/sqrt(0im+x + sqrt(0im+x-1.0).*sqrt(0im+x+1.0) ); 		-1im/sqrt(0im+x + sqrt(0im+x-1.0).*sqrt(0im+x+1.0) )    sqrt(0im+x + sqrt(0im+x-1.0).*sqrt(0im+x+1.0) ) ]*		[exp( 1im*(-1.0+0im)^(angle(x-1) <= 0)*(psi(x) +alpha*pi/2) )    0  ;   0    exp(-1im*(-1.0+0im)^(angle(x-1) <= 0)*(psi(x) +alpha*pi/2))]*		[ ((-1.0)^k)/k*(alpha^2+k/2-1/4)    -1im*(k-1/2)   ;   (-1.0)^k*1im*(k-1/2)    (alpha^2+k/2-1/4)/k]*		[exp(-1im*(-1.0+0im)^(angle(x-1) <= 0)*(psi(x) +alpha*pi/2))   0 ;   0   exp(1im*(-1.0+0im)^(angle(x-1) <= 0)*(psi(x) +alpha*pi/2) )]*		[sqrt(0im+x + sqrt(0im+x-1.0).*sqrt(0im+x+1.0) )     -1im/sqrt(0im+x + sqrt(0im+x-1.0).*sqrt(0im+x+1.0) )   ;   	1im/sqrt(0im+x + sqrt(0im+x-1.0).*sqrt(0im+x+1.0) )   sqrt(0im+x + sqrt(0im+x-1.0).*sqrt(0im+x+1.0) ) ]*[Dinf^(-1+0.0) 0; 0 Dinf]; end
	
    s = (1+1im)*zeros(2,2,nrT-1);
    R = (1+1im)*zeros(2,2,nrT-1);
    RR = eye(2) + (1+1im)*zeros(2,2);
    for k = 1:nrT-1
	for m = 1:ceil(Integer, k/2)
	    R[:,:,k] = R[:,:,k] + Uright[:,:,k,m]./(z-1).^m + Uleft[:,:,k,m]./(z+1).^m;
	end
	s[:,:,k] = DeltaR(k,z);
	if mod(k,2)==0 # k is even: add extra matrix
	    s[:,:,k] = s[:,:,k] -brac(k-1)*(4*alpha^2+2*k-1)/2.0^(k+1)/k/(log(z+sqrt(0im+z-1).*sqrt(0im+z+1) +0im))^k*eye(2);
	end
	RR = RR + (R[:,:,k]-s[:,:,k] )/n^k; # Avoids defining R_0^{right} = I
	for j = 1:(k-1)
	    RR = RR - R[:,:,k-j]*s[:,:,j]/n^k;
	end
    end
    b1 = cos(psi(z) +alpha*pi/2 + acos(z)/2)*besselj(alpha,n*acos(z) ) + sin(psi(z) +alpha*pi/2 + acos(z)/2)*(besselj(alpha-1,n*acos(z) ) -       besselj(alpha+1,n*acos(z) ) )/2;
    b2 = cos(psi(z) +alpha*pi/2 -acos(z)/2)*   besselj(alpha,n*acos(z) ) + sin(psi(z) +alpha*pi/2 -acos(z)/2)*   (besselj(alpha-1,n*acos(z) ) - besselj(alpha+1,n*acos(z) ) )/2;
end

if method == "m"
    bpiR = ([sqrt(0im+pi*n*acos(z))./(2.0^n.*sqrt(0im+w(z) ).*(1.0+0.0im+z).^(1/4).*(1.0+0.0im-z).^(1/4) )     0;]*RR*[Dinf*b1; -1im/Dinf*b2])[1,1];
elseif method == "o"
    bpiR = ([sqrt(0im+n*acos(z))./(sqrt(0im+w(z) ).*(1.0+0.0im+z).^(1/4).*(1.0+0.0im-z).^(1/4) )    0;]*RR*[b1; (-1im/Dinf^2*b2)]*nmz)[1,1]; 
elseif method == "mQ"
    bpiR = ([sqrt(0im+pi*n/sthosth)./(2.0.^n.*sqrt(0im+(1.0+0im+z).^beta.*h(z) ) )      0;]*RR*[Dinf*(cos(0im+psi(0im+z) +alpha*pi/2 + theta/2)*sbesowei +splsindbes); -1im/Dinf*(cos(0im+psi(0im+z) +alpha*pi/2 -theta/2)*sbesowei + smisindbes)])[1,1];
elseif method == "oQ"
    bpiR = ([sqrt(0im+n/sthosth)./(sqrt(0im+(1.0+0im+z).^beta.*h(z) ) )     0;]*RR*[(cos(0im+psi(0im+z) +alpha*pi/2 + theta/2)*sbesowei +splsindbes); -1im/Dinf^2*(cos(0im+psi(0im+z) +alpha*pi/2 -theta/2)*sbesowei + smisindbes)]*nmz)[1,1];

else
    error("Wrong method")
end

if dh == false return bpiR end

dRR = (1.0+1.0im)*zeros(2,2);
if ismatch(r"Q",method)
    dsbes = 0+0im;	dsden = 0+0im;	dsthosth = 0+0im;	dssinpl = 0+0im;	dssinmi = 0+0im;
    dsdbes = (n/2)^(alpha+1.0)*(-2.0)/gamma(2.0+0im+alpha)/2;
    for nn = 1:mo
	dsbes = dsbes + (-1.0)^nn*(n/2)^(alpha+2.0*nn)*theta^(2*nn-2.0)*2.0*nn/gamma(1.0+nn+alpha)/gamma(1.0+nn);
	dsden = dsden + (-1.0)^nn*theta^(2*nn-2.0)*2.0*nn/gamma(3.0+2*nn);
	dsthosth = dsthosth + (-1.0)^nn*theta^(2*nn-2.0)*2.0*nn/gamma(2.0+2*nn);
	dssinpl = dssinpl + (-1.0)^nn*(psi(0.0im+z)+alpha*pi/2+theta/2).^(2*nn-2.0)*2.0*nn/gamma(2.0+2*nn)*plcoef^2;
	dssinmi = dssinmi + (-1.0)^nn*(psi(0.0im+z)+alpha*pi/2-theta/2).^(2*nn-2.0)*2.0*nn/gamma(2.0+2*nn)*micoef^2;
	dsdbes = dsdbes + (-1.0)^nn/gamma(1.0+nn)*(n/2.0)^(2*nn+alpha+1.0)*(4.0/n^2.0*theta^(2*nn-2.0)*2*nn/gamma(nn+alpha+0.0)-theta^(2.0*nn)*(2*nn+2.0)/gamma(2.0+nn+alpha))/2;
    end
    dsbesowei = dsbes*(sden)^(-alpha/2+0.0) -alpha/2*sbes*sden^(-alpha/2-1.0)*dsden;
    dscosbesoweipl = cos(psi(0.0im+z) +alpha*pi/2 +theta/2)*dsbesowei - sbesowei*ssinpl*plcoef^2;
    dscosbesoweimi = cos(psi(0.0im+z) +alpha*pi/2 -theta/2)*dsbesowei - sbesowei*ssinmi*micoef^2;
	
    common1 = -alpha/2*sdbes*(sden)^(-alpha/2-1.0)*dsden + (sden)^(-alpha/2+0.0)*dsdbes;
    common2 = (dsthosth*contpsi(z+0im,1.0)+sthosth^2*contpsi(z+0im,2.0) )/(4im*pi)*sdbes*(sden)^(-alpha/2+0.0);
	
    dssindbesoweipl = plcoef*(common1*ssinpl + sdbes*(sden)^(-alpha/2+0.0)*dssinpl) + common2*ssinpl;
    dssindbesoweimi = micoef*(common1*ssinmi + sdbes*(sden)^(-alpha/2+0.0)*dssinmi) + common2*ssinmi;
	
    for k = (nrT-1):-1:1
	tmp = (1.0+1.0im)*zeros(2,2);
	for i = 2:size(Qright,4) #Start from i=2 to avoid constant term
	    tmp = tmp + Qright[:,:,k,i]*(-theta^2*sden)^(i-2.0)*(i-1.0)*(-sthosth);
	end
	dRR = (dRR + tmp)/n;
    end
else
    function dw(x) return -alpha*(1-x).^(alpha-1).*(1+x).^beta.*h(x)  +(1-x).^alpha.*beta*(1+x).^(beta-1).*h(x)+(1-x).^alpha.*(1+x).^beta.*dh(x); end
    function dSqrtPhi(x) 1/2./sqrt(0im+x + sqrt(0im+x-1.0).*sqrt(0im+x+1.0) ).*(1+x./sqrt(0im+x-1.0)./sqrt(0im+x+1.0) ); end
	
    function dDeltaR(k,x) DeltaR(k,x).*(log(x + sqrt(0im+x-1.0).*sqrt(0im+x+1.0)+0im).^k).*2*sqrt(0im+x+1.0)*sqrt(0im+x-1.0).*( (-k)*log(x + sqrt(0im+x-1.0).*sqrt(0im+x+1.0)+0im).^(-k-1.0)./(x + sqrt(0im+x-1.0).*sqrt(0im+x+1.0)).*(1 +x./sqrt(0im+x-1.0)./sqrt(0im+x+1.0) )./(2*sqrt(0im+x+1.0)*sqrt(0im+x-1.0))	+1./(log(x + sqrt(0im+x-1.0).*sqrt(0im+x+1.0)+0im).^k).*(-x/2).*(x+1.0+0im)^(-3/2+0.0).*(x-1.0+0im)^(-3/2+0.0) )    +     # End of derivative of term in front		
		brac(k-1)/2.0^k./(log(x + sqrt(0im+x-1.0).*sqrt(0im+x+1.0)+0im).^k)./(2*sqrt(0im+x+1.0)*sqrt(0im+x-1.0))*[Dinf 0; 0 Dinf^(-1+0.0)]*( [dSqrtPhi(x)    -1im*dSqrtPhi(x)./(x + sqrt(0im+x-1.0).*sqrt(0im+x+1.0) )   ;    1im*dSqrtPhi(x)./(x + sqrt(0im+x-1.0).*sqrt(0im+x+1.0) )    dSqrtPhi(x)]*[exp( 1im*(-1.0+0im)^(angle(x-1) <= 0)*(psi(x) +alpha*pi/2) )   0  ;   0   exp(-1im*(-1.0+0im)^(angle(x-1) <= 0)*(psi(x) +alpha*pi/2))]  +  [sqrt(0im+x + sqrt(0im+x-1.0).*sqrt(0im+x+1.0) )     1im/sqrt(0im+x + sqrt(0im+x-1.0).*sqrt(0im+x+1.0) )   ;   -1im/sqrt(0im+x + sqrt(0im+x-1.0).*sqrt(0im+x+1.0) )     sqrt(0im+x + sqrt(0im+x-1.0).*sqrt(0im+x+1.0) ) ]*[exp( 1im*(-1.0+0im)^(angle(x-1) <= 0)*(psi(x) +alpha*pi/2) )*1im*(-1.0+0im)^(angle(x-1) <= 0)*dpsi(x)    0   ;   0   exp(-1im*(-1.0+0im)^(angle(x-1) <= 0)*(psi(x) +alpha*pi/2))*(-1im)*(-1.0+0im)^(angle(x-1) <= 0)*dpsi(x)])*	[ ((-1)^k)/k*(alpha^2+k/2-1/4)    -1im*(k-1/2)    ;   (-1.0)^k*1im*(k-1/2)     (alpha^2+k/2-1/4)/k]*	[exp(-1im*(-1.0+0im)^(angle(x-1) <= 0)*(psi(x) +alpha*pi/2))   0  ;   0   exp(1im*(-1.0+0im)^(angle(x-1) <= 0)*(psi(x) +alpha*pi/2) )]*[sqrt(0im+x + sqrt(0im+x-1.0).*sqrt(0im+x+1.0) )     -1im/sqrt(0im+x + sqrt(0im+x-1.0).*sqrt(0im+x+1.0) )   ;   1im/sqrt(0im+x + sqrt(0im+x-1.0).*sqrt(0im+x+1.0) )   sqrt(0im+x + sqrt(0im+x-1.0).*sqrt(0im+x+1.0) ) ]*[Dinf^(-1+0.0) 0; 0 Dinf]   +  # End of derivative of Lambda, now derivative of Lambda^(-1)
		brac(k-1)/2.0^k./(log(x + sqrt(0im+x-1.0).*sqrt(0im+x+1.0)+0im).^k)./(2*sqrt(0im+x+1.0)*sqrt(0im+x-1.0))*[Dinf 0; 0 Dinf^(-1+0.0)]*	[sqrt(0im+x + sqrt(0im+x-1.0).*sqrt(0im+x+1.0) )     1im/sqrt(0im+x + sqrt(0im+x-1.0).*sqrt(0im+x+1.0) )  ;   -1im/sqrt(0im+x + sqrt(0im+x-1.0).*sqrt(0im+x+1.0) )    sqrt(0im+x + sqrt(0im+x-1.0).*sqrt(0im+x+1.0) ) ]*	[exp( 1im*(-1.0+0im)^(angle(x-1) <= 0)*(psi(x) +alpha*pi/2) )   0  ;  0   exp(-1im*(-1.0+0im)^(angle(x-1) <= 0)*(psi(x) +alpha*pi/2))]*[ ((-1.0)^k)/k*(alpha^2+k/2-1/4)    -1im*(k-1/2)   ;   (-1.0)^k*1im*(k-1/2)   (alpha^2+k/2-1/4)/k]*	(		[-1im*(-1.0+0im)^(angle(x-1) <= 0)*dpsi(x)*exp(-1im*(-1.0+0im)^(angle(x-1) <= 0)*(psi(x) +alpha*pi/2))   0  ;   0   1im*(-1.0+0im)^(angle(x-1) <= 0)*dpsi(x)*exp(1im*(-1.0+0im)^(angle(x-1) <= 0)*(psi(x) +alpha*pi/2) )]*		[sqrt(0im+x + sqrt(0im+x-1.0).*sqrt(0im+x+1.0) )    -1im/sqrt(0im+x + sqrt(0im+x-1.0).*sqrt(0im+x+1.0) )   ;    	1im/sqrt(0im+x + sqrt(0im+x-1.0).*sqrt(0im+x+1.0) )    sqrt(0im+x + sqrt(0im+x-1.0).*sqrt(0im+x+1.0) ) ] + 	[exp(-1im*(-1.0+0im)^(angle(x-1) <= 0)*(psi(x) +alpha*pi/2)) 0  ;   0   exp(1im*(-1.0+0im)^(angle(x-1) <= 0)*(psi(x) +alpha*pi/2) )]*	[dSqrtPhi(x)    1im*dSqrtPhi(x)./(x + sqrt(0im+x-1.0).*sqrt(0im+x+1.0) )   ;    -1im*dSqrtPhi(x)./(x + sqrt(0im+x-1.0).*sqrt(0im+x+1.0) )    dSqrtPhi(x)])*[Dinf^(-1+0.0) 0; 0 Dinf]; end
	
    ds = (1.0+1.0im)*zeros(2,2,nrT-1);
    dR = (1.0+1.0im)*zeros(2,2,nrT-1);
	
    for k = 1:nrT-1
	for m = 1:ceil(Integer, k/2)
	    dR[:,:,k] = dR[:,:,k] + Uright[:,:,k,m]./(z-1).^(m+1)*(-m) + Uleft[:,:,k,m]./(z+1).^(m+1)*(-m);
	end
		
	ds[:,:,k] = dDeltaR(k,z);
	if mod(k,2)==0 # k is even: add extra matrix
	    ds[:,:,k] = ds[:,:,k] -brac(k-1)*(4*alpha^2+2*k-1)/2.0^(k+1)/k*(-k)*log(z + sqrt(0im+z-1.0).*sqrt(0im+z+1.0)+0im).^(-k-1.0)./(z + sqrt(0im+z-1.0).*sqrt(0im+z+1.0)).*(1.0 +z./sqrt(0im+z-1.0)./sqrt(0im+z+1.0) )*eye(2);
	end
	dRR = dRR + (dR[:,:,k]-ds[:,:,k] )/n^k; # Avoids defining R_0^{right} = I
	for j = 1:(k-1)
	    dRR = dRR - (dR[:,:,k-j]*s[:,:,j]+R[:,:,k-j]*ds[:,:,j])/n^k;
	end
    end
    db1 = -sin(psi(z) +alpha*pi/2 + acos(z)/2)*(dpsi(z) -1/2/sqrt(0im+1.0-z)/sqrt(0im+1.0+z))*besselj(alpha,n*acos(z) ) +cos(psi(z) +alpha*pi/2 + acos(z)/2)*(besselj(alpha-1,n*acos(z) )    -besselj(alpha+1,n*acos(z) ) )/2*n*(-1.0)/sqrt(0im+1.0-z)/sqrt(0im+1.0+z) + # Above derivative of cos-term, below sine term
        cos(psi(z) +alpha*pi/2 + acos(z)/2)*(dpsi(z) -1/2/sqrt(0im+1.0-z)/sqrt(0im+1.0+z))*   (besselj(alpha-1,n*acos(z) ) - besselj(alpha+1,n*acos(z) ) )/2  +sin(psi(z) +alpha*pi/2 + acos(z)/2)*(besselj(alpha-2,n*acos(z) ) -2*besselj(alpha,n*acos(z) ) + besselj(alpha+2,n*acos(z) ) )/4*   n*(-1)/sqrt(0im+1.0-z)/sqrt(0im+1.0+z);
    db2 = -sin(psi(z) +alpha*pi/2 -acos(z)/2)* (dpsi(z) +1/2/sqrt(0im+1.0-z)/sqrt(0im+1.0+z) )*besselj(alpha,n*acos(z) )+cos(psi(z) +alpha*pi/2 -acos(z)/2)* (besselj(alpha-1,n*acos(z) ) - besselj(alpha+1,n*acos(z) ) )/2.0*(-n)./sqrt(0im+1.0-z)./sqrt(0im+1.0+z)  + cos(psi(z) +alpha*pi/2 -acos(z)/2)*(dpsi(z) +1/2/sqrt(0im+1.0-z)/sqrt(0im+1.0+z) )* (besselj(alpha-1,n*acos(z) ) - besselj(alpha+1,n*acos(z) ) )/2  + sin(psi(z) +alpha*pi/2 -acos(z)/2)*(besselj(alpha-2,n*acos(z) ) -2*besselj(alpha,n*acos(z) ) + besselj(alpha+2,n*acos(z) ) )/4*n*(-1)/sqrt(0im+1.0-z)/sqrt(0im+1.0+z) ;
end

if method == "m"
    dbpiR = (([sqrt(0im+pi*n)/2.0^n*(1/2/sqrt(0im+acos(z) )*(-1.0)./sqrt(0im+1.0-z)./sqrt(0im+1.0+z)./(sqrt(0im+w(z) ).*(1.0+0.0im+z).^(1/4).*(1.0+0.0im-z).^(1/4) ) + sqrt(0im+acos(z) ).*(dw(z)*(-1/2)./(sqrt(0im+w(z) ).^3.*(1.0+0.0im+z).^(1/4).*(1.0+0.0im-z).^(1/4) ) -1/4*(-2*z)./(sqrt(0im+w(z) ).*(1.0+0.0im+z).^(5/4).*(1.0+0.0im-z).^(5/4) ) ) )     0;]*RR + [sqrt(0im+pi*n*acos(z))./(2.0^n.*sqrt(0im+w(z) ).*(1.0+0.0im+z).^(1/4).*(1.0+0.0im-z).^(1/4) )     0;]*dRR        )*[Dinf*b1  ;   -1im/Dinf*b2]  + # Above derivative of row vector and matrix, below of column vector
        [sqrt(0im+pi*n*acos(z))./(2.0^n.*sqrt(0im+w(z) ).*(1.0+0.0im+z).^(1/4).*(1.0+0.0im-z).^(1/4) )    0;]*RR*[Dinf*db1;-1im/Dinf*db2])[1,1];
elseif method == "o"
    dbpiR = (( ([sqrt(0im+n)*(1/2/sqrt(0im+acos(z) )*(-1.0)./sqrt(0im+1.0-z)./sqrt(0im+1.0+z)./(sqrt(0im+w(z) ).*(1.0+0.0im+z).^(1/4).*(1.0+0.0im-z).^(1/4) ) + sqrt(0im+acos(z) ).*(dw(z)*(-1/2)./(sqrt(0im+w(z) ).^3.*(1.0+0.0im+z).^(1/4).*(1.0+0.0im-z).^(1/4) ) -1/4*(-2*z)./(sqrt(0im+w(z) ).*(1.0+0.0im+z).^(5/4).*(1.0+0.0im-z).^(5/4) ) ) )     0;]*RR + [sqrt(0im+n*acos(z))./(sqrt(0im+w(z) ).*(1.0+0.0im+z).^(1/4).*(1.0+0.0im-z).^(1/4) )     0;]*dRR)*[b1; -1im/Dinf^2*b2]  + #Above derivative of row vector and matrix, below of column vector
        [sqrt(0im+n*acos(z))./(sqrt(0im+w(z) ).*(1.0+0.0im+z).^(1/4).*(1.0+0.0im-z).^(1/4) )    0;]*RR*[db1; -1im/Dinf^2*db2] )*nmz )[1,1];
elseif method == "mQ"
    dbpiR = (([sqrt(0im+pi*n)/2.0^n*(-1/2*sthosth^(-3/2+0.0)*dsthosth./(sqrt((0.0im+1+z).^beta.*h(0.0im+z) ) ) -  (dh(0.0im+z)*(-1/2)./(sqrt(h(0.0im+z).^3.*(0.0im+1+z).^beta) ) -beta/2./(sqrt(h(0.0im+z).*(0.0im+1+z).^beta).*(0.0im+1+z) ) ) )     0;]*RR + [sqrt(0im+pi*n/sthosth)./(2.0.^n.*sqrt(h(0im+z).*(0im+1+z).^beta) )    0;]*dRR)*[Dinf*(cos(psi(0im+z) +alpha*pi/2 + theta/2)*sbesowei + splsindbes ); -1im/Dinf*(cos(psi(0im+z) +alpha*pi/2 -theta/2)* sbesowei + smisindbes)]  +
#Above derivative of row vector and matrix, below of column vector
        [sqrt(0im+pi*n/sthosth)./(2.0.^n.*sqrt(h(0im+z).*(0im+1+z).^beta) )      0;]*RR*[Dinf*(dscosbesoweipl +dssindbesoweipl);-1im/Dinf*(dscosbesoweimi+dssindbesoweimi)])[1,1];
	dbpiR = -dbpiR/sthosth;
elseif method == "oQ"
    dbpiR = (([sqrt(0im+n)*(-1/2*sthosth^(-3/2+0.0)*dsthosth./(sqrt((0.0im+1+z).^beta.*h(0.0im+z) ) ) -  (dh(0.0im+z)*(-1/2)./(sqrt(h(0.0im+z).^3.*(0.0im+1+z).^beta) ) -beta/2./(sqrt(h(0.0im+z).*(0.0im+1+z).^beta).*(0.0im+1+z) ) ) )     0;]*RR + [sqrt(0im+n/sthosth)./(sqrt(h(0im+z).*(0im+1+z).^beta) )    0;]*dRR)*[(cos(psi(0im+z) +alpha*pi/2 + theta/2)*sbesowei + splsindbes ); -1im/Dinf^2*(cos(psi(0im+z) +alpha*pi/2 -theta/2)* sbesowei + smisindbes)]  +
#Above derivative of row vector and matrix, below of column vector
        [sqrt(0im+n/sthosth)./(sqrt(h(0im+z).*(0im+1+z).^beta) )    0;]*RR*[(dscosbesoweipl +dssindbesoweipl);-1im/Dinf^2*(dscosbesoweimi+dssindbesoweimi)])[1,1];
	dbpiR = -dbpiR/sthosth*nmz;
end
return (bpiR, dbpiR)
end


# Evaluate the asymptotic expansion for the generalised Jacobi polynomial in the left boundary region with a certain method. Derivative is only computed when dw is specified, default is monic normalisation without Q-s or derivatives. When using the Q-s, a series of the cosine transform of the asymptotic expansion is also used and theta=arccos(-z).
# Input
#   n,z            - Degree and point at which to evaluate
#   alpha,beta,h   - Parts of the weight 
#   psi(x)         - Function for analytic phase function psi(x)
#   nrT            - Number of terms. 1 gives leading order term
#   Dinf           - Limit of the Szego function
#   Uright,Uleft   - Right and left U-matrices
#   [method        - 'm' for monic(default), 'o' for orthonormal, 'mQ' for monic using the Q-s, 'oQ' for orthonormal with Q-s]
#   [dh,dpsi       - Derivatives of part of weight and phase function]
#   [Qleft         - Left Q-matrices for when regexp(method,'Q') == 1]
#   [theta	   - theta=arccos(-z) when computing with Q-s to save digits]
#   [contpsi	   - Contour integral in psi, needed for series expansions]
# Output
#   bpiL           - Asymptotic expansion of the left boundary polynomial
#   [dbpiL         - Derivative of left boundary polynomial]
function asy_left(n,z, alpha,beta,h,psi,nrT,Dinf,Uright,Uleft, method="m",dh=false, dpsi=false,Qleft=false,theta=false,contpsi=false)

if ismatch(r"o",method)
    if (nrT == 1) || (length(Uright) == 0) nmz = 1; else
	nmz = sqrt(real(1+2*1im*Dinf^2*sum(reshape(Uright[2,1,1:(nrT-1),1] + Uleft[2,1,1:(nrT-1),1], (nrT-1), 1)./(n+1).^(1:(nrT-1) ) ) ) );
    end
end


if ismatch(r"Q",method)
    if (abs(z+1) >1e-2) || (abs(theta) > 1e-2)
	warn("Theta is not small so using Q-s will probably give an inaccurate result")
    end
    mo = size(Qleft,4)-1;
    sbes = 0im+0;	sden = 0im+0;	sthosth = 0im+0;	ssinpl = 0im+0;	ssinmi = 0im+0;	sdbes = 0im+0;
    for nn = 0:mo
	sbes = sbes + (-1)^nn*(n/2)^(beta+2*nn)*theta^(2*nn)/gamma(1+nn+beta)/gamma(1.0+nn);
	sden = sden + (-1)^nn*theta^(2*nn)/gamma(3+2*nn);
	sthosth = sthosth + (-1)^nn*theta^(2*nn)/gamma(2.0+2*nn);
				
	ssinpl = ssinpl + (-1)^(nn+1)*(psi(0im+z) -beta*pi/2 -theta/2).^(2*nn)/gamma(2.0+2*nn);
	ssinmi = ssinmi + (-1)^nn*(psi(0im+z) -beta*pi/2 +theta/2).^(2*nn)/gamma(2.0+2*nn);
	sdbes = sdbes + (-1)^nn/gamma(1.0+nn)*(n/2)^(2*nn+beta+1)*(4/n^2*theta^(2*nn)/gamma(nn+beta)-theta^(2*nn+2)/gamma(2+nn+beta))/2;
    end
    plcoef = (-alpha-beta-1)/2 +sthosth/(4im*pi)*contpsi(0im+z,1);
    micoef = (-alpha-beta+1)/2 +sthosth/(4im*pi)*contpsi(0im+z,1);
    sbesowei = sbes*(sden)^(-beta/2+0.0);
    splsindbes = plcoef*ssinpl*sdbes*sden^(-beta/2+0.0);
    smisindbes = micoef*ssinmi*sdbes*sden^(-beta/2+0.0);
    RL = (1+0im)*zeros(2,2);
    for k = (nrT-1):-1:1 # Counting down to avoid roundoff error
	tmp = zeros(2,2);
	for i = 1:size(Qleft,4)
	    tmp = tmp + Qleft[:,:,k,i]*(-theta^2*sden)^(i-1);		
	end
	RL = (RL + tmp)/n;
    end
    RL = RL + eye(2);
else
    function w(x) return (1-x).^alpha.*(1+x).^beta.*h(x); end
    if abs(z+1) < eps()^(1/3)
	warn("z is close to -1 so use Q-s and series expansions in 'method'")
    end
    function brac(k) convert(FloatingPoint, (k == 0)) + convert(FloatingPoint, (k != 0))*prod(4*beta^2 -(2*(1:k)-1).^2)/(2.0^(2*k)*gamma(1.0+k) ); end # = (beta,k)
	
    function DeltaL(k,x) brac(k-1)/2.0^k./(log(-x -sqrt(0im+x-1.0).*sqrt(0im+x+1.0)+0im).^k)./(2*sqrt(0im+x+1.0)*sqrt(0im+x-1.0))*[Dinf 0; 0 (Dinf+0.0)^(-1+0.0)]*[sqrt(0im+x + sqrt(0im+x-1.0).*sqrt(0im+x+1.0) )      1im/sqrt(0im+x + sqrt(0im+x-1.0).*sqrt(0im+x+1.0) ); -1im/sqrt(0im+x + sqrt(0im+x-1.0).*sqrt(0im+x+1.0) )      sqrt(0im+x + sqrt(0im+x-1.0).*sqrt(0im+x+1.0) ) ]*	[exp( 1im*(-1.0+0im)^(angle(x-1) <= 0)*(psi(x) -beta*pi/2) ) 0; 0 exp(-1im*(-1.0+0im)^(angle(x-1) <= 0)*(psi(x) -beta*pi/2))]*[ ((-1.0)^k)/k*(beta^2+k/2-1/4)     1im*(k-1/2)    ;    (-1.0+0im)^(k+1)*1im*(k-1/2)       (beta^2+k/2-1/4)/k]*	[exp(-1im*(-1.0+0im)^(angle(x-1) <= 0)*(psi(x) -beta*pi/2)) 0; 0 exp(1im*(-1.0+0im)^(angle(x-1) <= 0)*(psi(x) -beta*pi/2) )]*[sqrt(0im+x + sqrt(0im+x-1.0).*sqrt(0im+x+1.0) )     -1im/sqrt(0im+x + sqrt(0im+x-1.0).*sqrt(0im+x+1.0) ); 1im/sqrt(0im+x + sqrt(0im+x-1.0).*sqrt(0im+x+1.0) )     sqrt(0im+x + sqrt(0im+x-1.0).*sqrt(0im+x+1.0) ) ]*[(Dinf+0.0)^(-1+0.0) 0; 0 Dinf]; end
	
    s = (1.0+1.0im)*zeros(2,2,nrT-1);
    R = (1.0+1.0im)*zeros(2,2,nrT-1);
    RL = (1+0im)*eye(2);
    for k = 1:nrT-1
	s[:,:,k] = DeltaL(k,z);
	for m = 1:ceil(Integer, k/2)
	    R[:,:,k] = R[:,:,k] + Uright[:,:,k,m]./(z-1).^m + Uleft[:,:,k,m]./(z+1).^m;
	end
	if mod(k,2)==0 # k is even: add extra matrix
	    s[:,:,k] = s[:,:,k] -brac(k-1)*(4*beta^2+2*k-1)/2.0^(k+1)/k/(log(-z-sqrt(0im+z-1.0).*sqrt(0im+z+1.0) +0im) )^k*eye(2);
	end
	RL = RL + (R[:,:,k]-s[:,:,k] )/n^k;
	for j = 1:(k-1)
	    RL = RL - R[:,:,k-j]*s[:,:,j]/n^k;
	end
    end
    b1 = sin(psi(z) -beta*pi/2 + acos(z)/2)*besselj(beta,n*acos(-z) ) +   cos(psi(z) -beta*pi/2 + acos(z)/2)*(besselj(beta-1,n*acos(-z) ) - besselj(beta+1,n*acos(-z) ) )/2;
    b2 = sin(psi(z) -beta*pi/2 -acos(z)/2)*  besselj(beta,n*acos(-z) ) + cos(psi(z) -beta*pi/2 -acos(z)/2)*    (besselj(beta-1,n*acos(-z) ) - besselj(beta+1,n*acos(-z) ) )/2;
end

if method == "m"
    bpiL = ([sqrt(0im+pi*n*acos(-z))./((-2.0)^n.*sqrt(0im+w(z)).*(1.0+z+0.0im).^(1/4).*(1.0-z+0.0im).^(1/4) )      0;]*RL*[Dinf*b1; -1im/Dinf*b2])[1,1]; # *(-1.0+0.0im)^(abs(angle(z+1))>pi/2);
elseif method == "o"
    bpiL = ([sqrt(0im+n*acos(-z))./((-1.0)^n.*sqrt(0im+w(z)).*(1.0+z+0.0im).^(1/4).*(1.0-z+0.0im).^(1/4) )    0;]*RL*[b1; -1im/Dinf^2*b2]*nmz)[1,1]; # *(-1.0+0.0im)^(abs(angle(z+1))>pi/2);
elseif method == "mQ"
    bpiL = ([sqrt(0im+pi*n/sthosth)./((-2.0).^n*sqrt((0im+1-z).^alpha.*h(0im+z) ) )     0;]*RL*[Dinf*(sin(psi(0im+z) -beta*pi/2 + (pi-theta)/2)*sbesowei +splsindbes); -1im/Dinf*(sin(psi(z) -beta*pi/2 -(pi-theta)/2)*sbesowei + smisindbes)])[1,1]; # *(-1.0+0.0im)^(abs(angle(z+1))>pi/2);
elseif method == "oQ"
    bpiL = ([sqrt(0im+n/sthosth)./((-1).^n*sqrt((0im+1-z).^alpha.*h(0im+z) ) )     0;]*RL*[(sin(psi(0im+z) -beta*pi/2 + (pi-theta)/2)*sbesowei +splsindbes); -1im/Dinf^2*(sin(psi(z) -beta*pi/2 -(pi-theta)/2)*sbesowei + smisindbes)]*nmz)[1,1]; # *(-1.0+0.0im)^(abs(angle(z+1))>pi/2);
else
    error("Wrong method")
end

if dh == false return bpiL end

dRL = (1.0+1im)*zeros(2,2);
if ismatch(r"Q",method)
    dsbes = 0im+0;	dsden = 0im+0;	dsthosth = 0im+0;	dssinpl = 0im+0;	dssinmi = 0im+0;
    dsdbes = (n/2)^(beta+1)*(-2)/gamma(2.0+beta)/2;
    for nn = 1:mo
	dsbes = dsbes + (-1)^nn*(n/2)^(beta+2*nn)*theta^(2*nn-2)*2*nn/gamma(1.0+nn+beta)/gamma(1.0+nn);
	dsden = dsden + (-1.0)^nn*theta^(2.0*nn-2.0)*2*nn/gamma(3.0+2*nn);
	dsthosth = dsthosth + (-1)^nn*theta^(2.0*nn-2.0)*2.0*nn/gamma(2.0+2*nn);	
	dssinpl = dssinpl + (-1)^(nn+1.0)*(psi(0im+z)-beta*pi/2 -theta/2).^(2*nn-2)*2*nn/gamma(2.0+2*nn)*plcoef^2;
	dssinmi = dssinmi + (-1.0)^nn*(psi(0im+z)-beta*pi/2 +theta/2).^(2*nn-2)*2*nn/gamma(2.0+2*nn)*micoef^2;
	dsdbes = dsdbes + (-1.0)^nn/gamma(1.0+nn)*(n/2.0)^(2*nn+beta+1.0)*(4.0/n^2.0*theta^(2*nn-2.0)*2*nn/gamma(nn+beta+0.0)-theta^(2.0*nn)*(2*nn+2.0)/gamma(2.0+nn+beta))/2;
    end
    dsbesowei = dsbes*(sden)^(-beta/2+0.0) -beta/2*sbes*sden^(-beta/2-1.0)*dsden;
    dscosbesoweipl = sin(psi(0im+z) -beta*pi/2 +(pi-theta)/2)*dsbesowei + sbesowei*ssinpl*plcoef^2;
    dscosbesoweimi = sin(psi(0im+z) -beta*pi/2 -(pi-theta)/2)*dsbesowei + sbesowei*ssinmi*micoef^2;

    common1 = -beta/2*sdbes*(sden)^(-beta/2-1.0)*dsden + (sden)^(-beta/2+0.0)*dsdbes;
    common2 = (dsthosth*contpsi(0im+z,1.0)+sthosth^2*contpsi(0im+z,2.0) )/(4im*pi)*sdbes*(sden)^(-beta/2+0.0);
	
    dssindbesoweipl = plcoef*(common1*ssinpl + sdbes*(sden)^(-beta/2+0.0)*dssinpl) + common2*ssinpl;
    dssindbesoweimi = micoef*(common1*ssinmi + sdbes*(sden)^(-beta/2+0.0)*dssinmi) + common2*ssinmi;	
    for k = (nrT-1):-1:1
	tmp = zeros(2,2);
	for i = 2:size(Qleft,4)
	    tmp = tmp + Qleft[:,:,k,i]*(theta^2*sden)^(i-2.0)*(i-1.0)*sthosth;
	end
	dRL = (dRL + tmp)/n;
    end
else
    function dw(x) return -alpha*(1-x).^(alpha-1).*(1+x).^beta.*h(x)  +(1-x).^alpha.*beta*(1+x).^(beta-1).*h(x)+(1-x).^alpha.*(1+x).^beta.*dh(x);end
    function dSqrtPhi(x) 1/2./sqrt(0im+x + sqrt(0im+x-1).*sqrt(0im+x+1) ).*(1+x./sqrt(0im+x-1)./sqrt(0im+x+1) ); end
	
    function dDeltaL(k,x) DeltaL(k,x).*(log(-x -sqrt(0im+x-1).*sqrt(0im+x+1)+0im).^k).*2*sqrt(0im+x+1)*sqrt(0im+x-1).*(  (-k)*log(-x -sqrt(0im+x-1).*sqrt(0im+x+1)+0im).^(-k-1.0)./(-x -sqrt(0im+x-1).*sqrt(0im+x+1)).*(-1 -x./sqrt(0im+x-1)./sqrt(0im+x+1) )./(2*sqrt(0im+x+1)*sqrt(0im+x-1))+1./(log(-x -sqrt(0im+x-1).*sqrt(0im+x+1)+0im).^k).*(-x/2).*(x+1+0im)^(-3/2+0.0).*(x-1+0im)^(-3/2+0.0) )   +   # End of derivative of term in front
		brac(k-1)/2.0^k./(log(-x -sqrt(0im+x-1).*sqrt(0im+x+1)+0im).^k)./(2*sqrt(0im+x+1)*sqrt(0im+x-1))*[Dinf 0; 0 Dinf^(-1+0.0)]*( [dSqrtPhi(x)    -1im*dSqrtPhi(x)./(x + sqrt(0im+x-1).*sqrt(0im+x+1) );  1im*dSqrtPhi(x)./(x + sqrt(0im+x-1).*sqrt(0im+x+1) )     dSqrtPhi(x)]*	[exp( 1im*(-1+0im)^(angle(x-1) <= 0)*(psi(x) -beta*pi/2) ) 0; 0 exp(-1im*(-1+0im)^(angle(x-1) <= 0)*(psi(x) -beta*pi/2))]	+[sqrt(0im+x + sqrt(0im+x-1).*sqrt(0im+x+1) )     1im/sqrt(0im+x + sqrt(0im+x-1).*sqrt(0im+x+1) ); 	-1im/sqrt(0im+x + sqrt(0im+x-1).*sqrt(0im+x+1) )     sqrt(0im+x + sqrt(0im+x-1).*sqrt(0im+x+1) ) ]*	[exp( 1im*(-1+0im)^(angle(x-1) <= 0)*(psi(x) -beta*pi/2) )*1im*(-1+0im)^(angle(x-1) <= 0)*dpsi(x) 0; 0 exp(-1im*(-1+0im)^(angle(x-1) <= 0)*(psi(x) -beta*pi/2))*(-1im)*(-1+0im)^(angle(x-1) <= 0)*dpsi(x)] )*[ ((-1)^k)/k*(beta^2+k/2-1/4)       1im*(k-1/2) ; (-1+0im)^(k+1)*1im*(k-1/2)        (beta^2+k/2-1/4)/k]* [exp(-1im*(-1+0im)^(angle(x-1) <= 0)*(psi(x) -beta*pi/2)) 0; 0 exp(1im*(-1+0im)^(angle(x-1) <= 0)*(psi(x) -beta*pi/2) )]*[sqrt(0im+x + sqrt(0im+x-1).*sqrt(0im+x+1) )       -1im/sqrt(0im+x + sqrt(0im+x-1).*sqrt(0im+x+1) ); 		1im/sqrt(0im+x + sqrt(0im+x-1).*sqrt(0im+x+1) )      sqrt(0im+x + sqrt(0im+x-1).*sqrt(0im+x+1) ) ]*[Dinf^(-1.0) 0; 0 Dinf]	+   # End of derivative of Lambda, now derivative of Lambda^(-1)
		brac(k-1)/2.0^k./(log(-x -sqrt(0im+x-1).*sqrt(0im+x+1)+0im).^k)./(2*sqrt(0im+x+1)*sqrt(0im+x-1))*[Dinf 0; 0 Dinf^(-1.0)]*	[sqrt(0im+x + sqrt(0im+x-1).*sqrt(0im+x+1) )       1im/sqrt(0im+x + sqrt(0im+x-1).*sqrt(0im+x+1) ); 	-1im/sqrt(0im+x + sqrt(0im+x-1).*sqrt(0im+x+1) )      sqrt(0im+x + sqrt(0im+x-1).*sqrt(0im+x+1) ) ]*[exp( 1im*(-1+0im)^(angle(x-1) <= 0)*(psi(x) -beta*pi/2) ) 0; 0 exp(-1im*(-1+0im)^(angle(x-1) <= 0)*(psi(x) -beta*pi/2))]*	[ ((-1)^k)/k*(beta^2+k/2-1/4)       1im*(k-1/2) ; (-1)^(k+1)*1im*(k-1/2)        (beta^2+k/2-1/4)/k]*	(	[-1im*(-1+0im)^(angle(x-1) <= 0)*dpsi(x)*exp(-1im*(-1+0im)^(angle(x-1) <= 0)*(psi(x) -beta*pi/2)) 0; 	0 1im*(-1+0im)^(angle(x-1) <= 0)*dpsi(x)*exp(1im*(-1+0im)^(angle(x-1) <= 0)*(psi(x) -beta*pi/2) )]*	[sqrt(0im+x + sqrt(0im+x-1).*sqrt(0im+x+1) )       -1im/sqrt(0im+x + sqrt(0im+x-1).*sqrt(0im+x+1) ); 	1im/sqrt(0im+x + sqrt(0im+x-1).*sqrt(0im+x+1) )      sqrt(0im+x + sqrt(0im+x-1).*sqrt(0im+x+1) ) ] + 	[exp(-1im*(-1+0im)^(angle(x-1) <= 0)*(psi(x) -beta*pi/2)) 0; 0 exp(1im*(-1+0im)^(angle(x-1) <= 0)*(psi(x) -beta*pi/2) )]*	[dSqrtPhi(x)      1im*dSqrtPhi(x)./(x + sqrt(0im+x-1).*sqrt(0im+x+1) )    ;     -1im*dSqrtPhi(x)./(x + sqrt(0im+x-1).*sqrt(0im+x+1) )     dSqrtPhi(x)])*[Dinf^(-1.0) 0; 0 Dinf]; end
	
    ds = (1.0+1im)*zeros(2,2,nrT-1);
    dR = (1.0+1im)*zeros(2,2,nrT-1);
    for k = 1:nrT-1
	for m = 1:ceil(Integer, k/2)
	    dR[:,:,k] = dR[:,:,k] + Uright[:,:,k,m]./(z-1).^(m+1)*(-m) + Uleft[:,:,k,m]./(z+1).^(m+1)*(-m);
	end
	ds[:,:,k] = dDeltaL(k,z);
	if mod(k,2)==0 # k is even: add extra matrix
	    ds[:,:,k] = ds[:,:,k] -brac(k-1)*(4*beta^2+2*k-1)/2.0^(k+1)/k		*(-k)*log(-z -sqrt(0im+z-1).*sqrt(0im+z+1)+0im).^(-k-1.0)./(-z -sqrt(0im+z-1).*sqrt(0im+z+1))	.*(-1 -z./sqrt(0im+z-1)./sqrt(0im+z+1) )*eye(2);
	end
	dRL = dRL + (dR[:,:,k]-ds[:,:,k] )/n^k; # Avoids defining R_0^{right} = I
	for j = 1:(k-1)
	    dRL = dRL - (dR[:,:,k-j]*s[:,:,j]+R[:,:,k-j]*ds[:,:,j])/n^k;
	end
    end
    db1 = cos(psi(z) -beta*pi/2 + acos(z)/2)*(dpsi(z)-1/2/sqrt(0im+1-z)/sqrt(0im+1+z))*    besselj(beta,n*acos(-z) ) +    sin(psi(z) -beta*pi/2 + acos(z)/2)*(besselj(beta-1,n*acos(-z) )    -besselj(beta+1,n*acos(-z) ))/2*n*(+1)/sqrt(0im+1.0-z)/sqrt(0im+1.0+z)   +  # Above derivative of sine-term, below cos term
        -sin(psi(z) -beta*pi/2 + acos(z)/2)*(dpsi(z)-1/2/sqrt(0im+1.0-z)/sqrt(0im+1+z) )*   (besselj(beta-1,n*acos(-z) ) - besselj(beta+1,n*acos(-z) ) )/2  +cos(psi(z) -beta*pi/2 + acos(z)/2)*(besselj(beta-2,n*acos(-z) )  -2*besselj(beta,n*acos(-z) ) +besselj(beta+2,n*acos(-z) ) )/4*   n*(+1)/sqrt(0im+1.0-z)/sqrt(0im+1.0+z);
    db2 = cos(psi(z) -beta*pi/2 -acos(z)/2)*(dpsi(z)-1/2*(-1)/sqrt(0im+1.0-z)/sqrt(0im+1.0+z) )*    besselj(beta,n*acos(-z) )  + sin(psi(z) -beta*pi/2 -acos(z)/2)*(besselj(beta-1,n*acos(-z) )   -besselj(beta+1,n*acos(-z) ) )/2*n*(+1)/sqrt(0im+1.0-z)/sqrt(0im+1.0+z)   -sin(psi(z) -beta*pi/2 -acos(z)/2)*(dpsi(z)+1/2/sqrt(0im+1.0-z)/sqrt(0im+1.0+z))*   (besselj(beta-1,n*acos(-z) ) - besselj(beta+1,n*acos(-z) ) )/2 + cos(psi(z) -beta*pi/2 -acos(z)/2)*   (besselj(beta-2,n*acos(-z) ) -2*besselj(beta,n*acos(-z) )   +besselj(beta+2,n*acos(-z) ) )/4*n/sqrt(0im+1.0-z)/sqrt(0im+1.0+z);
end

if method == "m"
	dbpiL = (([sqrt(0im+pi*n)*(1/2./sqrt(0im+acos(-z) )./(sqrt(0im+1-z).*sqrt(0im+1+z).*(-2.0)^n.*sqrt(0im+w(z)).*	(1.0+z+0.0im).^(1/4).*(1.0-z+0.0im).^(1/4) ) 	+sqrt(0im+acos(-z))./(-2.0)^n*(dw(z)*(-1/2)./(sqrt(0im+w(z) ).^3.*(1.0+z+0.0im).^(1/4).*(1.0-z+0.0im).^(1/4) ) 	-1/4*(-2*z)./(sqrt(0im+w(z) ).*(1.0+z+0.0im).^(5/4).*(1.0-z+0.0im).^(5/4) ) ) )        0;]*RL +       [sqrt(0im+pi*n*acos(-z))./((-2.0)^n.*sqrt(0im+w(z)).*(1.0+z+0.0im).^(1/4).*(1.0-z+0.0im).^(1/4) )       0;]*dRL)*      [Dinf*b1; -1im/Dinf*b2]  +  # Above derivative of row vector and matrix, below of column vector
        [sqrt(0im+pi*n*acos(-z))./((-2.0)^n.*sqrt(0im+w(z)).*(1.0+z+0.0im).^(1/4).*(1.0-z+0.0im).^(1/4) )      0;]*RL*[Dinf*db1 ; -1im/Dinf*db2])[1,1]; # *(-1.0+0.0im)^(abs(angle(z+1))>pi/2);
elseif method == "o"
	dbpiL = ((([sqrt(0im+n)*(1/2./sqrt(0im+acos(-z) )./(sqrt(0im+1.0-z).*sqrt(0im+1.0+z).*(-1.0)^n.*sqrt(0im+w(z)).*	(1.0+z+0.0im).^(1/4).*(1.0-z+0.0im).^(1/4) ) +sqrt(0im+acos(-z))./(-1.0)^n*(dw(z)*(-1/2)./(sqrt(0im+w(z) ).^3.*(1.0+z+0.0im).^(1/4).*(1.0+0.0im-z).^(1/4) ) -1/4*(-2*z)./(sqrt(0im+w(z) ).*(1.0+z+0.0im).^(5/4).*(1.0-z+0.0im).^(5/4) ) ) )      0;]*RL +     [sqrt(0im+n*acos(-z))./((-1.0)^n.*sqrt(0im+w(z)).*(1.0+z+0.0im).^(1/4).*(1.0-z+0.0im).^(1/4) )       0;]*dRL)*    [b1; -1im/Dinf^2*b2]      +   # Above derivative of row vector and matrix, below of column vector
        [sqrt(0im+n*acos(-z))./((-1.0)^n.*sqrt(0im+w(z)).*(1.0+z+0.0im).^(1/4).*(1.0-z+0.0im).^(1/4) )      0;]*RL*[db1;  -1im/Dinf^2*db2])*nmz)[1,1]; # *(-1.0+0.0im)^(abs(angle(z+1))>pi/2);
elseif method == "mQ"
    dbpiL = (([sqrt(0im+pi*n).*(0im-2.0)^(-n+0.0)*(-1/2*(0.0im+sthosth)^(-3/2+0.0)*dsthosth./(sqrt(0im+(0im+1-z).^alpha.*h(0im+z) ) ) + (dh(0im+z)*(-1/2)./(sqrt(h(0im+z).^3.*(0im+1-z).^alpha) ) -alpha/2./(sqrt(h(0im+z).*(0im+1-z).^alpha).*(0im+1-z).*(-1) ) ) )     0;]*RL + [sqrt(0im+n*pi/sthosth).*(0im-2.0)^(-n+0.0)./(sqrt(h(0im+z) .*(0im+1-z).^alpha ) )     0;]*dRL )*[Dinf*(sin(psi(0im+z) -beta*pi/2 + (pi-theta)/2)*sbesowei +splsindbes ); -1im/Dinf*(sin(psi(0im+z) -beta*pi/2 -(pi-theta)/2)*sbesowei + smisindbes)] + [sqrt(0im+n*pi/sthosth).*(0im-2.0)^(-n+0.0)./(sqrt(h(0im+z).*(0im+1-z).^alpha) )     0;]*RL*[Dinf*(dscosbesoweipl +dssindbesoweipl);-1im/Dinf*(dscosbesoweimi+dssindbesoweimi)])[1,1];
    dbpiL = +dbpiL/sthosth; # *(-1.0+0.0im)^(abs(angle(z+1))>pi/2);
elseif method == "oQ"
    dbpiL = (([sqrt(0im+n).*(0.0im-1)^(-n+0.0)*(-1/2*(0.0im+sthosth)^(-3/2+0.0)*dsthosth./(sqrt(0im+(0im+1-z).^alpha.*h(0im+z) ) ) + (dh(0im+z)*(-1/2)./(sqrt(h(0im+z).^3.*(0im+1-z).^alpha) ) -alpha/2./(sqrt(h(0im+z).*(0im+1-z).^alpha).*(0im+1-z).*(-1) ) ) )     0;]*RL + [sqrt(0im+n/sthosth).*(-1)^(-n+0.0)./(sqrt(h(0im+z) .*(0im+1-z).^alpha ) )     0;]*dRL )*[(sin(psi(0im+z) -beta*pi/2 + (pi-theta)/2)*sbesowei +splsindbes ); -1im/Dinf^2*(sin(psi(0im+z) -beta*pi/2 -(pi-theta)/2)*sbesowei + smisindbes)] + [sqrt(0im+n/sthosth).*(-1)^(-n+0.0)./(sqrt(h(0im+z).*(0im+1-z).^alpha) )     0;]*RL*[(dscosbesoweipl +dssindbesoweipl);-1im/Dinf^2*(dscosbesoweimi+dssindbesoweimi)])[1,1];
    dbpiL = +dbpiL/sthosth*nmz; # *(-1.0+0.0im)^(abs(angle(z+1))>pi/2);
end

return (bpiL, dbpiL)
end


# Evaluate the asymptotic expansion for the gamma_n of the generalised Jacobi polynomial.
# Input
#   n            - Degree at which to evaluate
#   nrT          - Number of terms. 1 gives leading order term
#   Dinf         - Limit of the Szego function
#   Uright,Uleft - Right and left U-matrices
# Output
#   gamma        - Asymptotic expansion of the gamma_n
function gamman(n,nrT,Dinf,Uright,Uleft)
    if (nrT == 1) || (length(Uright) == 0)
	return 2.0.^n/Dinf/sqrt(pi);
    else
	return 2.0^(n)/Dinf/sqrt(pi)*sqrt(real(1+2*1im*Dinf^2*sum(reshape(Uright[2,1,1:(nrT-1),1] + Uleft[2,1,1:(nrT-1),1], (nrT-1), 1)./(n+1).^(1:(nrT-1) ) ) ));
    end
end

# Evaluate the asymptotic expansion for the recursion coefficient alpha_n of the generalised Jacobi polynomial.
# Input
#   n            - Degree at which to evaluate
#   nrT          - Number of terms. 1 gives leading order term
#   Uright,Uleft - Right and left U-matrices
# Output
#   alphan       - Asymptotic expansion of the alpha_n
function alphan(n,nrT,Uright,Uleft)
    if (nrT == 1) || (length(Uright) == 0)
	return 0;
    else
	return -sum(reshape(Uright[1,1,1:(nrT-1),1] + Uleft[1,1,1:(nrT-1),1], (nrT-1), 1)./(n+1).^(1:(nrT-1))) - sum(reshape(Uright[2,2,1:(nrT-1),1] + Uleft[2,2,1:(nrT-1),1],(nrT-1), 1)./n.^(1:(nrT-1)));
    end
end

# Evaluate the asymptotic expansion for the recursion coefficient beta_n of the generalised Jacobi polynomial.
# Input
#   n            - Degree at which to evaluate
#   nrT          - Number of terms. 1 gives leading order term
#   Dinf         - Limit of the Szego function
#   Uright,Uleft - Right and left U-matrices
# Output
#   betan       - Asymptotic expansion of the beta_n
function betan(n,nrT,Dinf,Uright,Uleft)
    if (nrT == 1) || (length(Uright) == 0)
	return 1/4;
    else
	return (1/(2im*Dinf^2) + sum(reshape(Uright[2,1,1:(nrT-1),1] + Uleft[2,1,1:(nrT-1),1],(nrT-1), 1)./n.^(1:(nrT-1))) )*(-Dinf^2/(2im) + sum(reshape(Uright[1,2,1:(nrT-1),1] + Uleft[1,2,1:(nrT-1),1],(nrT-1), 1)./n.^(1:(nrT-1))) );
    end
end


