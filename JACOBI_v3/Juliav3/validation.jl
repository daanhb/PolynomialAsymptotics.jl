# validation.jl: contains all functions enabling a comparison with the exact polynomials
# About
#   Author       - Peter Opsomer (Peter.Opsomer@cs.kuleuven.be)
#   History      - Created October 2014, last edit October 2015


# R_JACOBI Recurrence coefficients for monic Jacobi polynomials.
#    ab=R_JACOBI(n,a,b) generates the first n recurrence 
#    coefficients for monic Jacobi polynomials with parameters 
#    a and b. These are orthogonal on [-1,1] relative to the
#    weight function w(t)=(1-t)^a(1+t)^b. The n alpha-coefficients
#    are stored in the first column, the n beta-coefficients in
#    the second column, of the nx2 array ab. The call ab=
#    R_JACOBI(n,a) is the same as ab=R_JACOBI(n,a,a) and
#    ab=R_JACOBI(n) the same as ab=R_JACOBI(n,0,0).
#    Supplied by Dirk Laurie, 6-22-1998; edited by Walter Gautschi, 4-4-2002.
# About
#   Translated from https://www.cs.purdue.edu/archives/2002/wxg/codes/OPQ.html
function r_jacobi(N,a = 0,b = a)
if((N<=0)|(a<=-1)|(b<=-1)) 
	error("parameter(s) out of range")
end
nu=(b-a)/(a+b+2);
mu=2^(a+b+1)*gamma(a+1)*gamma(b+1)/gamma(a+b+2);
if N==1 return [nu mu]; end 
N=N-1; n=linspace(1,N,N); nab=2*n+a+b;
A=[nu; (b^2-a^2)./(nab.*(nab+2))]; 
n=2:N; nab=nab[n];
B1=4*(a+1)*(b+1)/((a+b+2)^2*(a+b+3));
B=4*(n+a).*(n+b).*n.*(n+a+b)./((nab.^2).*(nab+1).*(nab-1));
return [A'; [mu B1 B']]';
end

# GAUSS Gauss quadrature rule.
#    Given a weight function w encoded by the nx2 array ab of the 
#    first n recurrence coefficients for the associated orthogonal
#    polynomials, the first column of ab containing the n alpha-
#    coefficients and the second column the n beta-coefficients, 
#    the call xw=GAUSS(n,ab) generates the nodes and weights xw of
#    the n-point Gauss quadrature rule for the weight function w.
#    The nodes, in increasing order, are stored in the first 
#    column, the n corresponding weights in the second column, of
#    the nx2 array xw.
#
# About
#   Translated from https://www.cs.purdue.edu/archives/2002/wxg/codes/OPQ.html
function gauss(N,ab)
N0 = size(ab,1);
if N0 < N
	error("input array ab too short")
end
J=zeros(N,N);
for n=1:N
	J[n,n] = ab[n,1];
end
for n=2:N
  J[n,n-1] = sqrt(ab[n,2]);
  J[n-1,n] = J[n,n-1];
end
(D,V)=eig(J);
I=sortperm(D); #Actually already sorted..
V=V[:,I];
return (D, ab[1,2]*V[1,:]'.^2)
end

# Real orthonormal polynomials using the recurrence relation and the OPQ-library 
# of Walter Gautschi. These polynomials exhibit degrading accuracy from n about 
# 128 to about a relative error of 1e-8 around n=2^12, which can be checked by 
# evaluating exactPolys(0,0,@(x) 1,4096): then P(1,n+1) should be sqrt((2*n+1)/2).
# Input
#   alpha,beta,h - Parts of the weight w(x) = (1-x).^alpha.*(1+x).^beta.*h(x)
#   N            - Maximum degree
# Output 
#   P            - P(x,n) giving orthonormal poly of degree n at point x 
#                  as a function, ||pi(x,n)||_w = 1
#   gammaP       - gammaP(n+1)*2^n = \gamma_n = is the leading order coefficient 
#		   of the orthonormal polynomial of degree n. Multiplication by
#		   2^(-n) to avoid overflow of saved value gammaP(n+1)
#   alphaP,betaP - alphaP(n+1) & betaP(n+1) are the recurrence coefficients of
#		   the monic orthogonal polynomial of degree n
# Note: This uses another algorithm than the Matlab version.
function exactPolys(alpha,beta,h,N)
# Use alpha and beta to construct quadrature rules for computing the
# real alphaP & betaP. ab(:,1) != alphaP because it doesn't use h(x)
ab = r_jacobi(3*N,alpha,beta);
(z, w) = gauss(3*N,ab);
z = (1+0im)*z
w = (1+0im)*w
function inwprod(f,g)
	sum(w.*h(z+0im).*f(z+0im).*conj(g(z+0im) ) );
end
alphaP = zeros(N+1,1);
betaP = zeros(N+1,1);
gammaP = zeros(N+2,1);
gammaP[1] = 1/sqrt(inwprod(x -> ones(x), x -> ones(x) ) );
alphaP[1] = inwprod(x -> x, x -> ones(x) )*gammaP[1]^2;
function mon(x,n,alphaP,betaP)
	pin = [(1+1im)*zeros(length(x),1) (1+0im)*ones(length(x),1)];
	for j=1:n
	    pin = [pin[:,2] (x-alphaP[j]+0im).*pin[:,2]-betaP[j].*pin[:,1]];
	end
	return pin[:,2]
end
gammaP[2] = 1/sqrt(inwprod(x -> mon(x,1,alphaP,betaP), x -> mon(x,1,alphaP,betaP) ) );
start = time();
prevToc = start;
for k = 1:N
        if time()-prevToc > 1.0
		prevToc = time()
		print("OPQ iteration ", k, ", ", round((k-1)/N,3), "% done, estim sec left = ", (prevToc-start)*(N-k+1.0)/(k-1.0), '\n');
        end
	alphaP[k+1] = inwprod(x -> x.*mon(x,k,alphaP,betaP), x -> mon(x,k,alphaP,betaP) )*gammaP[k+1]^2;
        betaP[k+1] = (gammaP[k]/gammaP[k+1] )^2;
        gammaP[k+2] = 1/sqrt(inwprod(x -> mon(x,k+1,alphaP,betaP), x -> mon(x,k+1,alphaP,betaP) ) );
end
function P(x,n) return (mon(x,n,alphaP,betaP)*gammaP[n+1])[1]; end
return (P,gammaP,alphaP,betaP)
end

# Plot the convergence results with LS convergence in loglog, in 1 figure.
# Input
#   r     - Actual values
#   a     - Approximation
#   t     - Title of the plot
#   shift - Degree shift with respect to 2.0.^(1:maxP2)
function plotConv(r,a,t,shift=0)
(maxP2,maxOrder) = size(a);
data =  abs(a - repmat(r,1,maxOrder) )./abs(repmat(r,1,maxOrder) );
colors = "bgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmyk";
ns = 2.0.^(1:maxP2)+shift;
for i = 1:maxOrder
	idx = maxP2;
	while data[idx,i] < 0.1^11
		idx = idx-1;
	end
	if i == 1
	        loglog(ns,data[:,i],string("o",colors[1]),label="1 term")
	else
	        loglog(ns,data[:,i],string("o",colors[i]),label=string(i," terms") )
	end
	loglog(ns,data[idx,i]./((ns./ns[idx]).^i), string(colors[i],"-") );
end
legend(loc=3);
xlabel("n"); 
ylabel("Relative error"); 
title(t);
show()
end

