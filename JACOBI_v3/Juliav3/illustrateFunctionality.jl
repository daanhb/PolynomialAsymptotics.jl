# Illustrate functionality of the code with a test suite, heuristics and Gaussian quadrature
# About
#   Author       - Peter Opsomer (Peter.Opsomer@cs.kuleuven.be)
#   History      - Created October 2014, last edit October 2015


## Test suite: test all functionality of the code.
using PyPlot

include("validation.jl")
include("precomputations.jl")
include("asy_expansions.jl")
include("UExplicit.jl")

alpha = 2*rand()-0.9; beta = 3*rand()-0.7; cc = rand(); dd = rand();
function h(x) return exp(cc*x) + sin(dd*x).^2+1; end
function dh(x) return cc*exp(cc*x) +2.*sin(dd*x).*dd*cos(dd*x); end

maxOrder = rand(2:10,1)[1];
# maxP2 lower-> maybe no convergence, higher-> roundoff errors in exactPolys
maxP2 = 7; 
shift = (rand(0:3,1))[1];
print("alpha=", alpha, ", beta=", beta, ", cc=", cc, ", dd=", dd, ", maxOrder=", maxOrder, ", shift=", shift, '\n')
(P,gammaP,alphaP,betaP) = exactPolys(alpha,beta,h,2^maxP2+shift);

rs = rand(1:27,1,3); # The indexes of the regions to be tested
zr = (0+0im)*zeros(27,1); # At most 3 of these will be used
th = (0+0im)*zeros(27,1);
#  Taking rand iso randn for abs(zr) since it make points not lie inside the Bernstein ellipse
for ri = rs
	if ri <= 3
	elseif ri <= 7 zr[ri] = 0.5*rand()*exp(2im*rand()*pi); # lens
	elseif ri <= 11 zr[ri] = 1.1*rand()*exp(2im*rand()*pi); # outer
	elseif ri <= 15 zr[ri] = 1 + 0.4*rand()*exp(2im*rand()*pi); # right disk
	elseif ri <= 19 zr[ri] = -1 + 0.4*rand()*exp(2im*rand()*pi); # left disk
	elseif ri <= 23 th[ri] = eps()^0.4*rand()*exp(2im*rand()*pi); zr[ri] = cos(th[ri]); # right disk Q
	else th[ri] = eps()^0.4*rand()*exp(2im*rand()*pi); zr[ri] = -cos(th[ri]); # left disk Q
	end
end

legs = ["Recurrence coefficients alpha_n","Recurrence coefficients beta_n","Leading order coefficients gamma_n of the orthonormal polynomials",# = 3
	"Lens monic","Derivative lens monic","Lens orthonormal","Derivative lens orthonormal", # = 7
	"Outer monic","Derivative outer monic","Outer orthonormal","Derivative outer orthonormal", # = 11
	"Right disk monic","Derivative right disk monic","Right disk orthonormal","Derivative right disk orthonormal", # = 15
	"Left disk monic","Derivative Left disk monic","Left disk orthonormal","Derivative left disk orthonormal", # = 19
	"Right monic with Q-s","Derivative right monic with Q-s","orthonormal Right with Q-s", "Derivative orthonormal Right with Q-s", # = 23
	"Left monic with Q-s","Derivative left monic with Q-s","orthonormal left with Q-s","derivative orthonormal left with Q-s"]; # = 27

tic()
nrQ = ceil(Integer, 2+2*maxOrder*ones(1,maxOrder-1));
rho = 3; # Trapezium rules may become inaccurate if rho is (much) higher
if (maximum(abs(imag(zr))) > (rho/2 - 1/rho/2)*0.9 ) || (maximum(real(zr)) > (rho/2 + 1/rho/2)*0.9 ) ||	(minimum(real(zr)) < -(rho/2 + 1/rho/2)*0.9 )
	warn("Possible errors when evaluating far from the interval")
end


# All zr should lie within the contour. For Q-s, size(c&d) = ceil(nrT/2) so get higher number of Q-s by
(c, d, Dinf, psi, dpsi, contpsi) = contour_integrals(alpha,beta,h,4*maximum(nrQ+1),rho,2000);
(Uright,Uleft,Qright,Qleft) = UQ(alpha,beta,Dinf,c,d,maxOrder,"UQW",nrQ);
timePrecompute = toc()

(Ure,Ule) = UExplicit(alpha,beta,Dinf,c,d);
rk = min(size(Uright,3), size(Ure,3), 6); rm = min(size(Uright,4), size(Ure,4), 3);
lk = min(size(Uleft,3), size(Ule,3), 6); lm = min(size(Uleft,4), size(Ule,4), 3);
errorUright = sum(sum(sum(sum((Ure[:,:,1:rk,1:rm]-Uright[:,:,1:rk,1:rm]).^2))))
errorUleft = sum(sum(sum(sum((Ule[:,:,1:lk,1:lm]-Uleft[:,:,1:lk,1:lm]).^2))))
print("Error in Uright = ", errorUright, " and in Uleft = ", errorUleft, '\n')

exact = (0+0im)*zeros(maxP2,1,27);
asy = (0+0im)*zeros(maxP2,maxOrder,27);
gammas = (0+0im)*zeros(maxP2,1);
hh = sqrt(eps()); 
# Using complex finite differences (imag(P(zr[ri]+1i*hh,nt))/hh;) as explained in "William Squire and George Trapp. Using complex variables to estimate derivatives fo real functions. SIAM Review, 40(1) pp. 110-112, 1998." could make that hh may be smaller and thus better estimate of derivative, but not if zr is not real
for tn = 1:maxP2
	nt = convert(Int64,2^tn+shift);
	gammas[tn] = gammaP[nt+1];
	    for ri = rs
		if ri == 1 exact[tn,1,ri] = alphaP[nt+1];
			for i = 1:maxOrder   asy[tn,i,ri] = alphan(nt,i,Uright,Uleft); end
		elseif ri == 2 exact[tn,1,ri] = betaP[nt+1];
			for i = 1:maxOrder	asy[tn,i,ri] = betan(nt,i,Dinf,Uright,Uleft); end
		elseif ri == 3 exact[tn,1,ri] = gammaP[nt+1];
			for i = 1:maxOrder   asy[tn,i,ri] = gamman(nt,i,Dinf,Uright,Uleft); end
		elseif ri == 4	exact[tn,1,ri] = P(zr[ri],nt)/gammas[tn];
			
			for i = 1:maxOrder   asy[tn,i,ri] = asy_lens(nt,zr[ri],alpha,beta,h,psi,i,Dinf,Uright,Uleft,"m"); end
		elseif ri == 5	exact[tn,1,ri] = (P(zr[ri],nt)-P(zr[ri]-hh,nt))/hh/gammas[tn];
			for i = 1:maxOrder   
				(_, t) = asy_lens(nt,zr[ri],alpha,beta,h,psi,i,Dinf,Uright,Uleft,"m",dh,dpsi); asy[tn,i,ri] = t; 
			end
		elseif ri == 6	exact[tn,1,ri] = P(zr[ri],nt);
			for i = 1:maxOrder   asy[tn,i,ri] = asy_lens(nt,zr[ri],alpha,beta,h,psi,i,Dinf,Uright,Uleft,"o"); end
		elseif ri == 7	exact[tn,1,ri] = (P(zr[ri],nt)-P(zr[ri]-hh,nt))/hh;
			for i = 1:maxOrder  	
				(_, t) = asy_lens(nt,zr[ri],alpha,beta,h,psi,i,Dinf,Uright,Uleft,"o",dh,dpsi); asy[tn,i,ri] = t;
			end			
		elseif ri == 8	exact[tn,1,ri] = P(zr[ri],nt)/gammas[tn];
			for i = 1:maxOrder   asy[tn,i,ri] = asy_outer(nt,zr[ri],alpha,beta,h,psi,i,Dinf,Uright,Uleft,"m"); end
		elseif ri == 9	exact[tn,1,ri] = (P(zr[ri],nt)-P(zr[ri]-hh,nt))/hh/gammas[tn];
			for i = 1:maxOrder 	
				(_, t) = asy_outer(nt,zr[ri],alpha,beta,h,psi,i,Dinf,Uright,Uleft,"m",dh,dpsi);	asy[tn,i,ri] = t;
			end
		elseif ri == 10	exact[tn,1,ri] = P(zr[ri],nt);
			for i = 1:maxOrder   asy[tn,i,ri] = asy_outer(nt,zr[ri],alpha,beta,h,psi,i,Dinf,Uright,Uleft,"o"); end
		elseif ri == 11	exact[tn,1,ri] = (P(zr[ri],nt)-P(zr[ri]-hh,nt))/hh;
			for i = 1:maxOrder  	
				(_, t) = asy_outer(nt,zr[ri],alpha,beta,h,psi,i,Dinf,Uright,Uleft,"o",dh,dpsi);	asy[tn,i,ri] = t; 
			end				
		elseif ri == 12	exact[tn,1,ri] = P(zr[ri],nt)/gammas[tn];
			for i = 1:maxOrder    asy[tn,i,ri] = asy_right(nt,zr[ri],alpha,beta,h,psi,i,Dinf,Uright,Uleft,"m"); end		
		elseif ri == 13	exact[tn,1,ri] = (P(zr[ri],nt)-P(zr[ri]-hh,nt))/hh/gammas[tn];
			for i = 1:maxOrder    
				(_, t) = asy_right(nt,zr[ri],alpha,beta,h,psi,i,Dinf,Uright,Uleft,"m",dh,dpsi); 
				asy[tn,i,ri] = t; 
			end	
		elseif ri == 14	exact[tn,1,ri] = P(zr[ri],nt);
			for i = 1:maxOrder    asy[tn,i,ri] = asy_right(nt,zr[ri],alpha,beta,h,psi,i,Dinf,Uright,Uleft,"o"); end	
		elseif ri == 15	exact[tn,1,ri] = (P(zr[ri],nt)-P(zr[ri]-hh,nt))/hh;
			for i = 1:maxOrder    
				(_, t) = asy_right(nt,zr[ri],alpha,beta,h,psi,i,Dinf,Uright,Uleft,"o",dh,dpsi); 
				asy[tn,i,ri] = t; 
			end				
		elseif ri == 16	exact[tn,1,ri] = P(zr[ri],nt)/gammas[tn];
			for i = 1:maxOrder    asy[tn,i,ri] = asy_left(nt,zr[ri],alpha,beta,h,psi,i,Dinf,Uright,Uleft,"m"); end
		elseif ri == 17
			exact[tn,1,ri] = (P(zr[ri],nt)-P(zr[ri]-hh,nt))/hh/gammas[tn];
			for i = 1:maxOrder    
				(_ , t) = asy_left(nt,zr[ri],alpha,beta,h,psi,i,Dinf,Uright,Uleft,"m",dh,dpsi); 
				asy[tn,i,ri] = t; 
			end
		elseif ri == 18	exact[tn,1,ri] = P(zr[ri],nt);
			for i = 1:maxOrder    asy[tn,i,ri] = asy_left(nt,zr[ri],alpha,beta,h,psi,i,Dinf,Uright,Uleft,"o"); end	
		elseif ri == 19	exact[tn,1,ri] = (P(zr[ri],nt)-P(zr[ri]-hh,nt))/hh;
			for i = 1:maxOrder     
				(_, t) = asy_left(nt,zr[ri],alpha,beta,h,psi,i,Dinf,Uright,Uleft,"o",dh,dpsi); 
				asy[tn,i,ri] = t; 
			end
		elseif ri == 20	exact[tn,1,ri] = P(zr[ri],nt)/gammas[tn];
			for i = 1:maxOrder  
				asy[tn,i,ri] = asy_right(nt,zr[ri],alpha,beta,h,psi,i,Dinf,Uright,Uleft,"mQ",false,false,Qright,th[ri],contpsi);
			end	
		elseif ri == 21	exact[tn,1,ri] = (P(zr[ri],nt)-P(zr[ri]-hh,nt))/hh/gammas[tn];
			for i = 1:maxOrder
				(_, t) = asy_right(nt,zr[ri],alpha,beta,h,psi,i,Dinf,Uright,Uleft,"mQ",dh,dpsi,Qright,th[ri],contpsi); 
				asy[tn,i,ri] = t; 
			end	
		elseif ri == 22	exact[tn,1,ri] = P(zr[ri],nt);
			for i = 1:maxOrder   
				asy[tn,i,ri] = asy_right(nt,zr[ri],alpha,beta,h,psi,i,Dinf,Uright,Uleft,"oQ",false,false,Qright,th[ri],contpsi); end	
		elseif ri == 23	exact[tn,1,ri] = (P(zr[ri],nt)-P(zr[ri]-hh,nt))/hh;
			for i = 1:maxOrder
				(_, t) = asy_right(nt,zr[ri],alpha,beta,h,psi,i,Dinf,Uright,Uleft,"oQ",dh,dpsi,Qright,th[ri],contpsi); 
				asy[tn,i,ri] = t; 
			end	
			
		elseif ri == 24	exact[tn,1,ri] = P(zr[ri],nt)/gammas[tn];
			for i = 1:maxOrder   
				asy[tn,i,ri] = asy_left(nt,zr[ri],alpha,beta,h,psi,i,Dinf,Uright,Uleft,"mQ",false,false,Qleft,th[ri],contpsi); end	
		elseif ri == 25	exact[tn,1,ri] = (P(zr[ri],nt)-P(zr[ri]-hh,nt))/hh/gammas[tn];
			for i = 1:maxOrder
				(_, t) = asy_left(nt,zr[ri],alpha,beta,h,psi,i,Dinf,Uright,Uleft,"mQ",dh,dpsi,Qleft,th[ri],contpsi); 
				asy[tn,i,ri] = t; 
			end	
		elseif ri == 26	exact[tn,1,ri] = P(zr[ri],nt);
			for i = 1:maxOrder
				asy[tn,i,ri] = asy_left(nt,zr[ri],alpha,beta,h,psi,i,Dinf,Uright,Uleft,"oQ",false,false,Qleft,th[ri],contpsi); end	
		elseif ri == 27	exact[tn,1,ri] = (P(zr[ri],nt)-P(zr[ri]-hh,nt))/hh;
			for i = 1:maxOrder
				(_, t) = asy_left(nt,zr[ri],alpha,beta,h,psi,i,Dinf,Uright,Uleft,"oQ",dh,dpsi,Qleft,th[ri],contpsi); 
				asy[tn,i,ri] = t; 
			end	
		end
	end
end
# Testing order of convergence by plotting the errors in the (3 randomly) chosen regions, take into account the remarks in the article when interpreting these
for ri = rs
	te = (0+0im)*zeros(maxP2,1);
	ta = (0+0im)*zeros(maxP2,maxOrder);
	for tn = 1:maxP2
		te[tn,1] = exact[tn,1,ri];
		for i = 1:maxOrder
			ta[tn,i] = asy[tn,i,ri];
		end
	end
	if ri <= 3 plotConv(te,ta, legs[ri]);
	elseif ri <= 19 plotConv(te,ta, string(legs[ri], " at z= ", zr[ri]) );
	else plotConv(te, ta, string(legs[ri], " at theta= ", th[ri]));
	end
end


## Heuristics: see which expansion is more accurate where.
workspace();
using PyPlot
include("validation.jl"); include("precomputations.jl"); include("asy_expansions.jl");

alpha = rand(); beta = -1+3*rand(); function h(x) 1; end
print("Heuristics: alpha=", alpha, ", beta=", beta, '\n')
n = 100;
maxOrder = 3;

(P,gammaP,alphaP,betaP) = exactPolys(alpha,beta,h,n);
# Computing relative errors on nrP points on this line
zs = linspace(-1,1,25) +0.02im;

(c, d, Dinf, psi, _, contpsi) = contour_integrals(alpha,beta,h,2*(maxOrder+3) );
(Uright,Uleft,Qright,Qleft) = UQ(alpha,beta,Dinf,c,d,maxOrder,"UQV");

relerrs = zeros(length(zs),maxOrder,6);
# Will give many warnings when evaluating the asymptotic expansions in the disks with Q-s far from endpoint. "Currently the only fix is to eliminate the warning", so do $ julia illustrateFunctionality.jl 2> junk
for zi = 1:length(zs)
	exact = P(zs[zi],n);
	for i = 1:maxOrder
		relerrs[zi,i,1] =  abs(asy_lens( n,zs[zi],alpha,beta,h,psi,i,Dinf,Uright,Uleft,"o")-exact)/abs(exact);
		relerrs[zi,i,2] =  abs(asy_outer(n,zs[zi],alpha,beta,h,psi,i,Dinf,Uright,Uleft,"o")-exact)/abs(exact);
		relerrs[zi,i,3] =  abs(asy_right(n,zs[zi],alpha,beta,h,psi,i,Dinf,Uright,Uleft,"o")-exact)/abs(exact);
		relerrs[zi,i,4] =  abs(asy_left( n,zs[zi],alpha,beta,h,psi,i,Dinf, Uright,Uleft,"o")-exact)/abs(exact);
		relerrs[zi,i,5] = abs(asy_right(n,zs[zi],alpha,beta,h,psi,i,Dinf, Uright,Uleft,"oQ",false,false,Qright,acos(zs[zi]),contpsi)-exact)/abs(exact);
		relerrs[zi,i,6] = abs(asy_left(n,zs[zi],alpha,beta,h,psi,i,Dinf, Uright,Uleft,"oQ",false,false,Qleft,acos(-zs[zi]),contpsi)-exact)/abs(exact);
	end
end

colors = "bgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmyk";
semilogy(real(zs),relerrs[:,1,4],string("o",colors[1]),label=string("1 term") )
for i = 2:maxOrder
	semilogy(real(zs),relerrs[:,i,4],string("o",colors[i]),label=string(i," terms") )
end
legend(loc=0);
xlabel(string("Real part of z, imaginary part is ", mean(imag(zs))) ); 
ylabel("Relative error"); 
title("Relative errors of expansion in the left disk for varying z and number of terms");
show()


T = 3;
legs =  ["lens","outer","right","left","right Q","left Q"];
for i = 1:6
	semilogy(real(zs),relerrs[:,T,i],string("o",colors[i]),label=legs[i])
end
legend(loc=10);
xlabel(string("Real part of z, imaginary part is ", mean(imag(zs))));
ylabel("Relative error");
title(string("Relative errors of different expansion for varying z and T=", T))
show()


## Gaussian quadrature: Newton procedure to compute nodes and weights
workspace(); include("precomputations.jl"); include("asy_expansions.jl");

cheb = rand(Bool);
if cheb
	alpha = -1/2; beta = -1/2; 
	function h(x) 1; end;  
	function dh(x) 0; end
else
	alpha = 1/sqrt(3); beta = -1/pi; 
	function h(x) exp(x); end;  
	function dh(x) exp(x); end
end
n = 200;
nrT = 3;

(c, d, Dinf, psi, dpsi, _) = contour_integrals(alpha,beta,h,nrT);
(Uright,Uleft) = UQ(alpha,beta,Dinf,c,d,nrT,"UW");

zw = zeros(n,3);
for ni = 1:n
	if ni <=3 # (Approximation of) zero of bessel function to approximate first zeros
		zw[ni,1] = -1 + (2*sqrt(beta+1) + (ni-1)*pi)^2/2/n^2;
		if zw[ni,1] +1 < eps()^(1/3)
			warn("Use Q-s");
		end
	else
		zw[ni,1] = 3*zw[ni-1] -3*zw[ni-2] +zw[ni-3]; # Quadratic extrapolation
	end
	dpoly = "saveFromLoop"; polyM1 = "sfl"; polyP1 = "sfl";
	# Could improve the following by stopping Newton iterates when already tried that zero and then going over all (max. 5) closest floats.
	for attempt = 1:10
		if zw[ni,1]+1 < 0.2 # This bound is not general for all weights!
			(poly, dpoly) = asy_left(n,zw[ni,1],alpha,beta,h,psi,nrT,Dinf,Uright,Uleft,"o",dh,dpsi);
			polyM1 = asy_left(n-1,zw[ni,1],alpha,beta,h,psi,nrT,Dinf,Uright,Uleft,"o");
			polyP1 = asy_left(n+1,zw[ni,1],alpha,beta,h,psi,nrT,Dinf,Uright,Uleft,"o");
		elseif 1-zw[ni,1] < 0.2
			(poly, dpoly) = asy_right(n,zw[ni,1],alpha,beta,h,psi,nrT,Dinf,Uright,Uleft,"o",dh,dpsi);
			polyM1 = asy_right(n-1,zw[ni,1],alpha,beta,h,psi,nrT,Dinf,Uright,Uleft,"o");
			polyP1 = asy_right(n+1,zw[ni,1],alpha,beta,h,psi,nrT,Dinf,Uright,Uleft,"o");
		else
			(poly, dpoly) = asy_lens(n,zw[ni,1],alpha,beta,h,psi,nrT,Dinf, Uright,Uleft,"o",dh,dpsi);
			polyM1 = asy_lens(n-1,zw[ni,1],alpha,beta,h,psi,nrT,Dinf, Uright,Uleft,"o");
			polyP1 = asy_lens(n+1,zw[ni,1],alpha,beta,h,psi,nrT,Dinf, Uright,Uleft,"o");
		end
		zw[ni,1] = zw[ni,1] - real(poly)/real(dpoly); # We know the polynomial is real on the interval
	end
	if (zw[ni,1] < -1) || (zw[ni,1] > 1)
		error("Outside of the interval");
	elseif (ni > 1) && (zw[ni,1] < zw[ni-1,1] )
		error("Not going right");
	end
	zw[ni,2] = gamman(n,nrT,Dinf,Uright,Uleft)/real(dpoly)/gamman(n-1,nrT,Dinf,Uright,Uleft)/real(polyM1);
	zw[ni,3] = -gamman(n+1,nrT,Dinf,Uright,Uleft)/real(dpoly)/gamman(n,nrT,Dinf,Uright,Uleft)/real(polyP1);
end
errzw = norm(zw[:,2]-zw[:,3])/norm(zw[:,2]) # Formulas should give the same result

# Test by integrating (x^2-1) and (x^3) over the interval
ae1 = "saveFromIf"; ae2 = "sfi";
if cheb # % Weight 1/sqrt(1-x^2)
	chebyroots = cos((2*(1:n)-1)/2/n*pi);
	errorOnRoots = norm(zw[:,1]+chebyroots)
	print("Error in roots = ", errorOnRoots, '\n')
	ae1 = sum(zw[:,2].*(zw[:,1].^2-1) ) - (-pi/2)
	ae2 = sum(zw[:,2].*(zw[:,1].^3) )
else # w(x) = (1-x)^( 1/sqrt(3))*(1+x)^(-1/pi)*exp(x)
	ae1 = sum(zw[:,2].*(zw[:,1].^2-1) ) -(-1.293225037)
	ae2 = sum(zw[:,2].*(zw[:,1].^3) ) - (-.1659585353)
end
# Integrate a random polynomial with a degree low enough to be integrated exactly
r = [randn(2,1); rand(0:(2*n-3),1); rand(0:5,1)];
function fct(x) return 7*r[1]*x.^r[3] -23*r[2]*x.^r[4]; end
function w(x) return (1-x).^alpha.*(1+x).^beta.*h(x); end
# Better use the Cubature package for approximating the integral
xs = linspace(1/20000-1,1-1/20000,10000)
ae3 = sum(zw[:,2].*fct(zw[:,1]) ) -sum(w(xs).*fct(xs))/5000
print("Error in weights = ", errzw, ", errors in first, second and third integral are ", ae1, ", ", ae2, " and ", ae3, '\n')


