# Main script, computes the plots of the article
# About
#   Author       - Peter Opsomer (Peter.Opsomer@cs.kuleuven.be)
#   History      - Created October 2014, last edit October 2015
#using PyCall 
#@pyimport matplotlib.pyplot as plt
using PyPlot
include("validation.jl")
include("precomputations.jl")
include("asy_expansions.jl")

for fig = 1:2
	if fig == 1
		alpha = 0; beta = 0;
		function h(x) 	exp(-7*x.^(2*2) )	end
		maxOrder = 7;maxP2 = 8;
	elseif fig == 2
		alpha = -1/2; beta = 0; 
		function h(x) 1./sqrt(x+3) end
		maxOrder = 6;
		maxP2 = 7;
	end
	tic();
	(P,gammaP,alphaP,betaP) = exactPolys(alpha,beta,h,2^maxP2);
	timePrec = toc()
	tic();
	(c, d, Dinf, psi, dpsi) = contour_integrals(alpha,beta,h,maxOrder);
	(Uright,Uleft) = UQ(alpha,beta,Dinf,c,d,maxOrder);
	timePrecompute = toc()

	zInt = 0.2+0.5*1im;
	zL = -0.97;
	zOut = zInt; # Should lie within contour !

	piEvi = (1+1im)*zeros(maxP2,1); piEvbL = (1+1im)*zeros(maxP2,1); piEvo = (1+1im)*zeros(maxP2,1);
	ipiEvi = (1+1im)*zeros(maxP2,maxOrder); bpiLEv = (1+1im)*zeros(maxP2,maxOrder); opiEv = (1+1im)*zeros(maxP2,maxOrder);
	gammas = (1+1im)*zeros(maxP2,1); 
	for tn = 1:maxP2
		gammas[tn] = gammaP[2^tn+1];
	end
	
	for tn = 1:maxP2
		piEvi[tn] = P([zInt],2^tn)/gammas[tn];
		piEvbL[tn] = P([zL],2^tn)/gammas[tn];
		piEvo[tn] = P([zOut],2^tn)/gammas[tn];
		for i = 1:maxOrder
			ipiEvi[tn,i] = asy_lens(2^tn,zInt,alpha,beta,h,psi,i,Dinf,Uright,Uleft);
			bpiLEv[tn,i] = asy_left(2^tn,zL,alpha,beta,h,psi,i,Dinf,Uright,Uleft);
			opiEv[tn,i] = asy_outer(2^tn,zOut,alpha,beta,h,psi,i,Dinf,Uright,Uleft);
		end
	end
	if fig == 1
		plotConv(piEvbL,bpiLEv,"Left disk");
	elseif fig == 2
		plotConv(piEvi,ipiEvi,"Lens");
		plotConv(piEvo,opiEv,"Outside region");
	end
end
## Toda measure
ntest = 101;
order = 1;
alpha = -1/2;
beta = -1/2;

ts = [-2;0;2;5];
xs = linspace(-1+eps()^(1/3),1-eps()^(1/3),2000); # Avoid singularities and large roundoff errors
pa = (1+1im)*zeros(length(xs), length(ts));
pq = (1+1im)*zeros(length(xs), length(ts));
errX = 1.0*zeros(1,length(ts));
for ti = 1:length(ts)
	tToda = ts[ti];
	print("tToda = ", tToda, ", ");
	function h(x) exp(-tToda*x); end
	(P,) = exactPolys(alpha,beta,h,ntest);
	(c, d, Dinf, psi) = contour_integrals(alpha,beta,h,order);
	(Uright,Uleft) = UQ(alpha,beta,Dinf,c,d,order);
	for j = 1:length(xs)
		pq[j,ti] = P([xs[j]],ntest);
		pa[j,ti] = asy_lens(ntest, xs[j],alpha,beta,h,psi,order,Dinf,Uright,Uleft,"o");
	end
	errX[ti] = norm(pq[:,ti]-pa[:,ti])/norm(pq[:,ti]);
end
print('\n', "Error over x = ", errX, " and maximal elementwise absolute relative error =", maximum(abs(pq-pa)./abs(pq)), '\n');

plot(xs,real(pa[:,1]),"b-",label=string("t = ", ts[1]) )
plot(xs,real(pa[:,2]),"g-",label=string("t = ", ts[2]) )
plot(xs,real(pa[:,3]),"r-",label=string("t = ", ts[3]) )
plot(xs,real(pa[:,4]),"c-",label=string("t = ", ts[4]) )

legend(loc=3);
xlabel("x"); 
ylabel(r"$p_{101}(x,t)$");
show()

