% Make the plots from the article.
% About
%   Author       - Peter Opsomer (Peter.Opsomer@cs.kuleuven.be)
%   History      - Created October 2013, last edit October 2015
%% Initialising
format longe; close all; clear variables;
set(0,'DefaultFigureWindowStyle','docked');

for fig = 1:2
    switch fig
        case 1
            alpha = 0; beta = 0; h = @(x) exp(-7*x.^(2*2) );
			maxOrder = 7;maxP2 = 8;
        case 2
            alpha = -1/2; beta = 0; h = @(x) 1./sqrt(x+3); 
			maxOrder = 6;maxP2 = 7;
    end
    tic;
    [P,gammaP,alphaP,betaP] = exactPolys(alpha,beta,h,2^maxP2);
    timeOPQ = toc
    
    %% Computing higher order terms
    tic;
    [c, d, Dinf, psi, dpsi] = contour_integrals(alpha,beta,h,maxOrder);
    [Uright,Uleft] = UQ(alpha,beta,Dinf,c,d,maxOrder);
    timePrecompute = toc
    
    %% Computing exact results and asymptotic expansions
    zInt = 0.2+0.5*1i;
    zL = -0.97;
    zOut = zInt; % Should lie within contour !
    
    piEvi = zeros(maxP2,1); piEvbL = zeros(maxP2,1); piEvo = zeros(maxP2,1);
	ipiEvi = zeros(maxP2,maxOrder); bpiLEv = zeros(maxP2,maxOrder); opiEv = zeros(maxP2,maxOrder);
	gammas = zeros(maxP2,1);
	for tn = 1:maxP2
		gammas(tn) = gammaP(2^tn+1)*2^(2^tn);
	end
	
    for tn = 1:maxP2
        piEvi(tn) = P(zInt,2^tn +1)/gammas(tn);
        piEvbL(tn) = P(zL,2^tn +1)/gammas(tn);
        piEvo(tn) = P(zOut,2^tn +1)/gammas(tn);
        for i = 1:maxOrder
            ipiEvi(tn,i) = asy_lens(2^tn,zInt,alpha,beta,h,psi,i,Dinf,Uright,Uleft);
            bpiLEv(tn,i) = asy_left(2^tn,zL,alpha,beta,h,psi,i,Dinf,Uright,Uleft);
            opiEv(tn,i) = asy_outer(2^tn,zOut,alpha,beta,h,psi,i,Dinf,Uright,Uleft);
        end
    end
    
    %% Testing order of convergence
    switch fig
        case 1
            plotConv(piEvbL,bpiLEv,'left disk');
        case 2
            plotConv(piEvi,ipiEvi,'lens');
            plotConv(piEvo,opiEv,'outside region');
    end
    shg
end
%% Toda measure
ntest = 101;
order = 1;
alpha = -1/2;
beta = -1/2;

ts = [-2,0,2,5];
xs = linspace(-1+eps^(1/3),1-eps^(1/3),2000); % Avoid singularities and large roundoff errors
pa = zeros(length(xs), length(ts));
pq = zeros(length(xs), length(ts));
errX = zeros(1,length(ts));
for ti = 1:length(ts)
	tToda = ts(ti)
	h = @(x) exp(-tToda*x);
	P = exactPolys(alpha,beta,h,ntest);
	[c, d, Dinf, psi] = contour_integrals(alpha,beta,h,order);
	[Uright,Uleft] = UQ(alpha,beta,Dinf,c,d,order);
	for j = 1:length(xs)
		pq(j,ti) = P(xs(j),ntest+1);
		pa(j,ti) = asy_lens(ntest, xs(j),alpha,beta,h,psi,order,Dinf,Uright,Uleft,'o');
	end
	errX(ti) = norm(pq(:,ti)-pa(:,ti))/norm(pq(:,ti));
end
errX % Error over x
maxAbsErr = max(abs(pq-pa)./abs(pq)) % Elementwise error
figure;
plot(xs, real(pa(:,1)),'-b'); hold on
plot(xs, real(pa(:,2)),'--g');
plot(xs, real(pa(:,3)),':r');
plot(xs, real(pa(:,4)),'-.c');
xlabel('x');
ylabel('p_{101}(x,t)');
legend({['t=' num2str(ts(1)) ], ['t=' num2str(ts(2))], ['t=' num2str(ts(3))], ['t=' num2str(ts(4))]})
