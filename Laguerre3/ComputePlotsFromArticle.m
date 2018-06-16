% Make the plots from the article: run compRecLagBF.jl first to get the txt files and adjust the boolean below.
% About
%   Author       - Peter Opsomer (Peter.Opsomer@cs.kuleuven.be)
%   History      - Created July 2015, refactored May 2016 and November 2016
%% Initialising 
format longe; 
close all;
clear variables;
set(0,'DefaultFigureWindowStyle','docked');

compRecLagBF = 0; % Set to a nonzero when that Julia script has been run
prefix = ''; % Change to the relative path to the .txt result files from the Julia script


%% Standard Laguerre polynomials
clearvars -except compRecLagBF prefix
alpha = 0; q = [0 1];
maxOrder = 7; maxP2 = 8;

[P,gammaP,aP,bP] = exactLaguerre(alpha,q,2^maxP2);
s = getAsy(alpha,q,maxOrder);
% The following can be used to check correctness of the MRS numbers:
% dQ = @(x) (q(2:end).*(1:length(q)-1))*repmat(x, length(q)-1, 1).^repmat((0:length(q)-2)', 1, length(x)); 
pEv = zeros(maxP2,1);
aeEv = zeros(maxP2,maxOrder);
zt = 0.001;
for tn = 1:maxP2
	pEv(tn) = P(zt*s.betan(2^tn,maxOrder),2^tn +1); % Can also do:
%     pEv(tn) = laguerreL(2^tn, alpha, zt*s.betan(2^tn,maxOrder))*(-1)^(2^tn)*sqrt(gamma(2^tn+alpha+1)/factorial(2^tn));
	for i = 1:maxOrder
		aeEv(tn,i) = s.pd(2^tn,zt*s.betan(2^tn,maxOrder),i,0); % If the last argument would not be zero, we would use the Q-matrices.
	end
end
plotConv(pEv, aeEv, '',0,s);


%% Monomial Q
clearvars -except compRecLagBF prefix
alpha = 2.8; q = [1.5 0 0 0.7];
maxOrder = 6; maxP2 = 7;

[P,gammaP,aP,bP] = exactLaguerre(alpha,q,2^maxP2);
s = getAsy(alpha,q,maxOrder);

pEv = zeros(maxP2,1); pEvBf = zeros(maxP2,1); 
aeEv = zeros(maxP2,maxOrder);
if compRecLagBF
    aPbf = load([prefix 'anPsMon0p7ThirdDegImag.txt']);
    bPbf = load([prefix 'bnm1PsMon0p7ThirdDegImag.txt']);
else
    aPbf = nan*zeros(2^maxP2+1,1);
    bPbf = nan*zeros(2^maxP2+1,1);
end
zt = (1-1i)/100;
for tn = 1:maxP2
	pEvBf(tn) = orthonorm(zt*s.betan(2^tn,maxOrder), 2^tn, aPbf, bPbf);
	pEv(tn) = P(zt*s.betan(2^tn,maxOrder),2^tn +1);
	for i = 1:maxOrder
		aeEv(tn,i) = s.pd(2^tn,zt*s.betan(2^tn,maxOrder),i,1);
	end
end
% The following shows the problem of using only 16 digits (in Matlab) for calculating the recurrence coefficients.
plotConv(pEv, aeEv, 'Monomial Q: Recurrence coefficients computed without BigFloats',0,s); 
if compRecLagBF
    plotConv(pEvBf, aeEv, '',0,s);
end


%% Standard normalised Hermite
clearvars -except compRecLagBF prefix
maxP2 = 8;
maxOrder = 5;
H = exactLaguerre(nan,'Hermite', 2^maxP2*2+1);
se = getAsy(-1/2,[0 1], maxOrder);
zt = 0.97;
% zt = 0.01; % To compare bessel and bulk
heEv = zeros(maxP2,1);
evEv = zeros(maxP2,maxOrder);
evEvbulk = zeros(maxP2,maxOrder);
evEvbes = zeros(maxP2,maxOrder);
for tn = 1:maxP2
    nt = 2^tn;
    xt = sqrt(4*nt*zt);
    heEv(tn) = H(xt, 2*nt+1,1);
    for i = 1:maxOrder
        evEv(tn,i) = se.pc(nt,xt^2,i,1);
        evEvbulk(tn,i) = se.pb(nt,xt^2,i);
        evEvbes(tn,i) = se.pd(nt,xt^2,i,1);
%         evEvbes(tn,i) = se.pd(nt,xt^2,i,0);
    end
end
plotConv(heEv, evEv,'',0,se);
% norm(evEvbulk-evEvbes)/norm(evEvbulk) % To compare bessel and bulk


%% Nonstandard Hermite
clearvars -except compRecLagBF prefix
maxP2 = 7;
maxOrder = 11;
alpha = 0; q = [0 -3 1];
[PH,gammaPH,aPH,bPH] = exactLaguerre(alpha, q,-2^maxP2*2-1);
so = getAsy(alpha+1/2,q, maxOrder);

zt = 0.31;
% zt = 0.02; % To compare bessel and bulk
if compRecLagBF
    aPbf = load([prefix 'anPsHermite_x^4.txt']);
    bPbf = load([prefix 'bnm1PsHermite_x^4.txt']);
else
    aPbf = nan*zeros(2^maxP2*2+1,1);
    bPbf = nan*zeros(2^maxP2*2+2,1);
end
heEv = zeros(maxP2,1);
pEvbBf = zeros(maxP2,1);
evEv = zeros(maxP2,maxOrder);
evEvbulk = zeros(maxP2,maxOrder);
evEvbes = zeros(maxP2,maxOrder);

for tn = 1:maxP2
    nt = 2^tn;
    xt = sqrt(so.betan(nt,maxOrder)*zt);
    heEv(tn) = PH(xt, 2*nt+1+1);
	pEvbBf(tn) = orthonorm(sqrt(zt*so.betan(nt,maxOrder)), 2*nt+1, aPbf, bPbf);
    for i = 1:maxOrder
        evEv(tn,i) = xt*so.pb(nt,xt^2,i);
        evEvbulk(tn,i) = xt*so.pb(nt,xt^2,i);
        evEvbes(tn,i) = xt*so.pd(nt,xt^2,i,1);
%         evEvbes(tn,i) = se.pd(nt,xt^2,i,0);
    end
end
if compRecLagBF
    plotConv(pEvbBf, evEv,'Nonstandard Hermite: convergence of poly with recurrence coefficients computed with Bigfloat',0,so);
end
plotConv(heEv, evEv,'',0,so);
real((evEvbes-evEvbulk)./evEvbes)


%% Test Q = exp
clearvars -except compRecLagBF prefix
alpha = -1/2; q = @(x) exp(x); dQ = @(x) exp(x); Qinv = @(y) log(y);
% First check whether the explicit expression we found holds numerically: 
betas = linspace(-1.1e2,6e2, 20)';
intbs = nan*betas;
besbs = nan*betas;
for bi = 1:length(betas)
    intbs(bi) = quadgk(@(x) dQ(x).*sqrt(x./(betas(bi)-x)), 0, betas(bi))/2/pi;
    besbs(bi) = betas(bi)*exp(betas(bi)/2)*(besseli(0,betas(bi)/2) + besseli(1,betas(bi)/2) )/4;
end
transpose((intbs-besbs)./besbs)
figure; 
semilogx(besbs(betas>10), log(besbs(betas>10)) -betas(betas>10) );
hold on;
semilogx(besbs(betas>10), log(besbs(betas>10)) -log(log(8*pi)+2*log(besbs(betas>10)))/2 - betas(betas>10) );
semilogx(besbs(betas>10), log(besbs(betas>10)) -log(log(8*pi)+2*log(besbs(betas>10)))/2 +log(8*pi)/2- betas(betas>10) );
xlabel('n');
ylabel('Error');
legend('log(n)-\beta_n','log(n)-log(log(8\pi n^2))/2-\beta_n','log(n)-log(log(8\pi n^2))/2+log(8\pi)/2-\beta_n');

% Make the convergence plots:
maxP2 = 7; maxOrder = 5;
[P,gammaP,aP,bP] = exactLaguerre(alpha,q,2^maxP2);
pEvb = zeros(maxP2,1); pEvbBf = zeros(maxP2,1); 
bEv = zeros(maxP2,maxOrder); bnm1 = zeros(maxP2,maxOrder); 
bs = zeros(maxP2,1); bsbf = zeros(maxP2,1);
zt = 0.6;
if compRecLagBF
    aPbf = load([prefix 'anPsExp.txt']);
    bPbf = load([prefix 'bnm1PsExp.txt']);
else
    aPbf = nan*zeros(2^maxP2+1,1);
    bPbf = nan*zeros(2^maxP2+1,1);
end

for tn = 1:maxP2
	bs(tn) = bP(2^tn+1);
	bsbf(tn) = bPbf(2^tn+1);
    s = getAsy(alpha,q,maxOrder, dQ,Qinv, 2^tn); % Recompute the expansions for each n
	pEvbBf(tn) = orthonorm(zt*s.betan(2^tn,maxOrder), 2^tn, aPbf, bPbf);
	pEvb(tn) = P(zt*s.betan(2^tn,maxOrder),2^tn +1);
	for i = 1:maxOrder
		bEv(tn,i) = s.pb(2^tn,zt*s.betan(2^tn,maxOrder),i);
		bnm1(tn,i) = s.bnm1(2^tn,i);
	end
end
plotConv(bs,bnm1,'Q(x) = exp(x): b_{n-1} without using BigFloats',0,s);
plotConv(pEvb,bEv,'Q(x) = exp(x): Bulk poly without BigFloats',0,s);
if compRecLagBF
    plotConv(bsbf,bnm1,'',0,s);
    plotConv(pEvbBf,bEv,'',0,s);
end


%% General polynomial Q(z) = z^6 -21/10*z^5 + 3*z^3 -6*z^2 +9 
clearvars -except compRecLagBF prefix
alpha = 1.1;
q = [9 0 -6 3 0 -2.1 1];
maxP2 = 7;
mop = 7; % Maximal order for procedure with general polynomials
maxOrder = 7; % " functions
% mop should be m=6 times higher than maxOrder, but that would take too much time.

qf = @(x) x.^6 -21/10*x.^5 + 3*x.^3 -6*x.^2 +9; 
dQ = @(x) 6*x.^5 -21/10*5*x.^4 + 3*3*x.^2 -6*2*x; 
Qinv = @(y) y.^(1/6);

[P,gammaP,aP,bP] = exactLaguerre(alpha,q,2^maxP2);
tic
s = getAsy(alpha,q,mop);
timePrecGenPoly = toc;

pEv = zeros(maxP2,1);  pEvBf = zeros(maxP2,1); 
timeGenFct = zeros(maxP2,1);  timeEvGenPoly = zeros(maxP2,mop);  timeEvGenFct = zeros(maxP2,maxOrder);
aeEv = zeros(maxP2,mop);  aeEvf = zeros(maxP2,maxOrder);
if compRecLagBF
    aPbf = load([prefix 'anPsgenPoly.txt']);
    bPbf = load([prefix 'bnm1PsgenPoly.txt']);
else
    aPbf = nan*zeros(2^maxP2+1,1);
    bPbf = nan*zeros(2^maxP2+1,1);
end

zt = 0.99+0.02i;
for tn = 1:maxP2
	pEvBf(tn) = orthonorm(zt*s.betan(2^tn,maxOrder), 2^tn, aPbf, bPbf);
	pEv(tn) = P(zt*s.betan(2^tn,maxOrder),2^tn +1);
    tic
    sf = getAsy(alpha,qf,maxOrder, dQ,Qinv, 2^tn); % Recompute the expansions for each n
    timeGenFct(tn) = toc;
    for i = 1:mop
        tic
        aeEv(tn,i) = s.pc(2^tn,zt*s.betan(2^tn,maxOrder),i,1);
        timeEvGenPoly(tn,i) = toc;
    end
    for i=1:maxOrder
        tic
        aeEvf(tn,i) = sf.pc(2^tn,zt*s.betan(2^tn,maxOrder),i,1); % Use s.betan iso sf.betan to get the same x-coordinate.
        timeEvGenFct(tn,i) = toc;
    end
end
plotConv(pEv, aeEv, 'General poly: proc gen poly without BigFloats',0,s);
plotConv(pEv, aeEvf, 'General poly: proc gen fct without BigFloats',0,s);
if compRecLagBF
    plotConv(pEvBf, aeEv, '', 0, s);
    plotConv(pEvBf, aeEvf, '', 0, sf);
end
