% Illustrate functionality of the code with a test suite, heuristics and Gaussian quadrature
% About
%   Author       - Peter Opsomer (Peter.Opsomer@cs.kuleuven.be)
%   History      - Created October 2013, last edit February 2015
%% Test suite: test all functionality of the code.
close all; clear variables;%all; 
format longe; clc;
set(0, 'DefaultFigureWindowStyle', 'docked');

alpha = 2*rand-0.9, beta = 3*rand-0.7, cc = rand, dd = rand
h = @(x) exp(cc*x) + sin(dd*x).^2 +1; dh = @(x) cc*exp(cc*x) +2.*sin(dd*x).*dd*cos(dd*x);

maxOrder = randi([2 10],1)
% maxP2 lower-> maybe no convergence, higher-> roundoff errors in exactPolys
maxP2 = 7; 
shift = randi([0 3],1)
[P,gammaP,alphaP,betaP] = exactPolys(alpha,beta,h,2^maxP2+shift);

rs = randi(27,1,3); % The indexes of the regions to be tested
zr = zeros(27,1); % At most 3 of these will be used
th = zeros(27,1);
% Taking rand iso randn for abs(zr) since it make points not lie inside the Bernstein ellipse
for ri = rs
	if ri <= 3,
	elseif ri <= 7, zr(ri) = 0.5*rand*exp(2i*rand*pi); % lens
	elseif ri <= 11, zr(ri) = 1.1*rand*exp(2i*rand*pi); % outer
	elseif ri <= 15, zr(ri) = 1 + 0.4*rand*exp(2i*rand*pi); % right disk
	elseif ri <= 19, zr(ri) = -1 + 0.4*rand*exp(2i*rand*pi); % left disk
	elseif ri <= 23, th(ri) = eps^0.4*rand*exp(2i*rand*pi); zr(ri) = cos(th(ri)); % right disk Q
	else th(ri) = eps^0.4*rand*exp(2i*rand*pi); zr(ri) = -cos(th(ri)); % left disk Q
end, end

legs = {'Recurrence coefficients alpha_n','Recurrence coefficients beta_n',...
	'Leading order coefficients gamma_n of the orthonormal polynomials',...% = 3
	'Lens monic','Derivative lens monic','Lens orthonormal','Derivative lens orthonormal',... % = 7
	'Outer monic','Derivative outer monic','Outer orthonormal','Derivative outer orthonormal',...% = 11
	'Right disk monic','Derivative right disk monic','Right disk orthonormal','Derivative right disk orthonormal',...
	'Left disk monic','Derivative Left disk monic','Left disk orthonormal','Derivative left disk orthonormal',... % = 19
	'Right monic with Q-s','Derivative right monic with Q-s','orthonormal Right with Q-s',... % = 22
	'Derivative orthonormal Right with Q-s',...
	'Left monic with Q-s','Derivative left monic with Q-s','orthonormal left with Q-s',... % = 26
	'Derivative orthonormal left with Q-s'};

tic
nrQ = 2+2*maxOrder*ones(1,maxOrder-1);
rho = 3; % Trapezium rules may become inaccurate if rho is (much) higher
if (max(abs(imag(zr))) > (rho/2 - 1/rho/2)*0.9 ) || (max(real(zr)) > (rho/2 + 1/rho/2)*0.9 ) ||...
		(min(real(zr)) < -(rho/2 + 1/rho/2)*0.9 )
	warning('Possible errors when evaluating far from the interval')
end

% All zr should lie within contour. For Q-s, size(c&d) = ceil(nrT/2) so get higher number of Q-s by
[c, d, Dinf, psi, dpsi, contpsi] = contour_integrals(alpha,beta,h,4*max(nrQ+1),rho,2000);
[Uright,Uleft,Qright,Qleft] = UQ(alpha,beta,Dinf,c,d,maxOrder,'UQW',nrQ);
timePrecompute = toc

[Ure,Ule] = UExplicit(alpha,beta,Dinf,c,d);
rk = min([size(Uright,3), size(Ure,3),6]); rm = min([size(Uright,4), size(Ure,4),3]);
lk = min([size(Uleft,3), size(Ule,3),6]); lm = min([size(Uleft,4), size(Ule,4),3]);
errorUright = sum(sum(sum(sum((Ure(:,:,1:rk,1:rm)-Uright(:,:,1:rk,1:rm)).^2))))
errorUleft = sum(sum(sum(sum((Ule(:,:,1:lk,1:lm)-Uleft(:,:,1:lk,1:lm)).^2))))

exact = zeros(maxP2,1,27);
asy = zeros(maxP2,maxOrder,27);
gammas = zeros(maxP2,1);
hh = sqrt(eps); 
% Using complex finite differences (imag(P(zr(ri)+1i*hh,nt+1))/hh;) as explained in "William Squire and George Trapp. 
% Using complex variables to estimate derivatives fo real functions. SIAM Review, 40(1) pp. 110-112, 1998."
% could make that hh may be smaller and thus better estimate of derivative, but not if zr is not real
for tn = 1:maxP2
	nt = 2^tn+shift;
	gammas(tn) = gammaP(nt+1)*2^nt;
	for ri = rs,	switch ri
		case 1, exact(tn,1,ri) = alphaP(nt+1);
			for i = 1:maxOrder, asy(tn,i,ri) = alphan(nt,i,Uright,Uleft); end
		case 2, exact(tn,1,ri) = betaP(nt+1);
			for i = 1:maxOrder,	asy(tn,i,ri) = betan(nt,i,Dinf,Uright,Uleft); end
		case 3, exact(tn,1,ri) = gammaP(nt+1)*2^nt;
			for i = 1:maxOrder, asy(tn,i,ri) = gamman(nt,i,Dinf,Uright,Uleft); end
		case 4,	exact(tn,1,ri) = P(zr(ri),nt+1)/gammas(tn);
			
			for i = 1:maxOrder, asy(tn,i,ri) = asy_lens(nt,zr(ri),alpha,beta,h,psi,i,Dinf,Uright,Uleft,'m'); end
		case 5,	exact(tn,1,ri) = (P(zr(ri),nt+1)-P(zr(ri)-hh,nt+1))/hh/gammas(tn);
			for i = 1:maxOrder, [~, t] = asy_lens(nt,zr(ri),alpha,beta,h,psi,i,Dinf,Uright,Uleft,'m',dh,dpsi); asy(tn,i,ri) = t; end
		case 6,	exact(tn,1,ri) = P(zr(ri),nt+1);
			for i = 1:maxOrder, asy(tn,i,ri) = asy_lens(nt,zr(ri),alpha,beta,h,psi,i,Dinf,Uright,Uleft,'o'); end
		case 7,	exact(tn,1,ri) = (P(zr(ri),nt+1)-P(zr(ri)-hh,nt+1))/hh;
			for i = 1:maxOrder,	[~, t] = asy_lens(nt,zr(ri),alpha,beta,h,psi,i,Dinf,Uright,Uleft,'o',dh,dpsi); asy(tn,i,ri) = t;end
			
		case 8,	exact(tn,1,ri) = P(zr(ri),nt+1)/gammas(tn);
			for i = 1:maxOrder, asy(tn,i,ri) = asy_outer(nt,zr(ri),alpha,beta,h,psi,i,Dinf,Uright,Uleft,'m'); end
		case 9,	exact(tn,1,ri) = (P(zr(ri),nt+1)-P(zr(ri)-hh,nt+1))/hh/gammas(tn);
			for i = 1:maxOrder,	[~, t] = asy_outer(nt,zr(ri),alpha,beta,h,psi,i,Dinf,Uright,Uleft,'m',dh,dpsi);	asy(tn,i,ri) = t;end
		case 10,	exact(tn,1,ri) = P(zr(ri),nt+1);
			for i = 1:maxOrder, asy(tn,i,ri) = asy_outer(nt,zr(ri),alpha,beta,h,psi,i,Dinf,Uright,Uleft,'o'); end
		case 11,	exact(tn,1,ri) = (P(zr(ri),nt+1)-P(zr(ri)-hh,nt+1))/hh;
			for i = 1:maxOrder,	[~, t] = asy_outer(nt,zr(ri),alpha,beta,h,psi,i,Dinf,Uright,Uleft,'o',dh,dpsi);	asy(tn,i,ri) = t; end
				
		case 12,	exact(tn,1,ri) = P(zr(ri),nt+1)/gammas(tn);
			for i = 1:maxOrder, asy(tn,i,ri) = asy_right(nt,zr(ri),alpha,beta,h,psi,i,Dinf,Uright,Uleft,'m'); end		
		case 13,	exact(tn,1,ri) = (P(zr(ri),nt+1)-P(zr(ri)-hh,nt+1))/hh/gammas(tn);
			for i = 1:maxOrder, [~, t] = asy_right(nt,zr(ri),alpha,beta,h,psi,i,Dinf,Uright,Uleft,'m',dh,dpsi); 
				asy(tn,i,ri) = t; end	
		case 14,	exact(tn,1,ri) = P(zr(ri),nt+1);
			for i = 1:maxOrder, asy(tn,i,ri) = asy_right(nt,zr(ri),alpha,beta,h,psi,i,Dinf,Uright,Uleft,'o'); end	
		case 15,	exact(tn,1,ri) = (P(zr(ri),nt+1)-P(zr(ri)-hh,nt+1))/hh;
			for i = 1:maxOrder, [~, t] = asy_right(nt,zr(ri),alpha,beta,h,psi,i,Dinf,Uright,Uleft,'o',dh,dpsi); 
				asy(tn,i,ri) = t; end	
			
		case 16,	exact(tn,1,ri) = P(zr(ri),nt+1)/gammas(tn);
			for i = 1:maxOrder, asy(tn,i,ri) = asy_left(nt,zr(ri),alpha,beta,h,psi,i,Dinf,Uright,Uleft,'m'); end		
		case 17,	exact(tn,1,ri) = (P(zr(ri),nt+1)-P(zr(ri)-hh,nt+1))/hh/gammas(tn);
			for i = 1:maxOrder, [~, t] = asy_left(nt,zr(ri),alpha,beta,h,psi,i,Dinf,Uright,Uleft,'m',dh,dpsi); 
				asy(tn,i,ri) = t; end	
		case 18,	exact(tn,1,ri) = P(zr(ri),nt+1);
			for i = 1:maxOrder, asy(tn,i,ri) = asy_left(nt,zr(ri),alpha,beta,h,psi,i,Dinf,Uright,Uleft,'o'); end	
		case 19,	exact(tn,1,ri) = (P(zr(ri),nt+1)-P(zr(ri)-hh,nt+1))/hh;
			for i = 1:maxOrder, [~, t] = asy_left(nt,zr(ri),alpha,beta,h,psi,i,Dinf,Uright,Uleft,'o',dh,dpsi); 
				asy(tn,i,ri) = t; end
				
		case 20,	exact(tn,1,ri) = P(zr(ri),nt+1)/gammas(tn);
			for i = 1:maxOrder,
				asy(tn,i,ri) = asy_right(nt,zr(ri),alpha,beta,h,psi,i,Dinf,Uright,Uleft,'mQ',dh,dpsi,Qright,th(ri),contpsi); end	
		case 21,	exact(tn,1,ri) = (P(zr(ri),nt+1)-P(zr(ri)-hh,nt+1))/hh/gammas(tn);
			for i = 1:maxOrder
				[~, t] = asy_right(nt,zr(ri),alpha,beta,h,psi,i,Dinf,Uright,Uleft,'mQ',dh,dpsi,Qright,th(ri),contpsi); 
				asy(tn,i,ri) = t; 
			end	
		case 22,	exact(tn,1,ri) = P(zr(ri),nt+1);
			for i = 1:maxOrder,
				asy(tn,i,ri) = asy_right(nt,zr(ri),alpha,beta,h,psi,i,Dinf,Uright,Uleft,'oQ',dh,dpsi,Qright,th(ri),contpsi); end	
		case 23,	exact(tn,1,ri) = (P(zr(ri),nt+1)-P(zr(ri)-hh,nt+1))/hh;
			for i = 1:maxOrder
				[~, t] = asy_right(nt,zr(ri),alpha,beta,h,psi,i,Dinf,Uright,Uleft,'oQ',dh,dpsi,Qright,th(ri),contpsi); 
				asy(tn,i,ri) = t; 
			end	
			
		case 24,	exact(tn,1,ri) = P(zr(ri),nt+1)/gammas(tn);
			for i = 1:maxOrder,
				asy(tn,i,ri) = asy_left(nt,zr(ri),alpha,beta,h,psi,i,Dinf,Uright,Uleft,'mQ',dh,dpsi,Qleft,th(ri),contpsi); end	
		case 25,	exact(tn,1,ri) = (P(zr(ri),nt+1)-P(zr(ri)-hh,nt+1))/hh/gammas(tn);
			for i = 1:maxOrder
				[~, t] = asy_left(nt,zr(ri),alpha,beta,h,psi,i,Dinf,Uright,Uleft,'mQ',dh,dpsi,Qleft,th(ri),contpsi); 
				asy(tn,i,ri) = t; 
			end	
		case 26,	exact(tn,1,ri) = P(zr(ri),nt+1);
			for i = 1:maxOrder,
				asy(tn,i,ri) = asy_left(nt,zr(ri),alpha,beta,h,psi,i,Dinf,Uright,Uleft,'oQ',dh,dpsi,Qleft,th(ri),contpsi); end	
		case 27,	exact(tn,1,ri) = (P(zr(ri),nt+1)-P(zr(ri)-hh,nt+1))/hh;
			for i = 1:maxOrder
				[~, t] = asy_left(nt,zr(ri),alpha,beta,h,psi,i,Dinf,Uright,Uleft,'oQ',dh,dpsi,Qleft,th(ri),contpsi); 
				asy(tn,i,ri) = t; 
			end	
	end, end
	
end
% Testing order of convergence by plotting the errors in the (3 randomly) chosen regions,
% take into account the remarks in the article when interpreting these
for ri = rs
	if ri <= 3, plotConv(exact(:,1,ri), asy(:,:,ri), legs{ri});
	elseif ri <= 19, plotConv(exact(:,1,ri), asy(:,:,ri), [legs{ri} ' at z= ' num2str(zr(ri))]);
	else plotConv(exact(:,1,ri), asy(:,:,ri), [legs{ri} ' at theta= ' num2str(th(ri))]);
end, end

%% Heuristics: see which expansion is more accurate where.
clear variables;
alpha = rand, beta = -1+3*rand, h = @(x) 1;
n = 100;
maxOrder = 3;

[P,gammaP,alphaP,betaP] = exactPolys(alpha,beta,h,n);
% Computing relative errors on nrP points on this line
zs = linspace(-1,1,25) +0.02i;

[c, d, Dinf, psi, ~, contpsi] = contour_integrals(alpha,beta,h,2*(maxOrder+3) );
[Uright,Uleft,Qright,Qleft] = UQ(alpha,beta,Dinf,c,d,maxOrder,'UQV');

relerrs = zeros(length(zs),maxOrder,6);
warning('off','all');
for zi = 1:length(zs)
	exact = P(zs(zi),n+1);
	for i = 1:maxOrder
		relerrs(zi,i,1) =  abs(asy_lens( n,zs(zi),alpha,beta,h,psi,i,Dinf,Uright,Uleft,'o')-exact)/abs(exact);
		relerrs(zi,i,2) =  abs(asy_outer(n,zs(zi),alpha,beta,h,psi,i,Dinf,Uright,Uleft,'o')-exact)/abs(exact);
		relerrs(zi,i,3) =  abs(asy_right(n,zs(zi),alpha,beta,h,psi,i,Dinf,Uright,Uleft,'o')-exact)/abs(exact);
		relerrs(zi,i,4) =  abs(asy_left( n,zs(zi),alpha,beta,h,psi,i,Dinf, Uright,Uleft,'o')-exact)/abs(exact);
		relerrs(zi,i,5) = abs(asy_right(n,zs(zi),alpha,beta,h,psi,i,Dinf, Uright,Uleft,'oQ',[],...
			[],Qright,acos(zs(zi)),contpsi)-exact)/abs(exact);
		relerrs(zi,i,6) = abs(asy_left(n,zs(zi),alpha,beta,h,psi,i,Dinf, Uright,Uleft,'oQ',[],...
			[],Qleft,acos(-zs(zi)),contpsi)-exact)/abs(exact);
	end
end
warning('on','all');

figure;
semilogy(real(zs),relerrs(:,:,4),'*');
legend('1 term','2 terms','3 terms');
ylabel('Relative error');
xlabel(['Real part of z, imaginary part is ' num2str(mean(imag(zs)))]);
title('Relative errors of expansion in the left disk for varying z and number of terms')

figure;
T = 3;
semilogy(real(zs),reshape(relerrs(:,T,:),[length(zs),6]),'*');
legend('lens','outer','right','left','right Q','left Q');
ylabel('Relative error');
xlabel(['Real part of z, imaginary part is ' num2str(mean(imag(zs)))]);
title(['Relative errors of different expansion for varying z and T=' num2str(T)])


%% Gaussian quadrature: Newton procedure to compute nodes and weights
clear variables;
cheb = randi([0 1],1);
if cheb
	alpha = -1/2; beta = -1/2; h = @(x) 1; dh = @(x) 0;
else
	alpha = 1/sqrt(3); beta = -1/pi; h = @(x) exp(x); dh = @(x) exp(x);
end
n = 200;
nrT = 3;

[c, d, Dinf, psi, dpsi, contpsi] = contour_integrals(alpha,beta,h,nrT);
[Uright,Uleft] = UQ(alpha,beta,Dinf,c,d,nrT,'UW');

zw = zeros(n,3);
for ni = 1:n
	if ni <=3 % (Approximation of) zero of bessel function to approximate first zeros
		zw(ni,1) = -1 + (2*sqrt(beta+1) + (ni-1)*pi)^2/2/n^2;
		if zw(ni,1) +1 < eps^(1/3)
			warning('Use Q-s');
		end
	else
		zw(ni,1) = 3*zw(ni-1) -3*zw(ni-2) +zw(ni-3); % Quadratic extrapolation
	end
	% Could improve the following by stopping Newton iterates when already
	% tried that zero and then going over all (max. 5) closest floats.
	for attempt = 1:10
		if zw(ni,1)+1 < 0.2 % This bound is not general for all weights!
			[poly, dpoly] = asy_left(n,zw(ni,1),alpha,beta,h,psi,nrT,Dinf,Uright,Uleft,'o',dh,dpsi);
			polyM1 = asy_left(n-1,zw(ni,1),alpha,beta,h,psi,nrT,Dinf,Uright,Uleft,'o');
			polyP1 = asy_left(n+1,zw(ni,1),alpha,beta,h,psi,nrT,Dinf,Uright,Uleft,'o');
		elseif 1-zw(ni,1) < 0.2
			[poly, dpoly] = asy_right(n,zw(ni,1),alpha,beta,h,psi,nrT,Dinf,Uright,Uleft,'o',dh,dpsi);
			polyM1 = asy_right(n-1,zw(ni,1),alpha,beta,h,psi,nrT,Dinf,Uright,Uleft,'o');
			polyP1 = asy_right(n+1,zw(ni,1),alpha,beta,h,psi,nrT,Dinf,Uright,Uleft,'o');
		else
			[poly, dpoly] = asy_lens(n,zw(ni,1),alpha,beta,h,psi,nrT,Dinf, Uright,Uleft,'o',dh,dpsi);
			polyM1 = asy_lens(n-1,zw(ni,1),alpha,beta,h,psi,nrT,Dinf, Uright,Uleft,'o');
			polyP1 = asy_lens(n+1,zw(ni,1),alpha,beta,h,psi,nrT,Dinf, Uright,Uleft,'o');
		end
		zw(ni,1) = zw(ni,1) - real(poly)/real(dpoly); % We know the polynomial is real on the interval
	end
	if (zw(ni,1) < -1) || (zw(ni,1) > 1)
		error('Outside of the interval');
	elseif (ni > 1) && (zw(ni,1) < zw(ni-1,1) )
		error('Not going right');
	end
	zw(ni,2) = gamman(n,nrT,Dinf,Uright,Uleft)/real(dpoly)/gamman(n-1,nrT,Dinf,Uright,Uleft)/real(polyM1);
	zw(ni,3) = -gamman(n+1,nrT,Dinf,Uright,Uleft)/real(dpoly)/gamman(n,nrT,Dinf,Uright,Uleft)/real(polyP1);
end
errzw = norm(zw(:,2)-zw(:,3))/norm(zw(:,2)) % Formulas should give the same result
% Test by integrating (x^2-1) and (x^3) over the interval
if cheb % Weight 1/sqrt(1-x^2)
	chebyroots = cos((2*(1:n).'-1)/2/n*pi);
	errorOnRoots = norm(zw(:,1)+chebyroots)
	abserror = sum(zw(:,2).*(zw(:,1).^2-1) ) - (-pi/2)
	abserror = sum(zw(:,2).*(zw(:,1).^3) )
else % w(x) = (1-x)^( 1/sqrt(3))*(1+x)^(-1/pi)*exp(x)
	abserror = sum(zw(:,2).*(zw(:,1).^2-1) ) -(-1.293225037)
	abserror = sum(zw(:,2).*(zw(:,1).^3) ) - (-.1659585353)
end
% Integrate a random polynomial with a degree low enough to be integrated exactly
r = [randn(2,1); randi([0 2*n-3],1); randi([0 5],1)];
fct = @(x) 7*r(1)*x.^r(3) -23*r(2)*x.^r(4);
w = @(x) (1-x).^alpha.*(1+x).^beta.*h(x);
abserror = sum(zw(:,2).*fct(zw(:,1)) ) -integral(@(x) w(x).*fct(x),-1,1)

