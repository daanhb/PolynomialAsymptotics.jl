% Compute the contour integrals c_n, d_n, D_\infty, psi(z) [, its derivative and its contour integral].
% Input
%   alpha, beta - Parts of the weight w(x) = (1-x).^alpha.*(1+x).^beta.*h(x)
%   h           - Anonymous function for analytic function h(x)
%   nrT         - Number of terms, where 1 gives leading order term
%   [rho        - Berstein ellipse parameter: small to not encircle where log(h) is analytic and large to encircle the point x where to 
%				    evaluate psi(x). Choose it to optimise the conditioning number as in [F. Bornemann, 2011] & [Wang & Huybrechs 2015], 
%					so it is a vector with rho(n+2) optimal for log(h)/sqrt(zeta-1)/sqrt(zeta+1)/(z-zeta)^(n+1) with rho+1/rho going to 
%					twice the distance of where log h is not analytic to the interval, or it is a scalar which is acceptable for all n.]
%   [M          - Number of points in trapezoidal rules, with M(n+2) optimal for c&d(n+1) and increasing with n, or M a scalar]
% Output
%   c, d        - The c- and d-coefficients, size(c&d) == ceil(nrT/2)
%   Dinf        - Limit of the Szego function
%   psi         - Anonymous function for the phase function
%   [dpsi       - Derivative of the phase function]
%   [contpsi	- Contour integral in psi, contpsi(x,p) has (zeta-x)^p]
% About
%   Author       - Peter Opsomer (Peter.Opsomer@cs.kuleuven.be)
%   History      - Created October 2013, last edit October 2015
function [c, d, Dinf, psi, dpsi,contpsi] = contour_integrals(alpha,beta,h,nrT,rho,M)
c = zeros(ceil(nrT/2),1);
d = zeros(ceil(nrT/2),1);

ra = 2*rand-1; 
mm = 2; cc = 7; 
% Checking for known exact solutions by comparing value at one random z =
% ra, could also do if sum(chebfun(@(x) abs(h(x) -ones(size(x) ) ), domain(xc) ) ) == 0
if h(ra) == 1
    display('Using exact result for non-generalised weight in contour_integrals')
    Dinf = 2^(-alpha/2-beta/2);
    psi = @(x) 1/2.*( (alpha + beta).*acos(x)-alpha*pi);
    dpsi = @(x) 1/2.*(alpha+beta).*(-1./(sqrt(1-x).*sqrt(1+x) ) );
	contpsi = @(x,p) 0;
elseif abs(h(ra) - exp(-cc*ra.^(2*mm)) ) <100*eps
    display('Using exact result for Gibbs weight in contour_integrals');
	binom = @(x,n) prod(x:-1:(x-n+1) )/prod(1:n);
	% Could go until length(c)-1, but need exact result for psi, so:
    for n = 0:(2*mm-1)
		c(n+1) = 0;
        for j = 0:floor( (2*mm-n-1)/2)
            c(n+1) = c(n+1) -cc*binom(j-1/2,j)*binom(2*mm-1-2*j,2*mm-n-1-2*j);
        end
        d(n+1) = (-1)^(n+1)*c(n+1);
    end
    Dinf = 2^(-alpha/2-beta/2)*exp(-cc/2*binom(mm-1/2,mm));
    degs = (0:(length(c)-1) ).';
    psi = @(x) 1/2.*( (alpha + beta).*acos(x)-alpha*pi) +sqrt(1-x).*sqrt(1+x)/2.*sum( c.*(x-1).^(degs) );
    dpsi = @(x) 1/2.*(alpha+beta).*(-1./(sqrt(1-x).*sqrt(1+x) ) ) + ...
        (-x)./sqrt(1-x)./sqrt(1+x)/2.*sum(c.*(x-1).^(degs) ) + ...
        sqrt(1-x).*sqrt(1+x)/2.*sum(degs(2:end).*c(2:end).*(x-1).^(degs(2:end)-1) );
	contpsi = @(x,p) 2i*pi*sum( (degs(p:end)).^(p-1).*c(p:end).*(x-1).^(degs(p:end)-p+1) );
else % Trapezoidal rules
	if abs(h(ra) - 1./sqrt(ra+3)) < 100*eps
		display('Using trapezoidal rules for Fourier weight in contour_integrals');
		rho = 4; % M = 80; % But determine M automatically later
	elseif ~exist('rho','var') % Take rho smaller to reduce risk of encircling pole log h although not optimal:
		rho = 1.25; % > 1 so positive direction
		% M will also be determined by successive doubling
	end
	lr = length(rho);
	if exist('M','var') % Avoid successive doubling
		lm = length(M);
		for nr = 1:ceil(nrT/2) % nr = n+1 in c_n = c(nr)
			c(nr) = trap_rule(@(z) log(h(z ) )./( sqrt(z -1).*sqrt(z +1).*(z -1).^nr ), rho(min(nr+1,lr)),M(min(nr+1, lm)));
			d(nr) = trap_rule(@(z) log(h(z ) )./( sqrt(z -1).*sqrt(z +1).*(z +1).^nr ), rho(min(nr+1,lr)),M(min(nr+1, lm)));
		end
		% For Dinf use optimal values for log(h)/sqrt(zeta-1)/sqrt(zeta+1)
		Dinf = 2.^( -(alpha+beta)./2).*exp(real(trap_rule(@(z) log(h(z) )./sqrt(z -1)./sqrt(z +1), rho(1),M(1))/2) );
		% trap_rule already devides by 2i*pi, use optimal value for log(h)/sqrt(zeta-1)/sqrt(zeta+1)/(z-zeta)^{1 resp. 2}
		dpsi = @(x) 1/2.*(alpha + beta).*(-1)./(sqrt(1-x).*sqrt(1 +x)) + ...
			(-x)/(sqrt(1-x).*sqrt(1+x)*2)*trap_rule( @(z) log(h(z) )./sqrt(z-1)./sqrt(z+1)./(z-x),rho(min(2,lr)),...
			M(min(3,lm))) + sqrt(1-x).*sqrt(1+x)/2*trap_rule(@(z) log(h(z) )./sqrt(z-1)./sqrt(z+1)./(z-x).^2,...
			rho(min(3,lr)),M(min(3,lm)));
		% For (cont)psi use optimal values for log(h)/sqrt(zeta-1)/sqrt(zeta+1)/(z-zeta)
		psi = @(x) 1/2.*( (alpha + beta).*acos(x)-alpha*pi) + sqrt(1-x).*sqrt(1+x) ...
			.*trap_rule(@(z) log(h(z) )./sqrt(z-1)./sqrt(z+1)./(z-x), rho(min(2,lr)), M(min(2,lm)))/2;
		contpsi = @(x,p) 2i*pi*trap_rule(@(z) log(h(z) )./sqrt(z-1)./sqrt(z+1)./(z-x).^p, rho(min(2,lr)), M(min(2,lm)));
	else
		for nr = 1:ceil(nrT/2) % nr = n+1 in c_n = c(nr)
			c(nr) = trap_rule(@(z) log(h(z ) )./( sqrt(z -1).*sqrt(z +1).*(z -1).^nr ), rho(min(nr+1,lr)) );
			d(nr) = trap_rule(@(z) log(h(z ) )./( sqrt(z -1).*sqrt(z +1).*(z +1).^nr ), rho(min(nr+1,lr)) );
		end
		Dinf = 2.^( -(alpha+beta)./2).*exp(real(trap_rule(@(z) log(h(z) )./sqrt(z -1)./sqrt(z +1), rho(1))/2) );
		dpsi = @(x) 1/2.*(alpha + beta).*(-1)./(sqrt(1-x).*sqrt(1 +x)) + ...
			(-x)/(sqrt(1-x).*sqrt(1+x)*2)*trap_rule( @(z) log(h(z) )./sqrt(z-1)./sqrt(z+1)./(z-x),rho(min(2,lr)) )...
			 + sqrt(1-x).*sqrt(1+x)/2*trap_rule(@(z) log(h(z) )./sqrt(z-1)./sqrt(z+1)./(z-x).^2,rho(min(3,lr)) );
		psi = @(x) 1/2.*( (alpha + beta).*acos(x)-alpha*pi) + sqrt(1-x).*sqrt(1+x) ...
			.*trap_rule(@(z) log(h(z) )./sqrt(z-1)./sqrt(z+1)./(z-x), rho(min(2,lr)))/2;
		contpsi = @(x,p) 2i*pi*trap_rule(@(z) log(h(z) )./sqrt(z-1)./sqrt(z+1)./(z-x).^p, rho(min(2,lr)) );
	end
end
