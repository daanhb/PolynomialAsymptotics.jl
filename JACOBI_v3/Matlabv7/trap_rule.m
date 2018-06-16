% Compute the contour integral of the given function using the trapezoidal rule around the interval.
% Input
%   f       - the integrand
%   rho     - parameter for the Bernstein ellipse
%   [M       - number of terms, use successive doubling from [F. Bornemann 2011] if not specified]
% Output
%   cont    - the contour integral
% About
%   Author  - Peter Opsomer (Peter.Opsomer@cs.kuleuven.be)
%   History - Created October 2013, last edit October 2015
function cont = trap_rule(f,rho,M)
if ~exist('M','var')
	M = 8; tol = 1e-15;
	thi = @(k,l) 2i*pi*transpose(1:l:k)/k;
	bernst = @(k,l) rho/2*exp(thi(k,l)) + 1/rho/2*exp(-thi(k,l));
	s = f(bernst(M,1)).*(rho*1i/2*exp(thi(M,1)) -1i/rho/2*exp(-thi(M,1)) )/1i; cont = mean(s); err1 = NaN;
	while M < 1e6
		M = 2*M;
		s = [s; f(bernst(M,2)).*(rho*1i/2*exp(thi(M,2)) -1i/rho/2*exp(-thi(M,2)) )/1i];
		val = mean(s); kappa = mean(abs(s))/abs(val);
		err0 = abs(val-cont)/abs(val); 
		if (err0^3/err1^2 <= kappa*tol) || (~isfinite(kappa)), break; end
		cont = val; err1 = err0;
	end
elseif 1
	ths = linspace(0,2*pi,M+1); % = theta_k, k = 0 to M-1
	ths(end) = [];
	hh = ths(2)-ths(1);
	zs = rho/2*exp(1i*ths) + 1/rho/2*exp(-1i*ths);
	cont = hh*sum(f(zs).*(rho*1i/2*exp(1i*ths) -1i/rho/2*exp(-1i*ths) ) )/(2*pi*1i );
else
% Below an ellipse, but could also need other parametrisations, for example
% when evaluating psi(x) for x further from the interval than the point
% where log(h) is not analytic
    r1 = 1.8;
    r2 = 0.8;
    zs = r1*cos(ths) +r2*1i*sin(ths);
    cont = hh*sum(f(zs).*(-r1*sin(ths) +r2*1i*cos(ths) ) )/(2*pi*1i );
end
if rho < 1 % Bernstein in negative direction but need positive direction
    cont = -cont;
end
