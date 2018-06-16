% Real orthonormal polynomials using the recurrence relation and the OPQ-library 
% of Walter Gautschi. These polynomials exhibit degrading accuracy from n about 
% 128 to about a relative error of 1e-8 around n=2^12, which can be checked by 
% evaluating exactPolys(0,0,@(x) 1,4096): then P(1,n+1) should be sqrt((2*n+1)/2).
% Input
%   alpha,beta, h- Parts of the weight w(x) = (1-x).^alpha.*(1+x).^beta.*h(x)
%   N            - Maximum degree
% Output 
%   P            - P(x,n+1) giving orthonormal poly of degree n at point x 
%                  as a function, ||pi(x,n)||_w = 1
%   gammaP       - gammaP(n+1)*2^n = \gamma_n = is the leading order coefficient of 
%				   the orthonormal polynomial of degree n. Multiplication by
%				   2^(-n) to avoid overflow of saved value gammaP(n+1)
%   alphaP,betaP - alphaP(n+1) & betaP(n+1) are the recurrence coefficients 
%                  of the monic orthogonal polynomial of degree n
% Note: This uses another algorithm than the Julia version.
% About
%   Author       - Peter Opsomer (Peter.Opsomer@cs.kuleuven.be)
%   History      - Created October 2013, last edit February 2015
function [P,gammaP,alphaP,betaP] = exactPolys(alpha,beta,h,N)
% Use alpha and beta to construct quadrature rules for computing the
% real alphaP & betaP. ab(:,1) != alphaP because doesn't use h(x)
ab = r_jacobi(3*N,alpha,beta);
zw = gauss(3*N,ab);
inwprod = @(f,g) sum(zw(:,2).*h(zw(:,1) ).*f(zw(:,1) ).*conj(g(zw(:,1) ) ) );
alphaP = zeros(N+1,1);
betaP = ones(N+1,1); % contains beta(k+1) = beta_k and not sqrt(beta_k)
gammaP = zeros(N+2,1); % contains gamma_k/2^k = gamma(k+1)/2^k
betaP(1) = inwprod(@(x) ones(size(x) ), @(x) ones(size(x) ) );
gammaP(1) = 1/sqrt(betaP(1)); % = \gamma_0/2^0 = 1/sqrt(\beta_0)
tic
prevToc = toc;
for k = 1:N
	if (toc-prevToc > 2)
		prevToc = toc;
		display(['OPQ iteration ', num2str(k), ', ' num2str( (k-1)/N,'%7.3f'), ' % done, estim sec left = ' num2str(toc*(N-k+1)/(k-1)) ])
	end
	alphaP(k) = inwprod(@(x) orthonorm(x,k-1,alphaP,sqrt(betaP)).*x, @(x) orthonorm(x,k-1,alphaP,sqrt(betaP)) );
	betaP(k+1) =  inwprod(@(x) orthonorm(x,k,alphaP,sqrt(betaP)), @(x) orthonorm(x,k,alphaP,sqrt(betaP)) );
	gammaP(k+1) = gammaP(k)/sqrt(betaP(k+1) )/2;
end
% For when testing convergence of the asymptotic expansion of alpha_n :
alphaP(N+1) = inwprod(@(x) orthonorm(x,N,alphaP,sqrt(betaP)).*x, @(x) orthonorm(x,N,alphaP,sqrt(betaP)) );
P = @(x,n1) orthonorm(x,n1-1,alphaP,sqrt(betaP) );
end

function p = orthonorm(x,n,alphaP,sqbetaP)
p = [zeros(length(x),1), ones(length(x),1)./sqbetaP(1)];
for j=1:n
	p = [p(:,2) ((x-alphaP(j) ).*p(:,2)-sqbetaP(j).*p(:,1))./sqbetaP(j+1)];
end
p(:,1) = [];
end
