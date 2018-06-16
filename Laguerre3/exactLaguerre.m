% Real orthonormal polynomials using the recurrence relation. For q=[0 1], P(x,n+1) is equal to 
% laguerreL(n,alpha,x)*(-1)^n*sqrt(gamma(n+alpha+1)/factorial(n)). For q = 'Hermite', P(x,n+1) is 
% equal to the normalised Hermite polynomial of degree n at x. When N < 0, we compute the normalised
% generalized Hermite polynomials with respect to x^alpha*exp(-sum q_k x^(2k) ) via its recurrence
% coefficients.
% Input
%   alpha, q     - Parts of the weight w(x) =  x.^alpha.*exp(-sum(x.^(0:length(q)-1)*q ) or -q(x) )
%   N            - Maximal degree or -maximal degree for generalized Hermite.
% Output 
%   P            - P(x,n+1) gives the orthonormal poly of degree n at point x 
%                  as a function, ||P(x,n)||_w = 1
%   gammaP       - gammaP(n+1) = \gamma_n is the leading order coefficient of 
%				   the orthonormal polynomial of degree n.
%   aP,bP        - aP(n+1) & bP(n+1) are the recurrence coefficients 
%                  of the orthonormal orthogonal polynomial of degree n
% About
%   Author       - Peter Opsomer (Peter.Opsomer@cs.kuleuven.be)
%   History      - Created July 2015
function [P,gammaP,aP,bP] = exactLaguerre(alpha,q,N)

aP = zeros(N+1,1);
bP = ones(N+1,1); % contains beta(k+1) = beta_k and not sqrt(beta_k)
sel = @(ix,y) y(ix);
replnan = @(x) sel(isnan(x)+(1:2:2*length(x)),[x;zeros(size(x))]);
printtoc = 5;

if N < 0 % Hermite
    N = -N;
    % Could also integrate to x where Q(x) >= 750 iso to Inf to avoid e^(-Q) = realmin*eps giving NaN
    if isnumeric(q)
        inwprod = @(f,g) quadgk(@(x) replnan(x.^alpha.*exp(-q*repmat(x, length(q), 1).^repmat(2*(0:length(q)-1)', 1, ...
            length(x))).*f(x).*g(x)), -Inf, Inf); % 'quadgk' is more accurate than 'integral'.
    else
        inwprod = @(f,g) quadgk(@(x) replnan(x.^alpha.*exp(-q(x) ).*f(x).*g(x)), -Inf, Inf);
    end
    bP(1) = sqrt(inwprod(@(x) ones(size(x) ), @(x) ones(size(x) ) ) );
    gammaP = ones(N+1,1);
    tic;
    prevToc = toc;
    % Use the recurrence relation
    for k = 1:N
        aP(k) = inwprod(@(x) orthonorm(x,k-1,aP,bP).*x, @(x) orthonorm(x,k-1,aP,bP) );
        bP(k+1) =  sqrt(inwprod(@(x) orthonorm(x,k-1,aP,bP).*x, @(x) (x-aP(k)).*orthonorm(x,k-1,aP,bP) -bP(k)*orthonorm(x,k-2,aP,bP) ) );
        gammaP(k+1) = 1/prod(bP(1:k+1));
        if (toc-prevToc > printtoc)
            prevToc = toc; % Iteration k requires O(k) operations through integrals
            esl = prevToc*(sum((1:N).^1)-sum((1:k).^1))/sum((1:k).^1);
            display(['Iteration ' num2str(k) ' out of ' num2str(N) ' computing Hermite-type rec. coeff., est. sec. left = ' num2str(esl) ]);
        end
    end
    P = @(x,n1) orthonorm(x,n1-1,aP,bP);
    return
elseif isnumeric(q) && (length(q) == 2) && (q(1) == 0) && (q(2) == 1) 
    % Use the explicit expression for standard associated Laguerre orthonormal polynomials.
    bP(1) = sqrt(gamma(alpha+1)/1);
    aP(1) = 1+alpha;
    for n=1:N+1
        bP(n+1) = sqrt(n)*sqrt(alpha+n);
        aP(n+1) = 2*n+alpha+1;
    end
    P = @(x,n1) orthonorm(x,n1-1,aP,bP);
	gammaP = cumprod(1./bP);
	return
elseif strcmp(q,'Hermite')
    aP(:) = 2;
    bP = 2*(0:N);
    bP(1) = 1;
    gammaP = pi^(1/4)*sqrt(2.^(0:N).*factorial(0:N) )';
    P = @(x,n1,norm) herm(x,n1-1,aP,bP, norm);
    return
end

% Could also integrate to x where Q(x) >= 750 iso to Inf to avoid e^(-Q) = realmin*eps giving NaN
if isnumeric(q)
    inwprod = @(f,g) quadgk(@(x) replnan(x.^alpha.*exp(-q*repmat(x, length(q), 1).^repmat((0:length(q)-1)', 1, ...
        length(x))).*f(x).*g(x)), 0, Inf); % 'quadgk' is more accurate than 'integral'.
else
    inwprod = @(f,g) quadgk(@(x) replnan(x.^alpha.*exp(-q(x) ).*f(x).*g(x)), 0, Inf);
end
bP(1) = sqrt(inwprod(@(x) ones(size(x) ), @(x) ones(size(x) ) ) );
gammaP = ones(N+1,1);
gammaP(1) = 1/bP(1);
tic;
prevToc = toc;
% Use the recurrence relation
for k = 1:N
    aP(k) = inwprod(@(x) orthonorm(x,k-1,aP,bP).*x, @(x) orthonorm(x,k-1,aP,bP) );
    bP(k+1) =  sqrt(inwprod(@(x) orthonorm(x,k-1,aP,bP).*x, @(x) (x-aP(k)).*orthonorm(x,k-1,aP,bP) -bP(k)*orthonorm(x,k-2,aP,bP) ) );
    gammaP(k+1) = 1/prod(bP(1:k+1));
    if (toc-prevToc > printtoc)
        prevToc = toc; % Iteration k requires O(k) operations through integrals
        esl = prevToc*(sum((1:N).^1)-sum((1:k).^1))/sum((1:k).^1);
        display(['Iteration ' num2str(k) ' out of ' num2str(N) ' computing Laguerre-type rec. coeff., est. sec. left = ' num2str(esl) ]);
    end
end
aP(N+1) = inwprod(@(x) orthonorm(x,N,aP,bP).*x, @(x) orthonorm(x,N,aP,bP) );
P = @(x,n1) orthonorm(x,n1-1,aP,bP);

end


function p = herm(x,n,aP,bP, norm)
if n < 0
    p = zeros(size(x));
    return
end
flip = size(x,1) == 1;
if flip
    x = transpose(x);
end
if norm
    p = [zeros(length(x),1), 1/pi^(1/4)*ones(length(x),1)];
    for j=1:n
        p = [p(:,2) (x*2/sqrt(2*j).*p(:,2) -(j-1)/sqrt(j*max(1,j-1)).*p(:,1))];
    end
else
    p = [zeros(length(x),1), ones(length(x),1)];
    for j=1:n
        p = [p(:,2) (x*aP(j).*p(:,2) -bP(j).*p(:,1))];
    end
end
p(:,1) = [];
p(isinf(p)) = realmax;
if flip
    p = transpose(p);
end
end
