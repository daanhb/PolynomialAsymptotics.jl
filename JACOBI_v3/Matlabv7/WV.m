% Compute the W- or V-matrices to construct the asymptotic expansion of R 
% using the procedure with the convolutions as explained in the paper.
% Input
%   q            - alpha (when computing Wright) or beta (when computing Wleft)
%   t            - beta or alpha
%   Dinf         - Limit of the Szego function
%   cd           - c-or d- Coefficients needed for the U-matrices
%   maxOrder     - The maximum order of the error
%   r            - 1 when computing Wright, -1 when computing Wleft
%   isW          - 1 if the result are the Ws-, 0 if the V-s
%   [mos         - Maximum orders for each k, default round((maxOrder-2)/2) for all k]
% Output
%   WV           - Coefficient matrices for (z \pm 1)^m of s_k(z) or Delta_k(z)
% About
%   Author       - Peter Opsomer (Peter.Opsomer@cs.kuleuven.be)
%   History      - Created May 2014, last edit January 2015
function WV = WV(q,t,Dinf,cd,maxOrder,r,isW,mos)

if exist('mos')
    mo = max(mos);
else
    mo = round( (maxOrder-2)/2);
end
ns = 0:mo;

poch = @(x,n) prod(x+(0:(n-1)) ); 
% 'pochhammer' needs symbolic toolbox, could also use separate 1-line
% scripts or global functions for these
binom = @(x,n) prod(x:-1:(x-n+1) )/prod(1:n);
brac = @(k) ((k == 0) + (k~= 0)*prod(4*q^2 -(2*(1:k)-1).^2)/(2^(2*k)*(factorial(k) ) ) );

f = zeros(mo+1,1); % Coefficients in the expansion of \varphi(z)
for n = ns
    f(n+1) = poch(1/2,n)./(-r*2).^n./(1+2*n)./factorial(n);
end
% Or: f = poch(1/2,ns)./(-r*2).^ns./(1+2*ns)./factorial(ns);
g = zeros(maxOrder-1,mo+1);
g(1,1) = 1; % Coefficients in the expansion of \varphi(z)^{-i}
for n = 1:mo
    g(1,n+1) = -g(1,1:(n+1) )*f((n+1):-1:1);
end

st = zeros(maxOrder-1,mo+1); % = \rho_{i,n,alpha+beta}
stp = zeros(maxOrder-1,mo+1); % = \rho_{i,n,alpha+beta+1}
stm = zeros(maxOrder-1,mo+1); % = \rho_{i,n,alpha+beta-1}
for n = ns
    tmp = 0;
    for j = 0:n
       tmp = tmp + binom(1/2,j)/(r*2)^j*cd(n-j+1);
    end
    st(1,n+1) = r*tmp + f(n+1)*(q+t);
    stp(1,n+1) = r*tmp + f(n+1)*(q+t+1);
    stm(1,n+1) = r*tmp + f(n+1)*(q+t-1);
end

for i = 2:max(mo*2+1,(maxOrder-1))
    for n = ns
        g(i,n+1) = sum(g(i-1,1:(n+1) ).*g(1,(n+1):-1:1) );
        st(i,n+1) = sum(st(i-1,1:(n+1) ).*st(1,(n+1):-1:1) );
        stp(i,n+1) = sum(stp(i-1,1:(n+1) ).*stp(1,(n+1):-1:1) );
        stm(i,n+1) = sum(stm(i-1,1:(n+1) ).*stm(1,(n+1):-1:1) );
    end
end

OmOdd = zeros(mo+1,1); % = H_{n,alpha+beta}^{odd}
OmEven = zeros(mo+1,1); % = H_{n,alpha+beta}^{even}
XiOdd = zeros(mo+1,1); % = H_{n,alpha+beta+1}^{odd}
XiEven = zeros(mo+1,1); % = H_{n,alpha+beta+1}^{even}
ThOdd = zeros(mo+1,1); % = H_{n,alpha+beta-1}^{odd}
ThEven = zeros(mo+1,1); % = H_{n,alpha+beta-1}^{even}
for n = ns
    if ( (mod(maxOrder,2)==0) || (n ~= mo) )
        OmEven(n+1) = OmEven(n+1) + st(1,n+1);
        XiEven(n+1) = XiEven(n+1) + stp(1,n+1);
        ThEven(n+1) = ThEven(n+1) + stm(1,n+1);
        for j = 1:n
            OmOdd(n+1) = OmOdd(n+1) + (r*2)^j*st(2*j,n-j+1)/factorial(2*j);
            XiOdd(n+1) = XiOdd(n+1) + (r*2)^j*stp(2*j,n-j+1)/factorial(2*j);
            ThOdd(n+1) = ThOdd(n+1) + (r*2)^j*stm(2*j,n-j+1)/factorial(2*j);
            OmEven(n+1) = OmEven(n+1) + (r*2)^j*st(2*j+1,n-j+1)/factorial(2*j+1);
            XiEven(n+1) = XiEven(n+1) + (r*2)^j*stp(2*j+1,n-j+1)/factorial(2*j+1);
            ThEven(n+1) = ThEven(n+1) + (r*2)^j*stm(2*j+1,n-j+1)/factorial(2*j+1);
        end
	else % This splitting is needed to avoid wrong OmEven(n+1) etc when n=mo and maxOrder is odd
        for j = 1:n
            OmOdd(n+1) = OmOdd(n+1) + (r*2)^j*st(2*j,n-j+1)/factorial(2*j);
            XiOdd(n+1) = XiOdd(n+1) + (r*2)^j*stp(2*j,n-j+1)/factorial(2*j);
            ThOdd(n+1) = ThOdd(n+1) + (r*2)^j*stm(2*j,n-j+1)/factorial(2*j);
        end
    end
end
OmOdd(1) = 1;
XiOdd(1) = 1; % Overwrite because else wrong
ThOdd(1) = 1;

OmO = zeros(mo+1,1); % = X_{n,alpha+beta}^{odd}/sqrt(\pm 2)
OmE = zeros(mo+1,1); % = X_{n,alpha+beta}^{even}/sqrt(\pm 2)
XiO = zeros(mo+1,1); % = X_{n,alpha+beta+1}^{odd}/sqrt(\pm 2)
XiE = zeros(mo+1,1); % = X_{n,alpha+beta+1}^{even}
ThO = zeros(mo+1,1); % = X_{n,alpha+beta-1}^{odd}
ThE = zeros(mo+1,1); % = X_{n,alpha+beta-1}^{even}
for n = ns
    for j = 0:n
        OmO(n+1) = OmO(n+1) + binom(-1/2,j)/(r*2)^j*OmOdd(n-j+1)/sqrt(r*2);
        XiO(n+1) = XiO(n+1) + binom(-1/2,j)/(r*2)^j*XiOdd(n-j+1)/sqrt(r*2);
        ThO(n+1) = ThO(n+1) + binom(-1/2,j)/(r*2)^j*ThOdd(n-j+1)/sqrt(r*2);
        if ( (mod(maxOrder,2)==0) || (n ~= mo) )
            OmE(n+1) = OmE(n+1) + binom(-1/2,j)/(r*2)^j*OmEven(n-j+1);
            XiE(n+1) = XiE(n+1) + binom(-1/2,j)/(r*2)^j*XiEven(n-j+1);
            ThE(n+1) = ThE(n+1) + binom(-1/2,j)/(r*2)^j*ThEven(n-j+1);
        end
    end
end

Ts = zeros(2,2,mo+1);  % = T_{k,n}^{odd} or T_{k,n}^{even} depending on k, overwritten each new k
WV = zeros(2,2,maxOrder-1,mo+1);
for k = 1:(maxOrder-1)
    a = (q^2 + k/2-1/4)/k;
    b = -1i*r*(k-1/2);
    if exist('mos','var')
        mo = mos(k);
    else
        mo = round( (maxOrder-1-k)/2)-1+round(k/2);
    end
    Ts(:,:,:) = 0;
    if mod(k,2)
% Actually, for m = (-round(k/2)):(ceil((maxOrder-1-k)/2)-1)
% but Matlab does not allow negative indices so shift by round(k/2)+1 and
% easier to let n go to "mo" always iso ceil((maxOrder-1-k)/2)+round(k/2)
        for n = 0:mo
			% Or a*binom(-1/2,n)*r^(n+1/2)/2^(n+1/2)*(2*n+1)/(2*n-1)
			% Or a*(-r*binom(-1/2,n)/(r*2)^(n+1/2) -binom(-1/2,n-1)/(r*2)^(n-1/2))			
			Ts(:,:,n+1) = [ (-a*binom(n-3/2,n)*(-r*2)^(-n)*(2*n+1)*sqrt(r/2) +1i*b*OmO(n+1) ), ...
                Dinf^2*(1i*a*binom(-1/2,n)*(r*2)^(-n)/sqrt(r*2) +r*b*XiO(n+1) );...
                (1i*a*binom(-1/2,n)*(r*2)^(-n)/sqrt(r*2) +r*b*ThO(n+1) )/Dinf^2,...
                (a*binom(n-3/2,n)*(-r*2)^(-n)*(2*n+1)*sqrt(r/2) -1i*b*OmO(n+1) )];

			WV(:,:,k,n+1) = brac(k-1)/(r*2)^(3*k/2)*sum(repmat(...
                reshape(g(k,1:(n+1) ),[1,1,n+1]),[2,2,1]).*Ts(:,:,(n+1):-1:1),3);
        end
    else
        for n = 0:mo
            Ts(:,:,n+1) = [(a*(n==0) +1i*r*b*OmE(n+1)), Dinf^2*b*XiE(n+1) ;  ...
				b*ThE(n+1)/Dinf^2, (a*(n==0) -1i*b*r*OmE(n+1) )];
            WV(:,:,k,n+1) = brac(k-1)/(r*2)^(3*k/2)*sum(repmat(reshape(...
                g(k,1:(n+1) ),[1,1,n+1]),[2,2,1]).*Ts(:,:,(n+1):-1:1),3);
            if isW % The extra matrix of the theorem:
                WV(:,:,k,n+1) = WV(:,:,k,n+1) + brac(k-1)/(r*2)^(3*k/2)*( ...
                    -(4*q^2+2*k-1)*g(k,n+1)/2/k*eye(2) );
            end
        end
    end
end
