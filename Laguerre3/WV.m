% Compute the W- or V-matrices to construct the asymptotic expansion of R,
% using the procedure with the convolutions as explained in the paper.
% Input
%   s            - Structure containing alpha, m etc
%   maxOrder     - The maximum order of the error
%   r            - 1 when computing Wright, -1 when computing Wleft
%   isW          - 1 if the result is the W-s, 0 for the V-s
% Output
%   WV           - Coefficient matrices for (z - 1/2 \pm 1/2)^m of s_k(z) or Delta_k(z)
% About
%   Author       - Peter Opsomer (Peter.Opsomer@cs.kuleuven.be)
%   History      - Created July 2015 from Jacobi version, corrected October 2016
function WV = WV(s,maxOrder,r,isW)

mo = ceil(3*maxOrder/2)+4;
ns = 0:mo;
monomial = isnumeric(s.q) && (norm(s.q(2:end-1)) == 0);
if monomial || ~isnumeric(s.q)
    f = nan*zeros(mo+1,1); % Coefficients in the expansion of \bar{phi}_n(z) or \xi_n(z)
    g = nan*zeros(maxOrder-1,mo+1); % We multiply by nan to see which entries are used.
elseif 1
    f = nan*zeros(mo+1, 2*s.m*maxOrder-1); % first idx is power of z, second n^(1/s.m)
    g = nan*zeros(maxOrder-1,mo+1, 2*s.m*maxOrder-1);
end

if isnumeric(s.q)
    A = zeros(s.m+1,1);
    for k =0:s.m
        A(k+1) = prod((2*(1:k)-1)/2./(1:k));
    end
end

if (r == 1) % Right disk: near z=1
    ns = [ns, ns(end)+1]; % Extend by one because f(1) = 0 while not for left case
    u = zeros(ns(end)+1,ns(end)+1);
    v = zeros(ns(end)+1,ns(end)+1); 
    u(1,1) = 1;
    v(1,1) = 1;
    for n = [ns, ns(end)+1]
        u(2,n+1) = s.binom(1/2,n+1);
        v(2,n+1) = s.binom(1/2,n+2);
    end
    for kt = 2:ns(end)
        for n = ns
            u(kt+1,n+1) = sum(u(kt,(0:n)+1).*u(2,n-(0:n)+1));
            v(kt+1,n+1) = sum(v(kt,(0:n)+1).*v(2,n-(0:n)+1));
        end
    end
    q = zeros(ns(end)+1,1);
    rr = zeros(ns(end)+1,1); % Coefficients in the expansion of sqrt(2-2*sqrt(1-w))
    for kt = ns
        for l = 0:kt
            q(kt+1) = q(kt+1) + s.poch(1/2,kt-l)*u(kt-l+1,l+1)/(-2)^(kt-l)/factorial(kt-l)/(1+2*(kt-l) );
            rr(kt+1) = rr(kt+1) + s.binom(1/2,kt-l)*v(kt-l+1,l+1)*2^(kt-l);
        end
    end
    if isnumeric(s.q) && (length(s.q) == 2)
        for n = ns
            f(n+1,1) = -2*s.binom(1/2,n);
            for l = 0:n
                f(n+1,1) = f(n+1,1) + 2*q(l+1)*rr(n-l+1);
            end
        end
    elseif isnumeric(s.q)
        for j = ns
            f(j+1,1) = 0;
            for i=0:min(j,s.m-1)
                f(j+1,1) = f(j+1,1) + s.binom(1/2,j-i)*(-1)^(s.m-i-1)*gamma(-1/2-i)/gamma(1/2-s.m)/gamma(s.m-i); 
                % Simplification to binomial is not done here because s.binom is meant for integer second arguments:
%                f(j+1,1) = f(j+1,1) + s.binom(1/2,j-i)*(-1)^(s.m-i-1)*s.binom(-3/2-i,-1/2-s.m); 
            end
            f(j+1,1) = -f(j+1)/s.m/A(s.m+1);
            for l = 0:j
                f(j+1,1) = f(j+1,1) + 2*q(l+1)*rr(j-l+1);
            end
        end
        if ~monomial % Contributions in higher orders of n
            for l = 1:size(f,2)-1
                for j = ns
                    f(j+1,l+1) = 0;
                    for k = 0:s.m-1
                        tmp1 = 0;
                        for z = max(k+1,s.m-l):s.m
                            tmp1 = tmp1 + s.q(z+1)*A(z-k)*s.betac(z,z+l-s.m+1);
                        end
                        tmp2 = 0;
                        for i=0:min(k,j)
                            tmp2 = tmp2 + s.binom(1/2,j-i)*s.binom(k,i);
                        end
                        f(j+1,l+1) = f(j+1,l+1) + tmp1*tmp2;
                    end
                end
            end
            f(:,2:end) = -f(:,2:end)/2; 
        end
    else
        c = arrayfun(@(nn) s.trap_rule(@(y) sqrt(y).*s.dVn(y)./( sqrt(y -1).*(y-1).^nn), 1), [ns+1, ns(end)+2]);
        f = transpose([0 arrayfun(@(j) -1/(2*j+1)*sum(arrayfun(@(m) s.binom(-1/2,m)*c(j-1-m+1), 0:(j-1))), 1+ns)]);
    end
    if(abs(f(1)) > 10*eps)
        l = 2; fl = 0; for z=max(1,s.m-l):s.m, fl = fl + s.q(z+1)*s.betac(z,z+l-s.m+1)*2*z*A(z+1); end; fl
        error('xi_n should be O( (z-1)^(3/2) ): Expected f(1) to be zero');
    end
    ns = ns(1:end-1); % Reset ns to its value before computing f's
    
    g(1,1,1) = -1/f(2,1);
    for n = 1:mo
        g(1,n+1,1) = -g(1,1:n,1)*f((n+2):-1:3,1)/f(2,1);
    end
    if isnumeric(s.q) && ~monomial
        for l = 1:size(f,2)-1
            for i = ns
                g(1,i+1,l+1) = 0;
                for y=0:i
                    for p=0:(i-y)
                        tmp = 0;
                        for q=0:l-1
                            tmp = tmp + g(1,y+1,q+1)*f(p+2,l-q+1);
                        end
                        g(1,i+1,l+1) = g(1,i+1,l+1) + g(1,i-y-p+1,1)*tmp;
                        if isnan(g(1,i+1,l+1))
                            error('Wrong indices were implemented')
                        end
                    end
                end
            end
        end
    end
elseif 1 % Left disk: near z=0
    if isnumeric(s.q) && (length(s.q) == 2)
        for n = ns
            f(n+1,1) = -(s.binom(1/2,n)*(-1)^n + s.poch(1/2,n)./(1+2*n)./factorial(n));
        end
    elseif isnumeric(s.q)
        for n = ns
            f(n+1,1) = 0;
            for k = 0:min(s.m-1,n)
                f(n+1,1) = f(n+1,1) + s.binom(1/2,n-k)*(-1)^(n-k)*A(s.m-k);
            end
            f(n+1,1) = -f(n+1,1)/2/s.m/A(s.m+1)-s.poch(1/2,n)./(1+2*n)./factorial(n);
        end
        if ~monomial % Contributions in higher orders of n
            for l = 1:size(f,2)-1
                for j = ns
                    f(j+1,l+1) = 0;
                    for i = 0:min(j,s.m-1)
                        tmp = 0;
                        for n = max(i+1,s.m-l):s.m
                            tmp = tmp + s.q(n+1)*A(n-i)*s.betac(n,n+1+l-s.m);
                        end
                        f(j+1,l+1) = f(j+1,l+1) + (-1)^(j-i)*s.binom(1/2,j-i)*tmp;
                    end
                end
            end
            f(:,2:end) = -f(:,2:end)/4;
        end
    else
        d = arrayfun(@(nn) s.trap_rule(@(y) sqrt(y).*s.dVn(y)./( sqrt(y -1).*y.^nn), 0), [ns+1, ns(end)+2]);
        f = transpose(-1/2*arrayfun(@(j) 1/(2*j+1)*sum(arrayfun(@(m) s.binom(1/2,m)*(-1)^m*d(j-m+1), 0:j)), ns));
    end
    
    g(1,1,1) = 1/f(1,1);
    for n = 1:mo
        g(1,n+1,1) = -g(1,1:n,1)*f((n+1):-1:2,1)/f(1,1);
    end
    if isnumeric(s.q) && ~monomial
        for l = 1:size(f,2)-1
            for i = ns
                g(1,i+1,l+1) = 0;
                for y=0:i
                    for p=0:(i-y)
                        tmp = 0;
                        for q=0:l-1
                            tmp = tmp + g(1,y+1,q+1)*f(p+1,l-q+1);
                        end
                        g(1,i+1,l+1) = g(1,i+1,l+1) + g(1,i-y-p+1,1)*tmp;
                        if isnan(g(1,i+1,l+1))
                            error('Wrong indices were implemented')
                        end
                    end
                end
            end
        end
    end
end

rho = zeros(maxOrder,mo+1); 
for n = ns
    rho(2,n+1) = s.poch(1/2,n)/factorial(n)/(1+2*n)*(-r)^n;
end
rho(1,1) = 1;

for i = 2:(maxOrder-1)
    for n = ns
        g(i,n+1) = sum(g(i-1,1:(n+1) ).*g(1,(n+1):-1:1) );
    end
end
for i = 2:(mo*2+2)
    for n = ns     
        rho(i+1,n+1) = sum(rho(i,1:(n+1) ).*rho(2,(n+1):-1:1) );       
    end
end
if isnumeric(s.q) && ~monomial
    for k = 2:max(mo*2+2,(maxOrder-1))
        for l = 0:size(f,2)-1
            for n = ns
                g(k,n+1,l+1) = 0;
                for j=0:l
                    for i=0:n
                        g(k,n+1,l+1) = g(k,n+1,l+1) + g(k-1,i+1,j+1)*g(1,n-i+1,l-j+1);
                    end
                end
            end
        end
    end
end

OmOdd = zeros(mo+1,1); OmEven = zeros(mo+1,1);
XiOdd = zeros(mo+1,1); XiEven = zeros(mo+1,1); 
ThOdd = zeros(mo+1,1); ThEven = zeros(mo+1,1);

OmO = zeros(mo+1,1); OmE = zeros(mo+1,1);
XiO = zeros(mo+1,1); XiE = zeros(mo+1,1);
ThO = zeros(mo+1,1); ThE = zeros(mo+1,1);
for n = ns
    js = 0:n;
    for j = js
        OmOdd(n+1) = OmOdd(n+1) + (-1)^j/factorial(2*j)*(-2*s.alpha/sqrt(-r))^(2*j)*rho(2*j+1,n-j+1);
        XiOdd(n+1) = XiOdd(n+1) + (-1)^j/factorial(2*j)*(-2*(s.alpha+1)/sqrt(-r))^(2*j)*rho(2*j+1,n-j+1);
        ThOdd(n+1) = ThOdd(n+1) + (-1)^j/factorial(2*j)*(-2*(s.alpha-1)/sqrt(-r))^(2*j)*rho(2*j+1,n-j+1);
        
        OmEven(n+1) = OmEven(n+1) + (-1)^j/factorial(2*j+1)*(-2*s.alpha/sqrt(-r))^(2*j+1)*rho(2*j+2,n-j+1);
        XiEven(n+1) = XiEven(n+1) + (-1)^j/factorial(2*j+1)*(-2*(s.alpha+1)/sqrt(-r))^(2*j+1)*rho(2*j+2,n-j+1);
        ThEven(n+1) = ThEven(n+1) + (-1)^j/factorial(2*j+1)*(-2*(s.alpha-1)/sqrt(-r))^(2*j+1)*rho(2*j+2,n-j+1);
    end
    for j = js
        OmO(n+1) = OmO(n+1) + s.binom(-1/2,j)*(r)^j*OmOdd(n-j+1);
        XiO(n+1) = XiO(n+1) + s.binom(-1/2,j)*(r)^j*XiOdd(n-j+1);
        ThO(n+1) = ThO(n+1) + s.binom(-1/2,j)*(r)^j*ThOdd(n-j+1);
        
        OmE(n+1) = OmE(n+1) + s.binom(-1/2,j)*(r)^j*OmEven(n-j+1);
        XiE(n+1) = XiE(n+1) + s.binom(-1/2,j)*(r)^j*XiEven(n-j+1);
        ThE(n+1) = ThE(n+1) + s.binom(-1/2,j)*(r)^j*ThEven(n-j+1);
    end
end

Ts = zeros(2,2,mo+1);  % Proportional to G_{k,n}^{odd/even} or \Omega_{k,n}^{odd/even}, depends on k, overwritten on each new k.
WV = zeros(2,2,maxOrder-1,mo+1);
for k = 1:(maxOrder-1)
    Ts(:,:,:) = 0;
    if r == 1
        if mod(k,2)
            for n = 0:mo
                Ts(:,:,n+1) = s.nuk(k)*[-2*(2*s.binom(-1/2,n-1)*(n>0)+s.binom(-1/2,n)), 2i*4^(-s.alpha)*s.binom(-1/2,n) ; ...
                    2j*4^(s.alpha)*s.binom(-1/2,n), 2*(2*s.binom(-1/2,n-1)*(n>0) +s.binom(-1/2,n))] + ...
                    s.muk(k)*[-2*OmO(n+1), 4^(-s.alpha)*2j*XiO(n+1); 4^(s.alpha)*2j*ThO(n+1), 2*OmO(n+1)];
                WV(:,:,k,n+1) = sum(repmat(reshape(g(k,1:(n+1) ),[1,1,n+1]),[2,2,1]).*Ts(:,:,(n+1):-1:1),3)/8;
            end
        else
            for n = 0:mo
                 Ts(:,:,n+1) = s.nuk(k)*4*(n==0)*eye(2) -s.muk(k)*[-2i*OmE(n+1), -2*4^(-s.alpha)*XiE(n+1); -2*4^s.alpha*ThE(n+1), 2i*OmE(n+1)];
                 WV(:,:,k,n+1) = sum(repmat(reshape(g(k,1:(n+1) ),[1,1,n+1]),[2,2,1]).*Ts(:,:,(n+1):-1:1),3)/8;
            end
        end
    else
      if mod(k,2)
        for n = 0:mo	
            Ts(:,:,n+1) = -(s.alpha^2+k/2-1/4)/k*[-(-1)^n*(2*s.binom(-1/2,n-1)*(n>0)+s.binom(-1/2,n))*2, -1i*4^(-s.alpha)*2*(-1)^n*...
                s.binom(-1/2,n) ; -1j*4^(s.alpha)*2*(-1)^n*s.binom(-1/2,n), (-1)^n*(2*s.binom(-1/2,n-1)*(n>0) +s.binom(-1/2,n))*2] - ...
               (k-1/2)*[2*OmO(n+1), 4^(-s.alpha)*2j*XiO(n+1); 4^(s.alpha)*2j*ThO(n+1), -2*OmO(n+1)]; % binom(-1/2,-1) should be zero  
            WV(:,:,k,n+1) = -(-1)^(ceil(k/2)+1)*(1j*sqrt(2))^k*(-2)^(-k/2)/4^(k+1)*s.brac(k-1)*sum(repmat(reshape(...
                g(k,1:(n+1) ),[1,1,n+1]),[2,2,1]).*Ts(:,:,(n+1):-1:1),3);
        end
      else
        for n = 0:mo
            Ts(:,:,n+1) = (s.alpha^2+k/2-1/4)/k*4*(n==0)*eye(2) ...
                -2*(k-1/2)*[OmE(n+1), 4^(-s.alpha)*1j*XiE(n+1); 4^s.alpha*1j*ThE(n+1), -OmE(n+1)];
            WV(:,:,k,n+1) = -(-1)^(ceil(k/2)+1)*(1j*sqrt(2))^k*(-2)^(-k/2)/4^(k+1)*s.brac(k-1)*sum(repmat(reshape(...
                g(k,1:(n+1) ),[1,1,n+1]),[2,2,1]).*Ts(:,:,(n+1):-1:1),3);
        end
      end
    end 
    if isW && (mod(k,2) == 0)
        if r == 1
            WV(:,:,k,:) = WV(:,:,k,:) - s.vuk(k)*repmat(reshape(g(k,:,1),[1,1,1,mo+1]),[2,2,1,1]).*repmat(eye(2),[1,1,1,mo+1]);
        else
            WV(:,:,k,:) = WV(:,:,k,:) - s.brac(k-1)/2*(-1)^(k/2)/k*(4*s.alpha^2+2*k-1)/4^k*repmat(reshape(g(k,:,1),[1,1,1,mo+1]),[2,2,1,1]).*repmat(eye(2),[1,1,1,mo+1]);
        end
    end
end
% The following is for general polynomial Q(x)
WV = repmat(WV,[1,1,1,1,size(f,2)]);
for l=2:size(f,2)
    for k = 1:(maxOrder-1)
        Ts(:,:,:) = 0;
        if r == 1
            if mod(k,2)
                for n = 0:mo
                    Ts(:,:,n+1) = s.nuk(k)*[-2*(2*s.binom(-1/2,n-1)*(n>0)+s.binom(-1/2,n)), 2i*4^(-s.alpha)*s.binom(-1/2,n) ; ...
                        2j*4^(s.alpha)*s.binom(-1/2,n), 2*(2*s.binom(-1/2,n-1)*(n>0) +s.binom(-1/2,n))] + ...
                        s.muk(k)*[-2*OmO(n+1), 4^(-s.alpha)*2j*XiO(n+1); 4^(s.alpha)*2j*ThO(n+1), 2*OmO(n+1)];
                    WV(:,:,k,n+1,l) = sum(repmat(reshape(g(k,1:(n+1),l),[1,1,n+1]),[2,2,1]).*Ts(:,:,(n+1):-1:1),3)/8;
                end
            else
                for n = 0:mo
                    Ts(:,:,n+1) = s.nuk(k)*4*(n==0)*eye(2) ...
                        -s.muk(k)*[-2i*OmE(n+1), -2*4^(-s.alpha)*XiE(n+1); -2*4^s.alpha*ThE(n+1), 2i*OmE(n+1)];
                    WV(:,:,k,n+1,l) = sum(repmat(reshape(g(k,1:(n+1),l),[1,1,n+1]),[2,2,1]).*Ts(:,:,(n+1):-1:1),3)/8;
                end
            end
        else
            if mod(k,2)
                for n = 0:mo
                    Ts(:,:,n+1) = -(s.alpha^2+k/2-1/4)/k*[-(-1)^n*(2*s.binom(-1/2,n-1)*(n>0)+s.binom(-1/2,n))*2, -1i*4^(-s.alpha)*2*(-1)^n*...
                        s.binom(-1/2,n) ; -1j*4^(s.alpha)*2*(-1)^n*s.binom(-1/2,n), (-1)^n*(2*s.binom(-1/2,n-1)*(n>0) +s.binom(-1/2,n))*2] - ...
                        (k-1/2)*[2*OmO(n+1), 4^(-s.alpha)*2j*XiO(n+1); 4^(s.alpha)*2j*ThO(n+1), -2*OmO(n+1)]; % binom(-1/2,-1) should be zero
                    WV(:,:,k,n+1,l) = -(-1)^(ceil(k/2)+1)*(1j*sqrt(2))^k*(-2)^(-k/2)/4^(k+1)*s.brac(k-1)*sum(repmat(reshape(...
                        g(k,1:(n+1),l),[1,1,n+1]),[2,2,1]).*Ts(:,:,(n+1):-1:1),3);
                end
            else
                for n = 0:mo
                    Ts(:,:,n+1) = (s.alpha^2+k/2-1/4)/k*4*(n==0)*eye(2) ...
                        -2*(k-1/2)*[OmE(n+1), 4^(-s.alpha)*1j*XiE(n+1); 4^s.alpha*1j*ThE(n+1), -OmE(n+1)];
                    WV(:,:,k,n+1,l) = -(-1)^(ceil(k/2)+1)*(1j*sqrt(2))^k*(-2)^(-k/2)/4^(k+1)*s.brac(k-1)*sum(repmat(reshape(...
                        g(k,1:(n+1),l),[1,1,n+1]),[2,2,1]).*Ts(:,:,(n+1):-1:1),3);
                end
            end
        end
        if isW && (mod(k,2) == 0)
            if r == 1
                WV(:,:,k,:,l) = WV(:,:,k,:,l) + s.muk(k)/6/k*repmat(reshape(g(k,:,l),[1,1,1,mo+1]),[2,2,1,1]).*repmat(eye(2),[1,1,1,mo+1]);
            else
                WV(:,:,k,:,l) = WV(:,:,k,:,l) - s.brac(k-1)/2/k*(s.alpha^2+2*k-1)/4^k*repmat(reshape(g(k,:,l),[1,1,1,mo+1]),[2,2,1,1]).*repmat(eye(2),[1,1,1,mo+1]);
            end
        end
    end
end
