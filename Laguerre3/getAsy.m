% Get the asymptotics of the Laguerre-type orthonormal polynomials. 
% Input
%   alpha        - Part of the weight w(x) =  x.^alpha.*exp(-Q(x) ), default zero
%   q            - Q(x) = sum(x.^(0:length(q)-1)*q ) if q (default [0 1]) is numeric, else Q(x) = q(x)
%   maxOrder     - The maximal order of the error, default 5
%   [dQ          - Derivative of Q(x) to calculate \beta_n using integral of dQ if Q(x) is not a polynomial]
%   [Qinv        - Inverse of Q(x) to provide initial guess \beta_n = Qinv(n) if Q(x) is not a polynomial]
%   [np          - The degree for the polynomial if Q(x) is not a polynomial]
% Output 
%   s            - Structure containing the asymptotic information: s.alpha and s.q are the inputs above, s.betan(n,T) is the expansion of the MRS number.
%                    s.pa(n,x,z) is the expansion of the orthonormal polynomial in the outer region, s.pb(n,x,z) the expansion in the inner
%                    region (0,\beta_n), s.pc near the right disk and s.pd near the left disk: when appended with 'woq', e^(n V_n(z)/2 is not included. 
%                    s.bnm1 and s.an are the expansions of the recurrence coefficients, s.gamman is the expansion of the normalizing constants.
%                    When isnumeric(q), s.m is length(q)-1; else, s.np is the input above. The other fields of s contain temporary variables.
% About
%   Author       - Peter Opsomer (Peter.Opsomer@cs.kuleuven.be)
%   History      - Created July 2015, betas April 2016, refactored May 2016 and December 2016
function s = getAsy(alpha,q,maxOrder, dQ, Qinv,np)
if ~exist('alpha','var')
    alpha = 0; % Assume non-associated Laguerre
end
if ~exist('q','var')
    q = [0 1]; % Assume std associated Laguerre
end
if ~exist('maxOrder','var')
    maxOrder = 5;
end
if ~exist('dQ','var')
    hh = sqrt(eps); % Use numeric differentiation in the complex plane for real input.
    dQ = @(x) imag(q(x+1i*hh*abs(x)))/hh./abs(x);
end
if ~exist('Qinv','var')
    Qinv = @(y) fzero(@(x) q(x)-y, 1); % Error when isnumeric(q) but will not be called then
end

s = struct;
s.alpha = alpha;
s.q = q;
if ~isnumeric(q)
    s.np = np;
    bet = fzero(@(bet) quadgk(@(x) dQ(x).*sqrt(x./(bet-x)), 0, bet)/2/pi/np-1, Qinv(np) );
    s.betan = @(n,T) bet/(n == np); % Safety for when not evaluating at correct n, but maybe f_n is not that dependent on n
    Qx = @(bz) q(bz);
    
    s.dVn = @(z) bet*dQ(bet*z)/np;
    rho = @(z) 1.5*max([0.5,abs(z/2+1/4 -[0,1])]); % Make sure z lies within the contour
    M = 200; % should probably increase with (n and) the desired accuracy
    ths = linspace(0,2*pi,M+1); % = theta_k, k = 0 to M-1
    ths(end) = [];
    hh = ths(2)-ths(1);
    s.trap_rule = @(f,z) hh*sum(f(rho(z).*exp(1i*ths)+z/2+1/4).*rho(z).*1i.*exp(1i*ths) )/(2*pi*1i );
    hn = @(z) arrayfun(@(zz) s.trap_rule(@(y) sqrt(y).*s.dVn(y)./( sqrt(y -1).*(y - zz) ), zz), z);
    
    zt = 0.5; % Test-z to calculate l_n (3.24)
    lln = @(n,T) 2*quadgk(@(y) log(abs(zt-y)).*sqrt(1-y)./sqrt(y)/2/pi.*hn(y), 0, 1) -Qx(bet*zt)/np/(n==np);
    xin = @(n,z,T) -1/2*quadgk(@(x) sqrt(x-1)./sqrt(x).*hn(x), 1, z)/(n==np);
    
elseif(norm(q(2:end-1)) == 0) % 'Monomial' Q(x) including linear Q(x)
    s.m = length(q)-1;
    coeff = (1/2*s.m*q(end)*prod((2*(1:s.m)-1)/2./(1:s.m)))^(-1/s.m);
	s.betan = @(n,T) n.^(1/s.m)*coeff;
	lln = @(n,T) -2/s.m -4*log(2) -q(1)/n;
    if (length(q) == 2)
        Hn = @(n,z,T) 4*ones(size(z));
    else
        A = ones(s.m+1,1);
        for k =1:s.m
            A(k+1) = A(k)*(2*k-1)/2/k;
        end
        %  Also correct: A = zeros(s.m+1,1); for k =0:s.m, A(k+1) = prod((2*(1:k)-1)/2./(1:k)); end
        Hn = @(n,z,T) arrayfun(@(zz) 2/A(s.m+1)*zz.^(0:s.m-1)*A(s.m-(0:s.m-1))/s.m, z);
        % Also correct but slower: Hn = @(n,z,T) 4*s.m/(2*s.m-1)*double(hypergeom([1, 1-s.m], 3/2-s.m, z))/s.m;
    end
else % Non-monomial Q(x)
    s.m = length(q)-1;
    A = ones(s.m+1,1);
    for k =1:s.m
        A(k+1) = A(k)*(2*k-1)/2/k;
    end
    s.A = A;
    s.betac = nan*zeros(s.m, 2*s.m*maxOrder-1); % Only need first row for s.betan, but also need powers in WV so save all in s
    j = 0;
    s.betac(1,j+1) = (s.m/2*q(s.m+1)*A(s.m+1))^(-1/s.m);
    for kt = 2:s.m
        s.betac(kt,j+1) = sum(s.betac(kt-1,(0:j)+1).*s.betac(1,(j-(0:j))+1));
    end
    for j = 1:(size(s.betac,2)-1)
        s.betac(1,j+1) = 0;
        for i = 1:(j-1) % We could use s.betac(s.m-2,1) instead of s.betac(1,1)^(s.m-2) but that gives errors if s.m == 2.
            s.betac(1,j+1) = s.betac(1,j+1) + s.betac(1,1)^(s.m-2)*s.betac(1,j-i+1)*s.betac(1,i+1);
        end
        for i = 0:(s.m-3)
            for k = 1:(j-1)
                s.betac(1,j+1) = s.betac(1,j+1) + s.betac(1,1)^i*s.betac(1,j-k+1)*s.betac(s.m-1-i,k+1);
            end
        end
        s.betac(1,j+1) = s.betac(1,j+1)*s.m*q(s.m+1)*A(s.m+1);
        for kk = max(s.m-j, 1):(s.m -1)
            s.betac(1,j+1) = s.betac(1,j+1) +kk*q(kk+1)*A(kk+1)*s.betac(kk,j-s.m+kk+1);
        end
        s.betac(1,j+1) = -s.betac(1,j+1)/(s.m^2*q(s.m+1)*A(s.m+1)*s.betac(s.m-1,1));
        
        for kt = 2:s.m
            s.betac(kt,j+1) = sum(s.betac(kt-1,(0:j)+1).*s.betac(1,(j-(0:j))+1));
        end
    end
    s.betan = @(n,T) n.^(1/s.m)*n.^(-(0:T-1)/s.m)*transpose(s.betac(1,1:T)); % Make it possible to get O(n^(-1-1/m))
    lln = @(n,T) -4*log(2)-sum(arrayfun( @(k) q(k+1)/n*s.betan(n,T)^k*A(k+1), 0:s.m)); % Or more explicitly in n:
% lln = @(n,T) -4*log(2)-q(1+0)*n^(0-1)*1*A(1+0)-sum(arrayfun(@(k) q(k+1)*n^(k/s.m-1)*(s.betac(k,1:s.m*T)*n.^transpose(-(0:s.m*T-1)/s.m))*A(k+1), 1:s.m)); 

    Hn = @(n,z,T) arrayfun(@(zz) sum(arrayfun(@(k) sum(arrayfun(@(j) q(j+1)/n*s.betan(n,T)^j*A(j-k), k+1:s.m)).*zz.^k, 0:s.m-1)), z); % Or explicitly:
%     Hn = @(n,z,T) arrayfun(@(zz) sum(arrayfun(@(k) sum(arrayfun(@(j) q(j+1)*n^(j/s.m-1)*A(j-k)*...
%         sum(arrayfun(@(i) s.betac(j,i)*n^(-(i-1)/s.m), 1:T)), k+1:s.m)).*zz.^k, 0:s.m-1)), z);
end
s.lln = @(n,T) lln(n,T);
if isnumeric(q)
    xin = @(n,z,T) -1i.*(sqrt(z).*sqrt(1-z).*Hn(n,z,T)/2 -2*acos(sqrt(z)));
%     s.xin = @(n,z,T) xin(n,z,T);
    Qx = @(bz) q*repmat(bz,length(q),1).^repmat(transpose(0:length(s.q)-1),1,length(bz));
end


phi = @(z) 2*z-1+2*sqrt(z)*sqrt(z-1);
fn = @(n,z,T) n^(2/3)*(z-1)*(1/(z-1)^(3/2)*(-3/2)*(-1).^(angle(z-1)<=0).*xin(n,z,T))^(2/3); % Do not simplify to keep analytic continuation
sqphitn = @(n,z,T) -pi*1i/2 + xin(n,z,T)/2;
% s.sqphitn = @(n,z,T) sqphitn(n,z,T);

s3 = @(z) [z, 0; 0, 1/z];
Pinf = @(z) s3(2^(-alpha))*[sqrt(phi(z)), 1i/sqrt(phi(z)); -1i/sqrt(phi(z)), sqrt(phi(z))]/2/z^(1/4)/(z-1)^(1/4)*...
	s3(phi(z)^(alpha/2)/z^(alpha/2) );
PinfInv = @(z) s3(phi(z)^(-alpha/2)/z^(-alpha/2) )*[sqrt(phi(z)), -1i/sqrt(phi(z)); 1i/sqrt(phi(z)), sqrt(phi(z))]...
    /2/z^(1/4)/(z-1)^(1/4)*s3(2^(alpha));

s.poch = @(x,j) prod(x+(0:(j-1)) ); % not to be used vectorized, x may be non-integer
% 'pochhammer' needs symbolic toolbox, could also use separate 1-line function files or global functions for these functions
s.binom = @(x,n) prod(x:-1:(x-n+1) )/prod(1:n)*0^((n<0) + abs(rem(n,1)) ); % Not to be used vectorized (and n should be a positive integer)
s.nuk = @(n) -gamma(3*n-1/2)*2^n/27^n/2/n/sqrt(pi)/gamma(n*2);
s.muk = @(n) 3*gamma(3*n-1/2)*2^n/27^n/sqrt(pi)/gamma(n*2); % s.muk(k) = -6*k*s.nuk(k)
s.brac = @(n) prod(4*alpha^2-(2*(1:n)-1).^2 )/(2^(2*n)*factorial(n));

[s.UR,s.UL,Qright,Qleft] = UQ(s,maxOrder,'UQV');
% [s.UR,s.UL,Qright,Qleft] = UQ(s,maxOrder,'UQW');
% [s.UR,s.UL] = UQ(s,maxOrder,'UW');
s.QR = Qright; s.QL = Qleft;
zpow = repmat(reshape(1:1:size(s.UL,4),[1,1,1,size(s.UL,4)]),[2,2,maxOrder-1,1]);

Rk = @(z,T) sum(s.UL(:,:,1:T-1,:)./z.^zpow(:,:,1:T-1,:)+s.UR(:,:,1:T-1,:)./(z-1).^zpow(:,:,1:T-1,:),4);
powm = isnumeric(s.q) && (norm(s.q(2:end-1)) ~= 0);
if powm
    R = @(n,z,T) eye(2) + sum(Rk(z,T)./n.^repmat( (reshape(1:T-1,[1,1,T-1])-1)./s.m+1,[2,2,1]),3);
    s.bnm1 = @(n,T) s.betan(n,T)/4*sqrt(1+4i*sum(((s.UL(2,1,1:(T-1),1)+s.UR(2,1,1:(T-1),1))*4^(-s.alpha) ...
        -4^(s.alpha)*(s.UL(1,2,1:(T-1),1)+s.UR(1,2,1:(T-1),1)) )./n.^reshape( (0:(T-2))/s.m+1,[1,1,T-1])) ...
        +16*sum((s.UL(2,1,1:(T-1),1)+s.UR(2,1,1:(T-1),1))./n.^reshape((0:(T-2))/s.m+1,[1,1,T-1]) )*...
        sum((s.UL(1,2,1:(T-1),1)+s.UR(1,2,1:(T-1),1))./n.^reshape((0:(T-2))/s.m+1,[1,1,T-1]) ));
    s.an = @(n,T) s.betan(n,T)*(-alpha/4+sum((s.UL(1,1,1:(T-1),1)+s.UR(1,1,1:(T-1),1))./n.^reshape((0:(T-2))/s.m+1,[1,1,T-1]) ) +...
        (4^(-s.alpha)*1i*(s.alpha+2)/16 +sum(s.UL(1,2,1:(T-1),2)./n.^reshape((0:(T-2))/s.m+1,[1,1,T-1]) ) +... %3:T-1 changed to 1:T-1 for ease
        sum((s.UR(1,2,1:(T-1),2)+s.UR(1,2,1:(T-1),1))./n.^reshape((0:(T-2))/s.m+1,[1,1,T-1]) ) + ...
        sum(((s.UL(1,1,1:(T-1),1)+s.UR(1,1,1:(T-1),1))*4^(-1-s.alpha)*1i +(s.UL(1,2,1:(T-1),1)+s.UR(1,2,1:(T-1),1))*s.alpha/4)./...
        n.^reshape((0:(T-2))/s.m+1,[1,1,T-1]) ) )...
        /(1i*4^(-s.alpha-1)+ sum((s.UL(1,2,1:(T-1),1)+s.UR(1,2,1:(T-1),1))./n.^reshape((0:(T-2))/s.m+1,[1,1,T-1]) ) ));
    s.gamman = @(n,T) s.betan(n,T)^(-n-s.alpha/2-1/2)*exp(-n*lln(n,T)/2)*sqrt(2/pi)*2^alpha...
        /sqrt(1-4i*4^alpha*sum((s.UL(1,2,1:(T-1),1)+s.UR(1,2,1:(T-1),1))./n.^reshape((0:(T-2))/s.m+1,[1,1,T-1]) ) );
    % Remvove the factors that appear inversely in the expansions of the polynomials to avoid rapid over/underflow:
    gammaDen = @(n,T) s.betan(n,T)^(-1/2)*sqrt(2/pi)*2^alpha...
        /sqrt(1-4i*4^alpha*sum((s.UL(1,2,1:(T-1),1)+s.UR(1,2,1:(T-1),1))./n.^reshape((0:(T-2))/s.m+1,[1,1,T-1]) ) );
    s.ratGam = @(n,T) (1-1/n)^(n+s.alpha/2)*gammaDen(n,T)/(4*n-4)/exp(lln(n,T)/2)/gammaDen(n-1,T); % gamma_n/gamma_{n-1} for Qx(x) = x
else
    R = @(n,z,T) eye(2) + sum(Rk(z,T)./n.^repmat(reshape(1:T-1,[1,1,T-1]),[2,2,1]),3);
    s.bnm1 = @(n,T) s.betan(n,T)/4*sqrt(1+4i*sum(((s.UL(2,1,1:(T-1),1)+s.UR(2,1,1:(T-1),1))*4^(-s.alpha) ...
        -4^(s.alpha)*(s.UL(1,2,1:(T-1),1)+s.UR(1,2,1:(T-1),1)) )./n.^reshape(1:(T-1),[1,1,T-1])) ...
        +16*sum((s.UL(2,1,1:(T-1),1)+s.UR(2,1,1:(T-1),1))./n.^reshape(1:(T-1),[1,1,T-1]) )*...
        sum((s.UL(1,2,1:(T-1),1)+s.UR(1,2,1:(T-1),1))./n.^reshape(1:(T-1),[1,1,T-1]) ));
    s.an = @(n,T) s.betan(n,T)*(-alpha/4+sum((s.UL(1,1,1:(T-1),1)+s.UR(1,1,1:(T-1),1))./n.^reshape(1:(T-1),[1,1,T-1]) ) +...
        (4^(-s.alpha)*1i*(s.alpha+2)/16 +sum(s.UL(1,2,3:(T-1),2)./n.^reshape(3:(T-1),[1,1,max(0,T-3)]) ) +...
        sum((s.UR(1,2,1:(T-1),2)+s.UR(1,2,1:(T-1),1))./n.^reshape(1:(T-1),[1,1,T-1]) ) + ...
        sum(((s.UL(1,1,1:(T-1),1)+s.UR(1,1,1:(T-1),1))*4^(-1-s.alpha)*1i +(s.UL(1,2,1:(T-1),1)+s.UR(1,2,1:(T-1),1))*s.alpha/4)./n.^reshape(1:(T-1),[1,1,T-1]) ) )...
        /(1i*4^(-s.alpha-1)+ sum((s.UL(1,2,1:(T-1),1)+s.UR(1,2,1:(T-1),1))./n.^reshape(1:(T-1),[1,1,T-1]) ) ));
    s.gamman = @(n,T) s.betan(n,T)^(-n-s.alpha/2-1/2)*exp(-n*lln(n,T)/2)*sqrt(2/pi)*2^alpha...
        /sqrt(1-4i*4^alpha*sum((s.UL(1,2,1:(T-1),1)+s.UR(1,2,1:(T-1),1))./n.^reshape(1:(T-1),[1,1,T-1]) ) )*exp(q(1)/2);
    gammaDen = @(n,T) s.betan(n,T)^(-1/2)*sqrt(2/pi)*2^alpha...
        /sqrt(1-4i*4^alpha*sum((s.UL(1,2,1:(T-1),1)+s.UR(1,2,1:(T-1),1))./n.^reshape(1:(T-1),[1,1,T-1]) ) );
    s.ratGam = @(n,T) (1-1/n)^(n+s.alpha/2)*gammaDen(n,T)/(4*n-4)/exp(lln(n,T)/2)/gammaDen(n-1,T); % gamma_n/gamma_{n-1} for Q = x
end

s.pa = @(n,bz,T) gammaDen(n,T)*[1, 0]*R(n,bz/s.betan(n,T),T)*[2^(-alpha)*exp(1i*(-1).^(angle(bz/s.betan(n,T)-1)<=0)*acos(2*bz/s.betan(n,T)-1)/2); ...
    -1i*2^(alpha)*exp(-1i*(-1).^(angle(bz/s.betan(n,T)-1)<=0)*acos(2*bz/s.betan(n,T)-1)/2) ]/2/(bz/s.betan(n,T))^(1/4)/(bz/s.betan(n,T)-1)^(1/4)*...
	exp(1i*(-1).^(angle(bz/s.betan(n,T)-1)<=0)*acos(2*bz/s.betan(n,T)-1)*alpha/2)/(bz/s.betan(n,T))^(alpha/2)*s.betan(n,T)^(-alpha/2)*...
	exp(n*(-1).^(angle(bz/s.betan(n,T)-1) <=0 ).*xin(n,bz/s.betan(n,T),T)+Qx(bz)/2);
% Or use the shorter expression:
% s.pa = @(n,bz,T) gammaDen(n,T)*[1, 0]*R(n,bz/s.betan(n,T),T)*Pinf(bz/s.betan(n,T) )*[s.betan(n,T)^(-alpha/2)*...
% 	exp(n*(-1).^(angle(bz-1) <=0 ).*xin(n,bz/s.betan(n,T),T)+Qx(bz)/2); 0];

s.pb = @(n,bz,T) gammaDen(n,T)*(s.betan(n,T)^(-alpha/2))*[1, 0]*exp(Qx(bz)/2)*R(n,bz/s.betan(n,T),T)/(bz/s.betan(n,T))^(1/4 + alpha/2)*[...
    2^(-alpha)*cos(acos(2*bz/s.betan(n,T)-1)*(1/2+alpha/2) +n*xin(n,bz/s.betan(n,T),T)/1i -pi/4); ...
    -1i*2^(alpha)*cos(acos(2*bz/s.betan(n,T)-1)*(-1/2+alpha/2)+n*xin(n,bz/s.betan(n,T),T)/1i -pi/4)]/(1 - bz/s.betan(n,T))^(1/4);
s.pbwoq = @(n,bz,T) gammaDen(n,T)*[1, 0]*R(n,bz/s.betan(n,T),T)*Pinf(bz/s.betan(n,T) )*...
	[1, 0; (bz/s.betan(n,T))^(-alpha)*exp(-2*n*xin(n,bz/s.betan(n,T),T)), 1]*[s.betan(n,T)^(-alpha/2)*...
	exp(n*xin(n,bz/s.betan(n,T),T) ); 0];

Delta = @(k,z,n,T) 1/2/(-xin(n,z,T))^k*Pinf(z)*s3(z^(alpha/2))*[(-1)^k*s.nuk(k), -6*k*1i*s.nuk(k); ...
    6*k*1i*(-1)^k*s.nuk(k), s.nuk(k)]*s3(z^(-alpha/2))*PinfInv(z);
Deltat = @(k,z,n,T) s.brac(k-1)/(4*sqphitn(n,z,T))^k*Pinf(z)*s3( (-z)^(alpha/2) )*[((-1)^k)/k*(alpha^2+k/2-1/4), (k-1/2)*1i; ...
    -((-1)^k)*(k-1/2)*1i,(alpha^2+k/2-1/4)/k]*s3((-z)^(-alpha/2))*(PinfInv(z));

sR = @(k,z,n,T) Delta(k,z,n,T) - mod(k+1,2)*s.nuk(k)/(-xin(n,z,T))^k*eye(2);
sL = @(k,z,n,T) Deltat(k,z,n,T) -mod(k+1,2)*s.brac(k-1)*(4*s.alpha^2+2*k-1)*eye(2)/2/k/(4*sqphitn(n,z,T))^k;
Rko = @(z,k) sum(s.UL(:,:,k,:)./z.^zpow(:,:,k,:)+s.UR(:,:,k,:)./(z-1).^zpow(:,:,k,:),4);
function RRk = RRkf(z,T,n,usq)
    RRk = zeros(2,2,length(T));
    for ti = 1:length(T)
        if usq || (abs(z-1) < eps^(1/3)) % Too close to the endpoints to employ pole cancellation so use series expansion
            for l = 0:size(Qright,4)-1
                RRk(:,:,ti) = RRk(:,:,ti) + Qright(:,:,ti,l+1)*(z-1)^l;
            end
        else
            RRk(:,:,ti) = Rko(z,T(ti))-sR(T(ti),z,n, max(T)+1); % Take max of T's to avoid using betan with 1 term in sR(1,...)
            for m = 1:T(ti)-1
                RRk(:,:,ti) = RRk(:,:,ti) -Rko(z,T(ti)-m)*sR(m,z,n, max(T)+1);
            end
        end
    end
end
function RLk = RLkf(z,T,n,usq)
    RLk = zeros(2,2,length(T));
    for ti = 1:length(T)
        if usq || (abs(z) < eps^(1/3)) % Too close to the endpoints to employ pole cancellation so use series expansion
            for l = 0:size(Qright,4)-1
                RLk(:,:,ti) = RLk(:,:,ti) + Qleft(:,:,ti,l+1)*z^l;
            end
        else
            RLk(:,:,ti) = Rko(z,T(ti))-sL(T(ti),z,n,max(T)+1);
            for m = 1:T(ti)-1
                RLk(:,:,ti) = RLk(:,:,ti) -Rko(z,T(ti)-m)*sL(m,z,n, max(T)+1);
            end
        end
    end
end
function Rr = RR(n,z,T,usq)
    if powm && usq
        Rr = eye(2);
        for kr = 1:(T-1)
            for l = 0:size(Qright,4)-1
                Rr = Rr + Qright(:,:,kr,l+1)*(z-1)^l*n^(-1-(kr-1)/s.m);
            end
        end
    elseif powm
        Ro = eye(2);
        for kr = 1:(T-1) % Make it possible to get O(n^(-1-1/m) ) by not going to s.m*T
            for m = 1:ceil(kr/s.m)
                Ro = Ro + (s.UL(:,:,kr,m)./z.^m + s.UR(:,:,kr,m)./(z-1).^m)/n^( (kr-1)/s.m+1);
            end
        end
        Rr = eye(2);
        for kr = 1:ceil(T/s.m)
            Rr = Rr + sR(kr,z,n,max(T))/n^kr;
        end
        Rr = Ro*Rr;
    else
        Rr = eye(2) + sum(RRkf(z,1:T-1,n,usq)./n.^repmat(reshape(1:T-1,[1,1,T-1]),[2,2,1]),3);
    end
end
function Rl = RL(n,z,T,usq)
    if powm && usq
        Rl = eye(2);
        for kr = 1:(T-1)
            for l = 0:size(Qleft,4)-1
                Rl = Rl + Qleft(:,:,kr,l+1)*z^l*n^(-1-(kr-1)/s.m);
            end
        end
    elseif powm
        Ro = eye(2);
        for kr = 1:(T-1)
            for m = 1:ceil(kr/s.m)
                Ro = Ro + (s.UL(:,:,kr,m)./z.^m + s.UR(:,:,kr,m)./(z-1).^m)/n^( (kr-1)/s.m+1);
            end
        end
        Rl = eye(2);
        for kr = 1:ceil(T/s.m)
            Rl = Rl + sL(kr,z,n,max(T))/n^kr;
        end
        Rl = Ro*Rl;
    else
        Rl = eye(2) + sum(RLkf(z,1:T-1,n,usq)./n.^repmat(reshape(1:T-1,[1,1,T-1]),[2,2,1]),3);
    end
end

Uc = @(T,n,z,usq) z^(-s.alpha/2)*sqrt(pi)/z^(1/4)/(z-1)^(1/4)*[1,0]*RR(n,z,T,usq)*s3(2^(-s.alpha))*...
    [cos((s.alpha+1)/2*acos(2*z-1) ), -1i*sin((s.alpha+1)/2*acos(2*z-1) ); -1i*cos((s.alpha-1)/2*acos(2*z-1) ), -sin((s.alpha-1)/2*acos(2*z-1) )]*...
    [fn(n,z,T)^(1/4), 0; 0, (-1).^(angle(z-1) <= 0).*fn(n,z,T)^(-1/4)]*[airy(0,fn(n,z,T)); airy(1,fn(n,z,T))];
s.pc = @(n,bz,T,usq) gammaDen(n,T)*s.betan(n,T)^(-s.alpha/2)*exp(Qx(bz)/2)*Uc(T,n,bz/s.betan(n,T),usq);
s.pcwoq = @(n,bz,T,usq) gammaDen(n,T)*s.betan(n,T)^(-s.alpha/2)*Uc(T,n,bz/s.betan(n,T),usq);

Ud = @(T,n,z,usq) (-1)^n*sqrt(pi*1i*n*sqphitn(n,z,T))/z^(1/4)/(1-z)^(1/4)*z^(-alpha/2)*[1,0]*RL(n,z,T,usq)*s3(2^(-s.alpha))*...
    [sin( (s.alpha+1)/2*acos(2*z-1)-pi*s.alpha/2), cos( (s.alpha+1)/2*acos(2*z-1)-pi*s.alpha/2); ... 
    -1i*sin( (s.alpha-1)/2*acos(2*z-1)-pi*s.alpha/2), -1i*cos( (s.alpha-1)/2*acos(2*z-1)-pi*s.alpha/2)]*[besselj(alpha,2i*n*sqphitn(n,z,T)); ...
    (besselj(alpha-1,2i*n*sqphitn(n,z,T)) -alpha/(2i*n*sqphitn(n,z,T) )*besselj(alpha, 2i*n*sqphitn(n,z,T)))];
s.pd = @(n,bz,T,usq) gammaDen(n,T)*s.betan(n,T)^(-s.alpha/2)*exp(Qx(bz)/2)*Ud(T,n,bz/s.betan(n,T),usq);
s.pdwoq = @(n,bz,T,usq) gammaDen(n,T)*s.betan(n,T)^(-s.alpha/2)*Ud(T,n,bz/s.betan(n,T),usq);
end % getAsy
