% Get the U-matrices to construct the asymptotic expansion of R using the
% procedure with the convolutions as explained in the paper. Optionally,
% specify the method.
% Input
%   s            - Structure containing alpha, m etc
%   maxOrder     - The maximal order of the error (/m if powm)
%   [method      - 'UW'  to get the U-s through the W-s, 
%				   'UQW' to get the U- & Q-s through the W-s and 
%				   'UQV' (default) to get the U- & Q-s through the V-s]
% Output
%   Uright       - Coefficient matrices of R_k(z) for (z-1)^(-m)
%   Uleft        - Coefficient matrices of R_k(z) for z^(-m)
%   [Qright      - Coefficient matrices of R_k^{right}(z) for (z-1)^n]
%   [Qleft       - Coefficient matrices of R_k^{left}(z) for z^n]
% About
%   Author       - Peter Opsomer (Peter.Opsomer@cs.kuleuven.be)
%   History      - Created August 2015 from Jacobi version, corrected October 2016
function [Uright,Uleft,Qright,Qleft] = UQ(s,maxOrder,method)

if ~exist('method','var'), method = 'UQV'; end

powm = isnumeric(s.q) && (norm(s.q(2:end-1)) ~= 0);
fact = 1;
if powm
    fact = s.m;
end
Uright = zeros(2,2,fact*(maxOrder-1), ceil(3*maxOrder/2) );
Uleft = zeros(2,2,fact*(maxOrder-1), ceil(3*maxOrder/2) );
% About half of these tensors will not be used. We could multiply these with nan to see
% which elements are not used, but then using the expressions afterwards will become
% more technical, as would choosing ceil(maxOrder/2) for size(Uleft,4).
if regexp(method,'Q')
    Qright = zeros(2,2,fact*(maxOrder-1),ceil(3*maxOrder/2)+2);
    Qleft = zeros(2,2,fact*(maxOrder-1),ceil(3*maxOrder/2)+2);
end

if regexp(method,'W')
	Wr = WV(s,maxOrder,1,1);
    Wl = WV(s,maxOrder,-1,1);
    if ~powm % 'Monomial' Q(x)
        for k = 1:(maxOrder-1)
            % Actually, WV(:,:,k,-m) but Matlab does not allow negative indices so shift by round(3*k/2)+1
            for m = 1:round(3*k/2)
                Uright(:,:,k,m) = Wr(:,:,k,round(3*k/2)+1-m);
                for j=1:(k-1)
                    for l = max(m-ceil(3*j/2),1):ceil(3*(k-j)/2)
                        Uright(:,:,k,m) = Uright(:,:,k,m) + Uright(:,:,k-j,l)*Wr(:,:,j,round(3*j/2)+1+l-m);
                    end
                end
                for j=1:(k-1)
                    for n = 0:(ceil(3*j/2)-m)
                        for i = 1:ceil((k-j)/2)
                            Uright(:,:,k,m) = Uright(:,:,k,m) + s.binom(-i,n)*Uleft(:,:,k-j,i)*Wr(:,:,j,round(3*j/2)+1-n-m);
                        end
                    end
                end
            end
            % Above right-matrices, below left without zero left-matrices
            for m = 1:round(k/2)
                Uleft(:,:,k,m) = Wl(:,:,k,round(k/2)+1-m);
                for j=1:(k-1)
                    for l = max(m-ceil(j/2),1):ceil((k-j)/2)
                        Uleft(:,:,k,m) = Uleft(:,:,k,m) + Uleft(:,:,k-j,l)*Wl(:,:,j,round(j/2)+1+l-m);
                    end
                end
                for j=1:(k-1)
                    for n = 0:(ceil(j/2)-m)
                        for i = 1:ceil(3*(k-j)/2)
                            Uleft(:,:,k,m) = Uleft(:,:,k,m) + s.binom(-i,n)/(-1)^(i+n)*Uright(:,:,k-j,i)*Wl(:,:,j,round(j/2)+1-n-m);
                        end
                    end
                end
            end
        end
    else % Non-monomial: powers of n^(1/m) present
        for k = 1:s.m*(maxOrder-1)
            for y = 1:ceil(3/2*ceil(k/s.m) )
                Uright(:,:,k,y) = 0; 
                for j = floor(2*y/3):floor((k-1)/s.m)
                    Uright(:,:,k,y) = Uright(:,:,k,y) + Wr(:,:,j+1,round(3*(j+1)/2)+1-y, k-1-j*s.m+1);
                end
                
                for l=0:k-1-s.m
                    for j=0:floor((k-1-l)/s.m-1)
                        for i = max(-ceil(3*j/2),1-y):ceil(3*l/2)-y
                            Uright(:,:,k,y) = Uright(:,:,k,y) + Uright(:,:,l,y+i)*Wr(:,:,j+1,round(3*(j+1)/2)+1+i, k-1-l-s.m*(1+j)+1);
                        end
                    end
                end
                for l=0:k-1-s.m
                    for p = 1:ceil(l/2)
                        tmp = 0;
                        for j=0:floor((k-1-l)/s.m-1)
                            for i = -ceil(3*j/2):-y
                                tmp = tmp + Wr(:,:,j+1,round(3*(j+1)/2)+1+i,k-1-l-s.m*(1+j)+1)*s.binom(-p,-y-i);
                            end
                        end
                        Uright(:,:,k,y) = Uright(:,:,k,y) + Uleft(:,:,l+1,p)*tmp;
                    end
                end
            end
            % Above right, below left without zero left-matrices
            for m = 1:round(ceil(k/s.m)/2)
                Uleft(:,:,k,y) = 0;
                for j = 2*y:floor((k-1)/s.m)
                    Uleft(:,:,k,y) = Uleft(:,:,k,y) + Wl(:,:,j+1,round((j+1)/2)+1-y, k-1-j*s.m+1);
                end
                
                for l=0:k-1-s.m
                    for j=0:floor((k-1-l)/s.m-1)
                        for i = max(-ceil(j/2),1-y):ceil(l/2)-y
                            Uleft(:,:,k,y) = Uleft(:,:,k,y) + Uleft(:,:,l,y+i)*Wl(:,:,j+1,round((j+1)/2)+1+i, k-1-l-s.m*(1+j)+1);
                        end
                    end
                end
                for l=0:k-1-s.m
                    for p = 1:ceil(3*l/2)
                        tmp = 0;
                        for j=0:floor((k-1-l)/s.m-1)
                            for i = -ceil(j/2):-y
                                tmp = tmp + Wl(:,:,j+1,round((j+1)/2)+1+i,k-1-l-s.m*(1+j)+1)*s.binom(-p,-y-i)*(-1)^(-y-i-p);
                            end
                        end
                        Uleft(:,:,k,y) = Uleft(:,:,k,y) + Uright(:,:,l+1,p)*tmp;
                    end
                end
            end
        end
    end
elseif regexp(method,'V')
    Vr = WV(s,maxOrder,1,0);
    Vl = WV(s,maxOrder,-1,0);
else
	error(['Wrong method format, got ' method]);
end

if strcmp(method,'UQW')
    if ~powm % 'Monomial' Q(x)
        % Already did the U's, calculate the Q-s from them:
        for k = 1:(maxOrder-1)
            for n = 0:(ceil(3*(maxOrder-k+2)/2)-1)
                Qright(:,:,k,n+1) = -Wr(:,:,k,round(3*k/2)+1+n);
                for i = 1:round(k/2)
                    Qright(:,:,k,n+1) = Qright(:,:,k,n+1) + s.binom(-i,n)*Uleft(:,:,k,i);
                end
                for j = 1:(k-1)
                    for i = -round(3*j/2):n
                        for l = 1:round((k-j)/2)
                            Qright(:,:,k,n+1) = Qright(:,:,k,n+1) - s.binom(-l,n-i)*Uleft(:,:,k-j,l)*Wr(:,:,j,i+1+round(3*j/2) );
                        end
                    end
                    for i = (n+1):(n+round(3*(k-j)/2) )
                        Qright(:,:,k,n+1) = Qright(:,:,k,n+1)-Uright(:,:,k-j,i-n)*Wr(:,:,j,i+1+round(3*j/2) );
                    end
                end
            end
            % Above right, below left
            for n = 0:(ceil(3*(maxOrder-k+2)/2)-1)
                Qleft(:,:,k,n+1) = -Wl(:,:,k,round(k/2)+1+n);
                for i = 1:round(3*k/2)
                    Qleft(:,:,k,n+1) = Qleft(:,:,k,n+1) + s.pbinom(-i,n)/(-1)^(i+n)*Uright(:,:,k,i);
                end
                for j = 1:(k-1)
                    for i = -round(j/2):n
                        for l = 1:round(3*(k-j)/2)
                            Qleft(:,:,k,n+1) = Qleft(:,:,k,n+1) - s.binom(-l,n-i)/(-1)^(n-i+l)*Uright(:,:,k-j,l)*Wl(:,:,j,i+1+round(j/2) );
                        end
                    end
                    for i = (n+1):(n+round((k-j)/2) )
                        Qleft(:,:,k,n+1) = Qleft(:,:,k,n+1) -Uleft(:,:,k-j,i-n)*Wl(:,:,j,i+1+round(j/2) );
                    end
                end
            end
        end
    else % General polynomial Q(x)
       error('Not implemented UQW yet for general polynomial Q(x).'); 
    end
elseif strcmp(method,'UQV') && ~powm % 'Monomial' Q(x)
    for kt = 0:(maxOrder-2) % Uright(:,:,(maxOrder-1)+1,:) will not be used later on because first term in expansions I is without U's
        for mt = 0:(ceil(3*(kt+1)/2)-1)
            Uright(:,:,kt+1,mt+1) = Vr(:,:,kt+1,ceil(3*(kt+1)/2)-mt);
            for j = 0:(kt-1)
                for l = 0:(ceil(3*(j+1)/2)-mt-1)
                    Uright(:,:,kt+1,mt+1) = Uright(:,:,kt+1,mt+1) + ...
                        Qright(:,:,kt-j,l+1)*Vr(:,:,j+1,ceil(3*(j+1)/2)-l-mt);
                end
            end
        end
        for mt = 0:(ceil((kt+1)/2)-1)
            Uleft(:,:,kt+1,mt+1) = Vl(:,:,kt+1,ceil((kt+1)/2)-mt);
            for j= 0:(kt-1)
                for l = 0:(ceil((j+1)/2)-mt-1)
                    Uleft(:,:,kt+1,mt+1) = Uleft(:,:,kt+1,mt+1) + ...
                        Qleft(:,:,kt-j,l+1)*Vl(:,:,j+1,ceil((j+1)/2)-l-mt);
                end
            end
        end
        for n = 0:(ceil(3*(maxOrder-kt+1)/2)-1)
            Qright(:,:,kt+1,n+1) = -Vr(:,:,kt+1,ceil(3*(kt+1)/2)+1+n);
            Qleft(:,:,kt+1,n+1) = -Vl(:,:,kt+1,ceil((kt+1)/2)+1+n);
            for i = 1:ceil((kt+1)/2)
                Qright(:,:,kt+1,n+1) = Qright(:,:,kt+1,n+1) + s.binom(-i,n)*Uleft(:,:,kt+1,i);
            end
            for i = 1:ceil(3*(kt+1)/2)
                Qleft(:,:,kt+1,n+1) = Qleft(:,:,kt+1,n+1) + s.binom(-i,n)*(-1)^(-i-n)*Uright(:,:,kt+1,i);
            end
            for j = 0:(kt-1)
                for l = 0:(round((j+1)/2)+n)
                    Qleft(:,:,kt+1,n+1) = Qleft(:,:,kt+1,n+1) ...
                        -Qleft(:,:,kt-j,l+1)*Vl(:,:,j+1,n-l+1+ceil((j+1)/2) );
                end
                for l = 0:(round(3*(j+1)/2)+n)
                    Qright(:,:,kt+1,n+1) = Qright(:,:,kt+1,n+1) ...
                        -Qright(:,:,kt-j,l+1)*Vr(:,:,j+1,n-l+1+round(3*(j+1)/2) );
                end
            end
        end
    end
elseif strcmp(method,'UQV') && powm % General polynomial Q(x) 
    for k = 1:s.m*(maxOrder-1)
        for p = 1:ceil(3/2*ceil(k/s.m) )
            Uright(:,:,k,p) = zeros(2,2);
            for q=ceil((2*p-1)/3):((k-1)/s.m+1)
                Uright(:,:,k,p) = Uright(:,:,k,p) + Vr(:,:,q,ceil(3*q/2)+1-p,k-1+s.m-q*s.m +1);
            end
            for q=1:((k-1)/s.m)
                for l = 0:(k-1-s.m*q)
                    for i = -ceil(3*q/2):-p
                        Uright(:,:,k,p) = Uright(:,:,k,p) + Qright(:,:,k-l-s.m*q,-i-p+1)*Vr(:,:,q,ceil(3*q/2)+1+i,l+1);
                    end
                end
            end
        end
        for p = 1:ceil(ceil(k/s.m)/2)
            Uleft(:,:,k,p) = zeros(2,2);
            for q=ceil(2*p-1):((k-1)/s.m+1)
                Uleft(:,:,k,p) = Uleft(:,:,k,p) + Vl(:,:,q,ceil(q/2)+1-p,k-1+s.m-q*s.m +1);
            end
            for q=1:((k-1)/s.m)
                for l = 0:(k-1-s.m*q)
                    for i = -ceil(q/2):-p
                        Uleft(:,:,k,p) = Uleft(:,:,k,p) + ...
                            Qleft(:,:,k-l-s.m*q,-i-p+1)*Vl(:,:,q,ceil(q/2)+1+i,l+1);
                    end
                end
            end
        end
        for n = 0:(ceil(3*(maxOrder-k+1)/2)-1)
            Qright(:,:,k,n+1) = zeros(2,2);
            Qleft(:,:,k,n+1) = zeros(2,2) ;
            for q = 1:((k-1)/s.m+1)
                Qright(:,:,k,n+1) = Qright(:,:,k,n+1) -Vr(:,:,q,ceil(3*q/2)+1+n, k-1+s.m+q*s.m+1);
                Qleft(:,:,k,n+1) = Qleft(:,:,k,n+1) -Vl(:,:,q,ceil(q/2)+1+n, k-1+s.m+q*s.m+1);
            end
            for p = 1:ceil(ceil(k/s.m)/2) 
                Qright(:,:,k,n+1) = Qright(:,:,k,n+1) + Uleft(:,:,k,p)*s.binom(-p,n);
            end
            for p = 1:ceil(ceil(k/s.m)*3/2) 
                Qleft(:,:,k,n+1) = Qleft(:,:,k,n+1) + Uright(:,:,k,p)*s.binom(-p,n)*(-1)^(n-p);
            end
            for q = 1:(k-1)/s.m
                for l = 0:k-1-q*s.m
                    for i = -ceil(3*q/2):n
                        Qright(:,:,k,n+1) = Qright(:,:,k,n+1) -Qright(:,:,k-l-q*s.m,n-i+1)*Vr(:,:,q,ceil(3*q/2)+i+1, l+1);
                    end
                    for i = -ceil(q/2):n
                        Qleft(:,:,k,n+1) = Qleft(:,:,k,n+1) ...
                            -Qleft(:,:,k-l-q*s.m,n-i+1)*Vl(:,:,q,ceil(q/2)+i+1, l+1);
                    end
                end
            end
        end
    end
end
