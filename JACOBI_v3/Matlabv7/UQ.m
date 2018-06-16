% Get the U-matrices to construct the asymptotic expansion of R using the
% procedure with the convolutions as explained in the paper. Optionally,
% specify method or number of matrices to compute or get the Q-s.
% Input
%   alpha, beta  - Parts of the weight w(x) = (1-x).^alpha.*(1+x).^beta.*h(x)
%   Dinf         - Limit of the Szego function
%   c,d          - Coefficients needed for the U-matrices
%   maxOrder     - The maximum order of the error
%   [method      - 'UW' (default) to get the U-s through the W-s, 
%				   'UQW' to get the U- & Q-s through the W-s and 
%				   'UQV' to get the U- & Q-s through the V-s]
%   [ns          - If method contains 'Q', specify how many Q-s to compute
%                  so compute Q_{k,0:(ns(k)-1)}^{right/left}, also for W/V]
% Output
%   Uright       - Coefficient matrices of R_k(z) for (z-1)^(-m)
%   Uleft        - Coefficient matrices of R_k(z) for (z+1)^(-m)
%   [Qright      - Coefficient matrices of R_k^{right}(z) for (z-1)^n]
%   [Qleft       - Coefficient matrices of R_k^{left}(z) for (z+1)^n]
% About
%   Author       - Peter Opsomer (Peter.Opsomer@cs.kuleuven.be)
%   History      - Created January 2014, last edit February 2015
function [Uright,Uleft,Qright,Qleft] = UQ(alpha,beta,Dinf,c,d,maxOrder,method,ns)

if ~exist('method') method = 'UW'; end

existed = exist('ns','var');
if ~existed
    if regexp(method,'Q')
        ns = max([ (maxOrder-2-(1:(maxOrder-1))+1+round((1:(maxOrder-1))/2) ); ...
            zeros(1,maxOrder-1)],[],1);
    else 
        % Could use the following to specify the necessary (exact) number 
        % of W- or V-matrices to compute, but not really needed...
%         ns = round( (maxOrder-1-(1:(maxOrder-1)))/2)-1+round((1:(maxOrder-1))/2);
    end
end
Uright = zeros(2,2,maxOrder-1,ceil( (maxOrder-1)/2) );
Uleft = zeros(2,2,maxOrder-1,ceil( (maxOrder-1)/2));
% About half of these tensors will not be used
if regexp(method,'Q')
    if existed
        Qright = zeros(2,2,maxOrder-1,max(ns) );
        Qleft = zeros(2,2,maxOrder-1,max(ns) );
    else
        Qright = zeros(2,2,maxOrder-1,ceil( (maxOrder-1)/2) );
        Qleft = zeros(2,2,maxOrder-1,ceil( (maxOrder-1)/2));
    end
end

poch = @(x,n) prod(x+(0:(n-1) ) ); % 'pochhammer' needs symbolic toolbox, 
% could also use separate 1-line script or global function

if regexp(method,'W')
	if exist('ns','var')
		Wr = WV(alpha,beta,Dinf,c,maxOrder,1,1,ns+2+round( (1:(maxOrder-1)) /2) );
		% Wr = W_{k,m}^{right}, k=1..maxOrder-1, m given by ns; Wl analogous.
		Wl = WV(beta,alpha,Dinf,d,maxOrder,-1,1,ns+2+round( (1:(maxOrder-1)) /2));
	else
		Wr = WV(alpha,beta,Dinf,c,maxOrder,1,1);
		Wl = WV(beta,alpha,Dinf,d,maxOrder,-1,1);
	end
	for k = 1:(maxOrder-1)
		% Actually, for m = (-round(k/2)):1
		% but Matlab does not allow negative indices so shift by round(k/2)+1
		for m = 1:round(k/2)
			Uright(:,:,k,m) = Wr(:,:,k,round(k/2)+1-m);
			Uleft(:,:,k,m) = Wl(:,:,k,round(k/2)+1-m);
			for j=1:(k-1)
				for l = max(m-ceil(j/2),1):ceil((k-j)/2)
					Uright(:,:,k,m) = Uright(:,:,k,m) + ...
						Uright(:,:,k-j,l)*Wr(:,:,j,round(j/2)+1+l-m);
					Uleft(:,:,k,m) = Uleft(:,:,k,m) + ...
						Uleft(:,:,k-j,l)*Wl(:,:,j,round(j/2)+1+l-m);
				end
			end
			for j=1:(k-1)
				for n = 0:(ceil(j/2)-m)
					for i = 1:ceil((k-j)/2)
						Uright(:,:,k,m) = Uright(:,:,k,m) + ...
							poch(1-i-n,n)/2^(i+n)/factorial(n)*...
							Uleft(:,:,k-j,i)*Wr(:,:,j,round(j/2)+1-n-m);
						Uleft(:,:,k,m) = Uleft(:,:,k,m) + ...
							poch(1-i-n,n)/(-2)^(i+n)/factorial(n)*...
							Uright(:,:,k-j,i)*Wl(:,:,j,round(j/2)+1-n-m);
					end
				end
			end
		end
	end
elseif regexp(method,'V')
	if exist('ns','var')
		Vr = WV(alpha,beta,Dinf,c,maxOrder,1,0,ns+2+round( (1:(maxOrder-1)) /2) );
		% Vr = V_{k,m}^{right}, k=1..maxOrder-1, m given by ns; Vl analogous.
		Vl = WV(beta,alpha,Dinf,d,maxOrder,-1,0,ns+2+round( (1:(maxOrder-1)) /2) );
	else
		Vr = WV(alpha,beta,Dinf,c,maxOrder,1,0);
		Vl = WV(beta,alpha,Dinf,d,maxOrder,-1,0);
	end
else
	error(['Wrong method format, got ' method]);
end

if strcmp(method,'UQW')
    % Already did the U's, calculate the Q-s from them:
    for k = 1:(maxOrder-1)
        for n = 0:ns(k)
            Qright(:,:,k,n+1) = -Wr(:,:,k,round(k/2)+1+n);
            Qleft(:,:,k,n+1) = -Wl(:,:,k,round(k/2)+1+n);
            for i = 1:round(k/2)
                Qright(:,:,k,n+1) = Qright(:,:,k,n+1) + poch(1-i-n,n)/2^(i+n)...
                    /factorial(n)*Uleft(:,:,k,i);
                Qleft(:,:,k,n+1) = Qleft(:,:,k,n+1) + poch(1-i-n,n)/(-2)^(i+n)...
                    /factorial(n)*Uright(:,:,k,i);
            end
            for j = 1:(k-1)
                for i = -round(j/2):n
                    for l = 1:round((k-j)/2)
                        Qright(:,:,k,n+1) = Qright(:,:,k,n+1) - poch(1-l-n+i,n-i)/2^(n-i+l)...
                            /factorial(n-i)*Uleft(:,:,k-j,l)*Wr(:,:,j,i+1+round(j/2) );
                        Qleft(:,:,k,n+1) = Qleft(:,:,k,n+1) - poch(1-l-n+i,n-i)/(-2)^(n-i+l)...
                            /factorial(n-i)*Uright(:,:,k-j,l)*Wl(:,:,j,i+1+round(j/2) );
                    end
                end
                for i = (n+1):(n+round((k-j)/2) )
                    Qright(:,:,k,n+1) = Qright(:,:,k,n+1) ...
                        -Uright(:,:,k-j,i-n)*Wr(:,:,j,i+1+round(j/2) );
                    Qleft(:,:,k,n+1) = Qleft(:,:,k,n+1) ...
                        -Uleft(:,:,k-j,i-n)*Wl(:,:,j,i+1+round(j/2) );
                end
            end
        end
    end
elseif strcmp(method,'UQV')
	for k = 1:(maxOrder-1)
		for m = 1:round(k/2)
			Uright(:,:,k,m) = Vr(:,:,k,round(k/2)+1-m);
			Uleft(:,:,k,m) = Vl(:,:,k,round(k/2)+1-m);
			for j=1:(k-1)
				for l = 0:(ceil(j/2)-m)
					Uright(:,:,k,m) = Uright(:,:,k,m) + ...
						Qright(:,:,k-j,l+1)*Vr(:,:,j,round(j/2)+1-l-m);
					Uleft(:,:,k,m) = Uleft(:,:,k,m) + ...
						Qleft(:,:,k-j,l+1)*Vl(:,:,j,round(j/2)+1-l-m);
				end
			end
		end
		for n = 0:ns(k)
			Qright(:,:,k,n+1) = -Vr(:,:,k,round(k/2)+1+n);
			Qleft(:,:,k,n+1) = -Vl(:,:,k,round(k/2)+1+n);
			for i = 1:round(k/2)
				Qright(:,:,k,n+1) = Qright(:,:,k,n+1) + poch(1-i-n,n)/2^(i+n)...
					/factorial(n)*Uleft(:,:,k,i);
				Qleft(:,:,k,n+1) = Qleft(:,:,k,n+1) + poch(1-i-n,n)/(-2)^(i+n)...
					/factorial(n)*Uright(:,:,k,i);
			end
			for j = 1:(k-1)
				for l = 0:(round(j/2)+n)
					Qright(:,:,k,n+1) = Qright(:,:,k,n+1) ...
						-Qright(:,:,k-j,l+1)*Vr(:,:,j,n-l+1+round(j/2) );
					Qleft(:,:,k,n+1) = Qleft(:,:,k,n+1) ...
						-Qleft(:,:,k-j,l+1)*Vl(:,:,j,n-l+1+round(j/2) );
				end
			end
		end
	end
end
