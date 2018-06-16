% Evaluate the asymptotic expansion for the generalised Jacobi polynomial 
% in the lens-shaped region with a certain normalisation. Derivative is only
% computed when assigned, default is monic normalisation without derivatives.
% Input
%   n,z          - Degree and point at which to evaluate
%   alpha,beta,h - Parts of the weight 
%   psi          - Anonymous function for analytic phase function psi(x)
%   nrT          - Number of terms. 1 gives leading order term
%   Dinf         - Limit of the Szego function
%   Uright,Uleft - Right and left U-matrices
%   [nor         - Normalisation: 'm' for monic(default), 'o' for orthonormal]
%   [dh,dpsi     - Derivatives of part of weight and phase function]
% Output
%   ipi          - Asymptotic expansion of the polynomial in the lens
%   [dipi        - Derivative of polynomial in the lens]
% About
%   Author       - Peter Opsomer (Peter.Opsomer@cs.kuleuven.be)
%   History      - Created May 2014, last edit February 2015
function [ipi, dipi] = asy_lens(n,z,alpha,beta,h,psi,nrT,Dinf,Uright,Uleft,nor,dh,dpsi)

if ~exist('nor','var') nor = 'm'; end

if regexp(nor,'o')
	if (nrT == 1) || (numel(Uright) == 0),		nmz = 1;	else
		nmz = sqrt(real(1+2*1i*Dinf^2*sum(reshape(Uright(2,1,1:(nrT-1),1) + ...
			Uleft(2,1,1:(nrT-1),1), [(nrT-1) 1])./(n+1).^(1:(nrT-1) )' ) ) );
	end
end
w = @(x) (1-x).^alpha.*(1+x).^beta.*h(x);
RI = eye(2);
for k = 1:nrT-1
    for m = 1:round(k/2)
        RI = RI + (Uright(:,:,k,m)./(z-1).^m + Uleft(:,:,k,m)./(z+1).^m)/n.^k;
    end
end

if strcmp(nor,'m')
    ipi = [2.^(1/2-n)./(sqrt(w(z) ).*(1+z).^(1/4).*(1-z).^(1/4) ), 0]*RI*...
        [Dinf*cos((n+1/2)*acos(z) +psi(z) -pi/4); -1i/Dinf*cos((n-1/2)*acos(z) +psi(z) -pi/4)];
elseif strcmp(nor,'o')
    ipi = [(2/pi)^(1/2)./(sqrt(w(z) ).*(1+z).^(1/4).*(1-z).^(1/4) ), 0]*RI*...
        [cos((n+1/2)*acos(z) +psi(z) -pi/4); ...
        -1i/Dinf^2*cos((n-1/2)*acos(z) +psi(z) -pi/4)]*nmz;
else
    error('Wrong normalisation')
end

if nargout == 1 return; end

dw = @(x) -alpha*(1-x).^(alpha-1).*(1+x).^beta.*h(x)  +(1-x).^alpha.*beta*(1+x).^(beta-1).*h(x)...
    +(1-x).^alpha.*(1+x).^beta.*dh(x);
dRI = zeros(2);
for k = 1:nrT-1
    for m = 1:round(k/2)
        dRI = dRI + (Uright(:,:,k,m)./(z-1).^(m+1)*(-m) + ...
            Uleft(:,:,k,m)./(z+1).^(m+1)*(-m))/n.^k;
    end
end

if strcmp(nor,'m')
    dipi = ([2.^(1/2-n).*(dw(z)*(-1/2)./(sqrt(w(z) ).^3.*(1+z).^(1/4).*(1-z).^(1/4) ) ...
        -1/4*(-2*z)./(sqrt(w(z) ).*(1+z).^(5/4).*(1-z).^(5/4) ) ), 0]*RI ...
        + [2.^(1/2-n)./(sqrt(w(z) ).*(1+z).^(1/4).*(1-z).^(1/4) ), 0]*dRI)*...
        [Dinf*cos((n+1/2)*acos(z) +psi(z) -pi/4); -1i/Dinf*cos((n-1/2)*acos(z) +psi(z) -pi/4)]...
        + [2.^(1/2-n)./(sqrt(w(z) ).*(1+z).^(1/4).*(1-z).^(1/4) ), 0]*RI*...
        [-Dinf*sin((n+1/2)*acos(z) +psi(z) -pi/4)*(dpsi(z)-(n+1/2)./sqrt(1-z)./sqrt(1+z) ); ...
        1i/Dinf*sin((n-1/2)*acos(z) +psi(z) -pi/4)*(dpsi(z)-(n-1/2)./sqrt(1-z)./sqrt(1+z) )];
elseif strcmp(nor,'o')
    dipi = (([sqrt(2/pi).*(dw(z)*(-1/2)./(sqrt(w(z) ).^3.*(1+z).^(1/4).*(1-z).^(1/4) ) ...
        -1/4*(-2*z)./(sqrt(w(z) ).*(1+z).^(5/4).*(1-z).^(5/4) ) ), 0]*RI ...
        + [sqrt(2/pi)./(sqrt(w(z) ).*(1+z).^(1/4).*(1-z).^(1/4) ), 0]*dRI)*...
        [cos((n+1/2)*acos(z) +psi(z) -pi/4); -1i/Dinf^2*cos((n-1/2)*acos(z) +psi(z) -pi/4)]...
        + [sqrt(2/pi)./(sqrt(w(z) ).*(1+z).^(1/4).*(1-z).^(1/4) ), 0]*RI*...
        [-sin((n+1/2)*acos(z) +psi(z) -pi/4)*(dpsi(z)-(n+1/2)./sqrt(1-z)./sqrt(1+z) ); ...
        1i/Dinf^2*sin((n-1/2)*acos(z) +psi(z) -pi/4)*(dpsi(z)-(n-1/2)./sqrt(1-z)./sqrt(1+z) )])*nmz;
end