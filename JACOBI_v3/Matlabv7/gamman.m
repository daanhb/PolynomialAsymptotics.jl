% Evaluate the asymptotic expansion for the leading order coefficient gamma_n 
% of the generalised orthonormal Jacobi polynomial.
% Input
%   n            - Degree at which to evaluate
%   nrT          - Number of terms. 1 gives leading order term
%   Dinf         - Limit of the Szego function
%   Uright,Uleft - Right and left U-matrices
% Output
%   gamma        - Asymptotic expansion of gamma_n
% About
%   Author       - Peter Opsomer (Peter.Opsomer@cs.kuleuven.be)
%   History      - Created May 2014, last edit February 2015
function gamma = gamman(n,nrT,Dinf,Uright,Uleft)

if (nrT == 1) || (numel(Uright) == 0)
	gamma = 2.^(n)/Dinf/sqrt(pi);
else
	gamma = 2.^(n)/Dinf/sqrt(pi)*sqrt(real(1+2*1i*Dinf^2*sum(reshape(Uright(2,1,1:(nrT-1),1) + ...
		Uleft(2,1,1:(nrT-1),1), [(nrT-1) 1])./(n+1).^(1:(nrT-1) )' ) ) );
end