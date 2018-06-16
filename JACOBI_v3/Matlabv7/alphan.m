% Evaluate the asymptotic expansion for the recursion coefficient alpha_n 
% of the generalised Jacobi polynomial.
% Input
%   n            - Degree at which to evaluate
%   nrT          - Number of terms. 1 gives leading order term
%   Uright,Uleft - Right and left U-matrices
% Output
%   alphan       - Asymptotic expansion of alpha_n
% About
%   Author       - Peter Opsomer (Peter.Opsomer@cs.kuleuven.be)
%   History      - Created May 2014, last edit February 2015
function alphan = alphan(n,nrT,Uright,Uleft)

if (nrT == 1) || (numel(Uright) == 0)
	alphan = 0;
else
	alphan = -sum(reshape(Uright(1,1,1:(nrT-1),1) + Uleft(1,1,1:(nrT-1),1) , ...
		[(nrT-1) 1])./(n+1).^(1:(nrT-1))')   -  sum(reshape(Uright(2,2,1:(nrT-1),1) + ...
		Uleft(2,2,1:(nrT-1),1), [(nrT-1) 1])./n.^(1:(nrT-1) )');
end