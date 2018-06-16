% Evaluate the asymptotic expansion for the recursion coefficient beta_n 
% of the generalised Jacobi polynomial.
% Input
%   n            - Degree at which to evaluate
%   nrT          - Number of terms. 1 gives leading order term
%   Dinf         - Limit of the Szego function
%   Uright,Uleft - Right and left U-matrices
% Output
%   betan        - Asymptotic expansion of beta_n
% About
%   Author       - Peter Opsomer (Peter.Opsomer@cs.kuleuven.be)
%   History      - Created May 2014, last edit February 2015
function betan = betan(n,nrT,Dinf,Uright,Uleft)

if (nrT == 1) || (numel(Uright) == 0)
	betan = 1/4;
else
	betan = (1/(2i*Dinf^2) + sum(reshape(Uright(2,1,1:(nrT-1),1) + Uleft(2,1,1:(nrT-1),1),...
		[(nrT-1) 1])./n.^(1:(nrT-1))') )*(-Dinf^2/2i + sum(reshape(Uright(1,2,1:(nrT-1),1) + ...
		Uleft(1,2,1:(nrT-1),1), [(nrT-1) 1])./n.^(1:(nrT-1))') );
end