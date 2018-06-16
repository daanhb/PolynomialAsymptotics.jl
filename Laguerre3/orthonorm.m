% Compute the orthonormal polynomial from its recurrence coefficients.
% Input
%   x, n    - Point and degree at which to evaluate the polynomial
%   aP, bP  - Recurrence coefficients
% Output
%   p       - Value of the orthonormal polynomial
% About
%   Author  - Peter Opsomer (Peter.Opsomer@cs.kuleuven.be)
%   History - Refactored from exactLaguerre in December 2016
function p = orthonorm(x,n,aP,bP)
if n < 0
    p = zeros(size(x));
    return
end
flip = size(x,1) == 1;
if flip
	x = transpose(x);
end
p = [zeros(length(x),1), ones(length(x),1)./bP(1)];
for j=1:n
	p = [p(:,2) ((x-aP(j) ).*p(:,2)-bP(j).*p(:,1))./bP(j+1)];
end
p(:,1) = [];
p(isinf(p)) = realmax;
if flip
	p = transpose(p);
end
end