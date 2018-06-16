% Plot the convergence results with LS convergence in loglog, in 1 figure.
% Input
%   r     - Actual values
%   a     - Approximation
%   t     - Name of the region or coefficient
%   shift - Degree shift with respect to 2.^(1:maxP2)'
% About
%   Author       - Peter Opsomer (Peter.Opsomer@cs.kuleuven.be)
%   History      - Created October 2013, last edit October 2015
function plotConv(r,a,t,shift)
if ~exist('shift'), shift=0; end
[maxP2,maxOrder] = size(a);
data =  abs(a - repmat(r,1,maxOrder) )./abs(repmat(r,1,maxOrder) );
legends = cell(maxOrder,1);
ns = 2.^(1:maxP2)'+shift;
m = '*+o.xsd><phv^';
l =  {'-', '--', ':', '-.'};
c = 'bgrcmk';
figure;
for i = 1:maxOrder
    if i == 1
        legends{i} = [num2str(i) ' term'];
    else
        legends{i} = [num2str(i) ' terms'];
	end
	loglog(ns,data(:,i), [m(mod(i-1,length(m))+1) c(mod(i-1,length(c))+1)]); hold on;
end

inds = maxP2*ones(1,maxOrder);
tmp = sub2ind(size(data), inds, 1:maxOrder);
% Define `accuracy': which lowest relative error the lines have to interpolate
acc = 10^(-11);
ws = find(data(tmp) < acc );
[rowIdx,colIdx] = find(data(1:maxP2, ws) >= acc );
% with accumarray we take the maximum column index for every row
inds(ws) = accumarray(colIdx,rowIdx,[],@max)';
vals = data(sub2ind(size(data), inds,1:maxOrder) );

for i = 1:maxOrder
	loglog(ns, 2.^( (log2(ns(inds(i))' ) -(1:maxP2)' )*i)*vals(i), [l{mod(i-1,length(l))+1} c(mod(i-1,length(c))+1)] );
end

legend(legends); xlabel('n'); 
ylabel(['Relative error for the ' t]);
set(gca, 'XTick', ns);
mir = floor(log10(min(min(data) ) ) );
mar = ceil(log10(max(max(data) ) ) );
set(gca, 'YTick', 10.^(mir:2:mar));
axis([1 2*ns(end) min(min(data))/2 2*max(max(data))]);
