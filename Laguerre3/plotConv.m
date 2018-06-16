% Plot the convergence results with LS convergence in loglog, in one figure.
% Input
%   r       - Actual values
%   a       - Approximations
%   t       - Title of the plot
%   shift   - Degree shift with respect to 2.^(1:maxP2)'
%   s       - Structure containing info on the weight function to determine the slopes
% About
%   Author  - Peter Opsomer (Peter.Opsomer@cs.kuleuven.be)
%   History - Created October 2013, edited February 2015 and December 2016
function plotConv(r,a,t,shift,s)
if ~exist('shift', 'var'), shift=0; end
[maxP2,maxOrder] = size(a);
data =  abs(a - repmat(r,1,maxOrder) )./abs(repmat(r,1,maxOrder) );
legends = cell(maxOrder,1);
ns = 2.^(1:maxP2)'+shift;
[~, inds] = min(data,[], 1);
vals = data(sub2ind(size(data), inds,1:maxOrder) );
clrs = repmat('brgkcm',1,ceil(maxOrder/6));
mark = repmat({'+', 'x', 'o', 'v', '*', 'h', 'd'},1,ceil(maxOrder/7));
lin = repmat({'-','--',':', '-.', '-','--',':','-.'},1,ceil(maxOrder/8));
figure;
for i = 1:maxOrder
    if i == 1
        legends{i} = [num2str(i) ' term'];
    else
        legends{i} = [num2str(i) ' terms'];
    end
    loglog(ns,data(:,i),[clrs(i) mark{i}]); hold on;
end
powm = isnumeric(s.q) && (norm(s.q(2:end-1)) ~= 0);
if powm
    for i=1:maxOrder
        loglog(ns, 2.^( (repmat(log2(ns(inds(i))' ),maxP2,1)-(1:maxP2)' ).*...
            (repmat( (i-1)/s.m+1/2,maxP2,1) ) ).*repmat(vals(i),maxP2,1), [clrs(i) lin{i}] );
    end
else
    for i=1:maxOrder
        loglog(ns, 2.^( (repmat(log2(ns(inds(i))' ),maxP2,1)-(1:maxP2)' ).*...
            (repmat(i-1/2,maxP2,1) ) ).*repmat(vals(i),maxP2,1), [clrs(i) lin{i}] );
    end
end
legend(legends); xlabel('n'); ylabel('Relative error'); 
if ~isempty(t)
    title(t);
end
set(gca, 'XTick', ns);
mir = max(floor(log10(min(min(data) ) ) ),-16);
mar = min(16,ceil(log10(max(max(data) ) ) ));
if mar-mir > 2
    set(gca, 'YTick', 10.^(mir:2:mar));
    axis([1 2*ns(end) min(min(data))/2 2*max(max(data))]);
end
