% Allow the choice of heuristics for the region and the number of terms of the expansions.
% About
%   Author       - Peter Opsomer (Peter.Opsomer@cs.kuleuven.be)
%   History      - Created June 2016, refactored December 2016
%% Initialising
format longe; close all; clear variables;
set(0,'DefaultFigureWindowStyle','docked');

alpha = 0;
q = [0 1];
n = 200;
maxOrder = 5;

[P,gammaP,aP,bP] = exactLaguerre(alpha,q,n);
s = getAsy(alpha,q,maxOrder);


%% Get heuristics for the most accurate expansion for z in [0,1]
zs = linspace(0,1,100); % These can be shifted in the complex plane
ers = zeros(length(zs), 6, maxOrder); % Errors with respect to the exact polynomial
ersRec = zeros(length(zs), 1, maxOrder); % Error on the recurrence relation for standard Laguerre polynomials
[P,gammaP,aP,bP] = exactLaguerre(alpha,q,n);
% 1 is lens (b), 2 is outer (a), 3 is right with Q's (s.pc(...,1)), 4 is left usq (d), 5 right no Q, 6 left s.pd(...,0)
leg = cell(6,maxOrder);
for i = 1:maxOrder
    leg{1,i} = ['lens ' num2str(i)];
    leg{2,i} = ['outer ' num2str(i)];
    leg{3,i} = ['right usq ' num2str(i)];
    leg{4,i} = ['left usq ' num2str(i)];
    leg{5,i} = ['right woq ' num2str(i)];
    leg{6,i} = ['left woq ' num2str(i)];
end
for zi = 1:length(zs)
    pEx = P(zs(zi)*s.betan(n,maxOrder),n +1);
    for i = 1:maxOrder
        ers(zi,1,i) = s.pb(n,zs(zi)*s.betan(n,maxOrder),i);
        ers(zi,2,i) = s.pa(n,zs(zi)*s.betan(n,maxOrder),i);
        ers(zi,3,i) = s.pc(n,zs(zi)*s.betan(n,maxOrder),i,1);
        ers(zi,4,i) = s.pd(n,zs(zi)*s.betan(n,maxOrder),i,1);
        ers(zi,5,i) = s.pc(n,zs(zi)*s.betan(n,maxOrder),i,0);
        ers(zi,6,i) = s.pd(n,zs(zi)*s.betan(n,maxOrder),i,0);
        ersRec(zi,1,i) = abs( ( (zs(zi)*s.betan(n,maxOrder)-(2*n-1+alpha))*real(s.pb(n-1,zs(zi)*s.betan(n,maxOrder),i)) - ...
            sqrt((n-1)*(n-1+alpha))*real(s.pb(n-2,zs(zi)*s.betan(n,maxOrder),i)) )/...
            sqrt((n)*(n+alpha))/real(s.pb(n,zs(zi)*s.betan(n,maxOrder),i)) -1);
    end
    ers(zi,:,:) = abs((pEx-ers(zi,:,:))/pEx);
end
figure; semilogy(real(zs), squeeze(ers(:,1,:)) ); hold on; semilogy(real(zs), squeeze(ers(:,4,:)), '--'); 
ltmp = cell(2*maxOrder,1); ltmp(1:maxOrder) = leg(1,:); ltmp((maxOrder+1):(2*maxOrder)) = leg(4,:); legend(ltmp);
ord = maxOrder; figure; semilogy(real(zs),ers(:,:,ord) ); hold on; semilogy(real(zs), ersRec(:,1,maxOrder) );
legend([leg(:,maxOrder); ['Error on recurrence relation for lens ' num2str(maxOrder)]]);


%% Get heuristics for the number of terms as a function of n
ns = [32, 128, 256, 2.6e3, 1.8e5, 1e11]; % Approximate numbers got for maxP2=8, zc = 0.994; zd = 1e-4;
Ts = [6 5 4 3 2 1];
-Ts.*log(ns)


%% Get heuristics for switching from pd with Qs to pb with T=25/log(n)
ns = round(logspace(2.2, 4.5, 10))';
zsl = logspace(-22,-2,100)';
zsr = 1-zsl;
ers = zeros(length(zsl), 2);
zac = zeros(length(ns), 2);
for ni = 1:length(ns)
    n = ns(ni);
    T = ceil(25/log(n));
    % Find lowest z for which pb gets 11 digits for z close to zero.
    for zi = 1:length(zsl)
        ers(zi,1) = abs( ( (zsl(zi)*s.betan(n,maxOrder)-(2*n-1+alpha))*real(s.pbwoq(n-1,zsl(zi)*s.betan(n,maxOrder),T)) - ...
            sqrt((n-1)*(n-1+alpha))*real(s.pbwoq(n-2,zsl(zi)*s.betan(n,maxOrder),T)) )/...
            sqrt((n)*(n+alpha))/real(s.pbwoq(n,zsl(zi)*s.betan(n,maxOrder),T)) -1);
%         ers(zi,2) = abs( ( (zsr(zi)*s.betan(n,maxOrder)-(2*n-1+alpha))*real(s.pbwoq(n-1,zsr(zi)*s.betan(n,maxOrder),T)) - ...
%             sqrt((n-1)*(n-1+alpha))*real(s.pbwoq(n-2,zsr(zi)*s.betan(n,maxOrder),T)) )/...
%             sqrt((n)*(n+alpha))/real(s.pbwoq(n,zsr(zi)*s.betan(n,maxOrder),T)) -1);
        ers(zi,2) = abs( ( (zsr(zi)*s.betan(n,maxOrder)-(2*n-1+alpha))*real(s.pcwoq(n-1,zsr(zi)*s.betan(n,maxOrder),T,1)) - ...
            sqrt((n-1)*(n-1+alpha))*real(s.pcwoq(n-2,zsr(zi)*s.betan(n,maxOrder),T,1)) )/...
            sqrt((n)*(n+alpha))/real(s.pcwoq(n,zsr(zi)*s.betan(n,maxOrder),T,1)) -1);
    end
    zac(ni,1) = zsl(find(ers(:,1) < 1e-11, 1, 'first'));
%     zac(ni,2) = zsr(find(ers(:,2) < 1e-8, 1, 'last')); % For when using s.pbwoq in ers(zi,2)
    zac(ni,2) = zsr(find(ers(:,2) < 1e-8, 1, 'first'));
end
figure; loglog(ns, zac(:,1) ); hold on; plot(exp(25./(2:5)), exp(mean(log(zac(:,1))))*ones(1,4), 'r*');
idxs = find(zac(:,1) > eps);
s1i = [ones(size(ns(idxs)) ), log(ns(idxs))]\log(zac(idxs,1)) 
hold on; loglog(ns, exp(s1i(1))*ns.^s1i(2));

figure; loglog(ns, 1-zac(:,2) ); hold on; plot(exp(25./(2:5)), 0.05*ones(1,4), 'r*');
