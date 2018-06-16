% Perform tests and make plots related to Gauss-Jacobi quadrature
% About
%   Author       - Peter Opsomer (peteropsomer@gmail.com)
%   History      - Created 2017, finished April 2018

%% Initialising
close all; clear variables; 
format longe; clc;
set(0, 'DefaultFigureWindowStyle', 'docked');
addpath(genpath('/../JACOBI_v3/Matlabv7/')); % Add the path to Jacobi Matlabv7 
% Add the path to chebfun: addpath(genpath('/dir-to/chebfun/'))

%% Test accuracy and efficiency
% alpha=0; beta=0;
% alpha = -0.5; beta = 0.7;
alpha = 0.42; beta = -1/sqrt(5);
tims = nan(500, 8);
ers = nan(500, 16);
prev = cputime; 
printtoc = 3;
doRec = 1;
for n = 100:10^10 % Can interrupt this loop and use its saved results in the next execution cells
    if doRec
        tic; start = cputime;
        [xrec,wrec] = jacptsTest(n,alpha,beta, 'rec');
        tims(n,5) = toc; tims(n,6) = cputime-start;
    end
    tic; start = cputime;
    [xe,we] = jacptsTest(n,alpha,beta, [-1,1], 'exp');
    tims(n,3) = toc; tims(n,4) = cputime-start;
    
    tic; start = cputime;
    [x,w] = jacptsTest(n,alpha,beta, [-1,1], 'asy');
    tims(n,1) = toc; tims(n,2) = cputime-start;
    
    ers(n,1) = max(abs(x-xe)./x);
    ers(n,2) = max(abs((w-we)./w));
    ers(n,3) = max(abs(x-xe));
    ers(n,4) = max(abs(w-we));
    if doRec
        ers(n,5) = max(abs(xrec-xe)./xrec);
        ers(n,6) = max(abs((wrec-we)./wrec));
        ers(n,7) = max(abs(xrec-xe));
        ers(n,8) = max(abs(wrec-we));
        
        ers(n,9) = max(abs(x-xrec)./x);
        ers(n,10) = max(abs((w-wrec)./w));
        ers(n,11) = max(abs(x-xrec));
        ers(n,12) = max(abs(w-wrec));
    end
    if rem(log2(n),1) == 0 % Takes a lot of time so only at these n
        tic; start = cputime;
        [xrh,wrh] = jacptsTest(n,alpha,beta, 'rh');
        tims(n,7) = toc; tims(n,8) = cputime-start;
        
        ers(n,13) = max(abs(x-xrh)./x);
        ers(n,14) = max(abs((w-wrh)./w));
        ers(n,15) = max(abs(x-xrh));
        ers(n,16) = max(abs(w-wrh));
    end
    if cputime > prev+printtoc
        save timErrGJ.mat tims ers alpha beta
        disp([num2str(n) '=n']);
        prev = cputime;
        if n > 3e3
            doRec = 0;
        end
    end
end


%% Plot the errors from the previous test
% load timErrGJ.mat
% shift = 0; %EXP vs ASY
% shift = 4; %EXP vs REC
shift = 8; % REC vs ASY
figure;
loglog(100:size(ers,1), ers(100:end,shift+1), 'k', 'LineWidth', 2);
hold on;
loglog(100:size(ers,1), ers(100:end,shift+2), 'r-.', 'LineWidth', 2);
loglog(100:size(ers,1), ers(100:end,shift+3), 'b:', 'LineWidth', 2);
loglog(100:size(ers,1), ers(100:end,shift+4), 'g--', 'LineWidth', 2);
xlabel('n');
ylabel('Maximal error');
legend({'Rel. error $x_k$', 'Rel. error $w_k$', 'Abs. error $x_k$', 'Abs. error $w_k$'}, 'interpreter', 'latex');
set(gca, 'FontSize', 20);

% Plot the error of RH
figure;
loglog(100:size(ers,1), ers(100:end,13), 'k*', 'LineWidth', 2);
hold on;
loglog(100:size(ers,1), ers(100:end,14), 'r.', 'LineWidth', 2);
loglog(100:size(ers,1), ers(100:end,15), 'bo', 'LineWidth', 2);
loglog(100:size(ers,1), ers(100:end,16), 'gh', 'LineWidth', 2);
xlabel('n');
ylabel('Maximal error');
legend({'Rel. error $x_k$', 'Rel. error $w_k$', 'Abs. error $x_k$', 'Abs. error $w_k$'}, 'interpreter', 'latex');
set(gca, 'FontSize', 20);


% Plot the time
figure;
loglog(100:size(tims,1), tims(100:end,1), 'b--', 'LineWidth', 2);
hold on;
loglog(100:size(tims,1), tims(100:end,3), 'm-.', 'LineWidth', 2);
loglog(100:size(tims,1), tims(100:end,5), 'g-', 'LineWidth', 2);
loglog(100:size(tims,1), tims(100:end,7), 'k*', 'LineWidth', 2);
xlabel('n');
ylabel('Wallclock time (s)');
legend({'ASY', 'EXP', 'REC', 'RH'}, 'interpreter', 'latex');
set(gca, 'FontSize', 20);

%return



%% Make heuristics: choose the loop over alphas or over Ts
ns = round(logspace(1.65, 2.5, 20))';
T = 2; beta=0; alphas = [-0.92, -0.5, -0.23, 0, 0.5, 1, 2.1]; rests = nan(length(ns),length(alphas),12); for idx = 1:length(alphas), alpha=alphas(idx); P = exactPolys(alpha, beta, @(x) ones(size(x)), max(ns) );
% alpha = -0.4; beta=0.7; P = exactPolys(alpha, beta, @(x) ones(size(x)), max(ns) ); Ts = 2:2:8; rests = nan(length(ns),length(Ts),12); for idx=1:length(Ts), T=Ts(idx);
  for ni = 1:length(ns)
    n = ns(ni);
    [x,w] = jacptsTest(n,alpha,beta); % Hale and Townsend fast asy algo
    [bes,wbes] = jacptsTest(n,alpha,beta,[-1,1], 'exp', [T,n,n+1]);
    [bulk,wbulk] = jacptsTest(n,alpha,beta,[-1,1], 'exp', [T,0,n+1]);
    [air,wair] = jacptsTest(n,alpha,beta,[-1,1], 'exp', [T,0,1]);
    % Use bes, ..., wair to compute heuristics
    rests(ni,idx,5) = find(abs((bulk-x)./x) < abs((bes-x)./x), 1, 'first');
    rests(ni,idx,6) = find(abs((bulk-x)./x) < abs((air-x)./x), 1, 'last');
    
    rests(ni,idx,9) = find(abs((wbulk-w)./w) < abs((wbes-w)./w), 1, 'first');
    rests(ni,idx,10) = find(abs((wbulk-w)./w) < abs((wair-w)./w), 1, 'last');
    for k = 1:n
        bes(k) = abs(P(bes(k), n+1)); %Absolute value of the orthonormal poly should be small at the zero
        bulk(k) = abs(P(bulk(k), n+1));
        air(k) = abs(P(air(k), n+1));
    end
    rests(ni, idx, 1) = find(bulk < bes, 1, 'first');
    rests(ni, idx, 2) = find(bulk < air, 1, 'last');
    repr = [1,2, 5,6, 9,10];
    rests(ni, idx, repr+2) = x(rests(ni,idx,repr));
  end
end

% Plot the heuristics
plAll = 0;
for typ = repr %1:4
% for typ = (repr+2)
    fct = @(dat) ns-dat;
    if mod(typ,2)
        fct = @(dat) dat;
    end
    figure;
    loglog(ns, fct(rests(:, 1, typ)), 'k+');
    hold on;
    if plAll, loglog(ns, fct(rests(:, 2, typ)), 'bo', 'MarkerSize', 2); end
    loglog(ns, fct(rests(:, 3, typ)), 'gx');
    loglog(ns, fct(rests(:, 4, typ)), 'm.');
    if size(rests,2) == 4 % Iteration over number of terms
        legs = cell(5,1);
        for ai = 1:length(Ts)
            legs{ai} = ['$T$ = ' num2str(Ts(ai))];
        end
    else
        if plAll, loglog(ns, fct(rests(:, 5, typ)), 'yd', 'MarkerSize', 2); end
        loglog(ns, fct(rests(:, 6, typ)), 'cs');
        loglog(ns, fct(rests(:, 7, typ)), 'rp');
        
        legs = cell(8,1);
        for ai = 1:length(alphas)
            legs{ai} = ['$\alpha$ = ' num2str(alphas(ai))];
        end
    end
    
    if mod(typ,2)
        loglog(ns, sqrt(ns), 'r'); legs{end} = '$\sqrt{n}$'; heu = '$k$';
    else
        loglog(ns, fct(ns-sqrt(ns)), 'r'); legs{end} = '$\sqrt{n}$'; heu = '$n-k$';
    end
    if plAll
        legend(legs, 'interpreter', 'latex', 'location', 'best');
    else
        legend(legs{[1, 3:4, 6:8]}, 'interpreter', 'latex', 'location', 'best');
    end
    xlabel('n');
    if typ <= 4
        ylabel(['Transition ' heu ' for lowest $|p_n(x_k)|$'], 'interpreter', 'latex');
    elseif typ <= 8
        ylabel(['Transition ' heu ' for lowest relative error of $x_k$'], 'interpreter', 'latex');
    elseif typ <= 12
        ylabel(['Transition ' heu ' for lowest relative error of $w_k$'], 'interpreter', 'latex');
    else
        error('not defined');
    end
    set(gca, 'FontSize', 20);
end

% Choose one previous plot
if size(rests,2) == 4
    figure(1); typ = 9; fct = @(dat) dat; heu = '$k$'; legs{end} = '$\sqrt{n}$';
    ylabel(['Transition ' heu], 'interpreter', 'latex');

    loglog(ns, fct(rests(:, 1, typ)), 'k-');
    loglog(ns, fct(rests(:, 2, typ)), 'b:');
    loglog(ns, fct(rests(:, 3, typ)), 'g-.');
    loglog(ns, fct(rests(:, 4, typ)), 'm--');
    legend({legs{[1:5, 1:4]}});%, 'interpreter', 'latex', 'location', 'best');
    yticks([2,4,7,10,20]); xticks([50,70,100,140, 200, 300, 500]); %xticks(ns(1:4:length(ns)) );
    set(gca, 'FontSize', 20);
elseif 1
    figure(2); typ = 10; fct = @(dat) ns-dat; heu = '$n-k$'; legs{end} = '$\sqrt{n}$';
    ylabel(['Transition ' heu], 'interpreter', 'latex');
    
    loglog(ns, fct(rests(:, 1, typ)), 'k-');
    if plAll, loglog(ns, fct(rests(:, 2, typ)), 'b:'); end
    loglog(ns, fct(rests(:, 3, typ)), 'g-.');
%    loglog(ns, fct(rests(:, 4, typ)), 'm--');
    loglog(ns, fct(rests(:, 4, typ)), 'm--', 'LineWidth', 2);
    if plAll, loglog(ns, fct(rests(:, 5, typ)), 'y-'); end
    loglog(ns, fct(rests(:, 6, typ)), 'c-');
    loglog(ns, fct(rests(:, 7, typ)), 'r-.');
    if plAll 
        legend(legs{[1:8, 1:7]}, 'interpreter', 'latex', 'location', 'best');
    else
        legend({legs{[1, 3:4, 6:8, 1,3:4,6:7]}}, 'interpreter', 'latex', 'location', 'best');
    end
    yticks([2,4,7,10,20]); xticks([50,70,100,140, 200, 300, 500]); 
    set(gca, 'FontSize', 20);
end

%return



%% Gaussian quadrature: Newton procedure to compute nodes and weights, from Matlabv7/illustrateFunctionality.m 
alpha = 1/sqrt(3); beta = -1/pi; h = @(x) exp(x); dh = @(x) exp(x);
n = 200;
maxOrder = 5;
[c, d, Dinf, psi, dpsi, contpsi] = contour_integrals(alpha,beta,h,maxOrder);
[Uright,Uleft] = UQ(alpha,beta,Dinf,c,d,maxOrder,'UW');
[P,gammaP,alphaP,betaP] = exactPolys(alpha,beta,h,n);
printtoc = 2;

zw = zeros(n, maxOrder+1);
sqbetaP = sqrt(betaP);
legs = cell(maxOrder,1);

tic;
prevToc = toc;
for ni = 1:n
    if(toc - prevToc > printtoc)
        prevToc = toc;
        display([num2str(ni) ' = ni, est sec left = ' num2str(prevToc*(n-ni+1)/(ni-1))]);
    end
    if ni <=3 % (Approximation of) zero of bessel function to approximate first zeros
        zw(ni,1) = -1 + (2*sqrt(beta+1) + (ni-1)*pi)^2/2/n^2;
        if zw(ni,1) +1 < eps^(1/3)
            warning('Use Q-s');
        end
    else
        zw(ni,1) = 3*zw(ni-1) -3*zw(ni-2) +zw(ni-3); % Quadratic extrapolation
    end
    for attempt = 1:10
        [p, dp] = jacpnRecDer(zw(ni,1), n, alphaP, sqbetaP);
        zw(ni,1) = zw(ni,1) - p/dp;
        if abs(p/dp) < 1e-19*zw(ni,1)%Added
            break;
        end
    end
    
    if (zw(ni,1) < -1) || (zw(ni,1) > 1)
        error([num2str(ni) ' = ni, Exact polynomial gives zero outside of the interval']);
    elseif (ni > 1) && (zw(ni,1) < zw(ni-1,1) )
        error([num2str(ni) ' = ni, Exact zeros not going right']);
    end
    
    % Could improve the following by stopping Newton iterates when already
    % tried that zero and then going over all (max. 5) closest floats.
    for nrT = 1:maxOrder
        legs{nrT} = [num2str(nrT) ' terms'];
        zw(ni,nrT+1) = zw(ni,1);
        for attempt = 1:10
            if zw(ni,nrT+1)+1 < 0.2 % This bound is not general for all weights!
                [poly, dpoly] = asy_left(n,zw(ni,nrT+1),alpha,beta,h,psi,nrT,Dinf,Uright,Uleft,'o',dh,dpsi);
                polyM1 = asy_left(n-1,zw(ni,nrT+1),alpha,beta,h,psi,nrT,Dinf,Uright,Uleft,'o');
                polyP1 = asy_left(n+1,zw(ni,nrT+1),alpha,beta,h,psi,nrT,Dinf,Uright,Uleft,'o');
            elseif 1-zw(ni,nrT+1) < 0.2
                [poly, dpoly] = asy_right(n,zw(ni,nrT+1),alpha,beta,h,psi,nrT,Dinf,Uright,Uleft,'o',dh,dpsi);
                polyM1 = asy_right(n-1,zw(ni,nrT+1),alpha,beta,h,psi,nrT,Dinf,Uright,Uleft,'o');
                polyP1 = asy_right(n+1,zw(ni,nrT+1),alpha,beta,h,psi,nrT,Dinf,Uright,Uleft,'o');
            else
                [poly, dpoly] = asy_lens(n,zw(ni,nrT+1),alpha,beta,h,psi,nrT,Dinf, Uright,Uleft,'o',dh,dpsi);
                polyM1 = asy_lens(n-1,zw(ni,nrT+1),alpha,beta,h,psi,nrT,Dinf, Uright,Uleft,'o');
                polyP1 = asy_lens(n+1,zw(ni,nrT+1),alpha,beta,h,psi,nrT,Dinf, Uright,Uleft,'o');
            end
            zw(ni,nrT+1) = zw(ni,nrT+1) - real(poly)/real(dpoly); % We know the polynomial is real on the interval
            if abs(poly/dpoly) < 1e-19*zw(ni,nrT+1)%Added
                break;
            end
        end
        if (zw(ni,nrT+1) < -1) || (zw(ni,nrT+1) > 1)
            warning([num2str(ni) ' = ni, Outside of the interval, nrT = ' num2str(nrT)]);
        elseif (ni > 1) && (zw(ni,nrT+1) < zw(ni-1,nrT+1) )
            warning([num2str(ni) ' = ni, Not going right, nrT = ' num2str(nrT)]);
        end
    end
    debug = 1;
end
legs{1} = '1 term';
% save('GJsq3piexp.mat')

figure; semilogy(abs((zw(:,2:end) - repmat(zw(:,1),1,maxOrder) )./repmat(zw(:,1),1,maxOrder))); 
ylabel('Relative error on $x_k$', 'interpreter', 'latex');
xlabel('k'); set(gca, 'FontSize', 20); 
legend(legs);

% Integrate a polynomial with a degree low enough to be integrated exactly
fct = @(x) 7.3*x.^2 -23.5*cos(74*acos(x) );
w = @(x) (1-x).^alpha.*(1+x).^beta.*h(x);
exVal = integral(@(x) w(x).*fct(x),-1,1)
abserrorThree = sum(zw(:,2).*fct(zw(:,1)) ) -exVal
relerr = abserrorThree/exVal

% Test by integrating (x^2-1) and (x^3) over the interval
abserror = sum(zw(:,2).*(zw(:,1).^2-1) ) -(-1.293225037)
abserror = sum(zw(:,2).*(zw(:,1).^3) ) - (-.1659585353)
    
% Integrate a random polynomial with a degree low enough to be integrated exactly
r = [randn(2,1); randi([0 2*n-3],1); randi([0 5],1)];
fct = @(x) 7*r(1)*x.^r(3) -23*r(2)*x.^r(4);
w = @(x) (1-x).^alpha.*(1+x).^beta.*h(x);
abserror = sum(zw(:,2).*fct(zw(:,1)) ) -integral(@(x) w(x).*fct(x),-1,1)


%% Plot results Newton on recur and asy: test explicit expansions of nodes (and weights) for this generalised weight function
% load GJsq3piexp.mat;

ms = 9; %Max size
legExpls = cell(ms,1);
for nrT = 1:ms
    legExpls{nrT} = ['T = ' num2str(nrT)];
end
besR = nan(size(zw,1), ms); wbesR = nan(size(zw,1), ms);

[p, dp] = jacpnRecDer(zw(:,1), n, alphaP, sqbetaP);
weiRec = gammaP(n+1)*2^(1+n-n)/gammaP(n+0)./dp./P(zw(:,1),n-1+1);
relErrWei = (sum(weiRec)-integral(@(x) (1-x).^alpha.*(1+x).^beta.*h(x), -1, 1))/sum(weiRec)

% First both Bessel regions
figure;
fs = 14;
for endpt = [1, -1]
    bes = -ones(size(zw,1), ms); wbes = zeros(size(zw,1), ms);
    if endpt == 1
        [c, d, Dinf, psi, dpsi,contpsi] = contour_integrals(beta,alpha, @(x) h(-x), ms-1);
        alphat = alpha; betat = beta; % Temporarily switch these
        alpha = betat; beta = alphat;
    else
        [c, d, Dinf, psi, dpsi,contpsi] = contour_integrals(alpha,beta,h,ms-1);
    end
    c0 = c(1); c1 = c(2); d0 = d(1); d1 = d(2); d2 = d(3);
    z = 1/(2*n+alpha+beta+1-d0);
    % series node bessel
    jbk = besselroots(beta,n);
    bes(:,8:end) = bes(:,8:end) -1/90720*(576*jbk.^8 - 1152*(7*alpha^2 + 5*beta^2 + 490*d0^2 - 2940*d0*d1 + 4410*d1^2 ...
        - 3)*jbk.^6 + 16*(1312*beta^4 + 4*(1512*alpha^2 - 575)*beta^2 + 945*(4*alpha^2 ...
        -1)*c0^2 - 7560*(4*alpha^2 - 1)*c0*d0 - 945*(36*alpha^2 + 32*beta^2 - 17)*d0^2 ...
        - 68040*(4*beta^2 - 1)*d1^2 + 2268*alpha^2 + 5670*(5*(4*alpha^2 - 1)*c0 + (20*alpha^2 ...
        + 32*beta^2 - 13)*d0)*d1 - 452)*jbk.^4 - (181440*alpha^6 + 15808*beta^6 ...
        + 90048*(3*alpha^2 - 1)*beta^4 + 45360*(4*alpha^2 - 1)*c0^4 + 181440*(4*alpha^2 ...
        - 1)*c0*d0^3 + 45360*(4*alpha^2 - 1)*d0^4 - 544320*alpha^4 + 1344*(405*alpha^4 ...
        - 600*alpha^2 + 133)*beta^2 + 945*(336*alpha^4 + 112*(4*alpha^2 - 1)*beta^2 ...
        - 88*alpha^2 + 1)*c0^2 - 68040*(16*alpha^4 - 40*alpha^2 + 9)*c0*c1 + 315*(2160*alpha^4 ...
        + 704*beta^4 + 32*(72*alpha^2 - 53)*beta^2 + 864*(4*alpha^2 - 1)*c0^2 - 3384*alpha^2 ...
        + 947)*d0^2 + 11340*(176*beta^4 - 280*beta^2 + 59)*d1^2 + 536256*alpha^2 ...
        + 1890*(96*(4*alpha^2 - 1)*c0^3 + (528*alpha^4 + 152*(4*alpha^2 - 1)*beta^2 ...
        - 608*alpha^2 + 119)*c0 - 36*(16*alpha^4 - 40*alpha^2 + 9)*c1)*d0 ...
        - 3780*(15*(4*(4*alpha^2 - 1)*beta^2 - 4*alpha^2 + 1)*c0 + (352*beta^4 + 20*(12*alpha^2 ...
        - 31)*beta^2 - 60*alpha^2 + 133)*d0)*d1 - 104512)*jbk.^2)*z^8;
    
    bes(:,7:end) = bes(:,7:end) + 1/1440*(256*(13*d0 -75*d1 + 90*d2)*jbk.^6 ...
        + 48*(5*(4*alpha^2 - 1)*c0 - (180*alpha^2 + 32*beta^2 - 53)*d0 +150*(4*alpha^2 - 1)*d1 + 60*(4*beta^2 - 1)*d2)*jbk.^4 ...
        - (720*(4*alpha^2 - 1)*c0^3 + 2160*(4*alpha^2 - 1)*c0*d0^2 + 720*(4*alpha^2 - 1)*d0^3 + 15*(336*alpha^4 + ...
        88*(4*alpha^2 - 1)*beta^2 - 352*alpha^2 + 67)*c0 - 270*(16*alpha^4 - 40*alpha^2 + ...
        9)*c1 + (6480*alpha^4 + 1792*beta^4 + 40*(252*alpha^2 - 95)*beta^2 + 2160*(4*alpha^2 ...
        - 1)*c0^2 - 10080*alpha^2 + 2323)*d0 - 60*(32*beta^4 + 20*(12*alpha^2 + 1)*beta^2 ...
        -60*alpha^2 - 7)*d1 - 540*(16*beta^4 - 40*beta^2 + 9)*d2)*jbk.^2)*z^7;
    
    bes(:,6:end) = bes(:,6:end) + 1/90*(8*jbk.^6 - 12*(5*alpha^2 + 3*beta^2 - 2)*jbk.^4 ...
        + (180*alpha^4 + 28*beta^4 + 80*(3*alpha^2 - 1)*beta^2 + 45*(4*alpha^2 - 1)*c0^2 ...
        + 90*(4*alpha^2 - 1)*c0*d0 + 45*(4*alpha^2 - 1)*d0^2 - 240*alpha^2 + 52)*jbk.^2)*z^6;
    
    bes(:,5:end) = bes(:,5:end) - 1/6*(16*(d0 - 3*d1)*jbk.^4 + (3*(4*alpha^2 - 1)*c0 + ...
        (12*alpha^2 + 8*beta^2 - 5)*d0 - 6*(4*beta^2 - 1)*d1)*jbk.^2)*z^5; % Actually need repmat but Matlab seems to cope...
    bes(:,4:end) = bes(:,4:end) - 2/3*(jbk.^4 - (3*alpha^2 + beta^2 - 1)*jbk.^2)*z^4;
    bes(:,2:end) = bes(:,2:end) + 2*jbk.^2*z^2;
    bes = real(bes);
    % Weights in the left bessel region
    wbes(:,7:end) = wbes(:,7:end) +  1/360*(768*(13*d0 - 75*d1 + 90*d2)*jbk.^4 - 720*(4*alpha^2 - 1)*c0^3 - 2160*(4*alpha^2 ...
        - 1)*c0*d0^2 - 720*(4*alpha^2 - 1)*d0^3 + 96*(5*(4*alpha^2 - 1)*c0 - (180*alpha^2 + ...
        32*beta^2 - 53)*d0 + 150*(4*alpha^2 - 1)*d1 + 60*(4*beta^2 - 1)*d2)*jbk.^2 ...
        - 15*(336*alpha^4 + 88*(4*alpha^2 - 1)*beta^2 - 352*alpha^2 + 67)*c0 + 270*(16*alpha^4 ...
        - 40*alpha^2 + 9)*c1 - (6480*alpha^4 + 1792*beta^4 + 40*(252*alpha^2 - 95)*beta^2 ...
        + 2160*(4*alpha^2 - 1)*c0^2 - 10080*alpha^2 + 2323)*d0 + 60*(32*beta^4 + 20*(12*alpha^2 ...
        + 1)*beta^2 - 60*alpha^2 - 7)*d1 + 540*(16*beta^4 - 40*beta^2 + 9)*d2)*z^7; 
    wbes(:,6:end) = wbes(:,6:end) + 2/45*(180*alpha^4 + 28*beta^4 + 24*jbk.^4 + 80*(3*alpha^2 - 1)*beta^2 + 45*(4*alpha^2 ...
        - 1)*c0^2 + 90*(4*alpha^2 - 1)*c0*d0 + 45*(4*alpha^2 - 1)*d0^2 - 24*(5*alpha^2 ...
        + 3*beta^2 - 2)*jbk.^2 - 240*alpha^2 + 52)*z^6;
    wbes(:,5:end) = wbes(:,5:end) - 2/3*(32*(d0 - 3*d1)*jbk.^2 +3*(4*alpha^2 - 1)*c0 ...
        + (12*alpha^2 + 8*beta^2 - 5)*d0 - 6*(4*beta^2 -1)*d1)*z^5;
    wbes(:,4:end) = wbes(:,4:end) + 8/3*(3*alpha^2 + beta^2 - 2*jbk.^2 - 1)*z^4;
    wbes(:,2:end) = wbes(:,2:end) + 8*z^2; 
    if endpt == 1
        pos = 5;
        bes = -flipud(bes); 
        alpha = alphat; beta = betat;
        wbes = real(flipud(wbes).*(1-bes).^alpha.*(1+bes).^beta.*h(bes).*flipud(besselj(alpha-1,jbk).^(-2)) );
        wbesR = wbes; besR = bes;
    else
        pos = 1
        wbes = real(wbes.*(1-bes).^alpha.*(1+bes).^beta.*h(bes).*besselj(beta-1,jbk).^(-2));
    end
    subplot(3,2,pos)
    data = abs((bes(:,[2,4:8]) - zw(:,1) )./repmat(zw(:,1),1,6));
    semilogy(data, 'LineWidth', 2);
    axis([1,n,min(min(data)),max(max(data))]);
    ylabel(['Error on $x_k$ near ' num2str(endpt)], 'interpreter', 'latex');
    xlabel('k'); set(gca, 'FontSize', fs);
    legend(legExpls{[2,4:8]}, 'location','best');
    set(gca,'YTick', 10.^(round(linspace(ceil(log10(min(min(data+eps)))), floor(log10(max(max(data)))), 4))) );
    
    subplot(3,2,pos+1); 
    data = abs((wbes(:,[2,4:7]) - weiRec )./repmat(weiRec,1,5));
    semilogy(data, 'LineWidth', 2);
    axis([1,n,min(min(data)),max(max(data))]);
    ylabel(['Error on $w_k$ near ' num2str(endpt)], 'interpreter', 'latex'); 
    xlabel('k'); set(gca, 'FontSize', fs);
    legend(legExpls{[2,4:7]}, 'location','best');
    set(gca,'YTick', 10.^(round(linspace(ceil(log10(min(min(data+eps)))), floor(log10(max(max(data)))), 4))) );
end

% Now the bulk region
ms = 4;
bulk = zeros(size(zw,1), ms);
wbulk = zeros(size(zw,1), ms);

for k = 1:n
    t = cos((4*n -4*k +2*alpha +3)/(4*n +2*alpha +2*beta +2)*pi ); %initial guess from standard Jacobi
    % Could change this to Newton iterations with contpsi for computing t faster
    t = fzero(@(t) real(acos(t) -pi*(4*(n-k)+2*alpha+3)/(4*n+2*alpha+2*beta+2) + sqrt(1-t.^2)/(2*n+alpha+beta+1)*contpsi(t,1)/2i/pi), t); %pi*(4*k+2*alpha+3)/(4
    h0 = trap_rule(@(z) log(h(z ) )./( sqrt(z -1).*sqrt(z +1).*(z -t).^1 ), 1.25);
    h1 = trap_rule(@(z) log(h(z ) )./( sqrt(z -1).*sqrt(z +1).*(z -t).^2 ), 1.25);
    h2 = trap_rule(@(z) log(h(z ) )./( sqrt(z -1).*sqrt(z +1).*(z -t).^3 ), 1.25);
    h3 = nan;
    z = 1/(2*n+alpha+beta+1-h0);
    
    bulk(k,4:end) = bulk(k,4:end) + 1/24*(12*(2*alpha^2 + 2*beta^2 - 1)*h1^2*t^7 ...
        + 24*((2*alpha^2 + 2*beta^2 - 1)*h0*h1 + (alpha^2 - beta^2)*h1^2)*t^6 + 6*(2*(2*alpha^2 + 2*beta^2 - 1)*h0^2 - 6*(2*alpha^2 + 2*beta^2 - 1)*h1^2 ...
        + ((4*alpha^2 - 1)*c0 - (4*beta^2 - 1)*d0 + 2*(10*alpha^2 + 2*beta^2 - 3)*h0)*h1)*t^5 + 6*(2*(8*alpha^2 + 4*beta^2 - 3)*h0^2 ...
        - 12*(alpha^2 - beta^2)*h1^2 + ((4*alpha^2 - 1)*c0 - (4*beta^2 - 1)*d0)*h0 + ((4*alpha^2 - 1)*c0 + (4*beta^2 - 1)*d0 ...
        - 4*(alpha^2 + 7*beta^2 - 2)*h0)*h1)*t^4 - 32*alpha^4 + 32*beta^4 + (16*alpha^4 + 16*beta^4 + 4*(12*alpha^2 - 5)*beta^2 ...
        + 6*(4*alpha^2 - 1)*c0^2 + 6*(4*beta^2 - 1)*d0^2 + 24*(5*alpha^2 - beta^2 - 1)*h0^2 + 36*(2*alpha^2 + 2*beta^2 - 1)*h1^2 - 20*alpha^2 ...
        + 12*(2*(4*alpha^2 - 1)*c0 - (4*beta^2 - 1)*d0)*h0 - 12*((4*alpha^2 - 1)*c0 - (4*beta^2 - 1)*d0 + 2*(10*alpha^2 + 2*beta^2 - 3)*h0)*h1 + 5)*t^3 ...
        - 6*(4*alpha^2 - 1)*c0^2 + 6*(4*beta^2 - 1)*d0^2 - 72*(alpha^2 - beta^2)*h0^2 - 24*(alpha^2 - beta^2)*h1^2 + 6*((4*alpha^2 - 1)*c0^2 ...
        - (4*beta^2 - 1)*d0^2 - 2*(2*alpha^2 + 10*beta^2 - 3)*h0^2 + 12*(alpha^2 - beta^2)*h1^2 + 4*alpha^2 - 4*beta^2 + 2*((4*alpha^2 - 1)*c0 ...
        + 2*(4*beta^2 - 1)*d0)*h0 - 2*((4*alpha^2 - 1)*c0 + (4*beta^2 - 1)*d0 + 2*(4*alpha^2 - 8*beta^2 + 1)*h0)*h1)*t^2 + 40*alpha^2 - 40*beta^2 ...
        - 18*((4*alpha^2 - 1)*c0 + (4*beta^2 - 1)*d0)*h0 + 6*((4*alpha^2 - 1)*c0 + (4*beta^2 - 1)*d0 + 12*(alpha^2 - beta^2)*h0)*h1 ...
        - 3*(16*alpha^4 + 16*beta^4 + 4*(4*alpha^2 - 7)*beta^2 + 2*(4*alpha^2 - 1)*c0^2 + 2*(4*beta^2 - 1)*d0^2 + 12*(4*alpha^2 - 1)*h0^2 ...
        + 4*(2*alpha^2 + 2*beta^2 - 1)*h1^2 - 28*alpha^2 + 4*(2*(4*alpha^2 - 1)*c0 - (4*beta^2 - 1)*d0)*h0 - 2*((4*alpha^2 - 1)*c0 - (4*beta^2 - 1)*d0 ...
        + 2*(10*alpha^2 + 2*beta^2 - 3)*h0)*h1 + 11)*t)*z^4/(t^2 - 1);
    bulk(k,3:end) = bulk(k,3:end) - 1/4*(2*(2*alpha^2 + 2*beta^2 - 1)*h1*t^3 + 2*((2*alpha^2 + 2*beta^2 - 1)*h0 + 2*(alpha^2 - beta^2)*h1)*t^2 ...
        + (4*alpha^2 - 1)*c0 + (4*beta^2 - 1)*d0 + 8*(alpha^2 - beta^2)*h0 - 4*(alpha^2 - beta^2)*h1 + ((4*alpha^2 - 1)*c0 - (4*beta^2 - 1)*d0 ...
        + 4*(3*alpha^2 + beta^2 - 1)*h0 - 2*(2*alpha^2 + 2*beta^2 - 1)*h1)*t)*z^3;
    bulk(k,2:end) = bulk(k,2:end) + 1/2*(2*alpha^2 - 2*beta^2 + (2*alpha^2 + 2*beta^2 - 1)*t)*z^2;
    bulk(k,:) = bulk(k,:) + t; 
    
    z = 1/(2*n+alpha+beta+1);
    wbulk(k,4:end) = wbulk(k,4:end) +(-2*h1^3*t^6 - 6*h0*h1^2*t^5 - 6*(h0^2*h1 - h1^3)*t^4 - 2*(h0^3 - 6*h0*h1^2 + (2*alpha^2 + 2*beta^2 -1)*h2)*t^3 ...
        + 2*h1^3 - (6*h1^3 + (10*alpha^2 + 10*beta^2 - 6*h0^2 - 5)*h1 + 4*(alpha^2 - beta^2)*h2)*t^2 ...
        -1/2*(4*alpha^2 - 1)*c0 + 1/2*(4*beta^2 - 1)*d0 - 2*(alpha^2 - beta^2)*h0 + 2*(2*alpha^2 + 2*beta^2 - 1)*h1 ...
        +4*(alpha^2 - beta^2)*h2 - (6*h0*h1^2 + 3*(2*alpha^2 + 2*beta^2 - 1)*h0 + 6*(alpha^2 - beta^2)*h1 -2*(2*alpha^2 + 2*beta^2 - 1)*h2)*t)*z^3;
    wbulk(k,3:end) = wbulk(k,3:end) + (2*h1^2*t^4 + 4*h0*h1*t^3 - 4*h0*h1*t + 2*(h0^2 - 2*h1^2)*t^2 + 2*alpha^2 + 2*beta^2 + 2*h1^2 - 1)*z^2;
    wbulk(k,2:end) = wbulk(k,2:end) +(-2*h1*t^2 - 2*h0*t + 2*h1)*z;
    wbulk(k,:) = real((1-bulk(k,:)).^alpha.*(1+bulk(k,:)).^beta.*h(bulk(k,:)).*pi.*sqrt(1-t^2).*z.*(2+wbulk(k,:)) );
end

subplot(3,2,3);
data = abs((bulk(:,1:4) - zw(:,1) )./repmat(zw(:,1),1,4));
semilogy(data, 'LineWidth', 2);
axis([1,n,min(min(data)),max(max(data))]);
ylabel('Error on $x_k$ in bulk', 'interpreter', 'latex');
xlabel('k'); set(gca, 'FontSize', fs);
legend(legExpls{[2:5]}, 'location','best');
set(gca,'YTick', 10.^(round(linspace(ceil(log10(min(min(data+eps)))), floor(log10(max(max(data)))), 4))) );

subplot(3,2,4);
data = abs((wbulk(:,1:4) - weiRec)./repmat(weiRec,1,4));
semilogy(data, 'LineWidth', 2);
axis([1,n,min(min(data)),max(max(data))]);
ylabel('Error on $w_k$ in bulk', 'interpreter', 'latex');
xlabel('k'); set(gca, 'FontSize', fs);
legend(legExpls{1:4}, 'location','best');
set(gca,'YTick', 10.^(round(linspace(ceil(log10(min(min(data+eps)))), floor(log10(max(max(data)))), 4))) );



function [p dp] = jacpnRecDer(x,n,alphaP,sqbetaP)
pnprev = 0;
p = 1/sqbetaP(1);
pndprev = 0;
dp = 0;
for j=1:n
    pnold = p;
	p = ((x-alphaP(j) ).*pnold -sqbetaP(j).*pnprev)./sqbetaP(j+1);
    pnprev = pnold;
    pndold = dp;
    dp = ((x-alphaP(j) ).*pndold + pnold -sqbetaP(j).*pndprev)./sqbetaP(j+1);
    pndprev = pndold;
end

end

