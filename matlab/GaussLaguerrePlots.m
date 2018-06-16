% Perform tests and make plots related to Gauss-Laguerre quadrature
% About
%   Author       - Peter Opsomer (peteropsomer@gmail.com)
%   History      - Created June 2016, refactored December 2016, finished April 2018

%% Initialising
format longe; close all; clear variables;
set(0,'DefaultFigureWindowStyle','docked');

addpath(genpath('/../Laguerre3/')) % Add the path to Laguerre3: 
% Add the path to chebfun for the test at the end: addpath('dir-to/chebfun')


%% Laguerre-type
alpha = 0; 
q = @(x) exp(x); dQ = @(x) exp(x); Qinv = @(x) log(x);

maxAt = 11;
ns = 100; %ns = [0; round(logspace(1.65,2.5,13))']; 
Ts = 1:4; 
errs = nan(length(ns),max(ns),length(Ts),9);
[P,gammaP,aP,bP] = exactLaguerre(alpha,q,max(ns) );
legs = cell(length(Ts),1);

rho = @(z) 1.5*max([0.5,abs(z/2+1/4 -[0,1])]); % Make sure z lies within the contour
M = 200; % should probably increase with (n and) the desired accuracy
ths = linspace(0,2*pi,M+1); % = theta_k, k = 0 to M-1
ths(end) = [];
hh = ths(2)-ths(1);
trap_rule = @(f,z) hh*sum(f(rho(z).*exp(1i*ths)+z/2+1/4).*rho(z).*1i.*exp(1i*ths) )/(2*pi*1i );

for idx=1:length(Ts)
  T=Ts(idx);
  legs{idx} = [num2str(T) ' = T'];
  for ni = 1:length(ns)
    n = ns(ni);
    d = 1/(4*n+2*alpha+2);
    
    s = getAsy(alpha,q,max(Ts), dQ, Qinv, n);
    cs = real(arrayfun(@(nn) s.trap_rule(@(y) sqrt(y).*s.dVn(y)./( sqrt(y -1).*(y-1).^nn), 1), 1:(ceil(3*max(Ts)/2)+5) ));
    ds = real(arrayfun(@(nn) s.trap_rule(@(y) sqrt(y).*s.dVn(y)./( sqrt(y -1).*y.^nn), 0), 1:(ceil(max(Ts)/2)+5) )); 
    d0 = ds(1); d1 = ds(2);
    c0 = cs(1); c1 = cs(2); c2 = cs(3);
    
    bet = fzero(@(bet) quadgk(@(x) dQ(x).*sqrt(x./(bet-x)), 0, bet)/2/pi/n-1, Qinv(n) );
    bets = fzero(@(bet) quadgk(@(x) dQ(x).*sqrt(x./(bet-x)), 0, bet)/2/pi/(n-1)-1, bet);
    hn = @(z) arrayfun(@(zz) trap_rule(@(y) sqrt(y).*bet.*dQ(bet*y)/n./( sqrt(y -1).*(y - zz) ), zz), z); % Need to recompute for each n
    xin = @(z) arrayfun(@(zz) -1/2*quadgk(@(x) sqrt(x-1)./sqrt(x).*hn(x), 1, zz), z);
    zt = 0.5; % Test-z to calculate l_n (3.24)
    lln = 2*quadgk(@(y) log(abs(zt-y)).*sqrt(1-y)./sqrt(y)/2/pi.*hn(y), 0, 1) -q(bet*zt)/n;
    llns = lln+q(bet*zt)/n -q(bets*zt)/(n-1);
    
    pt = (4*n -4*transpose(1:n) +3)*d;
    t = pi^2/16*(pt -1).^2;
    for it = 1:6
        % Some iterations to cheaply get at least some initialisation for the actual Newton iterations: not close to the actual t at all
        t = t - (pt*pi +2*sqrt(t-t.^2) -acos(2*t-1) ).*sqrt(t./(1-t))/2; 
    end
    if 0 % Because these iterations to compute t take too long as mentioned in the article
    figure;
    rang = transpose(1:2); %1:n;
    for it = 1:6
%         f = -3*pi +4*(1:n)*pi-4*n*pi +acos(2*t-1)*(alpha+1)*2 -n/2i*xin(t);
        if 0
            f = real(-3*pi +4*transpose(1:n)*pi-4*n*pi +acos(2*t-1)*(alpha+1)*2 -n/2i*xin(t));
            df = real(-2*(alpha+1)./sqrt(t)./sqrt(1-t) -n/2i*sqrt(t-1)./sqrt(t).*hn(t));
        t = t -f./df;
        else
            f = real(-3*pi +4*rang*pi-4*n*pi +acos(2*t(rang)-1)*(alpha+1)*2 -n/2i*xin(t(rang)));
            df = real(-2*(alpha+1)./sqrt(t(rang))./sqrt(1-t(rang)) -n/2i*sqrt(t(rang)-1)./sqrt(t(rang)).*hn(t(rang)));
            t(rang) = t(rang) -f./df;
        end
        disp([num2str(it) ' iter, norm = ' num2str(norm(f))]);
        semilogy(abs(f) ); hold on;
    end
    legend(1:6)
    end
    
    jak = besselroots(alpha, n);
    bes = jak*0;
    wbes = 0*bes;
    bulk = t*0;
    wbulk = t*0;
    
    ak = [-13.69148903521072; -12.828776752865757; -11.93601556323626;    -11.00852430373326; -10.04017434155809; ...
        -9.02265085340981; -7.944133587120853;    -6.786708090071759; -5.520559828095551; -4.08794944413097; -2.338107410459767];
    tair = 3*pi/2*( transpose(n:-1:12)-0.25);
    ak = [-tair.^(2/3).*(1 + 5/48./tair.^2 - 5/36./tair.^4 + 77125/82944./tair.^6 -10856875/6967296./tair.^8); ak];
    
    if T >= 4
        bes = bes + 16/3*(6*(alpha^3 + 3*alpha^2 + 3*alpha + 1)*c0^4*d0^4 - 3*(28*alpha^3 + 60*alpha^2 +45*alpha + 13)*c0^4*d0 ...
            - 3*(12*(alpha^3 + 3*alpha^2 + 3*alpha + 1)*c0^4 -4*(2*alpha^3 + 6*alpha^2 + 5*alpha + 1)*c0^3 ...
            + 2*(2*alpha^3 + 6*alpha^2 + 3*alpha -1)*c0^2 - 6*(alpha + 1)*c1^2 + 5*(alpha + 1)*c0*c2 + (2*(alpha + 1)*c0^2 ...
            + (4*alpha^3+ 12*alpha^2 + 5*alpha - 3)*c0)*c1)*d0^3 + (4*(22*alpha^3 + 58*alpha^2 + 53*alpha +17)*c0^4 ...
            - 3*(16*alpha^3 + 40*alpha^2 + 29*alpha + 5)*c0^3 + 9*(alpha +1)*c0^2*c1)*d0^2 + 8*((alpha + 1)*c0^4*d0^2 ...
            - 3*(alpha + 1)*c0^4*d0 - (2*(alpha +1)*c0^4*d0 - 5*(alpha + 1)*c0^4)*d1)*jak.^2 - (2*(4*alpha^3 + 4*alpha^2 - alpha -1)*c0^4*d0 ...
            - 5*(4*alpha^3 + 4*alpha^2 - alpha - 1)*c0^4)*d1)/(c0^4*d0^4)*d^3;
    end
    if T >= 3
        bes = bes + 4/3*(9*(alpha^2 + 2*alpha + 1)*c0^2*d0^3 + 2*(22*alpha^2 + 36*alpha + 17)*c0^2*d0 -(4*alpha^2 - 1)*c0^2*d1 ...
            - 3*(12*(alpha^2 + 2*alpha + 1)*c0^2 - 2*(2*alpha^2 + 4*alpha+ 1)*c0 + c1)*d0^2 + 4*(c0^2*d0 - 2*c0^2*d1)*jak.^2)/(c0^2*d0^3)*d^2;
    end
    if T >= 2
        bulk = bulk + 2*(alpha + 1)*d*t;
        bes = bes + 4*((alpha + 1)*ds(1) - 2*alpha - 2)/ds(1)*d;   
        wbes = wbes -4*((alpha + 1)*c0 - (alpha + 1)*d0)/(c0*d0)*d; 
    end
    bes = 16*jak.^2./ds(1)^2.*d^2.*bet.*(1+bes);
    bulk = bulk + t*bet;
    air = bet*(1 + (-2/3/n/(-c0/3))^(2/3)*ak);    
%     wbes = d*bes.^alpha.*exp(-q(bes))./(besselj(alpha-1, jak)).^2.*( 512/((alpha + 3)*d0^3 - 4*(alpha + 1)*d0^2) + wbes);    
    wbes = bes.^alpha.*exp(-q(bes))./(besselj(alpha-1, jak)).^2*256./ds(1)^3.*d.*(1+wbes);
    
    wbulk = bulk.^alpha.*exp(-bulk)*2*pi.*sqrt(t./(1-t)).*(1 +wbulk);
%     wair = 4^(1/3)*air.^(alpha+1/3).*exp(-air)./(airy(1,ak)).^2;
    wair = real(1/(2*exp((n-1)*log(bet) +(-n-alpha)*log(bets) +n*lln/2 -(n-1)*llns/2))*((bet-bets)/bets./ak + (-2/3/n/(-c0/3))^(2/3)).^(+1/4).* ...
        exp(-q(air)).*air.^alpha.*(airy(1,ak)).^(-2).*(-3*(-c0/3)/2)^(-3/2)/sqrt(n)); 

    xcorr = zeros(n,1); wcorr = zeros(n,1);
    nstop = n;
    for k = 1:n
        if k <= 4
            xcorr(k) = bes(k);
        else
            xcorr(k) = 3*xcorr(k-1) -3*xcorr(k-2) +xcorr(k-3); % Quadratic extrapolation
        end
        ov = inf;
        for attempt = 1:maxAt
            [pe, dp] = jacpnRecDer(xcorr(k), n, aP, bP);
            step = pe/dp;
            xcorr(k) = real(xcorr(k) - step);
            if isnan(xcorr(k) ) || ((abs(pe) > abs(ov) +1e-12) && (k ~= nstop)) % Recur just one iter increase
                error('Converging outside interval or away from zero')
            elseif abs(step) < 1e-12*xcorr(k) % Converged
                % Lowering this bound would only be meaningful if the expansions are made to reach machine precision
                break;
            elseif attempt == maxAt
                % Poly in lens does not converge for first zeros
                error('No convergence')
            end
            ov = pe;
        end
        if ((k > 1) && (xcorr(k) <= xcorr(k-1))) || (xcorr(k) < 0)
            warning([num2str(ni) ' = ni, Node ' num2str(xcorr(k)) ' should be increasing at k = ' num2str(k)]);
        end
        wcorr(k) = real(gammaP(n+1)/gammaP(n)/dp/jacpnRecDer(xcorr(k),n-1,aP,bP));
        
        errs(ni,k,idx,1) = abs(P(bes(k), n+1)); %Absolute value of the orthonormal poly should be small at the zero
        errs(ni,k,idx,2) = abs(bes(k) -xcorr(k))/xcorr(k);
        errs(ni,k,idx,3) = abs(wbes(k) -wcorr(k))/wcorr(k);
        errs(ni,k,idx,4) = abs(P(bulk(k), n+1));
        errs(ni,k,idx,5) = abs(bulk(k) -xcorr(k))/xcorr(k);
        errs(ni,k,idx,6) = abs(wbulk(k) -wcorr(k))/wcorr(k);
        errs(ni,k,idx,7) = abs(P(air(k), n+1));
        errs(ni,k,idx,8) = abs(air(k) -xcorr(k))/xcorr(k);
        errs(ni,k,idx,9) = abs(wair(k) -wcorr(k))/wcorr(k);
    end
  end
    
end
errCorrWeights = [(sum(wcorr) - integral(@(x) x.^alpha.*exp(-q(x)), 0, inf))/sum(wcorr), ...
    (sum(xcorr.*wcorr) - integral(@(x) x.^(alpha+1).*exp(-q(x)), 0, inf))/sum(xcorr.*wcorr), ...
    (sum(xcorr.^5.*wcorr) - integral(@(x) x.^(alpha+5).*exp(-q(x)), 0, inf))/sum(xcorr.^5.*wcorr)]

fs = 14;
legExpls = cell(length(Ts),1);
legExpls{1} = 'T = 1';
for nrT = 2:length(Ts)
    legExpls{nrT} = ['T = ' num2str(Ts(nrT) )];
end
figure; subplot(1,2,1);
semilogy(squeeze(1:ns(end)), squeeze(errs(end,:,:,2)), 'LineWidth', 2);
axis([1,n,min(min(errs(end,:,:,2))),max(max(errs(end,:,:,2)))]);
xlabel('k'); set(gca, 'FontSize', fs);
ylabel(['Error on $x_k$ near ' num2str(0)], 'interpreter', 'latex'); 
legend(legExpls, 'location','best'); 
set(gca,'YTick', 10.^(round(linspace(ceil(log10(min(min(squeeze(errs(end,:,:,2))+eps)))), floor(log10(max(max(squeeze(errs(end,:,:,2))) ))), 4))) );

subplot(1,2,2);
semilogy(squeeze(1:ns(end)), squeeze(errs(end,:,1,8)), 'LineWidth', 2);
axis([1,n,min(min(errs(end,:,1,8))),max(max(errs(end,:,1,8)))]);
xlabel('k'); set(gca, 'FontSize', fs);
ylabel(['Error on $x_k$ near $\beta_n$'], 'interpreter', 'latex'); 
legend(legExpls{1}, 'location','best');
set(gca,'YTick', 10.^(round(linspace(ceil(log10(min(min(squeeze(errs(end,:,1,8))+eps)))), floor(log10(max(max(squeeze(errs(end,:,1,8))) ))), 4))) );

%return




%% Test accuracy and efficiency: This is not plotted in the paper although there is a similar result for Jacobi there
% alpha = -0.5;
alpha = 0; % Fast algo only for alpha=0
tims = nan(500, 6);
ers = nan(500, 4);
prev = cputime; 
printtoc = 3;
for n = 100:10^10 % Can interrupt this loop and plot its saved results
%     tic; start = cputime;
%     [x,w] = lagptsTest(n, [0, inf], 'RH', alpha);
%     tims(n,5) = toc; tims(n,6) = cputime-start;
    
    tic; start = cputime;
    [xe,we] = lagptsTest(n,[0,inf], 'exp', alpha);
    tims(n,3) = toc; tims(n,4) = cputime-start;
    
    tic; start = cputime;
%     [x,w] = lagptsTest(n,[0,inf], 'fast', alpha);
    [x,w] = lagptsTest(n,[0,inf], 'GW', alpha);
    tims(n,1) = toc; tims(n,2) = cputime-start;
    
%     ers(n,1) = max(abs((x-xe)./(x+3*eps)));
    ers(n,1) = max(abs(x-xe)./x);
    ers(n,2) = max(abs((w-we)./w));
    ers(n,3) = norm(x-xe)/norm(x);
    ers(n,4) = norm(w-we)/norm(w);
    if cputime > prev+printtoc
        save timErrGL.mat tims ers
        disp([num2str(n) '=n']);
        prev = cputime;
    end
end

%return




%% Heuristics for nodes using abs value of orthogonal polynomial computed by recurrence at the asy exp of the zero

ns = [0; round(logspace(1.65,2.5,33))'];

% T = 8; alpha = 0; rests = nan(length(ns),length(alphas),4); for idx = 1:1
% T = 2; alphas = [-0.92, -0.5, -0.23, 0, 0.5, 1, 2.1]; rests = nan(length(ns),length(alphas),4); for idx = 1:length(alphas), alpha=alphas(idx);
alpha = 0; Ts = 2:2:8; rests = nan(length(ns),length(Ts),4); for idx=1:length(Ts), T=Ts(idx);
  P = exactLaguerre(alpha,[0 1],max(ns) );
  for ni = 2:length(ns)
    n = ns(ni);
    d = 1/(4*n+2*alpha+2);

    pt = (4*n -4*transpose(1:n) +3)*d;
    t = pi^2/16*(pt -1).^2;

    for it = 1:6
        t = t - (pt*pi +2*sqrt(t-t.^2) -acos(2*t-1) ).*sqrt(t./(1-t))/2;
    end

    jak = besselroots(alpha, n);
    bes = jak*0;
    wbes = 0*bes;
    bulk = t*0;
    wbulk = t*0;

    ak = [-13.69148903521072; -12.828776752865757; -11.93601556323626;    -11.00852430373326; -10.04017434155809; ...
        -9.02265085340981; -7.944133587120853;    -6.786708090071759; -5.520559828095551; -4.08794944413097; -2.338107410459767];
    tair = 3*pi/2*( transpose(n:-1:12)-0.25);
    ak = [-tair.^(2/3).*(1 + 5/48./tair.^2 - 5/36./tair.^4 + 77125/82944./tair.^6 -10856875/6967296./tair.^8); ak];

    if (T >= 7)
        bes = bes + (657*jak.^6 +36*jak.^4*(73*alpha^2-181) +2*jak.^2*(2459*alpha^4 -10750*alpha^2 +14051) +...
            4*(1493*alpha^6 -9303*alpha^4 +19887*alpha^2 - 12077) )*d^6/2835;
        wbes = wbes + (11944*alpha^6 + 5256*jak.^6 - (5061*alpha^5 + 5085*alpha^4 + 4830*alpha^3 -22724*alpha^2 ...
            - 22932*alpha + 39164)*jak.^4 - 74424*alpha^4 + 8*(2459*alpha^4 -10750*alpha^2 + 14051)*jak.^2 + ...
            159096*alpha^2 - 96616)/2835/2*d^6;
        
        bulk = bulk -d^5/181440*(9216*(21*alpha^6 - 105*alpha^4 + 147*alpha^2 - 31)*t.^10 ...
            -69120*(21*alpha^6 - 105*alpha^4 + 147*alpha^2 - 31)*t.^9 + 384*(12285*alpha^6 -61320*alpha^4 + 85785*alpha^2 - 18086)*t.^8 ...
            - 64*(136080*alpha^6 - 675675*alpha^4 +943110*alpha^2 - 198743)*t.^7 +...
            144*(70560*alpha^6 - 345765*alpha^4 +479850*alpha^2 - 101293)*t.^6 + 72576*alpha^6 ...
            - (8128512*alpha^6 - 38656800*alpha^4+ 52928064*alpha^2 - 13067711)*t.^5 + ...
            5*(1016064*alpha^6 - 4581360*alpha^4 +6114528*alpha^2 + 113401)*t.^4 ...
            - 317520*alpha^4 - 10*(290304*alpha^6 -1245888*alpha^4 + 1620864*alpha^2 - 528065)*t.^3 + ...
            5*(290304*alpha^6 -1234800*alpha^4 + 1598688*alpha^2 - 327031)*t.^2 + 417312*alpha^2 ...
            -5*(96768*alpha^6 - 417312*alpha^4 + 544320*alpha^2 - 111509)*t -85616)./(t-1).^8./t.^2;
        wbulk = wbulk + d^6/362880*(9216*(21*alpha^6 - 105*alpha^4 + 147*alpha^2 - 31)*t.^10 ...
            -1536*(945*alpha^6 - 4830*alpha^4 + 6825*alpha^2 - 1444)*t.^9 + ...
            384*(11340*alpha^6 -60165*alpha^4 + 86310*alpha^2 - 18289)*t.^8 ...
            - 2*(2903040*alpha^6 - 17055360*alpha^4+ 25401600*alpha^2 - 5*alpha - 5252997)*t.^7 ...
            - (11753280*alpha^4 - 23506560*alpha^2+ 67*alpha - 13987519)*t.^6 - 290304*alpha^6 + ...
            12*(1016064*alpha^6 -3578400*alpha^4 + 4108608*alpha^2 + 16*alpha + 7100871)*t.^5 ...
            - 5*(4064256*alpha^6 -16559424*alpha^4 + 20926080*alpha^2 + 61*alpha - 15239393)*t.^4 +...
            1270080*alpha^4 +10*(1741824*alpha^6 - 7386624*alpha^4 + 9547776*alpha^2 + 29*alpha - 1560107)*t.^3 ...
            - 15*(580608*alpha^6 - 2503872*alpha^4 + 3265920*alpha^2 + 11*alpha - 669051)*t.^2 ...
            - 1669248*alpha^2 + 4*(604800*alpha^6 - 2630880*alpha^4 + 3447360*alpha^2 + 13*alpha- 706850)*t ...
            - 7*alpha + 342463)./(t-1).^9./t.^3;
    end
    if (T >= 5)
        bes = bes + (11*jak.^4 +3*jak.^2.*(11*alpha^2-19) +46*alpha^4 -140*alpha^2 +94)*d^4/45;
        wbes = wbes + (46*alpha^4 + 33*jak.^4 +6*jak.^2*(11*alpha^2 -19) -140*alpha^2 +94)/45*d^4;
        
        bulk = bulk - d^3/720*(32*(15*alpha^4 - 30*alpha^2 + 7)*t.^6 -144*(15*alpha^4 - 30*alpha^2 + 7)*t.^5 +...
            16*(225*alpha^4 - 450*alpha^2 +104)*t.^4 - 240*alpha^4 - 480*(5*alpha^4 - 10*alpha^2 + 1)*t.^3 + 480*alpha^2 ...
            +45*(16*alpha^4 - 32*alpha^2 + 7)*t + 990*t.^2 - 105)./(t-1).^5./t;
        wbulk = wbulk + d^4/720*(16*(15*alpha^4 - 30*alpha^2 + 7)*t.^6 - 32*(45*alpha^4 - 90*alpha^2 +22)*t.^5 + ...
            48*(75*alpha^4 - 150*alpha^2 + 74)*t.^4 + 240*alpha^4 - 600*(8*alpha^4- 16*alpha^2 - 5)*t.^3 + ...
            45*(80*alpha^4 - 160*alpha^2 + 57)*t.^2 - 480*alpha^2 -90*(16*alpha^4 - 32*alpha^2 + 7)*t + 105)./(t-1).^6./t.^2;
    end
    if (T >= 3)
        bes = bes + (jak.^2 + 2*alpha^2 - 2)*d^2/3;
        wbes = wbes + (alpha^2 + jak.^2 -1)*2/3*d^2;
        
        bulk = bulk - d/12*(4*(3*alpha^2 - 1)*t.^2 +12*alpha^2 - 12*(2*alpha^2 - 1)*t - 3)./(t-1).^2;
        wbulk = wbulk  + d^2/6*(2*t + 3)./(t-1).^3;
    end
    bes = jak.^2*d.*(1 + bes);
    air = 1/d +ak*(d/4)^(-1/3) +ak.^2*(d*16)^(1/3)/5 + (11/35-alpha^2-12/175*ak.^3)*d + ...
        (16/1575*ak+92/7875*ak.^4)*2^(2/3)*d^(5/3) -(15152/3031875*ak.^5+1088/121275*ak.^2)*2^(1/3)*d^(7/3);
    bulk = bulk + t/d;

    bulkOld = bulk;
    for k = 1:n
        bes(k) = abs(P(bes(k), n+1)); %Absolute value of the orthonormal poly should be small at the zero
        bulk(k) = abs(P(bulk(k), n+1));
        air(k) = abs(P(air(k), n+1));
    end
    if any(isnan(bes)) || any(isnan(bes)) || any(isnan(bes))
        error('Overflow');
    end
    rests(ni, idx, 1) = find(bulk < bes, 1, 'first');
    rests(ni, idx, 2) = find(bulk < air, 1, 'last');
    rests(ni, idx, 3) = bulkOld(rests(ni,idx,1));
    rests(ni, idx, 4) = bulkOld(rests(ni,idx,2));
  end
  idx
end
% Plot heuristics for nodes
for typ =1:4
    figure;
    loglog(ns, rests(:, 1, typ), 'k+');
    hold on;
    loglog(ns, rests(:, 2, typ), 'bo');
    loglog(ns, rests(:, 3, typ), 'gx');
    loglog(ns, rests(:, 4, typ), 'm.');
    if alpha == 0 % Then iteration was over number of terms
        legs = cell(5,1);
        for ai = 1:length(Ts)
            legs{ai} = ['$T$ = ' num2str(Ts(ai))];
        end
    else
        loglog(ns, rests(:, 5, typ), 'yd');
        loglog(ns, rests(:, 6, typ), 'cs');
        loglog(ns, rests(:, 7, typ), 'rp');
        
        legs = cell(8,1);
        for ai = 1:length(alphas)
            legs{ai} = ['$\alpha$ = ' num2str(alphas(ai))];
        end
    end
    
    if typ == 1
        loglog(ns, sqrt(ns), 'k'); legs{end} = '$\sqrt{n}$'; heu = 'low $k$';
        yticks([5,7,10,15]);
    elseif typ == 2
        loglog(ns, 0.9*ns, 'k'); legs{end} = '$0.9n$'; heu = 'high $k$';
        yticks([30, 50, 70, 100, 200, 300]);
    elseif typ == 3
        loglog(ns, pi^2/4*ones(size(ns)), 'k'); legs{end} = '$\pi^2/4$'; heu = 'low $x_k$';
        yticks([0.3, 0.6, 1, 2, 3]);
    elseif typ == 4
        loglog(ns, 3.6*ns, 'k'); legs{end} = '$3.6n$'; heu = 'high $x_k$';
    end
    legend(legs, 'interpreter', 'latex', 'location', 'best');
    xlabel('n');
    ylabel(['Optimal transition for ' heu], 'interpreter', 'latex');
    axis([0.7*min(ns), 1.1*max(ns), 0.9*min(min(rests(:,:,typ))), max(2.8, 1.1*max(max(rests(:,:,typ)) ))]);
    set(gca, 'FontSize', 20);
end

%return




%% Heuristics for orth polys themselves

alpha = 0;
q = [0 1];
n = 200;
maxOrder = 5;

[P,gammaP,aP,bP] = exactLaguerre(alpha,q,n);
s = getAsy(alpha,q,maxOrder);


%% Plot errors to see the most accurate expansion of the orthogonal polynomial (not the zero) 

zs = linspace(0,1,1000); % for z in [0,1] but these can be shifted in the complex plane
ers = zeros(length(zs), 6, maxOrder); % Errors with respect to the exact polynomial: 1 is lens (b), 2 is outer (a), 3 is right with Q's (s.pc(...,1)), 4 is left usq (d), 5 right no Q, 6 left s.pd(...,0)
[P,gammaP,aP,bP] = exactLaguerre(alpha,q,n);
ord = maxOrder;

for zi = 1:length(zs)
    pEx = P(zs(zi)*s.betan(n,maxOrder),n +1);
    for i = ord:ord % Could also plot multiple orders on same plot
        ers(zi,1,i) = s.pb(n,zs(zi)*s.betan(n,maxOrder),i);
        ers(zi,2,i) = s.pa(n,zs(zi)*s.betan(n,maxOrder),i);
        ers(zi,3,i) = s.pc(n,zs(zi)*s.betan(n,maxOrder),i,1);
        ers(zi,4,i) = s.pd(n,zs(zi)*s.betan(n,maxOrder),i,1);
        ers(zi,5,i) = s.pc(n,zs(zi)*s.betan(n,maxOrder),i,0);
        ers(zi,6,i) = s.pd(n,zs(zi)*s.betan(n,maxOrder),i,0);
% Could also compute the error on the recurrence relation itself as:
%         ersRec(zi,1,i) = abs( ( (zs(zi)*s.betan(n,maxOrder)-(2*n-1+alpha))*real(s.pb(n-1,zs(zi)*s.betan(n,maxOrder),i)) - ...
%             sqrt((n-1)*(n-1+alpha))*real(s.pb(n-2,zs(zi)*s.betan(n,maxOrder),i)) )/...
%             sqrt((n)*(n+alpha))/real(s.pb(n,zs(zi)*s.betan(n,maxOrder),i)) -1);
    end
    ers(zi,:,:) = abs((pEx-ers(zi,:,:))/pEx);
end

figure; 
leg = cell(6,maxOrder);
leg{1,ord} = 'lens'; semilogy(zs*s.betan(n,maxOrder), ers(:,1,ord), 'k.');
hold on;
leg{2,ord} = 'outer'; semilogy(zs*s.betan(n,maxOrder), ers(:,2,ord), 'c-.');
leg{3,ord} = 'right Q'; semilogy(zs*s.betan(n,maxOrder), ers(:,3,ord), 'g:');
leg{4,ord} = 'left Q'; semilogy(zs*s.betan(n,maxOrder), ers(:,4,ord), 'm:');
leg{5,ord} = 'right U'; semilogy(zs*s.betan(n,maxOrder), ers(:,5,ord), 'b--');
leg{6,ord} = 'left U'; semilogy(zs*s.betan(n,maxOrder), ers(:,6,ord), 'r--');

legend(leg(:,maxOrder));
xlabel('x');
ylabel('Relative error'); 
set(gca,'FontSize', 20);

%return



%% Accuracy of asymptotic expansions of zeros of Laguerre polynomials in the bulk, possible Chebfun example on Newton iterations on functions
% In the bulk, t is given by the solution of 
%
% $$ 2 \cos^{-1}(\sqrt{t}) -2 \sqrt{t-t^2} = \frac{4 n -4 k+3}{4 n +2}\pi. $$
%
% We set up a Newton iteration with an initial guess:

splitting on
p = chebfun('p',[0.01,0.99]);
N = @(t) 2*acos(sqrt(t)) -2*sqrt(t -t.^2) -p*pi;
dN = @(t) -2*sqrt(1./t-1);
g = cos(pi/4*(p+1)).^2;
% We take 5 Newton iterations
resids = norm(N(g));
for k = 1:5
  g = g - N(g)./dN(g);
  resids = [resids, norm(N(g))];
end
% The resulting chebfun consists of two smooth pieces of modest length. Consecutive
% iterations indeed converge:
resids
% We compute some reference nodes
n = 100;
k = (5:97)';
X = lagptsTest(n);
xnk = X(k);
% Compute the expansion of the nodes using our chebfun
% and compare the expansion with them
t = g((4*n -4*k +3)./(4*n +2));
x1 = (4*n +2)*t;
x2 = x1 -1/6/(2*n +1)*(5/4./(1-t).^2 -1./(1-t) -1);
x3 = x2  - 1/(4*n +2)^3/720*(224*t.^6 -1008*t.^5 + 1664*t.^4  - 480*t.^3 ...
    +315*t + 990*t.^2 - 105)./(t-1).^5./t;
x4 = x3  -1/(4*n +2)^5/181440*(-285696*t.^10 +2142720*t.^9 + -6945024*t.^8 ...
    +12719552*t.^7 -14586192*t.^6  +13067711*t.^5 + 567005*t.^4  +5280650*t.^3 +...
    -1635155*t.^2  +557545*t -85616)./(t-1).^8./t.^2;
figure; semilogy(xnk, [abs(x1-xnk)./xnk, abs(x2-xnk)./xnk, ...
    abs(x3-xnk)./xnk, abs(x4-xnk)./xnk], '*')
xlabel('x_k'); ylabel('Relative error'); 
legend('1 term', '2 terms', '3 terms', '4 terms', 'Location', 'best');

%return



%% Test Newton iterations on asymptotic expansions
alpha = 0.7;
n = 269;

[xr, wr] = lagptsTest(n, [0 inf], 'rec', alpha);
[xrh, wrh] = lagptsTest(n, [0 inf], 'RH', alpha);

[xbes, wbes] = lagptsTest(n, [0 inf], 'RHbes', alpha);
[xbulk, wbulk] = lagptsTest(n, [0 inf], 'RHbulk', alpha);
[xair, wair] = lagptsTest(n, [0 inf], 'RHairy', alpha);
% [xexpl, wexpl] = lagptsTest(n, [0 inf], 'expl', alpha);

figure; semilogy(abs(wrh-wr)./wr ); hold on; semilogy(abs(xrh-xr)./xr );

figure;
semilogy(abs(xbes-xr)./xr, 'b-');
hold on;
semilogy(abs(wbes-wr)./wr, 'bo');
semilogy(abs(xbulk-xr)./xr, 'r-.');
semilogy(abs(wbulk-wr)./wr, 'r.');
semilogy(abs(xair-xr)./xr, 'k--');
semilogy(abs(wair-wr)./wr, 'kx');
xlabel('k');
ylabel('Relative error');
legend({'Left x_k', 'Left w_k', 'Lens x_k', 'Lens w_k', 'Right x_k', 'Right w_k'}, 'Location', 'best')
set(gca, 'FontSize', 20);


