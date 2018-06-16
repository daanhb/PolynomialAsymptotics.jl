% Compute a Gauss-Laguerre quadrature rule from RH-asymptotics.
% About
%   Author       - Peter Opsomer (Peter.Opsomer@cs.kuleuven.be)
%   History      - Created June 2016
%% Initialising and choosing parameters
format longe; close all; clear variables;
set(0,'DefaultFigureWindowStyle','docked');

alpha = 0;
% q = [0 1]; % Standard associated Laguerre polynomials
% q = [8 3 0.9];
q = [0 0 0 1.3];

maxOrder = 4; % Maximum negative order of n calculated in the expansions
n = 1e2;
maxAt = 20; % Maximum number of Newton iterations.
usq = 1;
nstop = n; % n from which the weights will probably have underflowed

stdLag = isnumeric(q) && (length(q) == 2) && (norm(q-[0 1]) < eps);
Qx = @(x) q*repmat(x, length(q), 1).^repmat((0:length(q)-1)', 1, length(x));
if ~isnumeric(q)
    Qx = @(x) q(x);
end

    
%% Computing functions
s = getAsy(alpha,q,maxOrder);
[P,gammaP,aP,bP] = exactLaguerre(alpha,q,n);
if stdLag
    sd = getAsy(alpha+1,q,maxOrder);
    Pd = exactLaguerre(alpha+1,q,n);
end


%% Loop
zw = zeros(nstop,2);
atts = maxAt*ones(nstop,1);

ercnt = 0;
printtoc = 2;
tic;
prevToc = toc;
for ni = 1:nstop
    if(toc - prevToc > printtoc)
        prevToc = toc;
        display([num2str(ni) ' = ni, est sec left = ' num2str(prevToc*(n-ni+1)/(ni-1))]);
    end
    % Lower and Upper bounds:
    lb = nan; % Inequalities will return false if not stdLag.
    ube = nan;
    if stdLag
        lb = (fzero(@(x) besselj(alpha,x), pi*(ni+alpha/2-1/4)-(4*alpha^2-1)/8/pi/(ni+alpha/2-1/4)+[-0.5,0.5]) )^2 ...
            /(4*n+2*alpha+2); % DLMF (18.16.10) with initial guess (10.21.19) and dist between zeros is about pi > 2 -> besselzeros in chebfun
        ube = (4*ni+2*alpha+2)*(2*ni+alpha+1+sqrt((2*ni+alpha+1)^2+1/4-alpha^2))/(4*n+2*alpha+2); % DLMF (18.16.11)
    end
    if stdLag && (ni <= 3)
        zw(ni,1) = lb; % DLMF (18.16.10)
        if ni == 1
            lbfz = lb/2;
        else
            lbfz = (lb+zw(ni-1,1))/2;
        end
        ubfz = ((fzero(@(x) besselj(alpha,x), pi*(ni+alpha/2+3/4)-(4*alpha^2-1)/8/pi/(ni+alpha/2+3/4)+[-0.5,0.5]) )^2 ...
        /(4*n+2*alpha+2)+lb)/2;
    elseif ni <= 3
        zw(ni,1) = (fzero(@(x) besselj(alpha,x), pi*(ni+alpha/2-1/4)-(4*alpha^2-1)/8/pi/(ni+alpha/2-1/4)+[-0.5,0.5]) )^2 ...
            *(2*s.m-1)^2/16/s.m^2/n^2*s.betan(n,maxOrder);
        lbfz = zw(ni,1)/2;
        ubfz = ((fzero(@(x) besselj(alpha,x), pi*(ni+alpha/+3/4)-(4*alpha^2-1)/8/pi/(ni+alpha/2+3/4)+[-pi/2,pi/2]) )^2 ...
            *(2*s.m-1)^2/16/s.m^2/n^2*s.betan(n,maxOrder)+zw(ni,1))/2;
    else
        lbfz = 15*zw(ni-1,1)/8 -5*zw(ni-2,1)/4 + 3*zw(ni-3,1)/8;
        ubfz = 35*zw(ni-1,1)/8 -21*zw(ni-2,1)/4 +15*zw(ni-3,1)/8;
        zw(ni,1) = 3*zw(ni-1,1) -3*zw(ni-2,1) +zw(ni-3,1); % Quadratic extrapolation
%         zw(ni,1) = lb; % This bound is said to have the best possible constant for standard Laguerre.
    end
    if ( zw(ni,1) < 4*n*(1-max(1e-8, 1.8/sqrt(n) ) ) ) % Heuristic    
        po = @(nn,x) real(s.pdwoq(nn,x,maxOrder,usq)); 
        % [DLMF] d L^alpha_n/dx = -L^{alpha+1}_{n-1} with norm L^alpha_n = gamma(n+alpha+1)/n! and leading order coeff (-1)^n/n!
        % so d p^alpha_n/dx = + p^{alpha+1}_(n-1) sqrt(n) and take e^(Q(x)/2) out
        dpo = @(nn,x) real(sd.pdwoq(nn-1,x,maxOrder,usq))*sqrt(n); % dpo will not be called if not stdLag
    else
        po = @(nn,x) real(s.pcwoq(nn,x,maxOrder,usq));
        dpo = @(nn,x) real(sd.pcwoq(nn,x,maxOrder,usq))*sqrt(n+alpha+1); % Using sd.pcwoq(nn-1... would risk x > s.betan(n-1)
    end
    % zw(ni,1) = fzero(@(x) po(n,x), [lbfz, ubfz]); % Instead of Newton procedure below
    ov = inf;
    for attempt = 1:maxAt
        pe = po(n,zw(ni,1) );
        if stdLag
            step = pe/(dpo(n,zw(ni,1)) - po(n,zw(ni,1))/2); % (s.pc*exp(-Q/2) )' = exp(-Q/2)*(pc' -pc/2)
        else
            hh = 1e-7;
            step = pe/(po(n,zw(ni,1)*(1+hh)) -po(n,zw(ni,1)))*(hh*zw(ni,1));
        end
        zw(ni,1) = zw(ni,1) - step;
        if (zw(ni,1) < lbfz) || (zw(ni,1) > ubfz) || isnan(zw(ni,1)) || (abs(pe) > abs(ov)+1e-12)
            error('Converging outside interval or away from zero')
        elseif abs(step) < 400*eps*zw(ni,1) % Converged
            % Lowering this bound would only be meaningful if the expansions are made to reach machine precision
            atts(ni) = attempt;
            break;
        elseif attempt == maxAt
            error('No convergence')
        end
        ov = pe;
    end
    % Correctness checks for the nodes
    if (ni > 1) && (zw(ni,1) <= zw(ni-1,1))
        error('Nodes should be increasing');
    elseif zw(ni,1) < lb 
        warning('Lower bound exceeded');
        % Count number of times the lower bound is not satisfied because the asymptotics may be more 
        % inaccurate than the lower bound, for example when the number of terms is too low.
        ercnt = ercnt + 1; 
    elseif zw(ni,1) >= ube
        error('Upper bound exceeded');
    end
    if stdLag
        zw(ni,2) = -s.ratGam(n+1,maxOrder)/dpo(n,zw(ni,1))/po(n+1,zw(ni,1))/exp(zw(ni,1)); % e^(Q/2) not included in both polys and zw(ni,1) < s.betan(n+1,T)
    else
        hh = max(sqrt(eps)*zw(ni,1), sqrt(eps) );
        % This leaves out a constant factor, given by a ratio of leading order coefficients and normalising constants
        zw(ni,2) = hh/(po(n,zw(ni,1)+hh) - po(n,zw(ni,1) ) )/po(n-1,zw(ni,1) )/exp(Qx(zw(ni,1)));
    end
    
    if zw(ni,2) < 0
        error('Weights should be positive')
    elseif (zw(ni,2) == 0) && (nstop == ni) % Can neglect weights from now on
        display(['Weights are below realmin from n = ' num2str(nstop)]);
        break
    end
end
% We left out a constant factor while computing w in the non standard Laguerre case
if isnumeric(q) && ~stdLag
    zw(:,2) = zw(:,2)/sum(zw(:,2))*gamma((alpha+1)/s.m)*q(end)^(-(alpha+1)/s.m)/s.m;
    % This is still wrong if Qx(x) is not a monomial.
end
toc

%% Plotting
if isnumeric(q) && (norm(q(1:end-1)) == 0) 
    % Weights are only correct for monomial Qx(x) with an error about sqrt(eps) for nonstandard Laguerre.
    errSumW = sum(zw(:,2))-quadgk(@(x) x.^alpha.*exp(-Qx(x)), 0, inf)
    errExSumW = quadgk(@(x) x.^alpha.*exp(-Qx(x)), 0, inf) - gamma((alpha+1)/s.m)*q(end)^(-(alpha+1)/s.m)/s.m
    fct = @(x) cos(x)+log(x+3)+x.^8;
    quadFct = sum(zw(:,2).*fct(zw(:,1)))
    errQuadFct = sum(zw(:,2).*fct(zw(:,1)))-quadgk(@(x) fct(x).*x.^alpha.*exp(-Qx(x) ), 0, inf)
end
figure; plot(zw(:,1)); title('Nodes');
figure; plot(atts); title('Number of Newton iterations for each node')
figure; semilogy(zw(:,1), zw(:,2) ); title('Weights');


%% Testing asymptotics of the first zero for standard Laguerre
if stdLag
    fzbes = fzero(@(x) besselj(alpha,x), alpha+1);
    frAsy = (fzbes)^2/(4*n+2*alpha+2)
    firstRoot = fzero(@(x) real(s.pd(n,x,maxOrder, usq)), frAsy)
    dlmfRoot = fsolve(@(x) sqrt(x-x^2)+asin(sqrt(x))-2*fzbes/(4*n+2*alpha+2), sqrt(eps))*(4*n+2*alpha+2)
    % The latter only takes one term, while frAsy tends to the first node if n goes to infinity
end