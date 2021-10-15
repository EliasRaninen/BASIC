%% Reproduces Fig 2. from the paper "Bias adjusted sign covariance matrix",
% Elias Raninen and Esa Ollila, 2021.

clear; clc; close all;
rng('default')

%% Define simulation parameters
nmc = 2000; % number of Monte Carlo trials
p   = 100;  % dimension
n   = 50;   % number of samples

%% covariance matrix
rho = 0.5;

M = toeplitz(rho.^(0:p-1)); setup = 'AR1'; % AR(1)
% M = rho*ones(p) + (1-rho)*eye(p); setup = 'CS'; % CS

Msq  = sqrtm(M);

%% compute theoretical population SSCM
ts          = @(A) p*A/trace(A); % normalize to shape
I           = eye(p);
fcomplexmat = @(t) (ts(M)/((1-t)*I+t*ts(M))) / det((1-t)*I + t*ts(M)) * (1-t)^(p-1);

% population SSCM
ESSCM = integral(@(t) fcomplexmat(t), 0, 1, 'ArrayValued',true);
% population SSCM shape matrix
Vsgn = p*ESSCM;

% theoretical NMSE
gamsgn      = trace(Vsgn^2)/p;
a           = p/n + (n-1)/n * gamsgn;
MSEtheory   = @(al) al^2 * p * (a - 1) + (1-2*al) * p * (gamsgn - 1);
NMSEtheory  = @(al) MSEtheory(al)/(gamsgn*p);

% distance to population SSCM
distSSCM    = @(EST) norm(ts(EST)-ts(Vsgn),'F')^2 / norm(ts(Vsgn),'F')^2;

%% Main loop
points      = 21; % number of points to compute alpha
al_arr      = linspace(0,1,points); % array of alpha values

al_iter     = nan(nmc,points); % estimated alpha
NMSERSSCM   = nan(nmc,points);

err         = zeros(points,1); % for NMSE
NMSEth      = nan(points,1); % NMSE theory

for ii=1:points % for each alpha
    for mc=1:nmc % average over nmc Monte Carlo runs
        
        %% Generate data
        X0 = 1/sqrt(2)*complex(randn(n,p),randn(n,p));
        X = X0*conj(Msq);
        
        %% compute SSCM
        xnorm = sqrt(sum(abs(X).^2,2));
        U     = X./xnorm;
        SSCM  = (U.'*conj(U))/n;
        V     = p*SSCM;
        
        %% compute shrinkage parameter
        a       = trace(V^2)/p;
        num     = n/(n-1)*(a - p/n) - 1;
        den     = a - 1;
        alhat   = num/den;
        alhat   = min(1,max(0,alhat));
        al_iter(mc,ii) = alhat;
        
        % empirical NMSE error of RSSCM with estimated alpha
        NMSERSSCM(mc,ii) = distSSCM(alhat*V + (1-alhat)*eye(p));
        
        %% empirical NMSE of RSSCM for alpha in [0,1]
        al = al_arr(ii);
        RSSCM = al*V + (1-al)*eye(p);
        
        % empirical NMSE between RSCCM and Vsgn
        err(ii) = err(ii) + distSSCM(RSSCM)/nmc;
    end
    % theoretical NMSE of RSSCM for alpha in [0,1]
    NMSEth(ii) = NMSEtheory(al);
    fprintf('.')
end

%% mean of estimated regularization parameter alpha and RSSCM
almean      = mean(al_iter(:));     % mean of estimated alphas
alstd       = std(al_iter(:));      % standard error of estimated alphas
mNMSERSSCM  = mean(NMSERSSCM(:));   % NMSE of estimated RSSCM
sNMSERSSCM  = std(NMSERSSCM(:));    % standard error of estimated RSSCM

%% Plot

% NMSE for alpha in [0,1]
figure(1); clf; hold on;
plot(al_arr, err,'k','linewidth',1,'displayname','Empirical NMSE');
plot(al_arr, NMSEth,'b--','linewidth',1,'displayname','Theoretical NMSE');

% estimated alpha and RSSCM
plot([almean-alstd, almean+alstd],[mNMSERSSCM mNMSERSSCM], 'r', 'linewidth', 2, 'displayname', 'estimated \alpha +- std')
plot([almean almean], [mNMSERSSCM-sNMSERSSCM mNMSERSSCM+sNMSERSSCM], 'm', 'linewidth', 2, 'displayname', 'NMSE of \Lambda_{RSSCM}(estimated \alpha) +- std')
legend

if isequal(setup,'CS')
    axis([0.5 1 0.05 0.28]) % CS
elseif isequal(setup,'AR1')
    axis([0 0.5 0.285 0.4]); % AR
end

