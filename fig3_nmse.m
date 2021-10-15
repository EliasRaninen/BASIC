%% Reproduces Fig 3. from the paper "Bias adjusted sign covariance matrix",
% Elias Raninen and Esa Ollila, 2021.

clear; clc; close all;
rng('default')

%% Define number of Monte Carlo trials
nmc = 2000;

%% Covariance models and filenames for saving data
p = 100;
rho = 0.5;

% M = toeplitz(rho.^(0:p-1));
M = rho*ones(p) + (1-rho)*eye(p);

%% sample sizes to be simulated
n_arr = ((p/2):5:2*p).'; 

%% Distance function to be used in simulation
ts      = @(A) p*A/trace(A); % normalized to shape
dist    = @(EST) norm(ts(EST)-ts(M),'F')^2/norm(ts(M),'F')^2; % NMSE

%% normalize covariance matrix to shape and compute eigenvalues
M       = ts(M);
Msq     = sqrtm(M);
eigM    = svd(M);
I       = eye(p);

%% read table for BASICS method
[lambdas,deltas] = BASICtable(p,'complex');

%% indices for methods
idxBASICS     = 1;
idxRSSCM    = 2;
idxRFP      = 3;
nmethods    = 3;
%% Main loop
err = nan(nmc,numel(n_arr),nmethods);      % for error of methods 
comptime = nan(nmc,numel(n_arr),nmethods); % for computation time
fprintf('                                                                 |\n');
for jj=1:numel(n_arr)
    n = n_arr(jj)
    for mc=1:nmc
        if mod(mc,30)==0; fprintf('.'); end
        
        % generate complex t-distributed data with scatter matrix M
        X0 = (1/sqrt(2))*complex(randn(n,p),randn(n,p));
        v  = 2;
        s  = chi2rnd(v,n,1);
        X  = (X0 ./ sqrt(s/v)) * conj(Msq);
        
        % BASICS estimator
        tic;
        BEST = BASICS(X,lambdas,deltas);
        comptime(mc,jj,idxBASICS) = toc;
        
        % REGSSCM
        tic;
        RSSCM = REGSSCM(X);
        comptime(mc,jj,idxRSSCM) = toc;
                           
        % REGFP
        tic;
        RFP = REGFP(X);
        comptime(mc,jj,idxRFP) = toc;
        
        % compute and save errors (NMSE)
        err(mc,jj,idxBASICS)  = dist(BEST);
        err(mc,jj,idxRSSCM) = dist(RSSCM);
        err(mc,jj,idxRFP)    = dist(RFP);

    end
    fprintf('\n')
end

%% Plot NMSE

% compute mean over MC trials
meanerr = squeeze(mean(err,1));

% plot error
figure(1); clf; hold on;
n = n_arr;
plot(n,meanerr(:,idxBASICS),'b','linewidth',2,'displayname','BASICS');
plot(n,meanerr(:,idxRSSCM),'displayname','RSSCM');
plot(n,meanerr(:,idxRFP),'displayname','RFP');

xlabel('n')
ylabel('NMSE')
legend

%% plot computation time
CT = squeeze(median(comptime,1));

figure(2); clf; hold on;
plot(n,CT(:,idxBASICS),'displayname','BASICS');
plot(n,CT(:,idxRSSCM),'displayname','RSSCM');
plot(n,CT(:,idxRFP),'displayname','RFP');

legend
xlabel('n')
ylabel('median computation time');
ylabel('time')
