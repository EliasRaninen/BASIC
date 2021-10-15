function [RSSCM, al] = REGSSCM(X)
%% COMPUTES the regularized spatial sign covariance matrix (RSSCM)
% The function assumes that the data is centered and does not contain
% zero-valued samples or nan's.
%
% Usage:    [RSSCM,al] = REGSSCM(X)
%
% Inputs:   X       -   nxp data matrix of n samples and p dimensions
%
% Outputs:  RSSCM   -   regularized SSCM
%           al      -   estimated regularization parameter
%
% By Elias Raninen 2021
%
% version 1.0 (Sep. 29, 2021)

% n samples, p dimensions
[n,p] = size(X);

% compute SSCM
xnorm = sqrt(sum(abs(X).^2,2));
U     = X./xnorm;
SSCM  = (U.'*conj(U))/n;
V     = p*SSCM;

% compute shrinkage parameter
a   = trace(V^2)/p;
num = n/(n-1)*(a - p/n) - 1;
den = a - 1;
al  = num/den;
al  = min(1,max(0,al));

%% RSSCM (no bias correction)
RSSCM = al*V + (1-al)*eye(p);