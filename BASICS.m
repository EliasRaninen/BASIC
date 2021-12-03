function [BEST, al] = BASICS(X,lambdas,deltas)
%% COMPUTES the Bias Adjusted SIgn Covariance Shrinkage estimator (BASICS),
% which is a bias corrected regularized spatial sign covariance matrix. The
% function assumes that the data is centered (zero-mean) and does not contain
% zero-valued samples or nan's.
%
% Usage:    [BEST, al] = BASICS(X,lambdas,deltas).
%
% Inputs:   X           -   nxp data matrix of n samples and p dimensions
%           lambdas     -   need to be precomputed with BASICtable.m.
%           deltas      -   need to be precomputed with BASICtable.m.
%
% Outputs:  BEST        -   BASICS estimate
%           al          -   estimated regularization parameter
%
%
% Note: lambdas (eigenvalues of shape) and the corresponding deltas
%       (eigenvalues of SSCM) are obtained with BASICtable.m, e.g.,
%       [lambdas,deltas] = BASICtable(p,field,points)
%
% By Elias Raninen 2021
%
% version 1.01 (Dec. 2, 2021)


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

%% BASICS (bias correction)

% compute eigenvectors and eigenvalues
[E, deltahat] = svd(RSSCM/p);
deltahat = diag(deltahat);

% use supplied table to correct the eigenvalues
lambdahat = interp1(deltas,lambdas,deltahat,'linear');
if any(isnan(lambdahat)) % extrapolate in case deltahat value is out of range of given deltas
    fprintf('BASICS.m: eigenvalues of SSCM out of range of given table. Using extrapolation via splines.\n');
    lambdahat = interp1(deltas,lambdas,deltahat,'spline');
end

% normalize to shape
lambdaBASICS = lambdahat/sum(lambdahat)*p;

% BASICS estimate
BEST = E*diag(lambdaBASICS)*E';
