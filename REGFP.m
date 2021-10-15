function [RFP, al] = REGFP(X, maxiter, TOL)
%% COMPUTES the FP bias corrected regularized spatial sign covariance matrix (RFP).
% The function assumes that the data is centered and does not contain
% zero-valued samples or nan's.
%
% The fixed point algorithm was proposed for real-valued elliptical data in:
% A. Durre, R. Fried, and D. Vogel, "The Spatial Sign Covariance Matrix
% and Its Application for Robust Correlation Estimation," Austrian Journal
% of Statistics, vol. 46, no. 3-4, pp. 13–22, Apr. 2017.
%
% Usage:    [RFP, al] = REGFP(X)
%           [RFP, al] = REGFP(X, maxiter, TOL)      
%
% Inputs:   X           -   nxp data matrix of n samples and p dimensions
%           maxiter     -   maximum number of iterations (100 by default)
%           TOL         -   tolerance for stopping criterion:
%                           norm(lambda-lambda_old)/norm(lambda_old) < TOL.
% Outputs:  RFP         -   regularized SSCM with FP bias correction
%           al          -   estimated regularization parameter
%
% By Elias Raninen 2021
%
% version 1.0 (Sep. 29, 2021)

if ~exist('maxiter','var'); maxiter = 100; end
if ~exist('TOL','var'); TOL = 1e-5; end

if isreal(X)
    be = 1/2; % for real-valued data
else
    be = 1;   % for complex-valued data
end

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

%% RFP (bias correction using a fixed point algorithm)

% compute eigenvectors and eigenvalues
[E, delta] = svd(RSSCM/p);
delta = diag(delta);
lambda = p*delta;

% f = @(x,lambda) (1+lambda.*x).^(-1) .* prod(1+lambda.*x).^(-be); % alternative form
f = @(x,lambda) (((1-x)+x.*lambda).^(-1) .* prod((1-x)+x.*lambda)^(-be))  .* (1-x).^(be*p-1);

for iter = 1:maxiter

   lambda_old = lambda;

%    integ = integral(@(x) f(x,lambda),0,inf,'ArrayValued',true); % alternative form
   integ = integral(@(x) f(x,lambda),0,1,'ArrayValued',true);

   lambda = be^(-1) * delta .* integ.^(-1);
   lambda = p*lambda/sum(lambda);
   
   crit = norm(lambda-lambda_old)/norm(lambda_old);
   if  crit < TOL
       break
   end
end
if iter == maxiter
    fprintf('RFP: max iterations (%d) reached.\n', maxiter)
end

RFP = E*diag(lambda)*E';