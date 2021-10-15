function [lambdas,deltas] = BASICtable(p,field,points)
%% BASICtable reads or creates a lookup table for BASIC (Bias Adjusted SIgn
% Covariance). If a table already exists, it will be loaded. Otherwise a
% new table will be created.
%
% When creating a new table, the function computes approximate eigenvalues
% of the SSCM corresponding to possible eigenvalues of the shape matrix,
% which are defined via a uniform grid going from 0 to max(lambdas).
%
% Usage: [lambdas,deltas] = BASICtable(p,field,points)
%        [lambdas,deltas] = BASICtable(p,field)          (uses 5000 points)
%
% Inputs:
%           p       -   dimension of samples.
%           field   -   for real-valued data 'real' and for complex-valued
%                       'complex'.
%           points  -   length of table (default 5000).
%
% Outputs:  lambdas -   eigenvalues of shape matrix.
%           deltas  -   eigenvalues of population SSCM.
%
% After the table has been created, BASICS.m can be called using
% [BEST, al] = BASICS(X,lambdas,deltas), where X is the nxp data matrix.
% See more information by typing: help BASICS.
%
% By Elias Raninen 2021
%
% version 1.0 (Oct. 13, 2021)

if ~exist('points','var')
    points = 5000;
end
assert(mod(p,1)==0);
assert(mod(points,1)==0);

% real or complex data
if strcmp(field,'real')
    BASICint = @(t,lam) (1/2)*(lam ./ (1-t+t*lam)) .* (1-t).^(p/2-1); % BASIC integral
    fname = ['BASICtable' '-real' '-p' num2str(p) '-points' num2str(points) '.dat'];
elseif strcmp(field,'complex')
    BASICint = @(t,lam) (lam ./ (1-t+t*lam)) .* (1-t).^(p-1); % BASIC integral
    fname = ['BASICtable' '-complex' '-p' num2str(p) '-points' num2str(points) '.dat'];
else
    error('BASICtable.m: Something wrong with inputs. Write help BASICtable.')
end

try % read table
    tmp     = readtable(fname);
    tmp     = table2array(tmp);
    deltas  = tmp(:,2);
    lambdas = tmp(:,1);
    if strcmp(field,'real')
        fprintf(['The table: ' fname ' has been read.\nIt can be used in BASICS.m for real-valued data of dimension = %d.\n'], p)
    elseif strcmp(field,'complex')
        fprintf(['The table: ' fname ' has been read.\nIt can be used in BASICS.m for complex-valued data of dimension = %d.\n'], p)
    end
catch % write table
    fprintf(['Creating new table: ' fname '\nPlease wait...'])
    maxdelta  = 0;
    maxlambda = p;
    while maxdelta < 1 % ensure that max(lambdas) is at least 1
        maxdelta = integral(@(t) BASICint(t,maxlambda), 0, 1, 'ArrayValued', true);
        maxlambda = maxlambda + p/10;
    end
    lambdas = linspace(0,maxlambda,points).';
%     deltas = lambdas/p .* hypergeom([1,1],be*p+1,1-lambdas); % much slower
    deltas = integral(@(t) BASICint(t,lambdas), 0, 1, 'ArrayValued', true);

    T = table(lambdas,deltas);
    writetable(T,fname);
    fprintf('Done.\n');
end
