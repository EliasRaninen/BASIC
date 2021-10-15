%% Reproduces Fig 1. from the paper "Bias adjusted sign covariance matrix",
% Elias Raninen and Esa Ollila, 2021.

clear; close all; clc;
%% Covariance model

%% DOA setup 1
% p = 20;
% theta = 15;
% K = numel(theta);
% A = exp(-1j*pi*sind(theta).*(0:p-1).');
% P = 2;
% I = eye(p);
% s2 = 1;
% M = A*P*A' + s2*I;

%% DOA setup 2
p = 40;
theta = [-5 0 5];
K = numel(theta);
A = exp(-1j*pi*sind(theta).*(0:p-1).');
P = diag([2 1 0.5]);
I = eye(p);
s2 = 1;
M = A*P*A' + s2*I;

%% Normalize covariance matrix to shape and compute eigenvalues
M = p*M/trace(M);
eigshape = svd(M);

%% Theoretical eigenvalues of the SSCM
fcomplex = @(t) eigshape./(((1-t)+t*eigshape) .* prod( (1-t) + t*eigshape )) .* (1-t)^(p-1);
eigSSCM  = integral(fcomplex,0,1,'ArrayValued',true);
eigSSCMshape = eigSSCM*p;

%% read table for BASIC
[lambdas,deltas] = BASICtable(p,'complex');

%% Compute BASIC eigenvalues (bias corrected SSCM)
eigBASICnotnormalized = interp1(deltas,lambdas,eigSSCM,'linear');
eigBASIC = eigBASICnotnormalized/sum(eigBASICnotnormalized)*p;

%% Plot
figure(1); clf;
semilogy(1:p,eigshape,'ks','MarkerSize',13,'displayname','\lambda(\Lambda) (shape matrix)'); hold on;
semilogy(1:p,eigSSCMshape,'ro','MarkerFaceColor','r', 'displayname','\lambda(\Lambda_{sgn}) (SSCM shape matrix)')
semilogy(1:p,eigBASIC,'bo','MarkerFaceColor','b','displayname','\lambda(\Lambda_{BASIC}) (bias corrected)')

legend
axis([0.5 5.5 0 25])
xticks(1:5)
