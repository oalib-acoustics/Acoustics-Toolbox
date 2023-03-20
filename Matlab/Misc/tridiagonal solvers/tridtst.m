% MATH5315: High Performance Numerical Computing
% File:     tridtst.m
% Type:     Matlab script
% Requires: Matlab function trid.m
% Usage:    tridtst
% Purpose:  Test solution of tridiagonal linear systems

format compact

n = input('Size of tridiagonal system: n = ? ');

a = -ones(n,1);
b = 2*ones(n,1);
c = a;
A = diag(a(2:n),-1) + diag(b,0) + diag(c(1:n-1),1);
Asp = spdiags([a b a], [-1:1], n, n);
d = ones(n,1);

disp('Matlab dense matrix A, x = A\d');
tic;
x = A \ d;
toc

disp('Function M-file trid to solve A x = d');
tic
xm = trid(a, b, c, d);
toc
errm = norm(xm-x)

disp('Matlab sparse matrix A, x = A \ d');
tic
xs = Asp \ d;
toc
errt = norm(xs-x)

% Estimate condition number
Acond = condest(A)
