% MATH5315: High Performance Numerical Computing
% File:     trid.m
% Type:     Matlab script
% Requires: -
% Usage:    cycred
% Purpose:  Illustration of cyclic reduction to solve a tridaigonal system
% Example:  -

format compact

% Size of system n = 2^p - 1 for some p
n = 15

% Simple tridiagonal matrix which is symmetric positive definite
% Also try making diagonal elements b(i) = 2 (central difference approx u''(x))
a = -ones(n,1);
b = 4*ones(n,1);
c = -ones(n,1);

% Setup sparse diagonal matrix
A = spdiags([a b c], [-1:1], n, n);
figure(1)
spy(A)

% RHS corresponding to an exact solutions XS
XS = [ones(n,1), [1:n]'];
B = A*XS;
R = [A B];
% Number of equations plus number of RHS
nb = size(B,2)
nt = n + nb;
rhs = [n+1:nt];
e = ones(1, nt);

figure(1)
spy(R)
title('Full sparse data structure R, including RHS');
disp('Press any key to continue');
pause

% Cyclic reduction: Elimination steps, including all RHS
k = 2
while k < n

   I = [k:k:n]
   
   % Practical implementation would avoid temporary storage
   % of aa, bb, cc
   aa = [zeros(k/2,1); diag(R,-k/2)];
   bb = diag(R);
   cc = [diag(R,k/2); zeros(k/2,1)];
   mm = aa(I)./bb(I);
   mp = cc(I)./bb(I);
   R(I,:) = R(I,:) - mm(:,e).*R(I-k/2,:);
   R(I,:) = R(I,:) - mp(:,e).*R(I+k/2,:);
   
   figure(1)
   spy(R)
   title('Full sparse data structure R, including RHS');

   figure(2)
   spy(R(I,[I rhs]))
   title(['Reduced data structure, including RHS: k = ' num2str(k)]);
   
   disp('Press any key to continue');
   pause
     
   k = 2*k

end;
k = k/2

% Cyclic reduction: Substitution steps to find solution
X = zeros(n,nb);
ee = ones(1,nb);

X(I,:) = R(I,rhs) ./ (R(I,I)*ee)
while k > 1
   
   Ip = [k:k:n]
   k = k/2;
   I = [k:2*k:n]
   
   figure(2)
   spy(R([k:k:n],[k:k:n rhs]))
   title(['Data structure for substitution: k = ' num2str(k)]);
   
   disp('Press any key to continue');
   pause
   
   % Should exploit the fact that in any row of R(I, Ip) at most two
   % non-rhs elements are nonzero to avoid matrix multiply R(I,Ip)*X(Ip,:)
   X(I,:) = (R(I,rhs) - R(I,Ip)*X(Ip,:) ) ./ (diag(R(I,I))*ee)

end;

Errchk = norm(X-XS,inf)
