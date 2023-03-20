function x = tridiag( a, b, c, y )
%      x = tridiag( a, b, c, y )
%
%  Solve the  N x N  tridiagonal system for x:
%
%  [ b(1)  c(1)                                  ] [  x(1)  ]   [  y(1)  ] 
%  [ a(1)  b(2)  c(2)                            ] [  x(2)  ]   [  y(2)  ] 
%  [       a(2)  b(3)  c(3)                      ] [        ]   [        ] 
%  [            ...   ...   ...                  ] [  ...   ] = [  ...   ]
%  [                    ...    ...    ...        ] [        ]   [        ]
%  [                        a(N-2) b(N-1) c(N-1) ] [ x(N-1) ]   [ y(N-1) ]
%  [                               a(N-1)  b(N)  ] [  x(N)  ]   [  y(N)  ]
%
%  y must be a vector (row or column); N is determined from its length.
%  a, b, c must be vectors of lengths at least N-1, N, and N-1 respectively,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Check that the input arrays have acceptable sizes

if min(size(y))~=1 | min(size(a))~=1 | min(size(b))~=1 | min(size(c))~=1
   error('a, b, c, y must be vectors');
end
N = length(y);
if length(a)<N-1 | length(b)<N | length(c)<N-1
   error('a, b, c must be vectors of length at least N-1, N, N-1');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Solve the problem by back substitution

gam = zeros(1,N);   %  hold the LU decomposition
x = zeros(size(y));

%  phase 1:  LU decomposition
beta = b(1);
if beta==0
   error('beta = 0 at j=1:  matrix is singular');
end

x(1) = y(1) / beta;
for j=2:N
   gam(j) = c(j-1) / beta;
   beta = b(j) - a(j-1)*gam(j);
   if beta==0
      error( [ 'beta = 0 at j=' num2str(j) ':  matrix is singular' ]);
   end
   x(j) = ( y(j) - a(j-1)*x(j-1) )/beta;
end

%  phase 2:  back-substitution
for j=N-1:-1:1
   x(j) = x(j) - gam(j+1)*x(j+1);
end
